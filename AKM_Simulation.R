rm(list = ls())

# Load packages -----------------------------------------------------------

list.of.packages <- c("data.table", "tidyverse", "janitor", "lubridate",
                      "reshape", "lattice", "gridExtra", "mvtnorm", "ggplot2", 
                      "futile.logger", "gtools",
                      "lfe", "fixest", "reshape2", "igraph", "ivreg",
                      "microbenchmark"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

invisible(lapply(list.of.packages, library, character.only = TRUE))
rm(list.of.packages, new.packages)

# Parameters --------------------------------------------------------------

p <- list()
# p$nk = 30  # firm types
p$nk = 4  # firm types
# p$nl = 10  # worker types
p$nl = 2  # worker types

p$alpha_sd = 1
p$psi_sd   = 1

# let's draw some FE
p$psi   = with(p, qnorm(1:nk/(nk+1)) * psi_sd)
p$alpha = with(p, qnorm(1:nl/(nl+1)) * alpha_sd)

# let's assume moving probability is fixed
p$lambda = 0.05

p$csort = 0 # 0.5 # sorting effect
p$cnetw = 0 # 0.2 # network effect
p$csig  = 0.1 # unobserved variance of transition probability beyond FEs

p$w_sigma = 0.8 # wage variance

p$nt = 10 # Number of time periods
p$ni = 1000 # Number of individuals
# p$ni = 300000 # Number of individuals

p$auto_rho <- 0.1

set.seed(12345)

# Creating the transition matrix ------------------------------------------


# lets create type specific transition matrices
# we are going to use joint normal centered on different values
# G[i,j,k] = Pr[worker i, at firm j, moves to firm k]

getG <- function(p){
  G = with(p, array(0, c(nl, nk, nk)))
  for (l in 1:p$nl) for (j in 1:p$nk) {
    # prob of moving is highest if dnorm(0)
    G[l, j,] = with(p, dnorm(psi - cnetw*psi[j] - csort*alpha[l], sd = csig))
    # normalize to get transition matrix
    G[l, j,] = G[l, j,]/sum(G[l, j,])
  }
  return(G)
}
G <- getG(p)

getH <- function(p,G){
  # we then solve for the stationary distribution over psis for each alpha value
  H = with(p,array(1/nk,c(nl,nk)))
  for (l in 1:p$nl) {
    M = G[l,,]
    for (i in 1:100) {
      H[l,] = t(G[l,,]) %*% H[l,]
    }
  }
  return(H)
}
H <- getH(p,G)


Plot1 = wireframe(G[1,,],aspect = c(1,1),xlab = "previous firm",ylab="next firm")
Plot2 = wireframe(G[p$nl,,],aspect = c(1,1),xlab = "previous firm",ylab="next firm")
grid.arrange(Plot1, Plot2,nrow=1)

wireframe(H,aspect = c(1,1), xlab = "worker", ylab = "firm")


# Simulate the network ----------------------------------------------------

sim <- function(p,G,H){
  
  # we simulate a panel
  network    = array(0,c(p$ni,p$nt))
  spellcount = array(0,c(p$ni,p$nt))
  A = rep(0,p$ni)
  
  for(i in 1:p$ni) {
    # we draw the worker type
    l = sample.int(p$nl, 1)
    A[i] = l
    # at time 1, we draw from H
    network[i,1] = sample.int(p$nk, 1, prob = H[l,]) # We draw firm type for individual i from given the stat dist for her type
    for (t in 2:p$nt) {
      if (runif(1) < p$lambda) { # If there is a transition
        network[i,t] = sample.int(p$nk, 1, prob = G[l, network[i,t-1],]) # transition P given by worker type and previous firm
        spellcount[i,t] = spellcount[i,t-1] + 1
      } else { # If no transition
        network[i,t]    = network[i,t-1]
        spellcount[i,t] = spellcount[i,t-1]
      }
    }
  }
  
  data  = data.table(melt(network, c('i','t')))
  data2 = data.table(melt(spellcount, c('i','t')))
  setnames(data,"value","k")
  data[,spell := data2$value]
  data[,l := A[i],i]              # Assign for worker type to each individual
  data[,alpha := p$alpha[l],l]    # Assign worker type FE
  data[,psi := p$psi[k],k]        # Assign firm type FE
}

data <- sim(p,G,H)


# Assign identities to firms ----------------------------------------------

addSpells <- function(p, dat){
  firm_size = 10
  f_class_count = p$ni/(firm_size*p$nk*p$nt)
  
  dspell <- dat[, list(len = .N), list(i, spell, k)]
  dspell[, fid := sample(1:pmax(1, sum(len)/f_class_count), .N, replace = TRUE), k]
  dspell[, fid := .GRP, list(k, fid)]
  
  setkey(dat, i, spell)
  setkey(dspell, i, spell)
  
  dat[, fid:= dspell[dat, fid]]
}

addSpells(p,data)  # adds by reference to the same data.table object (no copy needed)


# Mean firm size in cross-section
# data %>%
#   group_by(t, fid) %>%
#   summarise(
#     f_size = n()
#   ) %>% pull(f_size) %>% mean()



# Simulating wages --------------------------------------------------------

sigma <- diag(p$nt)

for(i in 1:p$nt){
  for(j in 1:p$nt){
    sigma[i, j] <- ifelse(abs(i - j) == 1, p$auto_rho, ifelse(i == j, p$w_sigma, 0))
  }
}

addWage <- function(p,data){
  data[, lw := alpha + psi + p$w_sigma * rnorm(.N)]
  # data[, lw := alpha + psi + rmvnorm(1, rep(0, p$nt), sigma)[ ,1:p$nt], by = i]
}
addWage(p,data)


# Event study graph -------------------------------------------------------

# es_plot <- data %>%
#   group_by(fid) %>%
#   mutate(
#     mw = mean(lw)
#   ) %>%
#   ungroup() %>%
#   mutate(
#     mw_q = quantcut(mw, q = 4, labels = FALSE)
#   ) %>%
#   group_by(i) %>%
#   arrange(i, t, .by_group = TRUE) %>%
#   mutate(
#     move = (spell != lag(spell, order_by = t))
#   ) %>%
#   mutate(
#     n_move = sum(move, na.rm = T)
#   ) %>%
#   filter(n_move == 1) 
# test <- es_plot %>%
#   mutate(
#     move = ifelse(is.na(move), FALSE, move),
#     t_move = t[move == TRUE],
#     rel_t = t - t_move
#   ) %>% 
#   filter(length(unique(mw_q)) == 2) %>%
#   mutate(
#     fromto = factor(paste(unique(mw_q)[1], "to", unique(mw_q)[2]))
#   ) %>%
#   group_by(rel_t, fromto) %>%
#   summarise(
#     m_wage = mean(lw)
#   )
# 
# p <- test %>%
#   ggplot(aes(x = rel_t, y = m_wage, color = fromto)) +
#   geom_line() + geom_point()
# p


# Restricting to the connected set ----------------------------------------

concomp <- function(data){
  setkey(data,i,t)
  data[ ,fid.l1 := data[J(i,t-1),fid]]
  data <- data[complete.cases(data)]
  cf = lfe::compfactor(list(f1=data[,factor(fid)],f2=data[,factor(fid.l1)]))
  fr = data.frame(f1=data[,factor(fid)],f2=data[,factor(fid.l1)],cf)
  data = data[fr$cf==1]
  return(data)
}

concomp2 <- function(data){
  data <- as.data.table(data)
  setkey(data,i,t)
  data[ ,fid.l1 := data[J(i,t-1),fid]]
  foo <- data[fid != fid.l1][,.(fid, fid.l1)]
  
  # foo <- distinct(foo)
  foo <- graph.data.frame(foo, directed=FALSE)
  comp <- components(foo)$membership
  
  data <- data %>%
    left_join(
      data.frame(
        fid = as.numeric(names(comp)),
        cs = comp
      ), by = "fid"
    ) 
  data <- data %>%
    filter(cs == 1) %>%
    select(-cs)
  
  return(data)
}

# data2 <- concomp2(data) 
# 
# View(data2 %>% filter(is.na(cs)))
# View(data2 %>% filter(fid == 2383))

# Build the data ----------------------------------------------------------

buildData <- function(p){
  G <- getG(p)
  H <- getH(p, G)
  data <- sim(p, G, H)
  addSpells(p,data)
  addWage(p,data)
  # data <- concomp(data)
  return(data)
}

# Running the regression --------------------------------------------------

# LFE

# est <- felm(lw ~ 1 | i + fid, data = data2)
# summary(est)



data <- buildData(p)
data2 <- concomp2(data)

# Fixest ------------------------------------------------------------------

est2 <- feols(lw ~ 1 | i + fid, data = data2)

attach_fe <- function(data, reg){
  fe <- fixef(reg)
  
  f_i <- rownames_to_column(data.frame(fe$i), "i") %>% mutate(i = as.numeric(i))
  f_fid <- rownames_to_column(data.frame(fe$fid), "fid") %>% mutate(fid = as.numeric(fid))
  
  return(data %>%
    left_join(f_i, by = "i") %>%
    left_join(f_fid, by = "fid"))
}

data2 <- attach_fe(data2, est2)

var(data2$alpha)
var(data2$f_est2.i)

var(data2$psi)
var(data2$f_est2.fid)

cov(data2$alpha, data2$psi)
cov(data2$f_est2.i, data2$f_est2.fid)

cor(data2$alpha, data2$psi)
cor(data2$f_est2.i, data2$f_est2.fid)

# fevcov(est2)

# write.csv(data2, "data.csv")


# Lfe ---------------------------------------------------------------------

est1 <- felm(lw ~ 1 | i + fid, data = data2)

fe <- getfe(est1) %>% mutate(idx = as.numeric(as.character(idx)))

fe_i <- fe %>% filter(fe == "i") %>% select(fe_i_felm = effect, i = idx)
fe_fid <- fe %>% filter(fe == "fid") %>% select(fe_fid_felm = effect, fid = idx)

data2 <- data2 %>%
  left_join(
    fe_fid, by = "fid" 
  ) %>%
  left_join(
    fe_i, by = "i" 
  )


var(data2$psi)
var(data2$f_est2.fid)
var(data2$fe_fid_felm)

var(data2$alpha)
var(data2$f_est2.i)
var(data2$fe_i_felm)

homo_cor <- fevcov(est1)
homo_cor
cov(data2$alpha, data2$psi)
cov(data2$f_est2.i, data2$f_est2.fid)
cov(data2$fe_i_felm, data2$fe_fid_felm)

data2 <- data2 %>%
  mutate(
    dif_fe_fid = f_est2.fid - fe_fid_felm
  )

summary(data2$dif_fe_fid)

# Limited mobility bias ---------------------------------------------------

# Run simulation
p$nlambda = 10

pdat = data.frame(lambda=seq(from=0.03,0.5,length.out = p$nlambda),varFid = 0, cov_alpha_psi = 0,varFid_true = 0, cov_alpha_psi_true = 0)

for (i in 1:nrow(pdat)){
  p$lambda = pdat[i,"lambda"]
  flog.info("---> doing for lambda = %f",p$lambda)
  data = buildData(p)
  
  est <- feols(lw ~ 1 | i + fid, data = data)
  f_est <- fixef(est)
  
  f_i <- rownames_to_column(data.frame(f_est$i), "i") %>% mutate(i = as.numeric(i))
  f_fid <- rownames_to_column(data.frame(f_est$fid), "fid") %>% mutate(fid = as.numeric(fid))
  
  data <- data %>%
    left_join(f_i, by = "i") %>%
    left_join(f_fid, by = "fid")
  
  pdat[i,"varFid"] <- var(data$f_est.fid)
  pdat[i,"cov_alpha_psi"] <- cov(data$f_est.i, data$f_est.fid)
  pdat[i,"varFid_true"] <- var(data$psi)
  pdat[i,"cov_alpha_psi_true"] <- var(data$alpha)
}

p1=ggplot(data=pdat,mapping=aes(x=lambda,y=varFid)) + geom_line() + ggtitle("var(psi)")
p2=ggplot(data=pdat,mapping=aes(x=lambda,y=cov_alpha_psi)) + geom_line()+ ggtitle("cov(alpha,psi)")
grid.arrange(p1, p2,nrow=1)

p1=ggplot(data=pdat,mapping=aes(x=lambda,y=varFid_true)) + geom_line() + ggtitle("var(psi)")
p2=ggplot(data=pdat,mapping=aes(x=lambda,y=cov_alpha_psi_true)) + geom_line()+ ggtitle("cov(alpha,psi)")
grid.arrange(p1, p2,nrow=1)


p1=ggplot(data=pdat,mapping=aes(x=lambda,y=varFid-varFid_true)) + geom_line() + ggtitle("var(psi)")
p2=ggplot(data=pdat,mapping=aes(x=lambda,y=cov_alpha_psi-cov_alpha_psi_true)) + geom_line()+ ggtitle("cov(alpha,psi)")
grid.arrange(p1, p2,nrow=1)



# Alleviating the bias using Split Sample Jackknife -----------------------

n_MC <- 100

res <- data.frame()

for(i in 1:n_MC){
  print(i)
  data <- buildData(p)
  data2 <- concomp2(data)
  
  v_fid_true <- var(data2$psi)
  v_i_true <- var(data2$alpha)
  cov_true <- cov(data2$psi, data2$alpha)
  
  est_full <- feols(lw ~ 1 | i + fid, data = data2)
  
  data2 <- attach_fe(data2, est_full)
  
  v_fid_full <- var(data2$fe.fid)
  v_i_full <- var(data2$fe.i)
  cov_full <- cov(data2$fe.i, data2$fe.fid)
  
  # data2 <- data2 %>%
  #   group_by(i) %>%
  #   mutate(
  #     mover = (sum(spell) > 0)
  #   ) %>%
  #   ungroup()
  
  data2 <- data2 %>%
    group_by(i) %>%
    mutate(
      mover = (sum(spell) > 0)
    ) %>%
    ungroup()
  
  tmp <- data2 %>%
    filter(t == 1) %>%
    group_by(fid, mover) %>%
    mutate(
      z = rand_n_list2(n(), 2)
    ) %>%
    ungroup() %>%
    select(i, z)

  data2 <- data2 %>%
    left_join(
      tmp, by = "i"
    )
  
  # data2 <- data2 %>%
  #     group_by(fid, mover) %>%
  #     mutate(
  #       z = rand_n_list2(n(), 2)
  #     ) %>%
  #     ungroup() %>%
  #     select(-fe.i, -fe.fid)  
  
  # data <- data %>%
  #   filter(mover == 1) %>%
  #     group_by(fid) %>%
  #     mutate(
  #       z = rand_n_list2(n(), 2)
  #     ) %>%
  #     ungroup()
  
  data_A <- data2 %>% filter(z == 1)
  data_A <- concomp2(data_A)
  data_B <- data2 %>% filter(z == 2)
  data_B <- concomp2(data_B)
  
  est_A <- tc_NA(feols(lw ~ 1 | i + fid, data = data_A))
  est_B <- tc_NA(feols(lw ~ 1 | i + fid, data = data_B))
  
  # data_A <- tc_NA(attach_fe(data_A, est_A))
  # data_B <- tc_NA(attach_fe(data_B, est_B))
  
  fe_A <- fixef(est_A)
  
  f_i <- rownames_to_column(data.frame(fe_A$i), "i") %>% mutate(i = as.numeric(i))
  f_fid <- rownames_to_column(data.frame(fe_A$fid), "fid") %>% mutate(fid = as.numeric(fid))
  
  data2 <- data2 %>% 
    left_join(f_i, by = "i") %>%
    left_join(f_fid, by = "fid")
  
  fe_B <- fixef(est_B)
  
  f_i <- rownames_to_column(data.frame(fe_B$i), "i") %>% mutate(i = as.numeric(i))
  f_fid <- rownames_to_column(data.frame(fe_B$fid), "fid") %>% mutate(fid = as.numeric(fid))
  
  data2 <- data2 %>% 
    left_join(f_i, by = "i") %>%
    left_join(f_fid, by = "fid")
  
  # tmp <- data2 %>%
  #   filter(mover == 1) %>%
  #   group_by(fid) %>%
  #   summarise(
  #     Mj = length(unique(i))
  #   )
  # 
  # data2 <- data2 %>%
  #   left_join(
  #     tmp, by = "fid"
  #   ) %>%
  #   mutate(
  #     Mj_inv = 1/Mj
  #   )
  # 
  # sigma_nu <- 1/mean(data2$Mj_inv, na.rm = T)*(v_fid_full - cov(data2$fe_A.fid, data2$fe_B.fid, use = "complete.obs"))
  # 
  # data2 <- data2 %>%
  #   mutate(
  #     sigmaj = Mj_inv*sigma_nu,
  #     fe.fid.shrink = v_fid_full/(v_fid_full + sigmaj)*fe.fid
  #   )
  # 
  # mean((data2$fe.fid - data2$psi)^2)
  # mean((data2$fe.fid.shrink - data2$psi)^2)
  
  
  # v_fid_A <- tc_NA(var(data_A$fe.fid))
  # v_i_A <- tc_NA(var(data_A$fe.i))
  # cov_A <- tc_NA(cov(data_A$fe.i, data_A$fe.fid))
  # 
  # v_fid_B <- tc_NA(var(data_B$fe.fid))
  # v_i_B <- tc_NA(var(data_B$fe.i))
  # cov_B <- tc_NA(cov(data_B$fe.i, data_B$fe.fid))
  
  res <- rbind(res,
               data.frame(
                 v_fid_true = v_fid_true,
                 v_i_true = v_i_true,
                 cov_true = cov_true,                 
                 v_fid_full = v_fid_full,
                 v_i_full = v_i_full,
                 cov_full = cov_full,
                 v_fid_ss = cov(data2$fe_A.fid, data2$fe_B.fid, use = "complete.obs"),
                 cov_ss = cov(data2$fe_A.fid, data2$fe_B.i, use = "complete.obs"),
                 cov_ss_cf = 0.5*(cov(data2$fe_A.fid, data2$fe_B.i, use = "complete.obs") + cov(data2$fe_B.fid, data2$fe_A.i, use = "complete.obs"))
                 # v_fid_A = v_fid_A,
                 # v_i_A = v_i_A,
                 # cov_A = cov_A,
                 # v_fid_B = v_fid_B,
                 # v_i_B = v_i_B,
                 # cov_B = cov_B,
                 # v_fid_SSJK = 2*v_fid_full - 0.5*(v_fid_A + v_fid_B),
                 # v_i_SSJK = 2*v_i_full - 0.5*(v_i_A + v_i_B),
                 # cov_SSJK = 2*cov_full - 0.5*(cov_A + cov_B)
               ))
}

mean(res$v_fid_true)
mean(res$v_fid_full)
mean(res$v_fid_ss)

mean(res$cov_true)
mean(res$cov_full)
mean(res$cov_ss)
mean(res$cov_ss_cf)

with(res, mean(v_fid_true - v_fid_ss))
with(res, mean(cov_true - cov_ss))
with(res, mean(cov_true - cov_ss_cf))

sd(res$v_fid_true)
sd(res$v_fid_full)
sd(res$v_fid_ss)

sd(res$cov_true)
sd(res$cov_full)
sd(res$cov_ss)
sd(res$cov_ss_cf)


# FE as dependent variable ------------------------------------------------

n_MC <- 100

res <- data.frame()

rand_n_list2 <- function(n, n_groups){
  return(sample(rep(seq(1, n_groups), n %/% n_groups + 1, length.out = n), n))
}

tc_NA <- function(x){
  tryCatch(x, error = function(e) {NA})
}

for(i in 1:n_MC){
  print(i)
  data <- buildData(p)
  data2 <- concomp2(data) %>%
    mutate(
      rw = 1 + 1*alpha + 3*psi + rnorm(n(), 0, 2)
      # rw = 1 + 3*psi + rnorm(n(), 0, 2)
    ) %>%
    group_by(i) %>%
    mutate(
      mover = (sum(spell) > 0)
    ) %>%
    ungroup()
  
  # est2 <- tc_NA(feols(lw ~ 1 | i + fid, data = data2))
  # f_est2 <- tc_NA(fixef(est2))
  # 
  # f_i <- tc_NA(rownames_to_column(data.frame(f_est2$i), "i") %>% mutate(i = as.numeric(i)))
  # f_fid <- tc_NA(rownames_to_column(data.frame(f_est2$fid), "fid") %>% mutate(fid = as.numeric(fid)))
  # 
  # data2 <- data2 %>%
  #   left_join(f_i, by = "i") %>%
  #   left_join(f_fid, by = "fid")
  # 
  # cov(data2$psi, data2$alpha)
  
  tmp <- data2 %>%
    filter(t == 1) %>%
    group_by(fid, mover) %>%
    mutate(
      z = rand_n_list2(n(), 2)
    ) %>%
    ungroup() %>%
    select(i, z)
  
  data2 <- data2 %>%
    left_join(
      tmp, by = "i"
    )
  
  data_A <- data2 %>% filter(z == 1)
  data_A <- concomp2(data_A)
  data_B <- data2 %>% filter(z == 2)
  data_B <- concomp2(data_B)
  
  # est2 <- tc_NA(feols(lw ~ 1 | i + fid, data = data_B))
  # f_est2 <- tc_NA(fixef(est2))
  # 
  # f_i <- tc_NA(rownames_to_column(data.frame(f_est2$i), "i") %>% mutate(i = as.numeric(i)))
  # f_fid <- tc_NA(rownames_to_column(data.frame(f_est2$fid), "fid") %>% mutate(fid = as.numeric(fid)))
  # 
  # data_A <- data_A %>%
  #   left_join(f_i, by = "i") %>%
  #   left_join(f_fid, by = "fid")
  # 
  # est2 <- tc_NA(feols(lw ~ 1 | i + fid, data = data_A))
  # f_est2 <- tc_NA(fixef(est2))
  # 
  # f_i <- tc_NA(rownames_to_column(data.frame(f_est2$i), "i") %>% mutate(i = as.numeric(i)))
  # f_fid <- tc_NA(rownames_to_column(data.frame(f_est2$fid), "fid") %>% mutate(fid = as.numeric(fid)))
  # 
  # data_A <- data_A %>%
  #   left_join(f_i, by = "i") %>%
  #   left_join(f_fid, by = "fid")
  
  reg1 <- tc_NA(lm(dat = data_A, rw ~ alpha + psi))
  
  reg2 <- tc_NA(lm(dat = data_A, rw ~  f_est2.fid.y))
  
  reg3 <- tc_NA(lm(dat = data_A, rw ~  f_est2.fid.x))
  
  # reg4 <- tc_NA(ivreg(dat = data_A, rw ~  f_est2.fid.y + f_est2.i.y |  f_est2.fid.x + f_est2.i.y))
  reg4 <- tc_NA(ivreg(dat = data_A, rw ~  f_est2.fid.y |  f_est2.fid.x))
  
  res <- rbind(res, data.frame(
    true_coef_i = tc_NA(coef(summary(reg1))[2,1]),
    true_sd_i = tc_NA(coef(summary(reg1))[2,2]),
    # compA_coef_i = tc_NA(coef(summary(reg2))[2,1]),
    # compA_sd_i = tc_NA(coef(summary(reg2))[2,2]),
    # compB_coef_i = tc_NA(coef(summary(reg3))[2,1]),
    # compB_sd_i = tc_NA(coef(summary(reg3))[2,2]),
    # compIV_coef_i = tc_NA(coef(summary(reg4))[2,1]),
    # compIV_sd_i = tc_NA(coef(summary(reg4))[2,2]),
    true_coef_fid = tc_NA(coef(summary(reg1))[3,1]),
    true_sd_fid = tc_NA(coef(summary(reg1))[3,2]),
    compA_coef_fid = tc_NA(coef(summary(reg2))[2,1]),
    compA_sd_fid = tc_NA(coef(summary(reg2))[2,2]),
    compB_coef_fid = tc_NA(coef(summary(reg3))[2,1]),
    compB_sd_fid = tc_NA(coef(summary(reg3))[2,2]),
    compIV_coef_fid = tc_NA(coef(summary(reg4))[2,1]),
    compIV_sd_fid = tc_NA(coef(summary(reg4))[2,2])
  )
  )
}


mean(res$true_coef_i)
# mean(res$compA_coef_i)
# mean(res$compB_coef_i, na.rm = T)
# mean(res$compIV_coef_i, na.rm = T)
# 
# hist(res$compIV_coef_i, breaks = 20)

mean(res$true_coef_fid)
mean(res$compA_coef_fid)
mean(res$compB_coef_fid, na.rm = T)
mean(res$compIV_coef_fid, na.rm = T)

hist(res$compIV_coef_fid, breaks = 20)
# hist(res$compB_coef_fid, breaks = 20)
# hist(res$compB_coef_i, breaks = 20)

sd(res$true_coef_fid)
# sd(res$compA_coef_i)
# sd(res$compB_coef_i, na.rm = T)
# sd(res$compIV_coef_i, na.rm = T)

# mean(res$compIV_sd_i, na.rm = T)
# median(res$compIV_sd_i, na.rm = T)


sd(res$compIV_coef_fid, na.rm = T)
mean(res$compIV_sd_fid, na.rm = T)

res <- res %>%
  mutate(
    covered = 3 < compIV_coef_fid + 1.96*compIV_sd_fid & 3 > compIV_coef_fid - 1.96*compIV_sd_fid
  )

mean(res$covered, na.rm = T)


# FE as dependent variable using Shrinkage ------------------------------------------------

n_MC <- 100

res <- data.frame()

rand_n_list2 <- function(n, n_groups){
  return(sample(rep(seq(1, n_groups), n %/% n_groups + 1, length.out = n), n))
}

rand_n_list3 <- function(n, n_groups){
  rem <- n %% n_groups
  mult <- n %/% n_groups
  tmp <- sample(rep(seq(1, n_groups), mult), mult*n_groups, replace = FALSE)
  return(c(tmp, sample(seq(1, n_groups), rem)))
}

tc_NA <- function(x){
  tryCatch(x, error = function(e) {NA})
}

for(i in 1:n_MC){
  print(i)
  data <- buildData(p)
  # data2 <- concomp2(data) %>%
  data2 <- data %>%
    mutate(
      rw = 1 + 1*alpha + 3*psi + rnorm(n(), 0, 2)
      # rw = 1 + 3*psi + rnorm(n(), 0, 2)
    ) %>%
    group_by(i) %>%
    mutate(
      mover = (sum(spell) > 0)
    ) %>% ungroup()
  
    zaza <- data2 %>%
    group_by(fid, t) %>%
    summarise(
      size_t = n()
    ) %>%
      group_by(fid) %>%
    summarise(
      m_size = mean(size_t)
    ) %>%
      filter(m_size > 8)
    
  data2 <- data2 %>%
    filter(fid %in% zaza$fid) %>%
    concomp2()
    
  est_full <- feols(lw ~ 1 | i + fid, data = data2)
  
  data2 <- attach_fe(data2, est_full)
  
  tmp <- data2 %>%
    filter(t == 1) %>%
    group_by(fid, mover) %>%
    mutate(
      z = rand_n_list2(n(), 2)
    ) %>%
    ungroup() %>%
    select(i, z)

  data2 <- data2 %>%
    left_join(
      tmp, by = "i"
    )
  
  # data2 <- data2 %>%
  #   group_by(fid, mover) %>%
  #   mutate(
  #     z = rand_n_list2(n(), 2)
  #   )
  
  data_A <- data2 %>% filter(z == 1)
  data_A <- concomp2(data_A)
  data_B <- data2 %>% filter(z == 2)
  data_B <- concomp2(data_B)
  
  est_A <- tc_NA(feols(lw ~ 1 | i + fid, data = data_A))
  est_B <- tc_NA(feols(lw ~ 1 | i + fid, data = data_B))
  
  fe_A <- fixef(est_A)
  
  f_i <- rownames_to_column(data.frame(fe_A$i), "i") %>% mutate(i = as.numeric(i))
  f_fid <- rownames_to_column(data.frame(fe_A$fid), "fid") %>% mutate(fid = as.numeric(fid))
  
  data2 <- data2 %>% 
    left_join(f_i, by = "i") %>%
    left_join(f_fid, by = "fid")
  
  fe_B <- fixef(est_B)
  
  f_i <- rownames_to_column(data.frame(fe_B$i), "i") %>% mutate(i = as.numeric(i))
  f_fid <- rownames_to_column(data.frame(fe_B$fid), "fid") %>% mutate(fid = as.numeric(fid))
  
  data2 <- data2 %>% 
    left_join(f_i, by = "i") %>%
    left_join(f_fid, by = "fid")
  
  tmp <- data2 %>%
    filter(mover == 1) %>%
    group_by(fid) %>%
    summarise(
      Mj = length(unique(i))
    )

  data2 <- data2 %>%
    left_join(
      tmp, by = "fid"
    ) %>%
    mutate(
      Mj_inv = 1/Mj
    )

  v_fid_full <- var(data2$fe.fid)
  
  sigma_nu <- tc_NA(1/mean(data2$Mj_inv, na.rm = T)*(v_fid_full - cov(data2$fe_A.fid, data2$fe_B.fid, use = "complete.obs")))

  data2 <- data2 %>%
    mutate(
      sigmaj = tc_NA(Mj_inv*sigma_nu),
      fe.fid.shrink = tc_NA(v_fid_full/(v_fid_full + sigmaj)*fe.fid),
      fe_A.fid.shrink = tc_NA(v_fid_full/(v_fid_full + sigmaj)*fe_A.fid),
      fe_B.fid.shrink = tc_NA(v_fid_full/(v_fid_full + sigmaj)*fe_B.fid)
    )
  
  reg1 <- tc_NA(lm(dat = data2, rw ~ alpha + psi))
  
  reg2 <- tc_NA(lm(dat = data2, rw ~ fe.fid))
  # 
  # reg2b <- tc_NA(lm(dat = data2, rw ~ fe.fid + fe.i))
  # 
  # reg2c <- tc_NA(lm(dat = data2, rw ~ fe.fid + alpha))
  
  # reg4 <- tc_NA(ivreg(dat = data2 %>% filter(z == 1), rw ~  fe_A.fid | fe_B.fid  | fe_A.i))

  reg4b <- tc_NA(ivreg(dat = data2 %>% filter(z == 1), rw ~  fe_A.fid | fe_B.fid))
  
  reg4d <- tc_NA(ivreg(dat = data2 %>% filter(z == 1), rw ~  fe_A.fid + fe_A.i | fe_B.fid + fe_B.i))
  
  # reg4c <- tc_NA(ivreg(dat = data2 %>% filter(z == 1), rw ~  fe_A.fid.shrink | fe_B.fid.shrink))
  # 
  # reg7 <- tc_NA(ivreg(dat = data2 %>% filter(z == 1), rw ~  fe_A.fid |  fe_B.fid | alpha))
  
  reg5 <- tc_NA(lm(dat = data2, rw ~ fe.fid.shrink))
  
  reg6 <- tc_NA(lm(dat = data2, rw ~ fe.fid.shrink + fe.i))
  
  res <- rbind(res, data.frame(
    true_coef_i = tc_NA(coef(summary(reg1))[2,1]),
    true_sd_i = tc_NA(coef(summary(reg1))[2,2]),
    true_coef_fid = tc_NA(coef(summary(reg1))[3,1]),
    true_sd_fid = tc_NA(coef(summary(reg1))[3,2]),
    comp1_coef_fid = tc_NA(coef(summary(reg2))[2,1]),
    comp1_sd_fid = tc_NA(coef(summary(reg2))[2,2]),
    # comp1b_coef_fid = tc_NA(coef(summary(reg2b))[2,1]),
    # comp1b_sd_fid = tc_NA(coef(summary(reg2b))[2,2]),
    # comp1c_coef_fid = tc_NA(coef(summary(reg2c))[2,1]),
    # comp1c_sd_fid = tc_NA(coef(summary(reg2c))[2,2]),
    # compIV_coef_fid = tc_NA(coef(summary(reg4))[2,1]),
    # compIV_sd_fid = tc_NA(coef(summary(reg4))[2,2]),
    compIVnoC_coef_fid = tc_NA(coef(summary(reg4b))[2,1]),
    compIVnoC_sd_fid = tc_NA(coef(summary(reg4b))[2,2]),
    compIVnoC_dumbsplit_coef_fid = tc_NA(coef(summary(reg4d))[2,1]),
    compIVnoC_dumbsplit_sd_fid = tc_NA(coef(summary(reg4d))[2,2]),
    # compIVnoCshrink_coef_fid = tc_NA(coef(summary(reg4c))[2,1]),
    # compIVnoCshrink_sd_fid = tc_NA(coef(summary(reg4c))[2,2]),
    compSH_coef_fid = tc_NA(coef(summary(reg5))[2,1]),
    compSH_sd_fid = tc_NA(coef(summary(reg5))[2,2]),
    compSHFEi_coef_fid = tc_NA(coef(summary(reg6))[2,1]),
    compSHFEi_sd_fid = tc_NA(coef(summary(reg6))[2,2])
    # compIVTA_coef_fid = tc_NA(coef(summary(reg7))[2,1]),
    # compIVTA_sd_fid = tc_NA(coef(summary(reg7))[2,2]),
    # bias_FFE = mean(data2$fe.fid - data2$psi),
    # bias_FFE_shrink = mean(data2$fe.fid.shrink - data2$psi),
    # v_FFE = var(data2$fe.fid),
    # v_FFE_shrink = var(data2$fe.fid.shrink),
    # MSE_FFE = mean((data2$fe.fid - data2$psi)^2),
    # MSE_FFE_shrink = mean((data2$fe.fid.shrink - data2$psi)^2)
  )
  )
}


# mean(res$true_coef_i)
# mean(res$compA_coef_i)
# mean(res$compB_coef_i, na.rm = T)
# mean(res$compIV_coef_i, na.rm = T)
# 
# hist(res$compIV_coef_i, breaks = 20)

mean(res$true_coef_fid)
# mean(res$comp1_coef_fid)
# mean(res$comp1b_coef_fid, na.rm = T)
# mean(res$comp1c_coef_fid, na.rm = T)
mean(res$compIVnoC_coef_fid, na.rm = T)
mean(res$compIVnoC_dumbsplit_coef_fid, na.rm = T)
# mean(res$compIVnoCshrink_coef_fid, na.rm = T)
# mean(res$compIV_coef_fid, na.rm = T)
# mean(res$compIVTA_coef_fid)
mean(res$compSH_coef_fid, na.rm = T)
mean(res$compSHFEi_coef_fid, na.rm = T)

sd(res$true_coef_fid)
# sd(res$comp1_coef_fid)
# sd(res$comp1b_coef_fid, na.rm = T)
# sd(res$comp1c_coef_fid, na.rm = T)
sd(res$compIVnoC_coef_fid, na.rm = T)
sd(res$compIVnoC_dumbsplit_coef_fid, na.rm = T)
# sd(res$compIVnoCshrink_coef_fid, na.rm = T)
# sd(res$compIV_coef_fid, na.rm = T)
# sd(res$compIVTA_coef_fid)
sd(res$compSH_coef_fid, na.rm = T)
sd(res$compSHFEi_coef_fid, na.rm = T)


sd(res$compIV_coef_fid, na.rm = T)
mean(res$compIV_sd_fid, na.rm = T)
mean(res$compIVnoC_dumbsplit_sd_fid, na.rm = T)

sd(res$compSH_coef_fid)
mean(res$compSH_sd_fid, na.rm = T)
mean(res$compSHFEi_sd_fid, na.rm = T)

res <- res %>%
  mutate(
    covered_IVnoC = 3 < compIVnoC_coef_fid + 1.96*compIVnoC_sd_fid & 3 > compIVnoC_coef_fid - 1.96*compIVnoC_sd_fid,
    covered_IVnoCdumbsplit = 3 < compIVnoC_dumbsplit_coef_fid + 1.96*compIVnoC_dumbsplit_sd_fid & 3 > compIVnoC_dumbsplit_coef_fid - 1.96*compIVnoC_dumbsplit_sd_fid,
    covered_IVnoCdumbsplit2 = 3 < compIVnoC_dumbsplit_coef_fid + 1.96*sd(res$compIVnoC_dumbsplit_coef_fid, na.rm = T) & 3 > compIVnoC_dumbsplit_coef_fid - 1.96*sd(res$compIVnoC_dumbsplit_coef_fid, na.rm = T)
  )

mean(res$covered_IVnoC, na.rm = T)
mean(res$covered_IVnoCdumbsplit, na.rm = T)
mean(res$covered_IVnoCdumbsplit2, na.rm = T)

mean(res$bias_FFE)
mean(res$bias_FFE_shrink)

mean(res$v_FFE)
mean(res$v_FFE_shrink)

mean(res$MSE_FFE)
mean(res$MSE_FFE_shrink)





# MC for shrinkage --------------------------------------------------------



# Testing -----------------------------------------------------------------

f1 <- function(){
datalist = list()
for(i in 1:5){
  dat <- concomp2(buildData(p))
  dat$sim <- i  # maybe you want to keep track of which iteration produced it?
  datalist[[i]] <- dat # add it to your list
}

bdat <- data.table::rbindlist(datalist)
}

f2 <- function(){
  cdat <- function(aux){
    dat <- concomp2(buildData(p))
    dat$sim <- i  # maybe you want to keep track of which iteration produced it?
    return(dat)
  }
  bdat <- map_df(1:5, cdat)
}

f3 <- function(){
  cdat <- function(aux){
    dat <- concomp2(buildData(p))
    dat$sim <- i  # maybe you want to keep track of which iteration produced it?
    return(dat)
  }
  bdat <- map_df(1:5, cdat)
  bdat <- setDT(bdat)
}

microbenchmark(
  test1 = f1(),
  test2 = f2(),
  test3 = f3(),
  times = 2
)

foo <- copy(bdat)

foo[, rw := ]