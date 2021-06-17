rm(list = ls())

# Load packages -----------------------------------------------------------

list.of.packages <- c("data.table", "tidyverse", "janitor", "lubridate",
                      "reshape", "lattice", "gridExtra", "mvtnorm", "ggplot2", 
                      "futile.logger", "gtools",
                      "lfe", "fixest", "reshape2", "igraph", "ivreg"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

invisible(lapply(list.of.packages, library, character.only = TRUE))
rm(list.of.packages, new.packages)

# Parameters --------------------------------------------------------------

p <- list()
p$nk = 30  # firm types
p$nl = 10  # worker types

p$alpha_sd = 1
p$psi_sd   = 1

# let's draw some FE
p$psi   = with(p, qnorm(1:nk/(nk+1)) * psi_sd)
p$alpha = with(p, qnorm(1:nl/(nl+1)) * alpha_sd)

# let's assume moving probability is fixed
p$lambda = 0.05

p$csort = 0.5 # sorting effect
p$cnetw = 0.2 # network effect
p$csig  = 0.5 # unobserved variance of transition probability beyond FEs

p$w_sigma = 0.8 # wage variance

p$nt = 5 # Number of time periods
p$ni = 130000 # Number of individuals

set.seed(12345)

# Creating the transition matrix ------------------------------------------


# lets create type specific transition matrices
# we are going to use joint normal centered on different values
# G[i,j,k] = Pr[worker i, at firm j, moves to firm k]

getG <- function(p){
  G = with(p,array(0,c(nl,nk,nk)))
  for (l in 1:p$nl) for (k in 1:p$nk) {
    # prob of moving is highest if dnorm(0)
    G[l,k,] = with(p, dnorm(psi - cnetw*psi[k] - csort*alpha[l], sd = csig))
    # normalize to get transition matrix
    G[l,k,] = G[l,k,]/sum(G[l,k,])
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
H = getH(p,G)


Plot1 = wireframe(G[1,,],aspect = c(1,1),xlab = "previous firm",ylab="next firm")
Plot2 = wireframe(G[p$nl,,],aspect = c(1,1),xlab = "previous firm",ylab="next firm")
grid.arrange(Plot1, Plot2,nrow=1)

wireframe(H,aspect = c(1,1), xlab = "worker", ylab="firm")


# Simulate the network ----------------------------------------------------

sim <- function(p,G,H){
  
  # we simulate a panel
  network    = array(0,c(p$ni,p$nt))
  spellcount = array(0,c(p$ni,p$nt))
  A = rep(0,p$ni)
  
  for(i in 1:p$ni) {
    # we draw the worker type
    l = sample.int(p$nl,1)
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
data %>%
  group_by(t, fid) %>%
  summarise(
    f_size = n()
  ) %>% pull(f_size) %>% mean()



# Simulating wages --------------------------------------------------------

addWage <- function(p,data){
  data[, lw := alpha + psi + p$w_sigma * rnorm(.N)]
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
  data <- sim(p,G,H)
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
f_est2 <- fixef(est2)

f_i <- rownames_to_column(data.frame(f_est2$i), "i") %>% mutate(i = as.numeric(i))
f_fid <- rownames_to_column(data.frame(f_est2$fid), "fid") %>% mutate(fid = as.numeric(fid))

data2 <- data2 %>%
  left_join(f_i, by = "i") %>%
  left_join(f_fid, by = "fid")

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
    ) %>%
    group_by(i) %>%
    # mutate(
    #   mover = length(unique(fid)) > 1
    # )
    mutate(
      mover = (sum(spell) > 0)
    ) %>%
    group_by(fid, mover) %>%
    mutate(
      z = rand_n_list2(n(), 2)
    )
  data_A <- data2 %>% filter(z == 1)
  data_A <- concomp2(data_A)
  data_B <- data2 %>% filter(z == 2)
  data_B <- concomp2(data_B)
  
  est2 <- tc_NA(feols(lw ~ 1 | i + fid, data = data_B))
  f_est2 <- tc_NA(fixef(est2))
  
  f_i <- tc_NA(rownames_to_column(data.frame(f_est2$i), "i") %>% mutate(i = as.numeric(i)))
  f_fid <- tc_NA(rownames_to_column(data.frame(f_est2$fid), "fid") %>% mutate(fid = as.numeric(fid)))
  
  data_A <- data_A %>%
    left_join(f_i, by = "i") %>%
    left_join(f_fid, by = "fid")
  
  est2 <- tc_NA(feols(lw ~ 1 | i + fid, data = data_A))
  f_est2 <- tc_NA(fixef(est2))
  
  f_i <- tc_NA(rownames_to_column(data.frame(f_est2$i), "i") %>% mutate(i = as.numeric(i)))
  f_fid <- tc_NA(rownames_to_column(data.frame(f_est2$fid), "fid") %>% mutate(fid = as.numeric(fid)))
  
  data_A <- data_A %>%
    left_join(f_i, by = "i") %>%
    left_join(f_fid, by = "fid")
  
  reg1 <- tc_NA(lm(dat = data_A, rw ~ alpha + psi))
  
  reg2 <- tc_NA(lm(dat = data_A, rw ~ f_est2.i.y + f_est2.fid.y))
  
  reg3 <- tc_NA(lm(dat = data_A, rw ~ f_est2.i.x + f_est2.fid.x))
  
  reg4 <- tc_NA(ivreg(dat = data_A, rw ~ f_est2.i.y + f_est2.fid.y | f_est2.i.x + f_est2.fid.x))
  
  res <- rbind(res, data.frame(
    true_coef_i = tc_NA(coef(summary(reg1))[2,1]),
    true_sd_i = tc_NA(coef(summary(reg1))[2,2]),
    compA_coef_i = tc_NA(coef(summary(reg2))[2,1]),
    compA_sd_i = tc_NA(coef(summary(reg2))[2,2]),
    compB_coef_i = tc_NA(coef(summary(reg3))[2,1]),
    compB_sd_i = tc_NA(coef(summary(reg3))[2,2]),
    compIV_coef_i = tc_NA(coef(summary(reg4))[2,1]),
    compIV_sd_i = tc_NA(coef(summary(reg4))[2,2]),
    true_coef_fid = tc_NA(coef(summary(reg1))[3,1]),
    true_sd_fid = tc_NA(coef(summary(reg1))[3,2]),
    compA_coef_fid = tc_NA(coef(summary(reg2))[3,1]),
    compA_sd_fid = tc_NA(coef(summary(reg2))[3,2]),
    compB_coef_fid = tc_NA(coef(summary(reg3))[3,1]),
    compB_sd_fid = tc_NA(coef(summary(reg3))[3,2]),
    compIV_coef_fid = tc_NA(coef(summary(reg4))[3,1]),
    compIV_sd_fid = tc_NA(coef(summary(reg4))[3,2])
  )
  )
}


mean(res$true_coef_i)
mean(res$compA_coef_i)
mean(res$compB_coef_i, na.rm = T)
mean(res$compIV_coef_i, na.rm = T)

hist(res$compIV_coef_i, breaks = 20)

mean(res$true_coef_fid)
mean(res$compA_coef_fid)
mean(res$compB_coef_fid, na.rm = T)
mean(res$compIV_coef_fid, na.rm = T)

hist(res$compIV_coef_fid, breaks = 20)
hist(res$compB_coef_fid, breaks = 20)
hist(res$compB_coef_i, breaks = 20)

sd(res$true_coef_fid)
sd(res$compA_coef_i)
sd(res$compB_coef_i, na.rm = T)
sd(res$compIV_coef_i, na.rm = T)

mean(res$compIV_sd_i, na.rm = T)
median(res$compIV_sd_i, na.rm = T)
