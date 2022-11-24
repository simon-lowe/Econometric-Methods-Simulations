rm(list = ls())

# Load packages -----------------------------------------------------------

list.of.packages <- c("data.table", "tidyverse", "janitor", "lubridate",
                      "reshape", "lattice", "gridExtra", "mvtnorm", "ggplot2", 
                      "futile.logger", "gtools",
                      "lfe", "fixest", "reshape2", "igraph", "ivreg",
                      "microbenchmark", "tictoc"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

invisible(lapply(list.of.packages, library, character.only = TRUE))
rm(list.of.packages, new.packages)

# Parameters --------------------------------------------------------------

p <- list()
p$nk = 30  # firm types
p$nl = 10  # worker types

p$alpha_sd = 1 # SD of worker FEs
p$psi_sd   = 1 # SD of firm FEs

# Drawing FEs
p$psi   = with(p, qnorm(1:nk/(nk+1)) * psi_sd)
p$alpha = with(p, qnorm(1:nl/(nl+1)) * alpha_sd)

# let's assume moving probability is fixed
p$lambda = 0.05

p$csort = 0.5 # sorting effect
p$cnetw = 0.0 # network effect
p$csig  = 0.5 # unobserved variance of transition probability beyond FEs

p$w_sigma = 0.8 # wage variance

p$nt = 10 # Number of time periods
p$ni = 50000 # Number of individuals

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
# G <- getG(p)

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
# H <- getH(p,G)

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
      # network[i,t] = sample.int(p$nk, 1, prob = G[l, network[i,t-1],]) # transition P given by worker type and previous firm
      # spellcount[i,t] = ifelse(network[i,t] != network[i,t-1], spellcount[i,t-1] + 1,spellcount[i,t-1])
    }
  }
  
  data  = data.table(melt(network, c('i','t')))
  data2 = data.table(melt(spellcount, c('i','t')))
  setnames(data, "value", "k")
  data[, spell := data2$value]
  data[, l := A[i], i]              # Assign for worker type to each individual
  data[, alpha := p$alpha[l], l]    # Assign worker type FE
  data[, psi := p$psi[k], k]        # Assign firm type FE
}

# data <- sim(p,G,H)

# Assign identities to firms ----------------------------------------------

addSpells <- function(p, dat){
  firm_size = 20
  
  dspell <- dat[,list(len=.N), list(i, spell, k)]
  dspell <- dspell[,fid := sample(1:pmax(1, sum(len)/(firm_size*p$nt)), .N, replace=TRUE) , k]
  dspell <- dspell[,fid := .GRP, list(k,fid)]
  
  setkey(dat,i,spell)
  setkey(dspell,i,spell)
  
  
  dat[, fid:= dspell[dat,fid]]
}

# addSpells(p,data)  # adds by reference to the same data.table object (no copy needed)

# Simulating wages --------------------------------------------------------

# sigma <- diag(p$nt)
# 
# for(i in 1:p$nt){
#   for(j in 1:p$nt){
#     sigma[i, j] <- ifelse(abs(i - j) == 1, p$auto_rho, ifelse(i == j, p$w_sigma, 0))
#   }
# }

addWage <- function(p,data){
  data[, lw := alpha + psi + p$w_sigma * rnorm(.N)]
  # data[, lw := alpha + psi + rmvnorm(1, rep(0, p$nt), sigma)[ ,1:p$nt], by = i]
}
# addWage(p,data)


# Restricting to the connected set ----------------------------------------

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

# There is better way from CASD code

# Build the data ----------------------------------------------------------

buildData <- function(p){
  G <- getG(p)
  H <- getH(p, G)
  data <- sim(p, G, H)
  addSpells(p,data)
  addWage(p,data)
  return(data)
}


# Creating the data -------------------------------------------------------

data <- buildData(p)

# Event study graph -------------------------------------------------------

data[, mcww := sapply(i, function(x) mean(lw[-match(x,i)])) , by= .(fid, t)]
# data[, mcww3 := (sum(lw) - lw) / (.N - 1), by = .(fid, t)]

data[ , mcww_q := cut(mcww, quantile(mcww, probs = 0:4/4),
                      labels = FALSE, include.lowest = TRUE), by = t]

data[, event := spell != shift(spell, 1, fill = 0), by = i][, event_sum := cumsum(event), by = i][, event_max := max(event_sum), by = i]

tmp <- data[event == TRUE, .(i, t, t_d = t-2, t_u = t+2)][, event2 := (t_d >= 0 & t_d >= shift(t_u)) &
                                                            (t_u <= 10 & t_u <= shift(t_d, -1)), by = i]

tmp <- data[event == TRUE, .(i, t, t_d = t-2, t_u = t+2)][
    , event2 := (t_d >= shift(t_u, fill = 1)) & (t_u <= shift(t_d, -1, fill = 11)), by = i][
      event2 == TRUE][
        , eve_counter := seq_len(.N), by = i]

# data <- merge(data, tmp, by = c("i", "t"), all = TRUE)

data[tmp[, .(i, t_d, t_u = t_u - 1, eve_counter)], on = .(i, t >= t_d, t <= t_u), bla := eve_counter]

data[, rel_t := seq_len(.N) - (2 + 1), by = .(i, bla)][is.na(bla), rel_t := NA]

to_plot <- data[!is.na(rel_t)][, keep := mcww_q[rel_t == -1] == 4 | mcww_q[rel_t == -1] == 1, by = .(i, bla)][
  keep == 1][
  , tg := paste0(mcww_q[rel_t == -1], " to ", mcww_q[rel_t == 0]), by = .(i, bla)]

fin <- to_plot[, .(m_wage = mean(lw)), by = .(rel_t, tg)]

p <- fin %>%
  ggplot(aes(x = rel_t, y = m_wage, color = tg)) +
  geom_line() + geom_point()
p


# Run regression ----------------------------------------------------------

data2 <- concomp2(data)

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

setDT(data2)
data2[, residuals := est2$residuals]

# Variance analysis -------------------------------------------------------

var(data2$lw)

var(data2$alpha)
var(data2$fe.i)

var(data2$psi)
var(data2$fe.fid)

2*cov(data2$alpha, data2$psi)
2*cov(data2$fe.i, data2$fe.fid)

var(data2$residuals)

var(data2$fe.i) + var(data2$fe.fid) + 2*cov(data2$fe.i, data2$fe.fid)  + var(data2$residuals) - var(data2$lw)

# Card et al residual graph -----------------------------------------------

data2[ , fe.i_q := cut(fe.i, quantile(fe.i, probs = 0:10/10), labels = FALSE, include.lowest = TRUE)]
data2[ , fe.fid_q := cut(fe.fid, quantile(fe.fid, probs = 0:10/10), labels = FALSE, include.lowest = TRUE)]

res_mean <- data2[, .(mean_res = mean(residuals)), by = .(fe.i_q, fe.fid_q)]

p <- ggplot(res_mean, aes(fe.i_q, fe.fid_q, fill = mean_res)) +
  geom_tile()
p
