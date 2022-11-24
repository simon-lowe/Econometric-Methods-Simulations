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


# Add gender --------------------------------------------------------------

# tmp <- data.table(i = 1:p$ni, woman = rbinom(p$ni, 1, 0.48))
# 
# data <- merge(data, tmp, by = "i")
# rm(tmp)

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


# Simulating wages --------------------------------------------------------


addWage <- function(p,data){
  data[, lw := alpha + 1*(1 - woman)*psi + 0.8*woman*psi + p$w_sigma * rnorm(.N)]
}
# addWage(p,data)

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

extract_big_comp <- function(data){
  tmp <- data[woman == TRUE, .(i, fid)] %>% unique()
  
  graph <- graph.data.frame(tmp, directed=FALSE)
  comp <- components(graph)$membership
  
  data <- data %>%
    left_join(
      data.frame(
        fid = as.numeric(names(comp)),
        cs = comp
      ), by = "fid"
    )
  setDT(data)
  d1 <- data[, comp_size := .N, by = cs][comp_size == max(comp_size)][, cs := NULL]
  data[, cs := NULL]
  
  tmp <- data[woman == FALSE, .(i, fid)] %>% unique()
  
  graph <- graph.data.frame(tmp, directed=FALSE)
  comp <- components(graph)$membership
  
  data <- data %>%
    left_join(
      data.frame(
        fid = as.numeric(names(comp)),
        cs = comp
      ), by = "fid"
    )
  setDT(data)
  d2 <- data[, comp_size := .N, by = cs][comp_size == max(comp_size)][, cs := NULL]
  return(rbind(d1, d2))
}

# Build the data ----------------------------------------------------------

buildData <- function(p){
  G <- getG(p)
  H <- getH(p, G)
  data <- sim(p, G, H)
  tmp <- data.table(i = 1:p$ni, woman = rbinom(p$ni, 1, 0.48))
  data <- merge(data, tmp, by = "i")
  addSpells(p,data)
  addWage(p,data)
  return(data)
}

# Running the regression --------------------------------------------------


data <- buildData(p)
data2 <- extract_big_comp(data)

data2[, fid.wom := paste(fid, woman, sep = ".")]

est2 <- feols(lw ~ 1 | i + fid.wom, data = data2)

fe <- fixef(est2)
summary(fe)

