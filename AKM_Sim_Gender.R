library(data.table)
library(reshape)
library(lattice)
library(gridExtra)
library(mvtnorm)
library(ggplot2)
library(futile.logger)

rm(list = ls())

p <- list()
p$nk = 30  # firm types
p$nl = 10  # worker types

p$alpha_sd = 1
p$psi_sd   = 1

# let's draw some FE
p$psi   = with(p,qnorm(1:nk/(nk+1)) * psi_sd)
p$alpha = with(p,qnorm(1:nl/(nl+1)) * alpha_sd)

# let's assume moving PR is fixed
p$lambda = 0.05

p$csort = 0.5 # sorting effect
p$cnetw = 0.2 # network effect
p$csig  = 0.5 # 

# lets create type specific transition matrices
# we are going to use joint normal centered on different values
# G[i,j,k] = Pr[worker i, at firm j, moves to firm k]
getG <- function(p){
  G = with(p,array(0,c(nl,nk,nk)))
  for (l in 1:p$nl) for (k in 1:p$nk) {
    # prob of moving is highest if dnorm(0)
    G[l,k,] = with(p,dnorm( psi - cnetw *psi[k] - csort * alpha[l],sd = csig ))
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


Plot1=wireframe(G[1,,],aspect = c(1,1),xlab = "previous firm",ylab="next firm")
Plot2=wireframe(G[p$nl,,],aspect = c(1,1),xlab = "previous firm",ylab="next firm")
grid.arrange(Plot1, Plot2,nrow=1)


p$nt = 6
p$ni = 1000000

sim <- function(p,G,H){
  set.seed(1)
  
  # we simulate a panel
  network    = array(0,c(p$ni,p$nt))
  spellcount = array(0,c(p$ni,p$nt))
  A = rep(0,p$ni)
  
  for (i in 1:p$ni) {
    # we draw the worker type
    l = sample.int(p$nl,1)
    A[i]=l
    # at time 1, we draw from H
    network[i,1] = sample.int(p$nk,1,prob = H[l,])
    for (t in 2:p$nt) {
      if (runif(1)<p$lambda) {
        network[i,t] = sample.int(p$nk,1,prob = G[l,network[i,t-1],])
        spellcount[i,t] = spellcount[i,t-1] +1
      } else {
        network[i,t]    = network[i,t-1]
        spellcount[i,t] = spellcount[i,t-1]
      }
    }
  }
  
  data  = data.table(melt(network,c('i','t')))
  data2 = data.table(melt(spellcount,c('i','t')))
  setnames(data,"value","k")
  data[,spell := data2$value]
  data[,l := A[i],i]
  data[,alpha := p$alpha[l],l]
  data[,psi := p$psi[k],k]
}

data <- sim(p,G,H)

addSpells <- function(p,dat){
  firm_size = 10
  f_class_count = p$ni/(firm_size*p$nk*p$nt)
  
  dspell <- dat[,list(len=.N),list(i,spell,k)]
  dspell[,fid := sample( 1: pmax(1,sum(len)/f_class_count )   ,.N,replace=TRUE) , k]
  dspell[,fid := .GRP, list(k,fid)]
  
  setkey(dat,i,spell)
  setkey(dspell,i,spell)
  
  dat[, fid:= dspell[dat,fid]]
}

addSpells(p,data)  # adds by reference to the same data.table object (no copy needed)

data[, size := .N, by = .(t, fid)]



# Randomly woman
# tmp <- data.table(i = 1:p$ni, woman = rbinom(p$ni, 1, 0.48))
# data <- merge(data, tmp, by = "i")

# Perfectly correlated with worker type
data[, woman := l < p$nl/2, by = i]
# data[, woman := l + rnorm(1, 0, 2) < p$nk/2, by = i]

data[, h_wom := mean(woman), by = fid]
summary(data$h_wom)

p$w_sigma = 0.8
addWage <- function(p,data){
  # data[, lw := alpha + psi + p$w_sigma * rnorm(.N) ]
  data[, lw := alpha + 1*(1 - woman)*psi + 0.8*woman*psi + p$w_sigma * rnorm(.N)]
}
addWage(p,data)

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

dat1 <- concomp2(data[woman == TRUE])
dat2 <- concomp2(data[woman == FALSE])

est1 <- feols(lw ~ 1 | i + fid, data = dat1)

fe1 <- fixef(est1)


est2 <- feols(lw ~ 1 | i + fid, data = dat2)

fe2 <- fixef(est2)

f1 <- rownames_to_column(data.frame(fe1$fid), "fid") %>% mutate(fid = as.numeric(fid))
f2 <- rownames_to_column(data.frame(fe2$fid), "fid") %>% mutate(fid = as.numeric(fid))

f <- merge(f1, f2, by = "fid")

reg <- feols(data = f, fe1.fid ~ fe2.fid)
etable(reg)

reg <- feols(data = f, fe2.fid ~ fe1.fid)
etable(reg)


data[, split := rnorm(1) > 0, by = i]


dat1 <- concomp2(data[woman == TRUE & split == TRUE])
dat2 <- concomp2(data[woman == FALSE & split == TRUE])

est1 <- feols(lw ~ 1 | i + fid, data = dat1)

fe1 <- fixef(est1)


est2 <- feols(lw ~ 1 | i + fid, data = dat2)

fe2 <- fixef(est2)

f1 <- rownames_to_column(data.frame(fe1$fid), "fid") %>% mutate(fid = as.numeric(fid))
f2 <- rownames_to_column(data.frame(fe2$fid), "fid") %>% mutate(fid = as.numeric(fid))

f <- merge(f1, f2)

dat1 <- concomp2(data[woman == TRUE & split == FALSE])
dat2 <- concomp2(data[woman == FALSE & split == FALSE])

est1 <- feols(lw ~ 1 | i + fid, data = dat1)

fe1 <- fixef(est1)


est2 <- feols(lw ~ 1 | i + fid, data = dat2)

fe2 <- fixef(est2)

f1 <- rownames_to_column(data.frame(fe1$fid), "fid") %>% mutate(fid = as.numeric(fid))
f2 <- rownames_to_column(data.frame(fe2$fid), "fid") %>% mutate(fid = as.numeric(fid))

f2 <- merge(f1, f2)

f <- merge(f, f2, by = "fid")

reg <- feols(data = f, fe1.fid.x ~ 1 | fe2.fid.x ~ fe2.fid.y)
etable(reg)
