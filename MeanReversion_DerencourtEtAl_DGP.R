rm(list = ls())
gc()

library(data.table)
library(tidyverse)
library(tictoc)
library(fixest)


# Creating the data -------------------------------------------------------

ng <- 2.7e5

dat <- data.table(g = rep(1:ng, each = 25), t = rep(-12:12, ng))

dat[, lw := rnorm(.N, 2.5, 0.32)]

create_group_obs <- function(group, g_av){
  ts <- c(sample(-12:-1, 2), sample(0:12, 3))
  sizes <- rpois(5, 10)
  obs <- c()
  for(i in 1:5){
    obs <- c(obs, rnorm(sizes[i], g_av[g == group & t == ts[i]]$lw, 0.1))
  }
  foo <- data.table(g = rep(group, sum(sizes)), t = rep(ts, sizes), lwi = obs)
  return(foo)
}

tic()
dat2 <- map_df(1:ng, create_group_obs, g_av = dat)
toc()

dat3 <- merge(dat, dat2, by = c("g", "t"), all = TRUE)


# Compute exposure --------------------------------------------------------

dat3[, exposure := sum((t < 0)*(lwi < log(15)), na.rm = T)/sum((t < 0)*!is.na(lwi)), by = g]

dat3[, treat := exposure > 0.5]

reg <- feols(data = dat3, lwi ~ i(t, exposure, -1) | g)
iplot(reg)

reg <- feols(data = dat3, lwi ~ i(t, treat, -1) | g)
iplot(reg)
