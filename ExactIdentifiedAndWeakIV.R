rm(list = ls())

library(tidyverse)
library(MASS)
library(ivreg)
library(tictoc)
library(future)
library(future.apply)
library(pbapply)
library(tidyverse)
library(parallel)


plan(multisession)

## Set seed (for reproducibility)
set.seed(1234)

rhoxz <- 0.15
rhoxe <- 0.4
rhoze <- 0
n <- 100

IVsim <- function(aux, rhoxz = 0.15, rhoxe = 0.4, rhoze = 0, n = 100){
  cor_mat <- matrix(c(1, rhoxz, rhoxe, rhoxz, 1, rhoze, rhoxe, rhoze, 1), 3, 3)
  xze <- mvrnorm(n, rep(0, 3), cor_mat)
  dat <- data.frame(
    x = xze[,1],
    z = xze[,2],
    e = xze[,3]
  ) %>% 
    mutate(
      y = 1 + x + e
  )
  
  return(data.frame(
    tsls = coef(ivreg(data = dat, y ~ x | z))[[2]],
    ols = coef(lm(data = dat, y ~ x))[[2]]
    # fs = coef(lm(data = dat, x ~ z))[[2]]
  ))
}

n_MC <- 2.5e5

# tic()
# sim_serial = lapply(1:n_MC, IVsim) %>% bind_rows()
# toc(log = TRUE)

tic()
sim_future = future_lapply(1:n_MC, IVsim, future.seed=123L) %>% bind_rows()
toc()

# tic()
# sim_pblapply = pblapply(1:n_MC, IVsim, cl = parallel::detectCores()) %>% bind_rows()
# toc()

# Histograms
hist(sim_future$ols, breaks = 20)
hist(sim_future$tsls, breaks = 20)
# with(hist(sim_future$tsls, breaks = 20)

hgA <- hist(sim_future$ols, breaks = 20, plot = FALSE)
hgB <- hist(sim_future$tsls, breaks = 20, plot = FALSE)

c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")

plot(hgB, col = c1) # Plot 1st histogram using a transparent color
plot(hgA, col = c2, add = TRUE)

res <- sim_future %>%
  # filter(abs(tsls) < 2) %>%
  filter(tsls < 2 & tsls > 0) %>%
  pivot_longer(
    # cols = c("tsls", "ols", "fs"),
    cols = c("tsls", "ols"),
    names_to = "type",
    values_to = "value"
  )

library(hrbrthemes)
p <- res %>%
  ggplot( aes(x=value, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  theme_ipsum() +
  labs(fill="")
p
