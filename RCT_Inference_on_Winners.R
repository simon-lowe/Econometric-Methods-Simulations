library(tidyverse)

rand_n_list2 <- function(n, n_groups){
  return(sample(rep(seq(1, n_groups), n %/% n_groups + 1, length.out = n), n))
}

n_size <- 100
n_arms <- 10
n_size_arm <- ceiling(n_size/n_arms)

n_MC <- 1000

res_fin <- NULL
res_fin2 <- NULL
res_fin3 <- NULL



# Basic bias --------------------------------------------------------------

t_eff <- 1

for(i in 1:n_MC){
  
  res_arm <- c()
  for(k in 1:n_arms){
    dat <- data.frame(
      Y0 = rnorm(n_size_arm, 0, 1),
      Y1 = rnorm(n_size_arm, ifelse(k == 1, t_eff, 0), 1)
    ) %>%
      mutate(
        # z = rbinom(n(), 1, 0.5),
        z = rand_n_list2(n(), 2) - 1,
        y = z*Y1 + (1-z)*Y0
      ) %>%
      summarise(mean(y[z == 1] - y[z == 0]))
    res_arm <- c(res_arm, dat[[1]])
  }
  if(i == 1){
    r <- which.max(res_arm)
  }
  res_fin <- c(res_fin, max(res_arm))
  res_fin2 <- c(res_fin2, res_arm[r])
  if(which.max(res_arm) == 1){
    res_fin3 <- c(res_fin3, res_arm[1])
  }
}


mean(res_fin, na.rm = T)
mean(res_fin2, na.rm = T)
mean(res_fin3, na.rm = T)


# Inference on the chosen one ---------------------------------------------

t_eff <- 1

cov_res <- c()

for(i in 1:100){
  
  res_arm <- c()
  sd_arm <- c()
  for(i in 1:n_arms){
    dat <- data.frame(
      Y0 = rnorm(n_size_arm, 0, 1),
      Y1 = rnorm(n_size_arm, ifelse(i == 1, t_eff, 0), 1)
    ) %>%
      mutate(
        # z = rbinom(n(), 1, 0.5),
        z = rand_n_list2(n(), 2) - 1,
        y = z*Y1 + (1-z)*Y0
      ) 
    reg <- lm(data = dat, y ~ z)
    res_arm <- c(res_arm, coef(summary(reg))[2,1])
    sd_arm <- c(sd_arm, coef(summary(reg))[2,2])
  }
  
  arm_chosen <- which.max(res_arm)
  CI_u <- res_arm[arm_chosen] + 1.96*sd_arm[arm_chosen]
  CI_d <- res_arm[arm_chosen] - 1.96*sd_arm[arm_chosen]
  
  reps <- c()
  for(j in 1:5000){
    dat <- data.frame(
      Y0 = rnorm(n_size_arm, 0, 1),
      Y1 = rnorm(n_size_arm, ifelse(arm_chosen == 1, t_eff, 0), 1)
    ) %>%
      mutate(
        z = rand_n_list2(n(), 2) - 1,
        y = z*Y1 + (1-z)*Y0
      ) %>%
    summarise(mean(y[z == 1] - y[z == 0]))
    reps <- c(reps, dat[[1]])
  }
  
  cov_res <- c(cov_res, mean(reps > CI_d & reps < CI_u))
}
