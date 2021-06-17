library(tidyverse)

n_size <- 100
n_arms <- 2
n_size_arm <- ceiling(n_size/n_arms)

n_MC <- 1000

res_fin <- NULL
res_fin2 <- NULL

for(i in 1:n_MC){
  dat <- data.frame()
  
  t_eff <- 1
  
  for(i in 1:n_arms){
    foo <- rnorm(n_size_arm, ifelse(i == 1, t_eff, 0), 1)
    dat <- rbind(dat, foo)
  }
  dat <- data.frame(t(dat))
  
  names(dat) <- paste0("Y", 1:n_arms)
  
  res <- dat %>%
    summarise(across(everything(), mean))
  
  res_fin <- c(res_fin, max(t(res)))
  res_fin2 <- c(res_fin2, res$Y1)
}


mean(res_fin)
mean(res_fin2)
