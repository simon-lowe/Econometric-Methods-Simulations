library(tictoc)
library(future)
library(future.apply)
library(pbapply)
library(tidyverse)
library(parallel)


plan(multisession)

## Set seed (for reproducibility)
set.seed(1234)
# Set sample size
n = 1e6

## Generate a large data frame of fake data for a regression
our_data = 
  tibble(x = rnorm(n), e = rnorm(n)) %>%
  mutate(y = 3 + 2*x + e)

## Function that draws a sample of 10,000 observations, runs a regression and
## extracts the coefficient value on the x variable (should be around 2).
bootstrp = 
  function(i) {
    ## Sample the data
    sample_data = sample_n(our_data, size = 1e4, replace = TRUE)
    ## Run the regression on our sampled data and extract the extract the x
    ## coefficient.
    x_coef = lm(y ~ x, data = sample_data)$coef[2]
    ## Return value
    return(tibble(x_coef = x_coef))
  }

n_MC <- 1e4

tic()
sim_serial = lapply(1:n_MC, bootstrp) %>% bind_rows()
toc(log = TRUE)

tic()
sim_future = future_lapply(1:n_MC, bootstrp, future.seed=123L) %>% bind_rows()
toc()

tic()
sim_pblapply = pblapply(1:n_MC, bootstrp, cl = parallel::detectCores()) %>% bind_rows()
toc()