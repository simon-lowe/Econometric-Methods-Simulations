library(tidyverse)

N_MC <- 10000
N <- 100

res <- data.frame()

for(i in 1:N_MC){
  dat <- rnorm(N, 0, 1)
  
  tmp <- data.frame(
    m = mean(dat),
    ci.u = mean(dat) + qnorm(0.975)*sd(dat)/sqrt(N),
    ci.d = mean(dat) - qnorm(0.975)*sd(dat)/sqrt(N)
  )
  res <- rbind(res, tmp)
}

# Bias
mean(res$m)

# Wrong CI interpretation
res <- res %>% mutate(
  ci.u.wrong = ci.u[1],
  ci.d.wrong = ci.d[1]
)

with(res, mean(m < ci.u.wrong & m > ci.d.wrong))

# Correct CI interpretation
with(res, mean(0 < ci.u & m > ci.d))
