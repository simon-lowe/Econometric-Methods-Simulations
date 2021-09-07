library(MASS)
library(data.table)
library(sandwich)
library(lmtest)


n <- 6*365
tmp.r <- matrix(0.2, n, n)
tmp.r <- tmp.r^abs(row(tmp.r)-col(tmp.r))

# tmp.r[1:8, 1:8]



# x <- mvrnorm(1, rep(100,n), tmp.r)
x <- mvrnorm(1, seq(100, 200, length.out = n), tmp.r)
# e <- mvrnorm(1, rep(0,n), tmp.r/10)
# acf(x, plot=FALSE, lag.max=5)

dates <- sample(1:n, ceiling(n/100), replace = FALSE)

d <- data.table(
  t = 1:n,
  x = x,
  # x = rep(100, n),
  # e = e,
  y = 0,
  z = 0
)

# d[dates, ':=' (x = x + 1, y = 1)]
d[dates, ':=' (x = x + 1, z = 1)][dates - 1, y := 1]

# d[, ':=' (xm1 = shift(x, 1, type = "lead"), xp1 = shift(x, 1, type = "lag"))]

# cols <- paste0("`rel_t_", -5:5, "`")
# cols <- paste0("t", -5:5 +5)
cols <- paste0("z", -5:5 +5)

# d[, (cols) := shift(.SD, n = -5:5, type = "shift"), .SDcols = "x"]
d[, (cols) := shift(.SD, n = -5:5, type = "shift"), .SDcols = "z"]

# reg <- lm(dat = d, as.formula(paste("y ~ ", paste(cols, collapse = "+"))))
# reg <- lm(dat = d, as.formula(paste("x ~ ", paste(cols, collapse = "+"))))
reg <- lm(dat = d, as.formula(paste("x ~ ", paste(cols, collapse = "+"), " + t")))

summary(reg)

# reg2 <- NeweyWest(lm(dat = d, as.formula(paste("y ~ ", paste(cols, collapse = "+")))))
# reg2 <- NeweyWest(lm(dat = d, as.formula(paste("x ~ ", paste(cols, collapse = "+")))))
reg2 <- NeweyWest(lm(dat = d, as.formula(paste("x ~ ", paste(cols, collapse = "+"), " + t"))))

# sqrt(diag(reg2))

coeftest(reg, vcov = reg2)
