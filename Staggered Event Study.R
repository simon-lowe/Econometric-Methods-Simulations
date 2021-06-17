# Andrew Baker's code for event studies
# Taken from https://andrewcbaker.netlify.app/2020/06/27/how-to-create-relative-time-indicators/
# On 13/04/2021

rm(list = ls())

# Load packages -----------------------------------------------------------

library(tidyverse)
library(lfe)
library(ggthemes)
library(fastDummies)
library(broom)

theme_set(theme_clean() + theme(plot.background = element_blank()))
select <- dplyr::select

set.seed(20200403)

# Setting parameters ------------------------------------------------------

n_units <- 1000
n_state <- 40

# Generate data - 250 firms are treated every period, with the tre --------

make_data <- function(n_units, n_state) {
  
  # Fixed Effects ------------------------------------------------
  # unit fixed effects
  unit <- tibble(
    unit = 1:n_units, 
    unit_fe = rnorm(n_units, 0, 1),
    # generate state
    state = sample(1:n_state, n_units, replace = TRUE),
    # generate treatment effect
    mu = rnorm(n_units, 0.3, 0.2))
  
  # year fixed effects 
  year <- tibble(
    year = 1980:2010) %>%
    mutate(
      year_fe = rnorm(n(), 0, 1)
    )
  
  # Trend Break -------------------------------------------------------------
  # Put the states into treatment groups
  treat_taus <- tibble(
    # sample the states randomly
    state = sample(1:n_state, n_state, replace = FALSE),
    # place the randomly sampled states into five treatment groups G_g
    # cohort_year = sort(rep(c(1986, 1992, 1998, 2004), 10))
    cohort_year = sample(c(1986, 1992, 1998, 2004), n_state, replace = TRUE)
  )
  
  # make main dataset
  # full interaction of unit X year 
  expand_grid(unit = 1:n_units, year = 1980:2010) %>% 
    left_join(., unit) %>% 
    left_join(., year) %>% 
    left_join(., treat_taus) %>% 
    # make error term and get treatment indicators and treatment effects
    mutate(error = rnorm(n(), 0, 0.5),
           treat = ifelse(year >= cohort_year, 1, 0),
           tau = ifelse(treat == 1, mu, 0)) %>% 
    # calculate cumulative treatment effects
    group_by(unit) %>% 
    mutate(tau_cum = cumsum(tau)) %>% 
    ungroup() %>% 
    # calculate the dep variable
    mutate(dep_var = unit_fe + year_fe + tau_cum + error)
  
}

# make data
data <- make_data(n_units, n_state)

# plot
plot <- data %>% 
  ggplot(aes(x = year, y = dep_var, group = unit)) + 
  geom_line(alpha = 1/8, color = "grey") + 
  geom_line(data = data %>% 
              group_by(cohort_year, year) %>% 
              summarize(dep_var = mean(dep_var)),
            aes(x = year, y = dep_var, group = factor(cohort_year),
                color = factor(cohort_year)),
            size = 2) + 
  labs(x = "", y = "Value") + 
  geom_vline(xintercept = 1986, color = '#E41A1C', size = 2) + 
  geom_vline(xintercept = 1992, color = '#377EB8', size = 2) + 
  geom_vline(xintercept = 1998, color = '#4DAF4A', size = 2) + 
  geom_vline(xintercept = 2004, color = '#984EA3', size = 2) + 
  scale_color_brewer(palette = 'Set1') + 
  theme(legend.position = 'bottom',
        legend.title = element_blank(), 
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
plot



# Event-study with Post and Pre indicators --------------------------------

# variables we will use
keepvars <- c("`rel_year_-5`",  "`rel_year_-4`",  "`rel_year_-3`",  "`rel_year_-2`",
              "rel_year_0", "rel_year_1", "rel_year_2", "rel_year_3", "rel_year_4", "rel_year_5")

# function to run ES DID
run_ES_DiD <- function(aux, n_units = 1000, n_state = 40) {
  
  # resimulate the data
  data <- make_data(n_units, n_state)
  
  # make dummy columns
  data <- data %>% 
    # make dummies
    mutate(rel_year = year - cohort_year) %>% 
    dummy_cols(select_columns = "rel_year") %>% 
    # generate pre and post dummies
    mutate(Pre = ifelse(rel_year < -5, 1, 0),
           Post = ifelse(rel_year > 5, 1, 0))
  
  # estimate the model
  mod <- felm(dep_var ~ Pre + `rel_year_-5` + `rel_year_-4` + `rel_year_-3` + `rel_year_-2` + 
                `rel_year_0` + `rel_year_1` + `rel_year_2` + `rel_year_3` + `rel_year_4` + 
                `rel_year_5` + Post | unit + year | 0 | state, data = data, exactDOF = TRUE)
  
  # grab the obs we need
  # broom::tidy(mod) %>% 
  #   filter(term %in% keepvars) %>% 
  #   mutate(t = c(-5:-2, 0:5)) %>% 
  #   select(t, estimate)
  tibble::enframe(coef(mod)) %>%
    filter(name %in% keepvars) %>%
    mutate(t = c(-5:-2, 0:5)) %>%
    select(t, estimate = value) %>%
    as.data.frame()
}

# run it 1000 times
set.seed(20200403)
data <- map_dfr(1:10, run_ES_DiD, n_units = n_units, n_state = n_state)

data %>% 
  group_by(t) %>% 
  summarize(avg = mean(estimate),
            sd = sd(estimate),
            lower.ci = avg - 1.96*sd,
            upper.ci = avg + 1.96*sd) %>% 
  bind_rows(tibble(t = -1, avg = 0, sd = 0, lower.ci = 0, upper.ci = 0)) %>% 
  mutate(true_tau = ifelse(t >= 0, (t + 1)*.3, 0)) %>% 
  ggplot(aes(x = t, y = avg)) + 
  geom_linerange(aes(ymin = lower.ci, ymax = upper.ci), color = 'darkgrey', size = 2) + 
  geom_point(color = 'blue', size = 4) + 
  geom_line(aes(x = t, y = true_tau), color = 'red', linetype = "dashed", size = 2) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  scale_x_continuous(breaks = -5:5) + 
  labs(x = "Relative Time", y = "Estimate") + 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))


# Event study using only Post ---------------------------------------------

# function to run ES DID
run_ES_DiD <- function(aux, n_units = 1000, n_state = 40) {
  
  # resimulate the data
  data <- make_data(n_units, n_state)
  
  # make dummy columns
  data <- data %>% 
    # make dummies
    mutate(rel_year = year - cohort_year) %>% 
    dummy_cols(select_columns = "rel_year") %>% 
    # generate pre and post dummies
    mutate(Pre = ifelse(rel_year < -5, 1, 0),
           Post = ifelse(rel_year > 5, 1, 0))
  
  # estimate the model
  mod <- felm(dep_var ~ `rel_year_-5` + `rel_year_-4` + `rel_year_-3` + `rel_year_-2` + 
                `rel_year_0` + `rel_year_1` + `rel_year_2` + `rel_year_3` + `rel_year_4` + 
                `rel_year_5` + Post | unit + year | 0 | state, data = data, exactDOF = TRUE)
  
  # grab the obs we need
  # broom::tidy(mod) %>% 
  #   filter(term %in% keepvars) %>% 
  #   mutate(t = c(-5:-2, 0:5)) %>% 
  #   select(t, estimate)
  tibble::enframe(coef(mod)) %>%
    filter(name %in% keepvars) %>%
    mutate(t = c(-5:-2, 0:5)) %>%
    select(t, estimate = value) %>%
    as.data.frame()
}

# run it 1000 times
set.seed(20200403)
data <- map_dfr(1:10, run_ES_DiD, n_units = n_units, n_state = n_state)

data %>% 
  group_by(t) %>% 
  summarize(avg = mean(estimate),
            sd = sd(estimate),
            lower.ci = avg - 1.96*sd,
            upper.ci = avg + 1.96*sd) %>% 
  bind_rows(tibble(t = -1, avg = 0, sd = 0, lower.ci = 0, upper.ci = 0)) %>% 
  mutate(true_tau = ifelse(t >= 0, (t + 1)*.3, 0)) %>% 
  ggplot(aes(x = t, y = avg)) + 
  geom_linerange(aes(ymin = lower.ci, ymax = upper.ci), color = 'darkgrey', size = 2) + 
  geom_point(color = 'blue', size = 4) + 
  geom_line(aes(x = t, y = true_tau), color = 'red', linetype = "dashed", size = 2) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  scale_x_continuous(breaks = -5:5) + 
  labs(x = "Relative Time", y = "Estimate") + 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))


# Event study with fully saturated model ----------------------------------

# function to run ES DID
run_ES_DiD <- function(aux, n_units = 1000, n_state = 40) {
  
  # resimulate the data
  data <- make_data(n_units, n_state)
  
  # make dummy columns
  data <- data %>% 
    # make relative year indicator
    mutate(rel_year = year - cohort_year)
  
  # get the minimum relative year - we need this to reindex
  min_year <- min(data$rel_year)
  
  # reindex the relative years
  data <- data %>% 
    mutate(rel_year = rel_year - min_year) %>% 
    dummy_cols(select_columns = "rel_year")
  
  # make regression formula 
  indics <- paste("rel_year", (1:max(data$rel_year))[-(-1 - min_year)], sep = "_", collapse = " + ")
  keepvars <- paste("rel_year", c(-5:-2, 0:5) - min_year, sep = "_")  
  formula <- as.formula(paste("dep_var ~", indics, "| unit + year | 0 | state"))
  
  # run mod
  mod <- felm(formula, data = data, exactDOF = TRUE)
  
  # grab the obs we need
  tibble::enframe(coef(mod)) %>%
    filter(name %in% keepvars) %>%
    mutate(t = c(-5:-2, 0:5)) %>%
    select(t, estimate = value) %>%
    as.data.frame()
}

# run it 1000 times
set.seed(20200403)
data <- map_dfr(1:10, run_ES_DiD)

data %>% 
  group_by(t) %>% 
  summarize(avg = mean(estimate),
            sd = sd(estimate),
            lower.ci = avg - 1.96*sd,
            upper.ci = avg + 1.96*sd) %>% 
  bind_rows(tibble(t = -1, avg = 0, sd = 0, lower.ci = 0, upper.ci = 0)) %>% 
  mutate(true_tau = ifelse(t >= 0, (t + 1)*.3, 0)) %>% 
  ggplot(aes(x = t, y = avg)) + 
  geom_linerange(aes(ymin = lower.ci, ymax = upper.ci), color = 'darkgrey', size = 2) + 
  geom_point(color = 'blue', size = 4) + 
  geom_line(aes(x = t, y = true_tau), color = 'red', linetype = "dashed", size = 2) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  scale_x_continuous(breaks = -5:5) + 
  labs(x = "Relative Time", y = "Estimate") + 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))


# Event study with Pre and Post dropping the last treated cohort ----------

keepvars <- c("`rel_year_-5`",  "`rel_year_-4`",  "`rel_year_-3`",  "`rel_year_-2`",
              "rel_year_0", "rel_year_1", "rel_year_2", "rel_year_3", "rel_year_4", "rel_year_5")

# function to run ES DID
run_ES_DiD <- function(aux, n_units = 1000, n_state = 40) {
  
  # resimulate the data
  data <- make_data(n_units, n_state) %>%
    filter(year <= 2003) %>% 
    mutate(cohort_year = ifelse(cohort_year == 2004, 0, cohort_year))
  
  # make dummy columns
  data <- data %>% 
    # make dummies
    mutate(rel_year = year - cohort_year) %>% 
    dummy_cols(select_columns = "rel_year") %>% 
    # generate pre and post dummies
    mutate(Pre = ifelse(rel_year < -5, 1, 0),
           Post = ifelse(rel_year > 5 & cohort_year != 0, 1, 0))
  
  # estimate the model
  mod <- felm(dep_var ~ Pre + `rel_year_-5` + `rel_year_-4` + `rel_year_-3` + `rel_year_-2` + 
                `rel_year_0` + `rel_year_1` + `rel_year_2` + `rel_year_3` + `rel_year_4` + 
                `rel_year_5` + Post | unit + year | 0 | state, data = data, exactDOF = TRUE)
  
  # grab the obs we need
  tibble::enframe(coef(mod)) %>%
    cbind(se = unname(mod$se)) %>%
    mutate(
      conf.low = value - 1.96*se,
      conf.high = value + 1.96*se,
    ) %>%
    filter(name %in% keepvars) %>%
    mutate(t = c(-5:-2, 0:5)) %>%
    select(t, estimate = value, conf.low, conf.high) %>%
    as.data.frame()
}

# run it once 
set.seed(20200403)
data <- map_dfr(1, run_ES_DiD, n_units = n_units, n_state = n_state)

data %>% 
  bind_rows(tibble(t = -1, estimate = 0, conf.low = 0, conf.high = 0)) %>% 
  mutate(true_tau = ifelse(t >= 0, (t + 1)*.3, 0)) %>% 
  ggplot(aes(x = t, y = estimate)) + 
  geom_linerange(aes(ymin = conf.low, ymax = conf.high), color = 'darkgrey', size = 2) + 
  geom_point(color = 'blue', size = 4) + 
  geom_line(aes(x = t, y = true_tau), color = 'red', linetype = "dashed", size = 2) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  scale_x_continuous(breaks = -5:5) + 
  labs(x = "Relative Time", y = "Estimate") + 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))


# Event study with  Post dropping the last treated cohort -----------------

# function to run ES DID
run_ES_DiD <- function(aux, n_units = 1000, n_state = 40) {
  
  # resimulate the data
  data <- make_data(n_units, n_state) %>%
    filter(year <= 2003) %>% 
    mutate(cohort_year = ifelse(cohort_year == 2004, 0, cohort_year))
  
  # make dummy columns
  data <- data %>% 
    # make dummies
    mutate(rel_year = year - cohort_year) %>% 
    dummy_cols(select_columns = "rel_year") %>% 
    # generate pre and post dummies
    mutate(Pre = ifelse(rel_year < -5, 1, 0),
           Post = ifelse(rel_year > 5 & cohort_year != 0, 1, 0))
  
  # estimate the model
  mod <- felm(dep_var ~ `rel_year_-5` + `rel_year_-4` + `rel_year_-3` + `rel_year_-2` + 
                `rel_year_0` + `rel_year_1` + `rel_year_2` + `rel_year_3` + `rel_year_4` + 
                `rel_year_5` + Post | unit + year | 0 | state, data = data, exactDOF = TRUE)
  
  # grab the obs we need
  tibble::enframe(coef(mod)) %>%
    cbind(se = unname(mod$se)) %>%
    mutate(
      conf.low = value - 1.96*se,
      conf.high = value + 1.96*se,
    ) %>%
    filter(name %in% keepvars) %>%
    mutate(t = c(-5:-2, 0:5)) %>%
    select(t, estimate = value, conf.low, conf.high) %>%
    as.data.frame()
}

# run it once 
set.seed(20200403)
data <- map_dfr(1, run_ES_DiD, n_units = n_units, n_state = n_state)

data %>% 
  bind_rows(tibble(t = -1, estimate = 0, conf.low = 0, conf.high = 0)) %>% 
  mutate(true_tau = ifelse(t >= 0, (t + 1)*.3, 0)) %>% 
  ggplot(aes(x = t, y = estimate)) + 
  geom_linerange(aes(ymin = conf.low, ymax = conf.high), color = 'darkgrey', size = 2) + 
  geom_point(color = 'blue', size = 4) + 
  geom_line(aes(x = t, y = true_tau), color = 'red', linetype = "dashed", size = 2) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  scale_x_continuous(breaks = -5:5) + 
  labs(x = "Relative Time", y = "Estimate") + 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))


# Event study with saturated model dropping the last treated cohort -------

# function to run ES DID
run_ES_DiD <- function(aux, n_units = 1000, n_state = 40) {
  
  # resimulate the data
  data <- make_data(n_units, n_state) %>%
    filter(year <= 2003) %>% 
    mutate(cohort_year = ifelse(cohort_year == 2004, 0, cohort_year))
  
  # make dummy columns
  data <- data %>% 
    # make dummies
    mutate(rel_year = year - cohort_year)
  
  # get the minimum relative year - we need this to reindex
  min_year <- min(data %>% filter(cohort_year != 0) %>% pull(rel_year))
  
  # reindex the relative years
  data <- data %>% 
    mutate(rel_year = rel_year - min_year) %>% 
    dummy_cols(select_columns = "rel_year")
  
  # make regression formula 
  indics <- paste("rel_year", (1:max(data %>% filter(cohort_year != 0) %>% pull(rel_year)))[-(-1 - min_year)], sep = "_", collapse = " + ")
  keepvars <- paste("rel_year", c(-5:-2, 0:5) - min_year, sep = "_")  
  formula <- as.formula(paste("dep_var ~", indics, "| unit + year | 0 | state"))
  
  # run mod
  mod <- felm(formula, data = data, exactDOF = TRUE)
  
  # grab the obs we need
  tibble::enframe(coef(mod)) %>%
    cbind(se = unname(mod$se)) %>%
    mutate(
      conf.low = value - 1.96*se,
      conf.high = value + 1.96*se,
    ) %>%
    filter(name %in% keepvars) %>%
    mutate(t = c(-5:-2, 0:5)) %>%
    select(t, estimate = value, conf.low, conf.high) %>%
    as.data.frame()
}

# run it once 
set.seed(20200403)
data <- map_dfr(1, run_ES_DiD, n_units = n_units, n_state = n_state)

data %>% 
  bind_rows(tibble(t = -1, estimate = 0, conf.low = 0, conf.high = 0)) %>% 
  mutate(true_tau = ifelse(t >= 0, (t + 1)*.3, 0)) %>% 
  ggplot(aes(x = t, y = estimate)) + 
  geom_linerange(aes(ymin = conf.low, ymax = conf.high), color = 'darkgrey', size = 2) + 
  geom_point(color = 'blue', size = 4) + 
  geom_line(aes(x = t, y = true_tau), color = 'red', linetype = "dashed", size = 2) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  scale_x_continuous(breaks = -5:5) + 
  labs(x = "Relative Time", y = "Estimate") + 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))



# Using Callaway and Sant'Anna --------------------------------------------

library(did)

reg <- att_gt(yname = "dep_var",
              tname = "year",
              idname = "unit",
              gname = "cohort_year",
              xformla = NULL,
              data = data %>% filter(year <= 2003)
)

summary(reg)
ggdid(reg)

agg.es <- aggte(reg, type = "dynamic")
summary(agg.es)
ggdid(agg.es)

reg2 <- att_gt(yname = "dep_var",
               tname = "year",
               idname = "unit",
               gname = "cohort_year",
               xformla = NULL,
               data = data,
               control_group = "notyettreated"
)

agg.es2 <- aggte(reg2, type = "dynamic")
summary(agg.es2)
ggdid(agg.es2)
