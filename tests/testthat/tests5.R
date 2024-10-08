
library(survstan)
library(dplyr)

data(ipass)

ipass <- ipass %>%
  mutate(
    arm = as.factor(arm)
  )


fit <- ypreg(
  Surv(time, status) ~ arm,
  dist = bernstein(m=3),
  data = ipass,
)

newdata1 <- data.frame(
  arm = factor(0, levels = c("0", "1"))
)

newdata2 <- data.frame(
  arm = factor(1, levels = c("0", "1"))
)


newdata <- data.frame(
  arm = factor(0:1, levels = c("0", "1"))
)

surv <- survfit(fit, newdata)

tcross <- cross_time(fit, newdata1, newdata2, nboot = 10)
tcross




fit <- ehreg(
  Surv(time, status) ~ arm,
  dist = bernstein(m=3),
  data = ipass,
)


surv <- survfit(fit, newdata)

tcross <- cross_time(fit, newdata1, newdata2, nboot = 10)
tcross
