
library(survstan)
library(dplyr)

ovarian <- ovarian %>%
  mutate(
    rx = as.factor(rx)
  )

fit1 <- aftreg(Surv(futime, fustat) ~ 1, data = ovarian, baseline = "weibull", init = 0)
fit2 <- aftreg(Surv(futime, fustat) ~ age , data = ovarian, baseline = "weibull", init = 0)
fit3 <- aftreg(Surv(futime, fustat) ~ age + rx, data = ovarian, baseline = "weibull", init = 0)
anova.survstan(fit1, fit2, fit3)



fit1 <- aftreg(Surv(futime, fustat) ~ 1, data = ovarian, baseline = "weibull", init = 0)
fit2 <- aftreg(Surv(futime, fustat) ~ age , data = ovarian, baseline = "weibull", init = 0)
fit3 <- phreg(Surv(futime, fustat) ~ age + rx, data = ovarian, baseline = "weibull", init = 0)
anova.survstan(fit1, fit2, fit3)
