
library(survstan)

ovarian$rx <- as.factor(ovarian$rx)

newdata <- expand.grid(
  age = c(50, 60),
  rx = as.factor(1:2)
)

baselines <- c("exponential", "weibull", "lognormal", "loglogistic")

for(baseline in baselines){

  aft <- aftreg(Surv(futime, fustat) ~ age + rx, data = ovarian,
                baseline = baseline, init = 0)

  ph <- phreg(Surv(futime, fustat) ~ age + rx, data = ovarian,
              baseline = baseline, init = 0)

  ah <- ahreg(Surv(futime, fustat) ~ age + rx, data = ovarian,
              baseline = baseline, init = 0)

  po <- poreg(Surv(futime, fustat) ~ age + rx, data = ovarian,
              baseline = baseline, init = 0)

  yp <- ypreg(Surv(futime, fustat) ~ age + rx, data = ovarian,
              baseline = baseline, init = 0)

  tidy(aft)
  tidy(ph)
  tidy(po)
  tidy(ah)
  tidy(yp)

  summary(aft)
  summary(ph)
  summary(po)
  summary(ah)
  summary(yp)
  coef(aft)
  coef(ph)
  coef(po)
  coef(ah)
  coef(yp)

  confint(aft)
  confint(ph)
  confint(po)
  confint(ah)
  confint(yp)
  estimates(aft)
  estimates(ph)
  estimates(po)
  estimates(ah)
  estimates(yp)
  vcov(aft)
  vcov(ph)
  vcov(po)
  vcov(ah)
  vcov(yp)
  AIC(aft, ph, po, ah)


  ggresiduals(aft, "coxsnell")
  ggresiduals(aft, "martingale")
  ggresiduals(aft, "deviance")

  ggresiduals(ph, "coxsnell")
  ggresiduals(ph, "martingale")
  ggresiduals(ph, "deviance")

  ggresiduals(po, "coxsnell")
  ggresiduals(po, "martingale")
  ggresiduals(po, "deviance")

  ggresiduals(ah, "coxsnell")
  ggresiduals(ah, "martingale")
  ggresiduals(ah, "deviance")

  ggresiduals(yp, "coxsnell")
  ggresiduals(yp, "martingale")
  ggresiduals(yp, "deviance")

  model.matrix(aft)
  model.matrix(ph)
  model.matrix(po)
  model.matrix(ah)
  model.matrix(yp)

  surv_aft <- survfit(aft, newdata)
  surv_ah <- survfit(ah, newdata)
  surv_ph <- survfit(ph, newdata)
  surv_po <- survfit(po, newdata)
  surv_yp <- survfit(yp, newdata)
}
