
library(survstan)

aft <- aftreg(Surv(futime, fustat) ~ 1, data = ovarian,
              baseline = "loglogistic", init = 0)

ph <- phreg(Surv(futime, fustat) ~ 1, data = ovarian,
            baseline = "loglogistic", init = 0)

ah <- ahreg(Surv(futime, fustat) ~ 1, data = ovarian,
            baseline = "loglogistic", init = 0)

po <- poreg(Surv(futime, fustat) ~ 1, data = ovarian,
            baseline = "loglogistic", init = 0)

tidy(aft)
tidy(ph)
tidy(po)
tidy(ah)

summary(aft)
summary(ph)
summary(po)
summary(ah)

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


estimates(aft)
estimates(ph)
estimates(po)
estimates(ah)

vcov(aft)
vcov(ph)
vcov(po)
vcov(ah)

confint(aft)
confint(ph)
confint(po)
confint(ah)
