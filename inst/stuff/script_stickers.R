


library(hexSticker)
library(tidyverse)
library(GGally)
library(showtext)
library(KMsurv)
library(survival)
library(ggfortify)


ekm <- fortify(survfit(Surv(time, status)~1, data = lung))

df <- tibble(
  time = seq(0, 5, by = 0.01),
  survival = exp(-time)
)

p <- ggplot(df, aes(x=time, y=survival)) +
  geom_line()
p <- p+ theme_transparent()
p


## Loading Google fonts (http://www.google.com/fonts)
font_add_google("Gochi Hand", "gochi")
## Automatically use showtext to render text for future devices
showtext_auto()

## use the ggplot2 example
sticker(
  p, package="survstan",
  p_family = "gochi",
  p_size=37.5, s_x=1, s_y=.75, s_width=1.3, s_height=1,
  filename="man/figures/logo.png"
  )

pkgdown::build_favicons(pkg = ".", overwrite = TRUE)

