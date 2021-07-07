# libraries
library(rprojroot)
library(plyr)
library(ggplot2)
library(magrittr)
library(scales)
library(lme4)
library(sensR)
library(tidyverse)
library(MCMCglmm)
library(nnet)
library(broom)
library(boot)
library(ordinal)
library(ggvis)
library(pander)
library(stargazer)
library(texreg)
library(rprojroot)
library(ggmcmc)
library(Rmisc)
library(texreg)

if(!require(printr)) {
  install.packages(
    'printr',
    type = 'source',
    repos = c('http://yihui.name/xran', 'http://cran.rstudio.com')
  )
} else {
  library(printr)
}

# ggplot settings
theme_set(theme_classic())

# tibble settings
options(tibble.width=Inf)

# knitr settings
# knitr::opts_knit$set(root.dir = '../')
knitr::opts_chunk$set(cache=2, message=FALSE, warning=FALSE, autodep = TRUE, echo=FALSE,
                      cache.path = find_rstudio_root_file('output/crt-cache/'), 
                      fig.path=find_rstudio_root_file("figures/crt/"),
                      dev=c("svg","png"), device="cairo", dpi=96)
options(width=120)


