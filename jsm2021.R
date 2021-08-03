# Packages
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(lme4)
library(HLMdiag)
library(qqplotr)
library(cowplot)
library(lmeresampler)

# Example data set, dropping NAs
osteo <- read_delim("ost phenotypes.dat", delim = "\t") %>%
  select(-contains("X")) %>%
  rename(BMD = `CalL2-L4BMD`) %>%
  mutate(Sex = as_factor(Gender) %>% fct_recode(male = "1", female = "2"), 
         CentreNumber = factor(CentreNumber)) %>%
  filter(Proband == "No") %>%
  drop_na(BMD, Sex, Age, Height, Weight, CentreNumber) %>%
  mutate(Age = Age - mean(Age),
         Height = Height - mean(Height),
         Weight = Weight - mean(Weight))



## ----galwey model-----------------------------------------------------------------------------------------
fm <- lmer(BMD ~ Sex + Age + I(Age^2) + Height*Weight + 
             I(Height^2) + I(Weight^2) + (1|CentreNumber), 
           data = osteo)



## ----include=FALSE, cache=TRUE----------------------------------------------------------------------------
stime <- system.time({resid_boot <- bootstrap(
  fm,                 # lme4/nlme output 
  .f = fixef,         # user-specified function
  type = "residual",  # bootstrap algorithm
  B = 15000           # No. resamples
)})


## ----include=FALSE, cache=TRUE----------------------------------------------------------------------------
library(foreach)
library(doParallel)

no_cores <- 5
cl <- makeCluster(no_cores, type = "FORK")
registerDoParallel(cores = no_cores)

ptime <- system.time({boot_parallel <- foreach(
  B = rep(3000, 5), 
  .combine = combine_lmeresamp,
  .packages = c("lmeresampler", "lme4")) %dopar% {
    bootstrap(fm, .f = fixef, type = "residual", B = B)
  }
})


stopCluster(cl)

save(stime, ptime, resid_boot, file = "jsm_boot_results.RData")