library(readr)
library(dplyr)
library(MASS)
source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/medsim.R")

load("Brader_et_al2008.RData")

# function for demeaning
demean <- function(x) x - mean(x, na.rm = TRUE)

Brader2 <- Brader %>%
  dplyr::select(immigr, emo, p_harm, tone_eth, ppage, ppeducat, ppgender, ppincimp) %>% na.omit() %>%
  mutate(immigr = as.numeric(5 - immigr),
         hs = (ppeducat == "high school"),
         sc = (ppeducat == "some college"),
         ba = (ppeducat == "bachelor's degree or higher"),
         female = (ppgender == "female")) %>%
  mutate_at(vars(ppage, female, hs, sc, ba, ppincimp), demean)

mydata <- Brader2 %>%
  rename(y = immigr, a = tone_eth, m = emo, l = p_harm, c1 = ppage,c2 = female, c3 = hs, c4 = sc, c5 = ba, c6 = ppincimp)

sd(mydata$y)

mydata$y <- as.factor(mydata$y)
mydata$a <- as.factor(mydata$a)
mydata$m <- as.factor(mydata$m)
mydata$l <- as.factor(mydata$l)

# Define the model specifications
model_spec_1_1 <- list(
  list(func = "polr", formula = l ~ a + c1 + c2 + c3 + c4 + c5 + c6, args = list()),
  list(func = "polr", formula = m ~ a + I(as.numeric(as.character(l))) + c1 + c2 + c3 + c4 + c5 + c6, args = list()),
  list(func = "polr", formula = y ~ a + I(as.numeric(as.character(l))) + I(as.numeric(as.character(m))) + c1 + c2 + c3 + c4 + c5 + c6, args = list())
)

model_spec_1_2 <- list(
  list(func = "polr", formula = l ~ a + c1 + c2 + c3 + c4 + c5 + c6, args = list()),
  list(func = "polr", formula = m ~ a * I(as.numeric(as.character(l))) + c1 + c2 + c3 + c4 + c5 + c6, args = list()),
  list(func = "polr", formula = y ~ a * I(as.numeric(as.character(l))) + a * I(as.numeric(as.character(m))) + c1 + c2 + c3 + c4 + c5 + c6, args = list())
)

# Path-specific Effects
# Additive
medsim(data=mydata, num_sim = 2000,
       cat_list = c("0", "1"), treatment = "a", model_spec = model_spec_1_1,
       boot = TRUE, boot_reps = 2000, seed = 257760)

# Exposure-Mediators Interactions
medsim(data=mydata, num_sim = 2000,
       cat_list = c("0", "1"), treatment = "a", model_spec = model_spec_1_2,
       boot = TRUE, boot_reps = 2000, seed = 257760)

# Interventional Effects
# Additive
medsim(data=mydata, num_sim = 2000,
       cat_list = c("0", "1"), treatment = "a",
       intv_med = ("m"), model_spec = model_spec_1_1,
       boot = TRUE, boot_reps = 2000, seed = 257760)

# Exposure-Mediators Interactions
medsim(data=mydata, num_sim = 2000,
       cat_list = c("0", "1"), treatment = "a",
       intv_med = ("m"), model_spec = model_spec_1_2,
       boot = TRUE, boot_reps = 2000, seed = 257760)
