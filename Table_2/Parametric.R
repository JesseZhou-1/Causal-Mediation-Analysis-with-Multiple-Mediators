library(readr)
library(dplyr)
source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/medsim.R")

mydata <- read_csv("cleaned_natl2003_bin.csv")

mydata$y <- as.factor(mydata$y)
mydata$a <- as.factor(mydata$a)
mydata$l <- as.factor(mydata$l)
mydata$m <- as.factor(mydata$m)
mydata$c1 <- as.factor(mydata$c1)
mydata$c2 <- as.factor(mydata$c2)
mydata$c3 <- as.factor(mydata$c3)
mydata$c4 <- as.factor(mydata$c4)

# Define the model specifications 
model_spec <- list(
  list(func = "polr", formula = l ~ a * (c1 + c2 + c3 + c4), args = list()),
  list(func = "glm", formula = m ~ a * l + a * (c1 + c2 + c3 + c4) + l * (c1 + c2 + c3 + c4), args = list(family = "binomial")),
  list(func = "glm", formula = y ~ a * l + a * m + l * m + (a + l + m) * (c1 + c2 + c3 + c4), args = list(family = "binomial"))
)

# Path-specific Effects
medsim(data=mydata, num_sim = 2000,
        cat_list = c("0", "1"), treatment = "a",
        intv_med = NULL, model_spec = model_spec)

# Interventional Effects
medsim(data=mydata, num_sim = 2000,
        cat_list = c("0", "1"), treatment = "a",
        intv_med = ("m"), model_spec = model_spec)

