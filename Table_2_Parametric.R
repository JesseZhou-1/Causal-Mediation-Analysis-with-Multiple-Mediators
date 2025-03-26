library(readr)
library(dplyr)


cleaned_natl2003 <- read_csv("cleaned_natl2003.csv")
cleaned_natl2003 <- na.omit(cleaned_natl2003)

cleaned_natl2003_bin <- subset(cleaned_natl2003, !(precare == 1 | precare == 3))
cleaned_natl2003_bin$precare <- ifelse(cleaned_natl2003_bin$precare == 2, 1, cleaned_natl2003_bin$precare)

# Rename columns and create a new data frame
mydata <- cleaned_natl2003_bin %>%
  rename(y = prebirth, a = precare, m = urf_eclam, l = cigs, c1 = age,c2 = somecollege, c3 = mracerec, c4 = mar)

mydata$y <- as.factor(mydata$y)
mydata$a <- as.factor(mydata$a)
mydata$m <- as.factor(mydata$m)
mydata$c1 <- as.factor(mydata$c1)
mydata$c2 <- as.factor(mydata$c2)
mydata$c3 <- as.factor(mydata$c3)
mydata$c4 <- as.factor(mydata$c4)

# Define the model specifications 
model_spec <- list(
  list(func = "glm", formula = l ~ a * (c1 + c2 + c3 + c4), args = list(family = "poisson")),
  list(func = "glm", formula = m ~ a * l + a * (c1 + c2 + c3 + c4) + l * (c1 + c2 + c3 + c4), args = list(family = "binomial")),
  list(func = "glm", formula = y ~ a * l + a * m + l * m + (a + l + m) * (c1 + c2 + c3 + c4), args = list(family = "binomial"))
)

# Path-specific Effects
medsim_core(data=mydata, num_sim = 2000,
            cat_list = c("0", "1"), treatment = "a",
            intv_med = NULL, model_spec = model_spec)

# Interventional Effects
medsim_core(data=mydata, num_sim = 2000,
            cat_list = c("0", "1"), treatment = "a",
            intv_med = ("m"), model_spec = model_spec)

