library(readr)
library(dplyr)


cleaned_natl2003 <- read_csv("cleaned_natl2003.csv")
cleaned_natl2003 <- na.omit(cleaned_natl2003)


mydata2 <- cleaned_natl2003 %>%
  rename(y = prebirth, a = precare, m = urf_eclam, l = cigs, c1 = age, c2 = somecollege, c3 = mracerec, c4 = mar)


mydata2$y <- as.factor(mydata2$y)
mydata2$a <- as.factor(mydata2$a)
mydata2$m <- as.factor(mydata2$m)
mydata2$c1 <- as.factor(mydata2$c1)
mydata2$c2 <- as.factor(mydata2$c2)
mydata2$c3 <- as.factor(mydata2$c3)
mydata2$c4 <- as.factor(mydata2$c4)

model_spec_2_1 <- list(
  list(func = "glm", formula = l ~ a + c1 + c2 + c3 + c4, args = list(family = "poisson")),
  list(func = "glm", formula = m ~ a + l + c1 + c2 + c3 + c4, args = list(family = "binomial")),
  list(func = "glm", formula = y ~ a + l + m + c1 + c2 + c3 + c4, args = list(family = "binomial"))
)

model_spec_2_2 <- list(
  list(func = "glm", formula = l ~ a * (c1 + c2 + c3 + c4), args = list(family = "poisson")),
  list(func = "glm", formula = m ~ a * (c1 + c2 + c3 + c4) + l * (c1 + c2 + c3 + c4), args = list(family = "binomial")),
  list(func = "glm", formula = y ~ (a + l + m) * (c1 + c2 + c3 + c4), args = list(family = "binomial"))
)

# Path-specfic Effects
# Additive
medsim_core(data=mydata2, num_sim = 2000,
            cat_list = c("0", "1"), treatment = "a",
            intv_med = NULL, model_spec = model_spec_2_1)

# All Two-way Interactions
medsim_core(data=mydata2, num_sim = 2000,
            cat_list = c("0", "1"), treatment = "a",
            intv_med = NULL, model_spec = model_spec_2_2)

# Interventional Effects
# Additive
medsim_core(data=mydata2, num_sim = 2000,
            cat_list = c("0", "1"), treatment = "a",
            intv_med = ("m"), model_spec = model_spec_2_1)

# All Two-way Interactions
medsim_core(data=mydata2, num_sim = 2000,
            cat_list = c("0", "1"), treatment = "a",
            intv_med = ("m"), model_spec = model_spec_2_2)


