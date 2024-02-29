# REPLICATION

cleaned_natl2003 <- read_csv("cleaned_natl2003.csv")
cleaned_natl2003 <- na.omit(cleaned_natl2003)

# Rename columns and create a new data frame
mydata <- cleaned_natl2003 %>%
  rename(y = prebirth, a = precare, m = urf_eclam, l = tobuse, c1 = age,c2 = somecollege, c3 = mracerec, c4 = mar)

mydata2 <- cleaned_natl2003 %>%
  rename(y = prebirth, a = precare_4cat, m = urf_eclam, l = cigs, c1 = age, c2 = somecollege, c3 = mracerec, c4 = mar)

mydata$y <- as.factor(mydata$y)
mydata$a <- as.factor(mydata$a)
mydata$m <- as.factor(mydata$m)
mydata$l <- as.factor(mydata$l)
mydata$c1 <- as.factor(mydata$c1)
mydata$c2 <- as.factor(mydata$c2)
mydata$c3 <- as.factor(mydata$c3)
mydata$c4 <- as.factor(mydata$c4)

# Define the model specifications based on your SEM
model_spec <- list(
  list(func = "glm", formula = l ~ a + c1 + c2 + c3 + c4, args = list(family = "binomial")),
  list(func = "glm", formula = m ~ a + l + c1 + c2 + c3 + c4, args = list(family = "binomial")),
  list(func = "glm", formula = y ~ a + l + m + c1 + c2 + c3 + c4, args = list(family = "binomial"))
)

model_spec_2 <- list(
  list(func = "glm", formula = l ~ a * (c1 + c2 + c3 + c4), args = list(family = "binomial")),
  list(func = "glm", formula = m ~ (a + l) * (c1 + c2 + c3 + c4), args = list(family = "binomial")),
  list(func = "glm", formula = y ~ (a + l + m) * (c1 + c2 + c3 + c4), args = list(family = "binomial"))
)

# Call the medsim function
medsim_core(data=mydata, num_sim = 2000,
            cat_list = c("0", "1"), treatment = "a",
            intv_med = NULL, model_spec = model_spec)

