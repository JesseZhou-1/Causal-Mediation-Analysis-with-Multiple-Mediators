cleaned_natl2003 <- read_csv("cleaned_natl2003.csv")
cleaned_natl2003 <- na.omit(cleaned_natl2003)

cleaned_natl2003_bin <- subset(cleaned_natl2003, !(precare == 1 | precare == 3))
cleaned_natl2003_bin$precare <- ifelse(cleaned_natl2003_bin$precare == 2, 1, cleaned_natl2003_bin$precare)

# Rename columns and create a new data frame
mydata <- cleaned_natl2003_bin %>%
  rename(y = prebirth, a = precare, m = urf_eclam, l = tobuse, c1 = age,c2 = somecollege, c3 = mracerec, c4 = mar)

mydata$y <- as.factor(mydata$y)
mydata$a <- as.factor(mydata$a)
mydata$m <- as.factor(mydata$m)
mydata$l <- as.factor(mydata$l)
mydata$c1 <- as.factor(mydata$c1)
mydata$c2 <- as.factor(mydata$c2)
mydata$c3 <- as.factor(mydata$c3)
mydata$c4 <- as.factor(mydata$c4)

# Create modified datasets with different values of 'a' and 'l'
mydata0 <- mydata
mydata0$a <- 0
mydata0$a <- as.factor(mydata0$a)

mydata1 <- mydata
mydata1$a <- 1
mydata1$a <- as.factor(mydata1$a)

# Fit logistic model for 'a' on 'c' and predict 'pa1'
fit_a <- glm(a ~ c1 + c2 + c3 + c4, data=mydata, family=binomial(link="logit"))
mydata$pa1 <- predict(fit_a, newdata=mydata, type="response")

# Fit logistic model for 'l' on 'a' and 'c', then predict for mydata1 and mydata0
fit_l <- glm(l ~ a + c1 + c2 + c3 + c4, data=mydata, family=binomial(link="logit"))
mydata$pl1 <- predict(fit_l, newdata=mydata1, type="response")
mydata$pl10 <- predict(fit_l, newdata=mydata0, type="response")

# Modify mydata for all combinations of 'a' and 'l'
mydata00 <- mydata
mydata00$a <- 0
mydata00$l <- 0
mydata00$a <- as.factor(mydata00$a)
mydata00$l <- as.factor(mydata00$l)

mydata10 <- mydata
mydata10$a <- 1
mydata10$l <- 0
mydata10$a <- as.factor(mydata10$a)
mydata10$l <- as.factor(mydata10$l)

mydata01 <- mydata
mydata01$a <- 0
mydata01$l <- 1
mydata01$a <- as.factor(mydata01$a)
mydata01$l <- as.factor(mydata01$l)

mydata11 <- mydata
mydata11$a <- 1
mydata11$l <- 1
mydata11$a <- as.factor(mydata11$a)
mydata11$l <- as.factor(mydata11$l)

# Assuming logistic model for 'm' has already been fitted on 'mydata'
fit_m <- glm(m ~ a + l + c1 + c2 + c3 + c4, data=mydata, family=binomial(link="logit"))

mydataw <- mydata

# Score data for different scenarios and create predictions
mydataw$pm1 <- predict(fit_m, newdata=mydata1, type="response")
mydataw$pm10 <- predict(fit_m, newdata=mydata0, type="response")
mydataw$pm100 <- predict(fit_m, newdata=mydata00, type="response")
mydataw$pm101 <- predict(fit_m, newdata=mydata01, type="response")
mydataw$pm110 <- predict(fit_m, newdata=mydata10, type="response")
mydataw$pm111 <- predict(fit_m, newdata=mydata11, type="response")

# Duplicate the dataset for two scenarios
mydatanew1 <- mydataw
mydatanew2 <- mydataw

# Initialize w1 with NA
mydatanew1$w1 <- NA
mydatanew2$w1 <- NA

# Ensure 'a', 'm', and 'l' are numeric for calculations (if they're not already)
mydatanew1$a <- as.numeric(as.character(mydatanew1$a))
mydatanew1$m <- as.numeric(as.character(mydatanew1$m))
mydatanew1$l <- as.numeric(as.character(mydatanew1$l))
mydatanew2$a <- as.numeric(as.character(mydatanew2$a))
mydatanew2$m <- as.numeric(as.character(mydatanew2$m))
mydatanew2$l <- as.numeric(as.character(mydatanew2$l))

# Initial setup for astar equal to a in mydatanew1 and calculate initial w1
mydatanew1$astar <- mydatanew1$a
mydatanew1$w1 <- with(mydatanew1, a / pa1 + (1 - a) / (1 - pa1))

# In mydatanew2, invert a into astar, then recalculate w1 based on the value of a
mydatanew2$astar <- 1 - mydatanew2$a

# Conditional recalculation of w1 in mydatanew2 based on the inverted astar
mydatanew2$w1[mydatanew2$a == 0] <- with(mydatanew2[mydatanew2$a == 0, ],
                                         (1/(1-pa1)) * (l*pl1/pl10 + (1-l)*(1-pl1)/(1-pl10)) * ((1-m)*(1-pm1)/(1-pm10) + m*pm1/pm10))

mydatanew2$w1[mydatanew2$a == 1] <- with(mydatanew2[mydatanew2$a == 1, ],
                                         (1/pa1) * ((l*pl10/pl1 + (1-l)*(1-pl10)/(1-pl1)) * ((1-m)*(1-pm10)/(1-pm1) + m*pm10/pm1)))

mydatanew <- rbind(mydatanew1, mydatanew2)

# Initialize w3 with NA
mydatanew$w3 <- NA

# Apply conditions directly for w3 calculation

# Condition 1: if (a = 0) & (astar = 0) & (m = 1)
index <- which(mydatanew$a == 0 & mydatanew$astar == 0 & mydatanew$m == 1)
mydatanew$w3[index] <- (1 / (1 - mydatanew$pa1[index])) * (mydatanew$pm100[index] * (1 - mydatanew$pl10[index]) + mydatanew$pm101[index] * mydatanew$pl10[index]) / mydatanew$pm10[index]

# Condition 2: if (a = 0) & (astar = 0) & (m = 0)
index <- which(mydatanew$a == 0 & mydatanew$astar == 0 & mydatanew$m == 0)
mydatanew$w3[index] <- (1 / (1 - mydatanew$pa1[index])) * ((1 - mydatanew$pm100[index]) * (1 - mydatanew$pl10[index]) + (1 - mydatanew$pm101[index]) * mydatanew$pl10[index]) / (1 - mydatanew$pm10[index])

# Condition 3: if (a = 0) & (astar = 1) & (m = 1)
index <- which(mydatanew$a == 0 & mydatanew$astar == 1 & mydatanew$m == 1)
mydatanew$w3[index] <- (1 / (1 - mydatanew$pa1[index])) * (mydatanew$pm110[index] * (1 - mydatanew$pl1[index]) + mydatanew$pm111[index] * mydatanew$pl1[index]) / mydatanew$pm10[index]

# Condition 4: if (a = 0) & (astar = 1) & (m = 0)
index <- which(mydatanew$a == 0 & mydatanew$astar == 1 & mydatanew$m == 0)
mydatanew$w3[index] <- (1 / (1 - mydatanew$pa1[index])) * ((1 - mydatanew$pm110[index]) * (1 - mydatanew$pl1[index]) + (1 - mydatanew$pm111[index]) * mydatanew$pl1[index]) / (1 - mydatanew$pm10[index])

# Condition 5: if (a = 1) & (astar = 0) & (m = 1)
index <- which(mydatanew$a == 1 & mydatanew$astar == 0 & mydatanew$m == 1)
mydatanew$w3[index] <- (1 / mydatanew$pa1[index]) * (mydatanew$pm100[index] * (1 - mydatanew$pl10[index]) + mydatanew$pm101[index] * mydatanew$pl10[index]) / mydatanew$pm1[index]

# Condition 6: if (a = 1) & (astar = 0) & (m = 0)
index <- which(mydatanew$a == 1 & mydatanew$astar == 0 & mydatanew$m == 0)
mydatanew$w3[index] <- (1 / mydatanew$pa1[index]) * ((1 - mydatanew$pm100[index]) * (1 - mydatanew$pl10[index]) + (1 - mydatanew$pm101[index]) * mydatanew$pl10[index]) / (1 - mydatanew$pm1[index])

# Condition 7: if (a = 1) & (astar = 1) & (m = 1)
index <- which(mydatanew$a == 1 & mydatanew$astar == 1 & mydatanew$m == 1)
mydatanew$w3[index] <- (1 / mydatanew$pa1[index]) * (mydatanew$pm110[index] * (1 - mydatanew$pl1[index]) + mydatanew$pm111[index] * mydatanew$pl1[index]) / mydatanew$pm1[index]

# Condition 8: if (a = 1) & (astar = 1) & (m = 0)
index <- which(mydatanew$a == 1 & mydatanew$astar == 1 & mydatanew$m == 0)
mydatanew$w3[index] <- (1 / mydatanew$pa1[index]) * ((1 - mydatanew$pm110[index]) * (1 - mydatanew$pl1[index]) + (1 - mydatanew$pm111[index]) * mydatanew$pl1[index]) / (1 - mydatanew$pm1[index])


mydatanew$a <- as.factor(mydatanew$a)
mydatanew$astar <- as.factor(mydatanew$astar)
mydatanew$m <- as.factor(mydatanew$m)
mydatanew$l <- as.factor(mydatanew$l)

# Multivariate Natural Direct Effects
fit1 <- glm(y ~ a, data = mydatanew[mydatanew$astar == 0, ], family = "binomial", weights = mydatanew$w1[mydatanew$astar == 0])

summary(fit1)

# Multivariate Natural Indirect Effects
fit2 <- glm(y ~ astar, data = mydatanew[mydatanew$a == 1, ], family = "binomial", weights = mydatanew$w1[mydatanew$a == 1])

summary(fit2)

# Interventional Direct Effects
fit3 <- glm(y ~ a, data = mydatanew[mydatanew$astar == 0, ], family = "binomial", weights = mydatanew$w3[mydatanew$astar == 0])

summary(fit3)

# Interventional Indirect Effects
fit4 <- glm(y ~ astar, data = mydatanew[mydatanew$a == 1, ], family = "binomial", weights = mydatanew$w3[mydatanew$a == 1])

summary(fit4)



# Duplicate the dataset for three scenarios
mydatanew3 <- mydataw
mydatanew4 <- mydataw
mydatanew5 <- mydataw


# Initialize w1 with NA
mydatanew3$w2 <- NA
mydatanew4$w2 <- NA
mydatanew5$w2 <- NA

# Ensure 'a', 'm', and 'l' are numeric for calculations (if they're not already)
mydatanew3$a <- as.numeric(as.character(mydatanew3$a))
mydatanew3$m <- as.numeric(as.character(mydatanew3$m))
mydatanew3$l <- as.numeric(as.character(mydatanew3$l))
mydatanew4$a <- as.numeric(as.character(mydatanew4$a))
mydatanew4$m <- as.numeric(as.character(mydatanew4$m))
mydatanew4$l <- as.numeric(as.character(mydatanew4$l))
mydatanew5$a <- as.numeric(as.character(mydatanew5$a))
mydatanew5$m <- as.numeric(as.character(mydatanew5$m))
mydatanew5$l <- as.numeric(as.character(mydatanew5$l))

# Condition 1: if (astar = a) & (astarstar = a)
mydatanew3$astar <- mydatanew3$a
mydatanew3$astarstar <- mydatanew3$a
mydatanew3$w2 <- with(mydatanew1, a / pa1 + (1 - a) / (1 - pa1))

# Condition 2: if (astar = 1-a) & (astarstar = a)
mydatanew4$astar <- 1 - mydatanew4$a
mydatanew4$astarstar <- mydatanew4$a
mydatanew4$w2 <- with(mydatanew4, (a / pa1) * (l * pl10 / pl1 + (1 - l) * (1 - pl10) / (1 - pl1)) +
                        ((1 - a) / (1 - pa1)) * (l * pl1 / pl10 + (1 - l) * (1 - pl1) / (1 - pl10)))

# Condition 3: if (astar = 1-a) & (astarstar = a)
mydatanew5$astar <- 1 - mydatanew5$a
mydatanew5$astarstar <- 1 - mydatanew5$a
mydatanew5$w2 <- with(mydatanew5, (a / pa1) * (l * (pl10 / pl1) * (m * (pm101 / pm111) + (1 - m) * (1 - pm101) / (1 - pm111)) +
                                                 (1 - l) * ((1 - pl10) / (1 - pl1)) * (m * (pm100 / pm110) + (1 - m) * (1 - pm100) / (1 - pm110))) +
                        ((1 - a) / (1 - pa1)) * (l * (pl1 / pl10) * (m * (pm111 / pm101) + (1 - m) * (1 - pm111) / (1 - pm101)) +
                                                   (1 - l) * ((1 - pl1) / (1 - pl10)) * (m * (pm110 / pm100) + (1 - m) * (1 - pm110) / (1 - pm100))))

mydatanewnew <- rbind(mydatanew3, mydatanew4, mydatanew5)
mydatanewnew$a <- as.factor(mydatanewnew$a)
mydatanewnew$astar <- as.factor(mydatanewnew$astar)
mydatanewnew$astarstar <- as.factor(mydatanewnew$astarstar)
mydatanewnew$m <- as.factor(mydatanewnew$m)
mydatanewnew$l <- as.factor(mydatanewnew$l)

# PSE(A -> Y)
fit5 <- glm(y ~ a, data = mydatanewnew[mydatanewnew$astar == 0 & mydatanewnew$astarstar == 0, ], family = "binomial", weights = mydatanewnew$w2[mydatanewnew$astar == 0 & mydatanewnew$astarstar == 0])
summary(fit5)

# PSE(A -> L ~> Y)
fit6 <- glm(y ~ astar, data = mydatanewnew[mydatanewnew$a == 1 & mydatanewnew$astarstar == 1, ], family = "binomial", weights = mydatanewnew$w2[mydatanewnew$a == 1 & mydatanewnew$astarstar == 1])
summary(fit6)

# PSE(A -> M -> Y)
fit7 <- glm(y ~ astarstar, data = mydatanewnew[mydatanewnew$a == 1 & mydatanewnew$astar == 0, ], family = "binomial", weights = mydatanewnew$w2[mydatanewnew$a == 1 & mydatanewnew$astar == 0])
summary(fit7)


