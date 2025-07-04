# Compute APNCU-1M Index based on Kotelchuck (1994), VanderWeele et al. (2009) and VanderWeele et al. (2014)
computeAPNCU1M <- function(mpcb, uprevis, combgest) {
  if (is.na(mpcb) || is.na(uprevis) || is.na(combgest)) {
    return(NA)
  } else {
    # Calculate the number of visits expected by ACOG schedule up to the gestational age
    expectedVisits <- ifelse(combgest <= 28, combgest %/% 4,
                             ifelse(combgest <= 36, 7 + (combgest - 28) %/% 2,
                                    11 + (combgest - 36)))

    # Adjust for the month prenatal care began (missed visits are not made up)
    visitsMissedDueToLateStart <- ifelse(mpcb <= 7, (mpcb - 1),
                                         6 + (mpcb - 7) * 2)

    adjustedExpectedVisits <- max(0, expectedVisits - 1 - visitsMissedDueToLateStart)

    # Compute the actual to expected visit ratio
    actualToExpectedRatio <- if(adjustedExpectedVisits > 0) uprevis / adjustedExpectedVisits else 0


    # Determine the APNCU Index category
    if (mpcb <= 4 && actualToExpectedRatio >= 1.1 && (uprevis - adjustedExpectedVisits) >= 2) {
      return("Adequate Plus")
    } else if (mpcb <= 4 && actualToExpectedRatio >= 0.8 && (uprevis - adjustedExpectedVisits) < 2) {
      return("Adequate")
    } else if (mpcb <= 4 && actualToExpectedRatio >= 0.5 && actualToExpectedRatio < 0.8) {
      return("Intermediate")
    } else {
      return("Inadequate")
    }

    # If there are no prenatal visits
    if (uprevis == 0) {
      return("Inadequate")
    }
  }
}

# Load the necessary library
library(readr)
library(dplyr)

# Read the data
natl2003 <- read_csv("linkco2003us_den.csv")

# Keep only selected variables
natl2003 <- natl2003 %>%
  select(mager41, umeduc, mracerec, mar, mpcb, uprevis, urf_eclam, combgest, cigs, mracehisp)

count(natl2003)
natl2003 <- na.omit(natl2003)
count(natl2003)

# Treat certain values as missing
natl2003$cigs[natl2003$cigs == 99] <- NA
natl2003$mar[natl2003$mar == 9] <- NA
natl2003$umeduc[natl2003$umeduc == 99] <- NA
natl2003$mpcb[natl2003$mpcb == 99] <- NA
natl2003$urf_eclam[natl2003$urf_eclam %in% c(8, 9)] <- NA
natl2003$combgest[natl2003$combgest == 99] <- NA

# Compute APNCU-1M Index
natl2003$precare <- apply(natl2003, 1, function(x) computeAPNCU1M(x['mpcb'], x['uprevis'], x['combgest']))

print(sum(is.na(natl2003$precare)))

# Compute prebirth variable
natl2003$prebirth <- ifelse(natl2003$combgest < 37, 1, 0)

natl2003$somecollege <- ifelse(natl2003$umeduc < 13, 0, 1)

natl2003$mracerec <- ifelse(natl2003$mracehisp %in% 1:5, 5, natl2003$mracerec)

natl2003$age <- ifelse(natl2003$mager41 < 20, 1,
                              ifelse(natl2003$mager41 <= 35, 2,
                                     3))

# View the modified data frame
head(natl2003)

precare_proportions <- table(natl2003$precare) / nrow(natl2003) * 100
print(precare_proportions)
print(sum(is.na(natl2003$precare))/nrow(natl2003) * 100)

# RECODE
natl2003$cigs <- natl2003$cigs <- as.numeric(as.character(
  cut(natl2003$cigs,
      breaks = c(-Inf, 0, 5, 10, 20, 30, Inf),
      labels = c(0, 1, 2, 3, 4, 5),
      right = TRUE)
))
natl2003$urf_eclam <- ifelse(natl2003$urf_eclam == 2, 0, natl2003$urf_eclam)
natl2003$mar <- ifelse(natl2003$mar == 2, 0, natl2003$mar)
natl2003$precare <- ifelse(natl2003$precare == "Inadequate", 0,
                           ifelse(natl2003$precare == "Intermediate", 1,
                                  ifelse(natl2003$precare == "Adequate", 2,
                                         ifelse(natl2003$precare == "Adequate Plus", 3, NA))))

count(natl2003)

# Save the cleaned data
write_csv(natl2003, "cleaned_natl2003.csv")

# For reanalysis
natl2003_bin <- subset(natl2003, !(precare == 1 | precare == 3))
natl2003_bin$precare <- ifelse(natl2003_bin$precare == 2, 1, natl2003_bin$precare)

# Rename columns and subset the data to keep only the required variables
natl2003_bin <- natl2003_bin %>%
  rename(y = prebirth, a = precare, m = urf_eclam, l = cigs, c1 = age,c2 = somecollege, c3 = mracerec, c4 = mar)

cleaned_natl2003_bin <- natl2003_bin[, c("y", "a", "m", "l", "c1", "c2", "c3", "c4")]

# Drop rows with missing values
cleaned_natl2003_bin <- na.omit(cleaned_natl2003_bin)

# Save the data
write_csv(cleaned_natl2003_bin, "cleaned_natl2003_bin.csv")

# References
#1. Kotelchuck, Milton. "Overview of adequacy of prenatal care utilization index." Unpublished manuscript, University of North Carolina at Chapel Hill, Department of Maternal and Child Health, Chapel Hill, NC (1994).
#2. VanderWeele, Tyler J., John D. Lantos, Juned Siddique, and Diane S. Lauderdale. "A comparison of four prenatal care indices in birth outcome models: comparable results for predicting small-for-gestational-age outcome but different results for preterm birth or infant mortality." Journal of clinical epidemiology 62, no. 4 (2009): 438-445.
#3. VanderWeele, Tyler J., Stijn Vansteelandt, and James M. Robins. "Effect decomposition in the presence of an exposure-induced mediator-outcome confounder." Epidemiology 25, no. 2 (2014): 300-306.
