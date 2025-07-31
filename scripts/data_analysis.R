#---------------------------------------------------------------------
# Data import and workspace preparation
#---------------------------------------------------------------------

#Import the packages
library(mgcv)
library(FactoMineR)

#Import the data
source("scripts/data_import.R")


#---------------------------------------------------------------------
# Summary Statistics
#---------------------------------------------------------------------

#Range of concentrations in FF
range(total_mps$concentration)

#Mean concentration in FF
mean(total_mps$concentration)

#Standard deviation
sd(total_mps$concentration)

#Most prevalent polymers
aggregate(concentration ~ category, FUN = "mean", data = data[which(data$sample == "FF"),])



#---------------------------------------------------------------------
# MPs and pregnancy rate
#---------------------------------------------------------------------

#Remove two patients who opted to freeze their oocytes
total_mps_reduced <- total_mps[!total_mps$patient %in% c(55,59),]

#Build a "pregnancy boolean
total_mps_reduced$pregnancy <- 0
total_mps_reduced$pregnancy[total_mps_reduced$cycle_output == "pregnant"] <- 1

#Mean concentration in FF of patients who became pregnant
mean(total_mps_reduced$concentration[total_mps_reduced$cycle_output == "pregnant"])
sd(total_mps_reduced$concentration[total_mps_reduced$cycle_output == "pregnant"])


#Mean concentration in FF of patients who did not become pregnant
mean(total_mps_reduced$concentration[total_mps_reduced$cycle_output == "not pregnant"])
sd(total_mps_reduced$concentration[total_mps_reduced$cycle_output == "not pregnant"])


fit <- gam(concentration ~ cycle_output,
           data = total_mps_reduced,
           family = tw(),
           method = "REML")

summary(fit)




fit <- gam(pregnancy ~ concentration,
           data = total_mps_reduced,
           family = binomial(),
           method = "REML")

summary(fit)

#Show the odds
exp(coef(fit))

#---------------------------------------------------------------------
# Correlation between blood mps and FF mps
#---------------------------------------------------------------------

#Difference between blood and FF MP concentrations
t.test(total_mps_merged$concentration.blood,
       total_mps_merged$concentration.ff,
       paired = T)

#Not significantly different than normal
shapiro.test(total_mps_merged$concentration.blood - total_mps_merged$concentration.ff) 

fit <- gam(concentration.ff ~ concentration.blood,
           data = total_mps_merged,
           family = tw(),
           method = "REML")

summary(fit)




#Subset to only data from patients with paired blood and ff samples
paired_data <- data[data$ID %in% total_mps_merged$ID,]

t_test_results <- data.frame(category = unique(data$category),
                             V = NA,
                             #diff = NA,
                             p = NA,
                             significant = NA)

#changes in polymer concentrations between blood and ff
for(i in 1:length(unique(paired_data$category))){
  message("Results for: ", unique(paired_data$category)[i])
  
  data_subset <- paired_data[paired_data$category == unique(paired_data$category)[i],]
  
  data_subset <- reshape(data_subset[,c("ID","sample","concentration")],
                         idvar = "ID",
                         timevar = "sample",
                         direction = "wide")
  
  data_subset[is.na(data_subset)]<- 0
  
  
  #Difference between blood and FF MP concentrations
  t_test <- wilcox.test(data_subset$concentration.BL,
                        data_subset$concentration.FF,
                        paired = T)
  
  print(t_test)
  
  #Store results of glm in data frame
  t_test_results[t_test_results$category == unique(data$category)[i],"V"] <- unname(t_test$statistic)
  #t_test_results[t_test_results$category == unique(data$category)[i],"diff"] <- unname(t_test$estimate)
  t_test_results[t_test_results$category == unique(data$category)[i],"p"] <- unname(t_test$p.value)
  t_test_results[t_test_results$category == unique(data$category)[i],"significant"] <- ifelse(unname(t_test$p.value)<0.05,1,0)
  
}


#Merge with polymer information
polymer_data <- merge(x = poly_props,
                      y = t_test_results)

polymer_data$significant <- as.factor(polymer_data$significant)
polymer_data$hydrophobic <- as.factor(polymer_data$hydrophobic)

#Tests looking at polymer properties and a significant reduction in FF vs. blood
#Polymer density
fit <- gam(density ~ significant,
           data = polymer_data,
           family = tw(),
           method = "REML")

summary(fit)

#Polymer elasticity
fit <- gam(elasticity ~ significant,
           data = polymer_data,
           family = tw(),
           method = "REML")

summary(fit)

#Polymer hydrophobic/philic
fit <- gam(hydrophobic ~ significant,
           data = polymer_data,
           family = binomial(),
           method = "REML")

summary(fit)



#Particle size - length
fit <- gam(Length..µm. ~ Sample,
           data = mp_props,
           family = tw(),
           method = "REML")

summary(fit)

#Particle size - Width
fit <- gam(Width..µm. ~ Sample,
           data = mp_props,
           family = tw(),
           method = "REML")

summary(fit)

#Particle size - Volume
fit <- gam(SE.Volume..µm.. ~ Sample,
           data = mp_props,
           family = tw(),
           method = "REML")

summary(fit)

#---------------------------------------------------------------------
# Supporting analyses
#---------------------------------------------------------------------

#Patient age
fit <- gam(patient_age ~ concentration,
           data = total_mps,
           family = tw(),
           method = "REML")

summary(fit)


#Patient BMI
fit <- gam(patient_BMI ~ concentration,
           data = total_mps,
           family = tw(),
           method = "REML")

summary(fit)


#Patient AMH
fit <- gam(therapy_last_AMH_value ~ concentration,
           data = total_mps,
           family = tw(),
           method = "REML")

summary(fit)

#Patient E2
fit <- gam(therapy_last_E2_value ~ concentration,
           data = total_mps,
           family = tw(),
           method = "REML")

summary(fit)

#Patient oocytes obtained
fit <- gam(obtained_oocytes ~ concentration,
           data = total_mps,
           family = poisson(),
           method = "REML")

summary(fit)

#Patient oocytes that reached maturity
fit <- gam(MII_oocytes/obtained_oocytes ~ concentration,
           data = total_mps,
           family = betar(),
           method = "REML")

summary(fit)

#Patient oocytes that did not reach maturity
fit <- gam((GV_oocytes+MI_oocytes)/obtained_oocytes ~ concentration,
           data = total_mps,
           family = betar(),
           method = "REML")

summary(fit)

#Patient oocytes that were fertilized
fit <- gam(X2PNs/MII_oocytes ~ concentration,
           data = total_mps,
           family = betar(),
           method = "REML")

summary(fit)


#Patient rate of abnormal or unfertilised oocytes
fit <- gam((MII_oocytes - X2PNs)/MII_oocytes ~ concentration,
           data = total_mps,
           family = betar(),
           method = "REML")

summary(fit)


