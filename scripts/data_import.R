#Import the polymer concentration data
data <- read.csv("data/polymer_concentrations.csv")

#Some processing of the polymer names
data[which(data$category == "Polyvinyl Chloride"),"category"] <- "PVC"
data$ID <- as.factor(data$ID)
data$category <- as.factor(data$category)

#Import the demographic data
demographics <- read.csv("data/demographics.csv")
demographics$ID <- as.factor(demographics$ID)

#---------------------------------------------------------------------
# Total concentration per patient
#---------------------------------------------------------------------

#Count totals in each patient/sample
total_mps <- aggregate(cbind(concentration, count) ~ ID + sample, FUN = "sum", data = data[which(data$sample == "FF"),])

#Merge the polymer counts with the demographic information
total_mps <- merge(x = total_mps,
                   y = demographics,
                   by.x = "ID",
                   by.y = "ID")

#---------------------------------------------------------------------
# Total concentration per patient - Blood
#---------------------------------------------------------------------

#Count totals in each patient/sample
total_mps_blood <- aggregate(cbind(concentration, count) ~ ID + sample, FUN = "sum", data = data[which(data$sample == "BL"),])

#Merge the polymer counts with the demographic information
total_mps_blood <- merge(x = total_mps_blood,
                         y = demographics,
                         by.x = "ID",
                         by.y = "ID")


#---------------------------------------------------------------------
# Paired concentrations blood and FF
#---------------------------------------------------------------------

#Merge the polymer counts with the demographic information
total_mps_merged <- merge(x = total_mps[,c("ID", "concentration")],
                          y = total_mps_blood[,c("ID", "concentration")],
                          by.x = "ID",
                          by.y = "ID")

names(total_mps_merged) <- c("ID", "concentration.ff", "concentration.blood")

#---------------------------------------------------------------------
# Paired concentrations blood and FF - full polymers
#---------------------------------------------------------------------

ff_data <- data[which(data$sample == "FF"),c("ID", "category", "concentration")]
names(ff_data) <- c("ID", "category", "concentration.ff")

blood_data <- data[which(data$sample == "BL"),c("ID", "category", "concentration")]
names(blood_data) <- c("ID", "category", "concentration.blood")

#Merge the polymer counts with the demographic information
polymers_merged <- merge(x = ff_data,
                          y = blood_data,
                          by.x = c("ID", "category"),
                          by.y = c("ID", "category"))

rm(ff_data); rm(blood_data)


#---------------------------------------------------------------------
# Polymer properties
#---------------------------------------------------------------------

poly_props <- read.csv("data/polymers_properties.csv")
poly_props$hydrophobic <- as.numeric(poly_props$hydrophobic )
poly_props[which(poly_props$category == "Polyvinyl Chloride"),"category"] <- "PVC"

#---------------------------------------------------------------------
# Particle properties
#---------------------------------------------------------------------

mp_props <- read.csv("data/particle_details.csv")
