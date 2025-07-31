#---------------------------------------------------------------------
# Data import and workspace preparation
#---------------------------------------------------------------------

#Import the packages
library(ggplot2)
library(ggpubr)
library(mgcv)
library(FactoMineR)

#Import the data
source("scripts/data_import.R")



#Subset to only data from patients with paired blood and ff samples
paired_data <- data[data$ID %in% total_mps_merged$ID,]

t_test_results <- data.frame(category = unique(data$category),
                             t = NA,
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
  t_test_results[t_test_results$category == unique(data$category)[i],"t"] <- unname(t_test$statistic)
  #t_test_results[t_test_results$category == unique(data$category)[i],"diff"] <- unname(t_test$estimate)
  t_test_results[t_test_results$category == unique(data$category)[i],"p"] <- unname(t_test$p.value)
  t_test_results[t_test_results$category == unique(data$category)[i],"significant"] <- ifelse(unname(t_test$p.value)<0.05,1,0)
  
}



res.pca <- PCA(poly_props[1:13,2:4]) # Conduct a PCA on the polymer propertiy information
PC1 <- res.pca$ind$coord[,1] #Store individual coordinates of PC1 as a vector
PC2 <- res.pca$ind$coord[,2] #Store individual coordinates of PC2 as a vector
PCs.ID <- data.frame(cbind(PC1,PC2)) #Bind the coordinates together as a dataframe
PCs.ID$category <- poly_props[1:13,1] #Add in polymer ID to the df as well as regression results
PCs.ID <- merge(x = PCs.ID,
                y = t_test_results)
PCs.ID$significant <- as.factor(PCs.ID$significant)

#Define axis labels based on % of data explained across each dimension of the PCA
DIM_1 <- paste("Dim 1 (", round(res.pca$eig[1,2], 1), "%)")
DIM_2 <- paste("Dim 2 (", round(res.pca$eig[2,2], 1), "%)")


#Then make the figure
a <- 
  ggplot(PCs.ID, aes(x=PC1, y=PC2, color = significant, fill = significant)) +
  #ggtitle("a") +
  geom_hline(aes(yintercept=0), linetype="dashed", lwd = 0.1, col = "grey70", alpha = 0.8) +
  geom_vline(aes(xintercept=0), linetype="dashed", lwd = 0.1, col = "grey70", alpha = 0.8) +
  stat_ellipse(geom = "polygon", alpha = 0.2, segments = 200, linewidth = 0.1, show.legend = FALSE) +
  theme_bw() +
  ylab(DIM_2) +
  xlab(DIM_1) + 
  geom_point(size=0.5, aes(color = significant)) +
  scale_fill_manual(values = c("#e56b6f","#355070"),
                    guide = "none") +
  scale_colour_manual(values = c("#e56b6f","#355070"),
                      guide = "none") +
  geom_segment(aes(x = 0, y = 0,
                   xend = res.pca$var$coord["density","Dim.1"],
                   yend = res.pca$var$coord["density","Dim.2"]),
               arrow = arrow(length = unit(0.2, "cm")),
               col = "black") +
  annotate(geom = "text",
           x = res.pca$var$coord["density","Dim.1"]*1.2,
           y = res.pca$var$coord["density","Dim.2"]*1.2, label = "Density") +
  geom_segment(aes(x = 0, y = 0,
                   xend = res.pca$var$coord["hydrophobic","Dim.1"],
                   yend = res.pca$var$coord["hydrophobic","Dim.2"]),
               arrow = arrow(length = unit(0.2, "cm")),
               col = "black") +
  annotate(geom = "text",
           x = res.pca$var$coord["hydrophobic","Dim.1"]*1.2,
           y = res.pca$var$coord["hydrophobic","Dim.2"]*1.2, label = "Hydrophobic") +
  geom_segment(aes(x = 0, y = 0,
                   xend = res.pca$var$coord["elasticity","Dim.1"],
                   yend = res.pca$var$coord["elasticity","Dim.2"]),
               arrow = arrow(length = unit(0.2, "cm")),
               col = "black") +
  annotate(geom = "text",
           x = res.pca$var$coord["elasticity","Dim.1"]*1.5,
           y = res.pca$var$coord["elasticity","Dim.2"]*1.5, label = "Elasticity") +
  ylab(DIM_2) +
  xlab(DIM_1) + 
  theme(panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", size=1),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0, size = 12, family = "sans", face = "bold"),
        axis.title.x  = element_text(size=10, family = "sans", face = "bold"),
        axis.title.y  = element_text(size=10, family = "sans", face = "bold"),
        axis.text.y  = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=8, family = "sans"),
        legend.position=c(0.8,0.9),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=6, family = "sans"),
        legend.key.size = unit(0.3, "cm"),
        legend.key = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA))


ggsave(a,
       file="figures/figure_S3.png",
       width = 6.86,
       height=6,
       units = "in",
       dpi = 600)
