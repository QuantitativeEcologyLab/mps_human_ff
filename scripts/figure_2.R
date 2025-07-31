#---------------------------------------------------------------------
# Data import and workspace preparation
#---------------------------------------------------------------------

#Import the packages
library(ggplot2)
library(ggpubr)
library(mgcv)

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


a <- 
  ggplot(data = total_mps_merged,
         aes(x = concentration.blood,
             y = concentration.ff)) +
  ggtitle("a") +
  geom_smooth(method = "gam", formula = y ~ x, method.args = list(family = "tw"), alpha = 0.2, col = "black") +
  geom_point(size = 1, shape = 16) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=10, family = "sans", face = "bold"),
        axis.title.x = element_text(size=10, family = "sans", face = "bold"),
        axis.text.y = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=8, family = "sans"),
        plot.title = element_text(hjust = -0.03, size = 12, family = "sans", face = "bold"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=8, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  #scale_y_continuous(limits = c(0,580), expand = c(0,6)) +
  xlab(expression(bold(Blood~MP~Concentration~(Particles~mL^-1)))) +
  ylab(expression(bold(Follicular~Fluid~MP~Concentration~(Particles~mL^-1)))) +
  guides(col = guide_legend(ncol = 8, nrow = 2, byrow = TRUE))




#Generate the figure

test <- rbind(total_mps[total_mps$ID %in% total_mps_merged$ID,],
              total_mps_blood)

b <-
  ggplot(data = test,
         aes(x = sample,
             y = concentration,
             fill = sample)) +
  ggtitle("b") +
  geom_boxplot(size = 0.1, outlier.size = 0.2, outlier.shape = 16, outlier.alpha = 0, alpha = 0.7) +
  stat_summary(fun.y=mean, geom="point", shape=18, size=2, color="black", fill="black") +
  geom_jitter(aes(col = sample), size = 1.5, shape = 16, position=position_jitter(height=0, width=0)) +
  geom_line(aes(group=patient), linewidth = 0.3, col = "grey20", alpha = 0.6) +
  scale_fill_manual(values = c("#edae49", "#00798c")) +
  scale_colour_manual(values = c("#edae49", "#00798c")) +
  ylab(expression(bold(MP~Concentration~(Particles~mL^-1)))) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=10, family = "sans", face = "bold"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=10, family = "sans", face = "bold", color = "black"),
        plot.title = element_text(hjust = -0.05, size = 14, family = "sans", face = "bold"),
        strip.text.x = element_text(size=10, family = "sans", face = "bold", color = "black"),
        #strip.background = element_blank(),
        #strip.text.x = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=5, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.2, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  scale_x_discrete(breaks = c("BL", "FF"), labels = c("Blood","Follicular fluid"))


#Generate the figure
polymer_data <- merge(x = poly_props,
                      y = t_test_results)

polymer_data$significant <- as.factor(polymer_data$significant)
polymer_data$hydrophobic <- as.factor(polymer_data$hydrophobic)

c <- 
  ggplot(data = polymer_data, 
         aes(x = hydrophobic,
             y = significant,
             fill = significant)) +
  ggtitle("c") +
  geom_bar(stat = "identity", position="fill", alpha = 0.8) +
  scale_fill_manual(labels = c("Not reduced in FF", "Reduced in FF"),
                    values = c("#e56b6f","#355070")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=10, family = "sans", face = "bold"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=10, family = "sans", face = "bold", colour = "black"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=8, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  scale_y_discrete(limits = c(0,1), expand = c(0,0)) +
  scale_x_discrete(breaks = c("0", "1"), labels = c("Hydrophilic", "Hydrophobic")) +
  ylab(expression(bold(Proportion~of~polymers)))

d <-
  ggplot(data = polymer_data,
         aes(x = significant,
             y = density)) +
  ggtitle("d") +
  geom_boxplot(size = 0.1, outlier.size = 0.2, outlier.shape = 16, outlier.alpha = 0,
               fill = c("#e56b6f", "#355070"), alpha = 0.8) +
  stat_summary(fun.y=mean, geom="point", shape=18, size=2, color="black", fill="black") +
  geom_jitter(aes(col = significant), size = 0.5, shape = 16, position=position_jitter(height=0, width=0.1)) +
  scale_fill_manual(values = c("#e56b6f", "#355070"), labels = c("Not reduced", "Reduced")) +
  scale_colour_manual(values = c("#e56b6f", "#355070"), labels = c("Not reduced", "Reduced")) +
  ylab(expression(bold(Polymer~density~(g/cm^3)))) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=10, family = "sans", face = "bold"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=10, family = "sans", face = "bold", color = "black"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        #strip.text.x = element_text(size=6, family = "sans", face = "bold", color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=5, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.2, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  scale_x_discrete(breaks = c("0", "1"), labels = c("Not reduced in FF", "Reduced in FF"))


e <-
  ggplot(data = polymer_data,
         aes(x = significant,
             y = elasticity)) +
  ggtitle("e") +
  geom_boxplot(size = 0.1, outlier.size = 0.2, outlier.shape = 16, outlier.alpha = 0,
               fill = c("#e56b6f", "#355070"), alpha = 0.8) +
  stat_summary(fun.y=mean, geom="point", shape=18, size=2, color="black", fill="black") +
  geom_jitter(aes(col = significant), size = 0.5, shape = 16, position=position_jitter(height=0, width=0.1)) +
  scale_fill_manual(values = c("#e56b6f", "#355070"), labels = c("Not reduced", "Reduced")) +
  scale_colour_manual(values = c("#e56b6f", "#355070"), labels = c("Not reduced", "Reduced")) +
  ylab(expression(bold(Polymer~elasticity~(MPa)))) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=10, family = "sans", face = "bold"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=10, family = "sans", face = "bold", color = "black"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        #strip.text.x = element_text(size=6, family = "sans", face = "bold", color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=5, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.2, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  scale_y_log10() +
  scale_x_discrete(breaks = c("0", "1"), labels = c("Not reduced in FF", "Reduced in FF"))









f <- 
  ggplot() +
  ggtitle("f") +
  geom_histogram(data = mp_props, aes(Length..µm., fill = Sample),
                 #alpha = 0.8,
                 bins = 60,
                 col = "black",
                 linewidth = 0.05) +
  scale_fill_manual(values = c("#edae49", "#00798c"),
                    labels = c("Blood","Follicular Fluid")) +
  scale_x_log10(expand = c(0,0.01))+
  scale_y_continuous(limits = c(0,260), expand = c(0,1)) +
  ylab("Number of particles") +
  xlab(bquote(bold('Particle length '(µm)))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=10, family = "sans", face = "bold"),
        axis.title.x = element_text(size=10, family = "sans", face = "bold"),
        axis.text.y = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=8, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        legend.position = c(0.7, 0.8),
        legend.title = element_blank(),
        legend.text = element_text(size=8, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.2, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  guides(fill = guide_legend(byrow = TRUE))



g <- 
  ggplot() +
  ggtitle("g") +
  geom_histogram(data = mp_props, aes(Width..µm., fill = Sample),
                 #alpha = 0.8,
                 bins = 50,
                 col = "black",
                 linewidth = 0.05) +
  scale_fill_manual(values = c("#edae49", "#00798c"),
                    labels = c("Blood","Follicular Fluid")) +
  scale_x_log10(expand = c(0,0.01))+
  scale_y_continuous(limits = c(0,370), expand = c(0,1)) +
  ylab("Number of particles") +
  xlab(bquote(bold('Particle width '(µm)))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=10, family = "sans", face = "bold"),
        axis.title.x = element_text(size=10, family = "sans", face = "bold"),
        axis.text.y = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=8, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=8, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.2, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  guides(fill = guide_legend(byrow = TRUE))


top <- ggarrange(a,b,
                 ncol=2, nrow=1,
                 widths = c(0.6, 0.4))

mid <- ggarrange(c,d,e,
                   ncol=3, nrow=1)

bot <- ggarrange(f,g,
                 ncol=2, nrow=1)

FIG <- ggarrange(top, mid, bot,
                 ncol=1, nrow=3,
                 heights = c(0.6,0.4,0.4))

ggsave(FIG,
       file="figures/figure_2.png",
       width = 6.86*1.6,
       height=6*1.6,
       units = "in",
       dpi = 600)
