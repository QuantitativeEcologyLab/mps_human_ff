#---------------------------------------------------------------------
# Data import and workspace preparation
#---------------------------------------------------------------------

#Import the packages
library(ggplot2)
library(ggpubr)
library(mgcv)

#Import the data
source("scripts/data_import.R")



#---------------------------------------------------------------------
# Figure S1a - MP concentrations and oocytes obtained
#---------------------------------------------------------------------


#Generate the figure
a <-
  ggplot(data = total_mps,
         aes(x = concentration,
             y = obtained_oocytes)) +
  ggtitle("a") +
  geom_smooth(method = "gam", formula = y ~ x, method.args = list(family = "tw"), col = "grey30", linetype = "solid", alpha = 0.2) +
  geom_point(size = 1, shape = 16) +
  ylab(expression(bold(Number~of~Oocytes~Obtained))) + 
  xlab(expression(bold(Total~MP~Concentration~(Particles~mL^-1)))) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=10, family = "sans", face = "bold"),
        axis.title.x = element_text(size=10, family = "sans", face = "bold"),
        axis.text.y = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=10, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 14, family = "sans", face = "bold"),
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
  scale_x_continuous(limits = c(0,580), expand = c(0,6))


#---------------------------------------------------------------------
# Figure S1b - MP concentrations and AMH
#---------------------------------------------------------------------

#Generate the figure
b <-
  ggplot(data = total_mps,
         aes(x = concentration,
             y = therapy_last_AMH_value)) +
  ggtitle("b") +
  geom_smooth(method = "gam", formula = y ~ x, method.args = list(family = "tw"), col = "grey30", linetype = "dashed", alpha = 0.2) +
  geom_point(size = 1, shape = 16) +
  ylab(expression(bold(Therapy~last~AMH~value~(ng/mL)))) + 
  xlab(expression(bold(Total~MP~Concentration~(Particles~mL^-1)))) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=10, family = "sans", face = "bold"),
        axis.title.x = element_text(size=10, family = "sans", face = "bold"),
        axis.text.y = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=10, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 14, family = "sans", face = "bold"),
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
  scale_x_continuous(limits = c(0,580), expand = c(0,6))


#---------------------------------------------------------------------
# Figure S1c - MP concentrations and Age
#---------------------------------------------------------------------

#Generate the figure
c <-
  ggplot(data = total_mps,
         aes(x = concentration,
             y = patient_age)) +
  ggtitle("c") +
  geom_smooth(method = "gam", formula = y ~ x, method.args = list(family = "tw"), col = "grey30", linetype = "dashed", alpha = 0.2) +
  geom_point(size = 1, shape = 16) +
  ylab(expression(bold(Patient~Age~(Years)))) +
  xlab(expression(bold(Total~MP~Concentration~(Particles~mL^-1)))) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=10, family = "sans", face = "bold"),
        axis.title.x = element_text(size=10, family = "sans", face = "bold"),
        axis.text.y = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=10, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 14, family = "sans", face = "bold"),
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
  scale_x_continuous(limits = c(0,580), expand = c(0,6))



#---------------------------------------------------------------------
# Figure S1d - MP concentrations and BMI
#---------------------------------------------------------------------

#Generate the figure
d <-
  ggplot(data = total_mps,
         aes(x = concentration,
             y = patient_BMI)) +
  ggtitle("d") +
  geom_smooth(method = "gam", formula = y ~ x, method.args = list(family = "tw"), col = "grey30", linetype = "dashed", alpha = 0.2) +
  geom_point(size = 1, shape = 16) +
  ylab(expression(bold(Patient~BMI))) +
  xlab(expression(bold(Total~MP~Concentration~(Particles~mL^-1)))) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=10, family = "sans", face = "bold"),
        axis.title.x = element_text(size=10, family = "sans", face = "bold"),
        axis.text.y = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=10, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 14, family = "sans", face = "bold"),
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
  scale_x_continuous(limits = c(0,580), expand = c(0,6))



#---------------------------------------------------------------------
# Figure 1d - MP concentrations and E2
#---------------------------------------------------------------------

#Generate the figure
e <-
  ggplot(data = total_mps,
         aes(x = concentration,
             y = therapy_last_E2_value)) +
  ggtitle("e") +
  geom_smooth(method = "gam", formula = y ~ x, method.args = list(family = "tw"), col = "grey30", linetype = "dashed", alpha = 0.2) +
  geom_point(size = 1, shape = 16) +
  ylab(expression(bold(Therapy~last~E2~value~(pg/mL)))) + 
  xlab(expression(bold(Total~MP~Concentration~(Particles~mL^-1)))) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=10, family = "sans", face = "bold"),
        axis.title.x = element_text(size=10, family = "sans", face = "bold"),
        axis.text.y = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=10, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 14, family = "sans", face = "bold"),
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
  scale_x_continuous(limits = c(0,580), expand = c(0,6))




#---------------------------------------------------------------------
# Figure 1f - MP concentrations and oocytes the reached maturity
#---------------------------------------------------------------------


#Generate the figure
f <-
  ggplot(data = total_mps,
         aes(x = concentration,
             y = MII_oocytes/obtained_oocytes)) +
  ggtitle("f") +
  geom_smooth(method = "gam", formula = y ~ x, method.args = list(family = "betar"), col = "grey30", linetype = "dashed", alpha = 0.2) +
  geom_point(size = 1, shape = 16) +
  ylab(expression(bold(Proportion~of~metaphase~II~oocytes))) + 
  xlab(expression(bold(Total~MP~Concentration~(Particles~mL^-1)))) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=10, family = "sans", face = "bold"),
        axis.title.x = element_text(size=10, family = "sans", face = "bold"),
        axis.text.y = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=10, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 14, family = "sans", face = "bold"),
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
  scale_x_continuous(limits = c(0,580), expand = c(0,6))



#---------------------------------------------------------------------
# Figure 1h - MP concentrations and normally fertilised oocytes
#---------------------------------------------------------------------

#Generate the figure
g <-
  ggplot(data = total_mps,
         aes(x = concentration,
             y = X2PNs/MII_oocytes)) +
  ggtitle("g") +
  geom_smooth(method = "gam", formula = y ~ x, method.args = list(family = "betar"), col = "grey30", linetype = "dashed", alpha = 0.2) +
  geom_point(size = 1, shape = 16) +
  ylab(expression(bold(Proportion~of~normally~fertilised~oocytes~'(2PNs)'))) + 
  xlab(expression(bold(Total~MP~Concentration~(Particles~mL^-1)))) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=10, family = "sans", face = "bold"),
        axis.title.x = element_text(size=10, family = "sans", face = "bold"),
        axis.text.y = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=10, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 14, family = "sans", face = "bold"),
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
  scale_x_continuous(limits = c(0,580), expand = c(0,6))


#---------------------------------------------------------------------
# Compile and save
#---------------------------------------------------------------------

FIG <- ggarrange(a,b,c,
                 d,e,f,
                 g,
                 ncol=3, nrow=3)


#Save the figures
ggsave(FIG,
       width = 6.86*2, height = 6*2, units = "in",
       dpi = 600,
       bg = "transparent",
       file="figures/figure_S1.png")
