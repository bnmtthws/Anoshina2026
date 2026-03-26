###### R code for Mosquito texture paper, 2026 ######
###
###
###
###### load required libraries
library(ggplot2)
library(readr)
library(dplyr)
library(tidyverse)
library(grid)
library(RColorBrewer)
library(plotrix)
library(cowplot)
library(gg.gap)
library(ggbeeswarm)
library(gghalves)
library(ggridges)
library(multcompView)
library(ggridges)
library(dplyr)
library(ggplot2)
library(Hmisc)
library(drc)


# Figure 1

######
## Figure 1
######
data.fig1 <- read.csv(file="TextureFig1.csv", head=TRUE,sep=",", skipNul=TRUE, stringsAsFactors = FALSE)


## panel B
fig1.panelB <- ggplot(data=data.fig1, aes(y=preference,x=(microns),group=grit2)) + 
  geom_boxplot(outlier.shape = NA, fill="#00BFC4") + 
  geom_jitter(position=position_jitter(1),alpha=.65,stroke=NA,colour="#00BFC4") +
  coord_cartesian(ylim=c(-1,1),xlim=c(0,200)) +
  geom_hline(yintercept=0,linetype = "dashed")+
  xlab("surface texture (µm)")+
  ylab("preference index (PI)")+
  theme(axis.text.x = element_text(angle = -60, vjust = 0.5, hjust=0.75, size=10, face="bold", colour="black"), 
        axis.title = element_text(size=13,colour="black", face="bold"),
        axis.text.y = element_text(face="bold", colour="black"))+
  theme_classic()

## panel C


### fit ll4 with upper bound fixed at 1 and lower bound at 0
fit <- drm(preference ~ microns, data = data.fig1, fct = LL.4(fixed=(c(NA,0,1,NA))))
summary(fit)
# Extract the EC50 value from the fitted model
ed50 <- summary(fit)$coefficients["e:(Intercept)", "Estimate"]

# Create a sequence of new x values to predict
new_data <- data.frame(microns = seq(min(data.fig1$microns), max(data.fig1$microns), length.out = 100))
# Predict the values using the fitted model
new_data$pred <- predict(fit, newdata = new_data)

# Calculate confidence intervals
pred_ci <- predict(fit, newdata = new_data, interval = "confidence")
new_data$lower <- pred_ci[, "Lower"]
new_data$upper <- pred_ci[, "Upper"]

fig1.panelC <- ggplot() +
  stat_summary(fun.data = "mean_cl_normal",data=data.fig1, aes(x = microns, y = preference)) +
  geom_line(data = new_data, aes(x = microns, y = pred), color = "green") +
  geom_ribbon(data = new_data, aes(x = microns, ymin = lower, ymax = upper), fill = "green",alpha = 0.2) +
  labs(
    x = "Texture (µM)",
    y = "Egg-laying preference index") + ylim(-1,1) + xlim(0,200) + 
  theme_classic()  + geom_hline(yintercept=0,linetype='dashed')

## panel D
fig1.panelD <- 
  ggplot(data=data.fig1, aes(y=(egg1+egg2)/10,x=microns)) + 
  stat_summary(fun.data = "mean_cl_normal") +
  geom_boxplot(outlier.shape = NA,aes(group=as.factor(microns)))+
  geom_jitter(position=position_jitter(1),alpha=.65,stroke=NA) +
  theme(axis.text.x = element_text(angle = -60, vjust = 0.5, hjust=0.75, size=7))+
  xlab("surface texture (µm)")+
  ylab("egg count per female")+
  theme_classic()+ theme(legend.position = "none") +
  geom_smooth(method="lm",linewidth = .5,color = "grey",fill="grey") + 
  coord_cartesian(ylim=c(0,130),xlim=c(0,200))

## panel D stats
fig1.panelD.lm <- lm(data=filter(data.fig1),(egg1+egg2)/10 ~ microns)
summary(fig1.panelD.lm)  # grit2         9  221203   24578   1.872 0.0632 .
#TukeyHSD(fig1.panelD.aov)


#########
## Figure 2
#########
data.fig2 = read.csv(file="TextureFig2.csv",head=TRUE,sep=",")



# panel 2B
fig2.faceted <- ggplot(data=data.fig2,aes(y=as.numeric(PI_without_mid),x=microns_right,group=microns_left, 
                                          fill=-1*(microns_left-microns_right))) + 
  geom_boxplot(aes(group=microns_right),outlier.shape = NA) +
  geom_jitter(position=position_jitter(1),alpha=.65,stroke=NA) +
  # geom_smooth(method="lm") +
  facet_grid(. ~ microns_left) +
  theme_classic() +
  theme(legend.position ="right") +
  coord_cartesian(ylim=c(-1, 1)) +
  geom_hline(yintercept=0,linetype="dashed") +
  xlab('Particle size (μm)') +
  ylab('Preference index') +
  scale_x_continuous(breaks=c(0, 35, 46, 68, 92, 190), labels=c('0', '35', '46', '68', '92','190'), guide = guide_axis(angle = 45)) +
  scale_fill_distiller(palette = "YlOrBr") +
  labs(fill='Texture series')


# panel 2C
### fit ll4 with upper bound fixed at 1 and lower bound at 0
fit.fig2 <- drm(as.numeric(PI_without_mid) ~ microns_difference, data = data.fig2, fct = LL.4(fixed=(c(NA,0,1,NA))))
summary(fit.fig2)
# Extract the EC50 value from the fitted model
ed50.fig2 <- summary(fit)$coefficients["e:(Intercept)", "Estimate"]

# Create a sequence of new x values to predict
new_data2 <- data.frame(microns = seq(min(data.fig2$microns_difference), max(data.fig2$microns_difference), length.out = 100))
# Predict the values using the fitted model
new_data2$pred <- predict(fit.fig2, newdata = new_data2)

# Calculate confidence intervals
pred_ci2 <- predict(fit.fig2, newdata = new_data2, interval = "confidence")
new_data2$lower <- pred_ci2[, "Lower"]
new_data2$upper <- pred_ci2[, "Upper"]

fig2.panelC <- ggplot() +
  stat_summary(fun.data = "mean_cl_normal",data=data.fig2, aes(x = microns_difference, y = as.numeric(PI_without_mid))) +
  geom_line(data = new_data2, aes(x = microns, y = pred), color = "magenta") +
  geom_ribbon(data = new_data2, aes(x = microns, ymin = lower, ymax = upper), fill = "magenta",alpha = 0.2) +
  labs(
    x = "Texture (µM)",
    y = "Egg-laying preference index") + ylim(-1,1) + xlim(0,200) + 
  theme_classic()  + geom_hline(yintercept=0,linetype='dashed')


fig2.relative <- ggplot(data=data.fig2,aes(y=as.numeric(PI_without_mid),x=microns_difference,group)) + 
  geom_point(color= "black", shape=16) +
  geom_smooth() +
  theme_classic() +
  coord_cartesian(ylim=c(-1, 1)) +
  geom_hline(yintercept=0,linetype="dashed") +
  xlab('Difference in texture (μm)') +
  ylab('Preference index') +
  scale_fill_brewer(palette = "BuPu")
 




########
## Figure 3a
########


data.fig3a = read.csv(file="TextureFig3a.csv")

fig3.panelA <- ggplot(data=data.fig3a,aes(y=as.numeric(PI_without_mid),
                                           x=microns_right,group=microns_right, 
                                           fill=as.factor(microns_right))) + 
  geom_boxplot(color="black") +
  geom_jitter(color= "black", shape=16, position=position_jitter(0.2)) +
  facet_grid(. ~ microns_left) +
  theme_classic() +
  coord_cartesian(ylim=c(-1, 1)) +
  geom_hline(yintercept=0) +
  xlab('Particle size (μm)') +
  ylab('Preference index') +
  scale_x_continuous(breaks=c(0, 35, 68), labels=c('0', '35',  '68' ), guide = guide_axis(angle = 45)) +
  scale_fill_brewer(palette = "BuPu") +
  labs(fill='Texture series')


data.fig3a <- data.fig3a %>% 
  mutate(microns_dif = microns_left-microns_right) %>%
  slice(-1)

### fit ll4 with upper bound fixed at 1 and lower bound at 0
fit.fig3a <- drm(as.numeric(PI_without_mid) ~ microns_dif, data = data.fig3a, fct = LL.4(fixed=(c(NA,0,1,NA))))
summary(fit.fig3a)
# Extract the EC50 value from the fitted model
ed50.fig3a <- summary(fit.fig3a)$coefficients["e:(Intercept)", "Estimate"]

# Create a sequence of new x values to predict
new_data3a <- data.frame(microns = seq(min(data.fig3a$microns_dif), 
                                      max(data.fig3a$microns_dif), 
                                      length.out = 100))
# Predict the values using the fitted model
new_data3a$pred <- predict(fit.fig3a, newdata = new_data3a)

# Calculate confidence intervals
pred_ci3a <- predict(fit.fig3a, newdata = new_data3a, interval = "confidence")
new_data3a$lower <- pred_ci3a[, "Lower"]
new_data3a$upper <- pred_ci3a[, "Upper"]

fig3.panelA_models <- ggplot(data=data.fig3a,
                            aes(y=as.numeric(PI_without_mid),
                                x=microns_left-microns_right)) + 
  geom_line( inherit.aes = FALSE,data = new_data3a, aes(x = microns, y = pred), color = "green",alpha=0.5) +
  geom_ribbon( inherit.aes = FALSE,data = new_data3a, aes(x = microns,ymin = lower, ymax = upper), fill = "green",alpha = 0.3) +
  geom_line( inherit.aes = FALSE,data = new_data2, aes(x = microns, y = pred), color = "magenta",alpha=0.3) +
  geom_ribbon( inherit.aes = FALSE,data = new_data2, aes(x = microns,ymin = lower, ymax = upper), fill = "magenta",alpha = 0.1) +
  stat_summary(fun.data = "mean_cl_normal",
               data=data.fig3a, aes(x = microns_left-microns_right, 
                                    y = as.numeric(PI_without_mid),
                                    fill=(-1*(microns_left-microns_right))),
               geom = "pointrange",
               shape = 21,
               size = 0.8,
               color = "black") +
    theme_classic() +
  coord_cartesian(ylim=c(-1, 1)) +
  geom_hline(yintercept=0,linetype="dashed") +
  xlab('Difference in texture (μm)') +
  ylab('Preference index')+
  scale_fill_distiller(palette = "YlOrBr")


### Figure 3C

data.fig3c = read.csv(file="TextureFig3.csv",head=TRUE,sep=",") %>% 
  mutate(salt_tex = replace(salt_tex, salt_tex == 80, 190)) %>%
  mutate(salt_tex = replace(salt_tex, salt_tex == 220, 68)) %>%
  mutate(salt_tex = replace(salt_tex, salt_tex == 400, 35))

fig3.panelC <- ggplot(data=data.fig3c,aes(x=as.factor(salt_tex),
                                         y=PI,group=as.factor(salt_tex))) + 
  ylim(-1,1) +
  geom_boxplot(width=0.5,position=position_dodge2(2),outlier.shape = NA, 
               aes(group=as.factor(salt_tex))) + 
  geom_point(position=position_jitter(0.15),aes(group=as.factor(salt_tex),
                                                colour=as.factor(salt))) + 
  geom_hline(yintercept=0,linetype = 'dashed') +
  facet_grid(.~salt) + theme_classic() 
  

# panel 3D

df_diff_ci <- data.fig3c %>%
  filter(salt %in% c(100, 200)) %>%
  group_by(salt_tex, salt) %>%
  summarise(
    mean_PI = mean(PI),
    se = sd(PI) / sqrt(n()),
    n = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = salt,
    values_from = c(mean_PI, se, n),
    names_sep = "_"
  ) %>%
  mutate(
    diff_PI = mean_PI_200 - mean_PI_100,
    se_diff = sqrt(se_200^2 + se_100^2),
    ci_low = diff_PI - 1.96 * se_diff,
    ci_high = diff_PI + 1.96 * se_diff
  )

fig3.panelD.diff <- ggplot(df_diff_ci,
                           aes(x = factor(salt_tex), y = diff_PI)) +
  
  geom_pointrange(
    aes(ymin = ci_low, ymax = ci_high),
    size = 0.6
  ) +
  geom_hline(yintercept=0,linetype = 'dashed') +
  
  labs(
    x = "Texture (µm)",
    y = "Δ Preference Index (200 − 100)"
  ) +
  coord_cartesian(ylim=c(-1, 1)) +
  
  theme_classic() + aes(x = factor(salt_tex, levels = c(0, 35, 68, 190)))

model <- lm(
  PI ~ factor(salt_tex) * factor(salt),
  data = data.fig3c
)

anova(model)
summary(model)



########
## Figure 4 clumping
########


# change texture from conventional to numerical
replacement_values <- c("120" = 115, "150" = 92, "80" = 190, "P180" = 82, "P220" = 68, 
                        "P320" = 46.2, "P360" = 40.5, "P400" = 35, "P500" = 30.2, "smooth" = 0)

data.fig4 <- read_tsv(file="TextureFig4.tsv") %>%   
  mutate(texture = recode(texture, !!!replacement_values),
         texture = as.numeric(texture)) %>% 
  filter(genotype=="ORL")

sum.avg1 <- data.fig4 %>%
  group_by(texture) %>%
  dplyr::summarize(
    mean = mean(avg1),
    ci_lower = Hmisc::smean.cl.normal(avg1)[2],
    ci_upper = Hmisc::smean.cl.normal(avg1)[3]
  )

sum.avg6 <- data.fig4 %>%
  group_by(texture) %>%
  dplyr::summarize(
    mean = mean(avg6),
    ci_lower = Hmisc::smean.cl.normal(avg6)[2],
    ci_upper = Hmisc::smean.cl.normal(avg6)[3]
  )


# panel 4B 
fig4.avg1 <- ggplot(data=data.fig4, aes(x=avg1,y=as.factor(texture),group=as.factor(texture),fill=as.factor(texture))) + 
  stat_density_ridges(quantile_lines = TRUE, alpha=0.4, quantiles = 0.5) +
  geom_pointrange(data = sum.avg1, aes(x = mean, xmin = ci_lower, xmax = ci_upper, y = as.factor(texture))) + 
  theme_minimal() +
  ylab("surface texture (µm)") +
  xlab("nearest neighbor (µm)") + coord_cartesian(xlim = c(0,100)) + 
  theme(legend.position = "none")

# panel 4C
fig3.avg6 <- ggplot(data=data.fig4, aes(x=avg6,y=as.factor(texture),group=as.factor(texture),fill=as.factor(texture))) + 
  stat_density_ridges(quantile_lines = TRUE, alpha=0.4, quantiles = 0.5) +
  geom_pointrange(data = sum.avg6, aes(x = mean, xmin = ci_lower, xmax = ci_upper, y = as.factor(texture))) + 
  theme_minimal() +
  ylab("surface texture (µm)") +
  xlab("average of 6 nearest neighbors (µm)") + coord_cartesian(xlim = c(0,100)) +
  theme(legend.position = "none")



