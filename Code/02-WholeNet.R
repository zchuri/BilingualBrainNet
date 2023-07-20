#!/usr/bin/Rscript

# Clear workspace
rm(list = ls())

# Load libraries
if(!require(emmeans)) install.packages("emmeans"); library(emmeans)
if(!require(ggplot2)) install.packages("ggplot2"); library(ggplot2)
if(!require(ggpubr)) install.packages("ggpubr"); library(ggpubr)
if(!require(gridExtra)) install.packages("gridExtra"); library(gridExtra)
if(!require(modelsummary)) install.packages("modelsummary"); library(modelsummary)
if(!require(psych)) install.packages("psych"); library(psych)
if(!require(rstatix)) install.packages("rstatix"); library(rstatix)

# Read phenotypic and whole-brain data
info <- read.csv("https://raw.githubusercontent.com/zchuri/BilingualBrainNet/main/Data/Phenotypic/phenotypic.csv")

########################################################################
# Sample description by group
########################################################################
# N
table(info$group)
# Age
describeBy(info$Age, info$group)
summary(aov(Age~group,info))
# Sex
table(info$group,info$Sex_Gender)
chisq.test(table(info$group,info$Sex_Gender))
# AoA
describeBy(info$AOA, info$group)
summary(aov(AOA~group,info))
# YoE
describeBy(info$YOE, info$group)
summary(aov(AOA~group,info))
# FD
describeBy(info$avgFD, info$group)
summary(aov(avgFD~group,info))

########################################################################
# Whole-brain group inference
########################################################################

# Group inference
# Set results folder
#pltdir <- file.path(getwd(),"Results","01-WholeNet")
#if(!dir.exists(pltdir)) dir.create(path = pltdir, recursive = T)

# All groups - Modularity
ff <- anova(lm(Q~group, info))
# Plot boxplot
plt_df <- info
plt_df$group <- factor(plt_df$group, levels = c("Mono","Late.B","Early.B","Sim.B"))
g3 <- ggplot(data = plt_df,
             mapping = aes(y=Q, x=group, color=group)) +
  geom_violin() +
  geom_boxplot(width=0.4, color="grey", alpha=0.2) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  labs(title = "Modularity",
       subtitle = paste0("F=",round(ff$`F value`[1],1)," (p=",signif(ff$`Pr(>F)`[1],2),")"),
       y = "Q") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))
show(g3)

# All groups - Efficiency
ff <- anova(lm(E~group, info))
# Plot boxplot
g4 <- ggplot(data = plt_df,
             mapping = aes(y=E, x=group, color=group)) +
  geom_violin() +
  geom_boxplot(width=0.4, color="grey", alpha=0.2) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  labs(title = paste0("Efficiency (F=",round(ff$`F value`[1],1),", p=",signif(ff$`Pr(>F)`[1],2),")"),
       #subtitle = paste0("F=",round(ff$`F value`[1],1)," (p=",signif(ff$`Pr(>F)`[1],2),")"),
       y = "E") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))
g4 <- g4 + geom_signif(comparisons = list(
  #c("Mono", "Late-B"),
  c("Mono", "Early.B"),
  c("Mono", "Sim.B")),
  test = "t.test",
  step_increase = 0.075,
  map_signif_level = c("*"=0.05),
  color = "black")
show(g4)

# Save plots
#ggsave("E_bili_groups.png", g4, "png", pltdir,1,6,4)
#ggsave("E_bili_groups.svg", g4, "svg", pltdir,1,4,3)

# Post-hoc ala tidyverse
# Prepare the data and inspect a random sample of the data
t_test(info, E ~ group, p.adjust.method = "holm")

########################################################################
# Age of acquisition (AoA)
########################################################################

# Bilinguals excluding simultaneous ones
idx <- which(info$AOA>0)
cor.test(info$AOA[idx], info$Q[idx])
show(ct <- cor.test(info$AOA[idx], info$E[idx]))

# Plot
# Fit model
plt_df <- info[idx,]
fit <- lm(E~AOA, plt_df)
newx <- seq(min(plt_df$AOA), max(plt_df$AOA), length.out=100)
preds <- predict(fit, newdata = data.frame(AOA=newx), interval = "confidence")
#create plot of x vs. y, but don't display individual points (type='n') 
ttl <- paste0("Age of acquisition effect\n(rho = ",round(ct$estimate,2),
              "; p = ",signif(ct$p.value,2),")")
plot(E~AOA, data = plt_df, type = "n", axes = F,
     main= ttl, xlab = "AoA", ylab = "E")
points(E~AOA, data = plt_df, pch = 23, cex = 1,
       bg = "yellow2", col = "gray40", lwd = 0.5)
#add fitted regression line
minX <- which.min(plt_df$AOA)
maxX <- which.max(plt_df$AOA)
lines(plt_df$AOA[c(minX,maxX)],predict(fit)[c(minX,maxX)],
      lwd = 1.5)
#add dashed lines for confidence bands
lines(newx, preds[ ,3], lwd = 1.5, lty = 3)
lines(newx, preds[ ,2], lwd = 1.5, lty = 3)
# add axes
axis(1); axis(2, las = 1)

########################################################################
# Years of experience (YoE)
########################################################################

# Only bilinguals
idx <- which(info$AOA>=0)
cor.test(info$YOE[idx], info$Q[idx])
cor.test(info$YOE[idx], info$E[idx])

########################################################################
# L2 proficiency
########################################################################

# All sample
cor.test(rank(info$L2proficiency), info$AOA)
cor.test(rank(info$L2proficiency), info$YOE)
cor.test(rank(info$L2proficiency), info$Q)
cor.test(rank(info$L2proficiency), info$E)

########################################################################
# Music experience affects the AoA relationship
########################################################################

# Set data.frame
plt_df <- info
plt_df$group <- factor(plt_df$group, levels = c("Mono","Late.B","Early.B","Sim.B"))

# Contingency table
table(plt_df$music, plt_df$group)

# All groups - Modularity
fit <- lm(Q~group+music, plt_df)
summary(fit)$coefficients
emmeans(fit, specs = pairwise ~ group, adjust = "holm")

# All groups - Efficiency
fit <- lm(E~group+music, plt_df)
summary(fit)$coefficients
(emfit <- emmeans(fit, specs = pairwise ~ group, adjust = "holm"))
# Plot
plot(x = emfit, xlab = "Efficiency")

# Age of Acquisition
phen <- plt_df[which(plt_df$AOA>0),]
fit <- lm(E~AOA+music, phen)
summary(fit)$coefficients
# Plot
modelplot(fit, coef_omit = 'Interc')
