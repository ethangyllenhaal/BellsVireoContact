# By: Ethan Gyllenhaal
# Updated: 1Feb2026
#
# Script for plotting adn exploring simulation output

# Load packages, I don't think I used beeswarm in any final figures, but have some exploratory ones here
library(tidyverse)
library(gridExtra)
library(ggbeeswarm)

# set working directory and path to output
setwd("/path/to/workdir/")
path <- "/path/to/workdir/simulation_output/"

# initialize dataframe
# don't remember why I went this in depth with type safety, but here we are
# note frequency did not vary in final version
df <- data.frame(step=numeric(), coeff=numeric(), freq=numeric(), asym=numeric(), mig=numeric(), rep=numeric(), slope=numeric(), p=numeric(), rsq=numeric(), intercept=numeric())
df[1,] <- NA
# but for loop for all output files
for(file in list.files(path, pattern=".txt")){
  # read input, make columns for the ratio and transformed ratios
  input <- read.csv(paste(path,file, sep=""), header=T, sep="\t") %>%
    mutate(ratio = pi1/pi2) %>%
    filter(ratio != Inf & ratio != 0) %>%
    mutate(corr_ratio = ratio/mean(ratio)) %>%
    mutate(log_ratio = log(ratio)) %>%
    mutate(standard_log_ratio = ((log_ratio-mean(log_ratio, na.rm = TRUE))/sd(log_ratio, na.rm = TRUE)))
  # if insufficient sites, pass
  if(nrow(input)<10){next}
  # simple linear model for the ratio vs windowed fst
  model <- lm(log_ratio~fst, data=input)
  # split out vector of base name for populating parameter fields
  base <- strsplit(strsplit(file, split='.txt')[[1]], split="_")
  # pulling out model slope, p value, r-squared, and intercept
  slope <- model$coefficients[[2]]
  p <- summary(model)$coefficients[,4][[2]]
  rsq <- summary(model)$adj.r.squared
  intercept <- model$coefficients[[1]]
  # vector of parameters then model output
  vec <- t(data.frame(c(base[[1]][2], base[[1]][3],base[[1]][4], base[[1]][5],base[[1]][6],base[[1]][7], slope, p, rsq, intercept)))
  # add column names then combine
  colnames(vec) <- c("step", "coeff", "freq", "asym", "mig", "rep", "slope", "p", "rsq", "intercept")
  df <- rbind(vec, df)
}
# kill row names then remove NAs
rownames(df) <- NULL
final_df <- drop_na(df)

# awful looking but informative plot splitting simulations out by a range of factors (left to right per fill/line color combo is mig=0,0.5, and 1)
ggplot(filter(final_df,step==1), aes(y=as.numeric(slope), x=as.factor(coeff), fill=as.factor(asym), color=as.factor(freq), shape=as.factor(mig))) +
  geom_boxplot()
# less awful but still confusing plot of focal values over time
# note that asym variable is opposite the final figure (i.e., left is highest asymmetry)
ggplot(filter(final_df, freq<0.001 & mig==0.5), aes(y=as.numeric(slope), x=as.factor(coeff), fill=as.factor(step), shape=as.factor(asym))) +
  geom_boxplot()


# more focal plots ####

# plot of slope
# note that hline denotes parity here, rather than the mean value for a set of simulations
# only includes 10,000 generations after contact
ggplot(filter(final_df, freq<0.001, step==10), aes(y=as.numeric(slope), x=as.factor(coeff), fill=as.factor(1/as.numeric(asym)))) +
  geom_boxplot(position = position_dodge(0.9)) +
  geom_hline(yintercept = 0, linetype="dotted") +
  facet_wrap(vars(as.factor(mig)), nrow=1) + theme_bw() +
  theme(panel.spacing = unit(0.5, "lines")) +
  ylab("FST-Diversity Ratio Slope") + xlab("Selection Coefficient")
mean(as.numeric(filter(final_df, freq<0.001, step==10, asym==0.1, mig==1, coeff==0)$slope), na.rm=T)

# plot of intercept
ggplot(filter(final_df, freq<0.001, step==10), aes(y=as.numeric(intercept), x=as.factor(coeff), fill=as.factor(1/as.numeric(asym)))) +
  geom_boxplot(position = position_dodge(0.9)) +
  geom_hline(yintercept = 0, linetype="dotted") +
  facet_wrap(vars(as.factor(mig)), nrow=1) + theme_bw() +
  theme(panel.spacing = unit(0.5, "lines")) +
  ylab("FST-Diversity Ratio Intercept") + xlab("Selection Coefficient")

# plot of linear model p value, dotted line is log10 of 0.05, solid is bonf-corrected
ggplot(filter(final_df, freq<0.001, step==10), aes(y=log10(as.numeric(p)), x=as.factor(coeff), fill=as.factor(1/as.numeric(asym)))) +
  geom_boxplot() +
  geom_hline(yintercept = log10(0.05), linetype="dotted") +
  geom_hline(yintercept = log10(0.05/nrow(filter(final_df, freq<0.001, step==10)))) +
  facet_wrap(vars(as.factor(mig)), nrow=1) + theme_bw() +
  theme(panel.spacing = unit(0.5, "lines")) +
  ylab("FST-Diversity Ratio log10 p") + xlab("Selection Coefficient")


# Not used in final paper, but it's interesting to explore autosomal vs Z chromosome outputs! ####
# because it wasn't use, it isn't commented well, but mostly matches the first part

# Z ####

df_Z <- data.frame(step=numeric(), coeff=numeric(), freq=numeric(), asym=numeric(), mig=numeric(), rep=numeric(), slope=numeric(), p=numeric(), rsq=numeric(), intercept=numeric())
df_Z[1,] <- NA

for(file in list.files(path)){
  input <- read.csv(paste(path,file, sep=""), header=T, sep="\t") %>%
    mutate(ratio = pi2/pi1) %>%
    filter(ratio != Inf & ratio != 0 & chr==3) %>%
    mutate(corr_ratio = ratio/mean(ratio)) %>%
    mutate(log_ratio = log(ratio)) %>%
    mutate(standard_log_ratio = ((log_ratio-mean(log_ratio, na.rm = TRUE))/sd(log_ratio, na.rm = TRUE)))
  if(nrow(input)<10){next}
  #print(ggplot(data=input, aes(x=fst, y=log_ratio))+geom_point()+geom_smooth(method="lm"))
  model <- lm(log_ratio~fst, data=input)
  #print(summary(model))
  base <- strsplit(strsplit(file, split='.txt')[[1]], split="_")
  slope <- model$coefficients[[2]]
  p <- summary(model)$coefficients[,4][[2]]
  rsq <- summary(model)$adj.r.squared
  intercept <- model$coefficients[[1]]
  vec <- t(data.frame(c(base[[1]][2], base[[1]][3],base[[1]][4], base[[1]][5],base[[1]][6],base[[1]][7], slope, p, rsq, intercept)))
  colnames(vec) <- c("step", "coeff", "freq", "asym", "mig", "rep", "slope", "p", "rsq", "intercept")
  df_Z <- rbind(vec, df_Z)
}
rownames(df_Z) <- NULL
final_df_Z <- drop_na(df_Z)

ggplot(filter(final_df_Z,step==1), aes(y=as.numeric(slope), x=as.factor(coeff), fill=as.factor(asym), color=as.factor(freq), shape=as.factor(mig))) +
  geom_boxplot()
ggplot(filter(final_df_Z, freq<0.001 & mig==0.5), aes(y=as.numeric(slope), x=as.factor(coeff), fill=as.factor(step), shape=as.factor(asym))) +
  geom_boxplot()

# no Z ####
df_noZ <- data.frame(step=numeric(), coeff=numeric(), freq=numeric(), asym=numeric(), mig=numeric(), rep=numeric(), slope=numeric(), p=numeric(), rsq=numeric(), intercept=numeric())
df_noZ[1,] <- NA

for(file in list.files(path)){
  input <- read.csv(paste(path,file, sep=""), header=T, sep="\t") %>%
    mutate(ratio = pi2/pi1) %>%
    filter(ratio != Inf & ratio != 0 & chr!=3) %>%
    mutate(corr_ratio = ratio/mean(ratio)) %>%
    mutate(log_ratio = log(ratio)) %>%
    mutate(standard_log_ratio = ((log_ratio-mean(log_ratio, na.rm = TRUE))/sd(log_ratio, na.rm = TRUE)))
  if(nrow(input)<10){next}
  #print(ggplot(data=input, aes(x=fst, y=log_ratio))+geom_point()+geom_smooth(method="lm"))
  model <- lm(log_ratio~fst, data=input)
  #print(summary(model))
  base <- strsplit(strsplit(file, split='.txt')[[1]], split="_")
  slope <- model$coefficients[[2]]
  p <- summary(model)$coefficients[,4][[2]]
  rsq <- summary(model)$adj.r.squared
  vec <- t(data.frame(c(base[[1]][2], base[[1]][3],base[[1]][4], base[[1]][5],base[[1]][6],base[[1]][7], slope, p, rsq)))
  colnames(vec) <- c("step", "coeff", "freq", "asym", "mig", "rep", "slope", "p", "rsq")
  df_noZ <- rbind(vec, df_noZ)
}
rownames(df_noZ) <- NULL
final_df_noZ <- drop_na(df_noZ)

ggplot(filter(final_df_noZ,step==1), aes(y=as.numeric(slope), x=as.factor(coeff), fill=as.factor(asym), color=as.factor(freq), shape=as.factor(mig))) +
  geom_boxplot()
ggplot(filter(final_df_noZ, freq<0.001 & mig==0.5), aes(y=as.numeric(slope), x=as.factor(coeff), fill=as.factor(step), shape=as.factor(asym))) +
  geom_boxplot()

steps <- ggplot(filter(final_df, freq<0.001 & mig==0.5), aes(y=as.numeric(slope), x=as.factor(coeff), fill=as.factor(asym), shape=as.factor(step))) +
  geom_boxplot() + ylim(-1.5,1.5)
steps + xlab("Selection Coefficient") + ylab("Slope") + ggtitle("FST-Diversity Ratio Slopes")

focus <- ggplot(filter(final_df, freq<0.001 & step == 3 & coeff == 0), aes(y=as.numeric(slope), x=as.factor(mig), fill=as.factor(asym))) +
  geom_boxplot() + ylim(-2,2) + theme(legend.position = "none") 
focus_Z <- ggplot(filter(final_df_Z, freq<0.001 & step == 3 & coeff == 0), aes(y=as.numeric(slope), x=as.factor(mig), fill=as.factor(asym))) +
  geom_boxplot() + ylim(-2,2) + theme(legend.position = "none") 
focus_noZ <- ggplot(filter(final_df_noZ, freq<0.001 & step == 3 & coeff == 0), aes(y=as.numeric(slope), x=as.factor(mig), fill=as.factor(asym))) +
  geom_boxplot() + ylim(-2,2) + theme(legend.position = "none") 

grid.arrange(focus_Z+ggtitle("Z"), focus_noZ+ggtitle("no Z"), nrow=1)

t.test(as.numeric(filter(final_df_noZ, freq<0.001 & mig==0.5 & asym==0.1 & coeff==0)$slope),
       as.numeric(filter(final_df_Z, freq<0.001 & mig==0.5 & asym==0.1 & coeff==0)$slope))

