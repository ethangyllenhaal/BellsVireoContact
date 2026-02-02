# By: Ethan Gyllenhaal
# Updated: 1Feb2026
#
# Script for plotting and exploring the landscape of genomic divergence
# Primarily focuses on diversity ratio, but also other pixy-related outputs 
## (e.g., FST vs Dxy, genome-wide weighted mean pi)

# load libraries
library(gridExtra)
library(viridis)
library(tidyverse)

# set working directory
setwd('C:/Documents/Projects/BEVI/diversity_fst_plots')

# primary color scheme
color_values = c("#CC9B2B", "#03688C")

# load in pixy data frames
pi <- read.csv("pixy_output_pi.txt", sep = "\t")
dxy <- read.csv("pixy_output_dxy.txt", sep = "\t")
fst <- read.csv("pixy_output_fst.txt", sep = "\t")


# function for making diversity ratio input
# needs pop 1 and 2 to match order of input
make_input <- function(name1, name2){
  # combine pi, dxy, and FST dataframes, make ratio column (P1/P2) and other related variables
  input <- pi %>%
    filter(pop==name1) %>% rename(P1.Pi = avg_pi) %>%
    inner_join(y=filter(pi, pop==name2) %>% rename(P2.Pi = avg_pi),
               by = c("chromosome" = "chromosome", "window_pos_1" = "window_pos_1")) %>%
    inner_join(y=filter(dxy, pop1==name1 & pop2==name2) %>% rename(Dxy=avg_dxy),
               by = c("chromosome" = "chromosome", "window_pos_1" = "window_pos_1")) %>%
    inner_join(y=filter(fst, pop1==name1 & pop2==name2) %>% rename(Mean.FST=avg_wc_fst),
               by = c("chromosome" = "chromosome", "window_pos_1" = "window_pos_1")) %>%
    mutate(ratio = P1.Pi/P2.Pi) %>%
    filter_all(all_vars(!is.infinite(.))) %>%
    mutate(ratio_div_mean = ratio/mean(ratio, na.rm = TRUE)) %>%
    mutate(standard_ratio = ((ratio-mean(ratio, na.rm = TRUE))/sd(ratio, na.rm = TRUE))) %>%
    drop_na() %>%
    mutate(chr_type = ifelse(chromosome == "CM022210.2_RagTag", "Z", "Autosome"))
  return(input)
}

# function for making diversity ratio input with the reciprocal of diversity ratio
## only used to confirm pattern isn't reliant on order of populations
make_input_flipDiv <- function(name1, name2){
  # combine pi, dxy, and FST dataframes, make ratio column (P2/P1) and other related variables
  input <- pi %>%
    filter(pop==name1) %>% rename(P1.Pi = avg_pi) %>%
    inner_join(y=filter(pi, pop==name2) %>% rename(P2.Pi = avg_pi),
               by = c("chromosome" = "chromosome", "window_pos_1" = "window_pos_1")) %>%
    inner_join(y=filter(dxy, pop1==name1 & pop2==name2) %>% rename(Dxy=avg_dxy),
               by = c("chromosome" = "chromosome", "window_pos_1" = "window_pos_1")) %>%
    inner_join(y=filter(fst, pop1==name1 & pop2==name2) %>% rename(Mean.FST=avg_wc_fst),
               by = c("chromosome" = "chromosome", "window_pos_1" = "window_pos_1")) %>%
    mutate(ratio = P2.Pi/P1.Pi) %>%
    filter_all(all_vars(!is.infinite(.))) %>%
    mutate(ratio_div_mean = ratio/mean(ratio, na.rm = TRUE)) %>%
    mutate(standard_ratio = ((ratio-mean(ratio, na.rm = TRUE))/sd(ratio, na.rm = TRUE))) %>%
    drop_na() %>%
    mutate(chr_type = ifelse(chromosome == "CM022210.2_RagTag", "Z", "Autosome"))
  return(input)
}


# function for making diversity ratio plots
make_ratioplot <- function(df){
  # plots the standardized ratio vs windowed FST, separated by chromosome type
  plot <- ggplot(df, aes(Mean.FST, standard_ratio, colour=chr_type)) +
    theme_bw(base_size = 10) + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    geom_point(size=1, show.legend = FALSE) + ylim(c(-6,28)) + xlim(c(-0.1,0.65)) +
    geom_smooth(method='lm', show.legend = FALSE) + 
    geom_hline(yintercept = mean(df$standard_ratio, na.rm = TRUE), linetype="dashed", color = "black") +
    scale_colour_manual(values = color_values)
  return(plot)
}


# function for making diversity ratio plots with the log ratio
make_log_ratioplot <- function(df){
  # plots the log of the raw ratio vs windowed FST, separated by chromosome type
  plot <- ggplot(df, aes(Mean.FST, log(ratio), colour=chr_type)) +
    theme_bw(base_size = 10) + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    geom_point(size=1, show.legend = FALSE) + ylim(c(-2.5,2.5)) + xlim(c(-0.1,0.65)) +
    geom_smooth(method='lm', show.legend = FALSE) + 
    geom_hline(yintercept = mean(log(df$ratio), na.rm = TRUE), linetype="dashed", color = "black") +
    scale_colour_manual(values = color_values)
  return(plot)
}

# function for making Dxy plots
make_dxyplot <- function(df){
  # plots FST vs Dxy, separated by chromosome
  plot <- ggplot(df, aes(Mean.FST, Dxy, colour=chr_type)) +
    theme_bw(base_size = 10) + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    geom_point(size=1, show.legend = FALSE) +
    geom_smooth(method='lm', show.legend = FALSE) + 
    scale_colour_manual(values = color_values)
  return(plot)
}



## Main plots ####
# Only the first instance of a given line is described

# WTX vs AZ

# make input for two populations
wtx_az <- make_input("WTX", "AZ")
# make plot of standardized ratio
wtx_az_plot <- make_ratioplot(wtx_az) + ggtitle("West Texas : Arizona") + 
  theme(plot.title = element_text(size = 12))
# make plot of log ratio for focal comparison
wtx_az_log_plot <- make_log_ratioplot(wtx_az) + ggtitle("West Texas : Arizona") + 
  theme(plot.title = element_text(size = 12))
# print plot
wtx_az_plot
# output linear model for data with FST above zero
summary(lm(log(ratio)~Mean.FST, data=filter(wtx_az, Mean.FST>=0)))
# output the mean log ratio value
mean(log(wtx_az$ratio), na.rm=T)
# output the mean ratio value
mean(wtx_az$ratio, na.rm=T)
# make FST vs Dxy plot
make_dxyplot(wtx_az)

# Sev vs Butte

butte_sev <- make_input("Butte", "Sev")
butte_sev_plot <- make_ratioplot(butte_sev) + ggtitle("Sevilleta : Elephant Butte") + 
  theme(plot.title = element_text(size = 12))
butte_sev_plot
make_dxyplot(butte_sev)

# Butte vs AZ

butte_az <- make_input("Butte", "AZ")
butte_az_plot <- make_ratioplot(butte_az) + ggtitle("Elephant Butte : Arizona") + 
  theme(plot.title = element_text(size = 12))
butte_az_plot
make_dxyplot(butte_az)


# WTX vs Sev

wtx_sev <- make_input("WTX", "Sev")
wtx_sev_plot <- make_ratioplot(wtx_sev) + ggtitle("West Texas : Sevilleta") + 
  theme(plot.title = element_text(size = 12))
wtx_sev_plot
make_dxyplot(wtx_sev)


# WTX vs Butte

wtx_butte <- make_input("WTX", "Butte")
wtx_butte_plot <- make_ratioplot(wtx_butte) + ggtitle("West Texas : Elephant Butte") + 
  theme(plot.title = element_text(size = 12))
wtx_butte_plot
make_dxyplot(wtx_butte)


# Sev vs AZ

sev_az <- make_input("Sev", "AZ")
sev_az_plot <- make_ratioplot(sev_az) + ggtitle("Sevilleta : Arizona") + 
  theme(plot.title = element_text(size = 12))
sev_az_plot
make_dxyplot(sev_az)

# arrange plots
grid.arrange(wtx_az_plot, wtx_sev_plot, wtx_butte_plot, butte_sev_plot, sev_az_plot, butte_az_plot, nrow=2)
grid.arrange(wtx_az_plot, butte_sev_plot, wtx_sev_plot, butte_az_plot, 
             left="Standardized Diversity Ratio", bottom="FST (50kb windows)", nrow=2)


# log plots

butte_sev_log_plot <- make_log_ratioplot(butte_sev) + ggtitle("Sevilleta : Elephant Butte") + 
  theme(plot.title = element_text(size = 12))
wtx_sev_log_plot <- make_log_ratioplot(wtx_sev) + ggtitle("West Texas : Sevilleta") + 
  theme(plot.title = element_text(size = 12))
butte_az_log_plot <- make_log_ratioplot(butte_az) + ggtitle("Elephant Butte : Arizona") + 
  theme(plot.title = element_text(size = 12))

grid.arrange(wtx_az_log_plot, butte_sev_log_plot, wtx_sev_log_plot, butte_az_log_plot, 
             left="Log Diversity Ratio", bottom="FST (50kb windows)", nrow=2)

wtx_az_log_plot + ylab("Log(WTX:AZ) Diversity Ratio") + xlab("FST (50kb windows)") + 
  ylim(-1,2.2) + xlim(-0.05, 0.62)


# Stats of diversity ratio ####
## Assumes all run above for input dataframes

mean(c(mean(filter(wtx_sev, chr_type=="Autosome")$P1.Pi, na.rm = TRUE)/
mean(filter(wtx_sev, chr_type=="Z")$P1.Pi, na.rm = TRUE),
mean(filter(wtx_sev, chr_type=="Autosome")$P2.Pi, na.rm = TRUE)/
mean(filter(wtx_sev, chr_type=="Z")$P2.Pi, na.rm = TRUE)))
mean(filter(wtx_sev, chr_type=="Autosome")$Dxy, na.rm = TRUE)/
mean(filter(wtx_sev, chr_type=="Z")$Dxy, na.rm = TRUE)

# T test for differences in pi
t.test(filter(wtx_sev, chr_type=="Z")$P1.Pi, filter(wtx_sev, chr_type=="Z")$P2.Pi)
t.test(filter(wtx_sev, chr_type=="Autosome")$P1.Pi, filter(wtx_sev, chr_type=="Autosome")$P2.Pi)

t.test(filter(butte_az, chr_type=="Z")$P1.Pi, filter(wtx_sev, chr_type=="Z")$P2.Pi)
t.test(filter(butte_az, chr_type=="Autosome")$P1.Pi, filter(wtx_sev, chr_type=="Autosome")$P2.Pi)

# linear models of FST's correlates, first just for the ratio then breaking it down into components
summary(lm(Mean.FST ~ standard_ratio, data=filter(wtx_az, chr_type=="Autosome")))
summary(lm(Mean.FST ~ standard_ratio, data=filter(butte_sev, chr_type=="Autosome")))
summary(lm(Mean.FST ~ standard_ratio, data=filter(wtx_sev, chr_type=="Autosome")))
summary(lm(Mean.FST ~ standard_ratio, data=filter(butte_az, chr_type=="Autosome")))

summary(lm(Mean.FST ~ standard_ratio, data=filter(wtx_az, chr_type=="Z")))
summary(lm(Mean.FST ~ standard_ratio, data=filter(butte_sev, chr_type=="Z")))
summary(lm(Mean.FST ~ standard_ratio, data=filter(wtx_sev, chr_type=="Z")))
summary(lm(Mean.FST ~ standard_ratio, data=filter(butte_az, chr_type=="Z")))


summary(lm(Mean.FST ~ ratio, data=filter(wtx_az, chr_type=="Autosome")))
summary(lm(Mean.FST ~ P1.Pi, data=filter(wtx_az, chr_type=="Autosome")))
summary(lm(Mean.FST ~ P2.Pi, data=filter(wtx_az, chr_type=="Autosome")))
summary(lm(Mean.FST ~ Dxy, data=filter(wtx_az, chr_type=="Autosome")))


lm(Mean.FST ~ P1.Pi + P2.Pi + Dxy , data=filter(wtx_az, chr_type=="Autosome"))
lm(Mean.FST ~ P1.Pi + P2.Pi + Dxy , data=filter(wtx_az, chr_type=="Z"))

lm(Mean.FST ~ Butte.Pi + AZ.Pi + Dxy , data=filter(modified_butteaz, chr_type=="Autosome"))
lm(Mean.FST ~ Butte.Pi + AZ.Pi + Dxy , data=filter(modified_butteaz, chr_type=="Z"))

# Calculating mean autosome and Z pi per population ####

mean_pi_z <- function(name1){
  input <- pi %>%
    filter(pop==name1) %>%
    filter(chromosome == "CM022210.2_RagTag")
  sum_diffs <- sum(input$count_diffs)
  sum_comps <- sum(input$count_comparisons)
  return(sum_diffs/sum_comps)
}

mean_pi_auto <- function(name1){
  input <- pi %>%
    filter(pop==name1) %>%
    filter(chromosome != "CM022210.2_RagTag")
  sum_diffs <- sum(input$count_diffs)
  sum_comps <- sum(input$count_comparisons)
  return(sum_diffs/sum_comps)
}

mean_pi_auto("AZ")
mean_pi_auto("Butte")
mean_pi_auto("Sev")
mean_pi_auto("WTX")

mean_pi_z("AZ")
mean_pi_z("Butte")
mean_pi_z("Sev")
mean_pi_z("WTX")

# calculate mean pi for putative inversions

sum(filter(pi, pop=="WTX" & chromosome == "CM022185.2_RagTag" & window_pos_1 >= 16.55e6 & window_pos_1 <=17.5e6)$count_diffs)/
  sum(filter(pi, pop=="WTX" & chromosome == "CM022185.2_RagTag"& window_pos_1 >= 16.55e6 & window_pos_1 <=17.5e6)$count_comparisons)
sum(filter(pi, pop=="AZ" & chromosome == "CM022185.2_RagTag"& window_pos_1 >= 16.55e6 & window_pos_1 <=17.5e6)$count_diffs)/
  sum(filter(pi, pop=="AZ" & chromosome == "CM022185.2_RagTag"& window_pos_1 >= 16.55e6 & window_pos_1 <=17.5e6)$count_comparisons)

sum(filter(pi, pop=="WTX" & chromosome == "CM022210.2_RagTag" & window_pos_1 >= 13.14e6 & window_pos_1 <=13.20e6)$count_diffs)/
  sum(filter(pi, pop=="WTX" & chromosome == "CM022210.2_RagTag"& window_pos_1 >= 13.14e6 & window_pos_1 <=13.20e6)$count_comparisons)
sum(filter(pi, pop=="AZ" & chromosome == "CM022210.2_RagTag"& window_pos_1 >= 13.14e6 & window_pos_1 <=13.20e6)$count_diffs)/
  sum(filter(pi, pop=="AZ" & chromosome == "CM022210.2_RagTag"& window_pos_1 >= 13.14e6 & window_pos_1 <=13.20e6)$count_comparisons)

sum(filter(pi, pop=="WTX" & chromosome == "CM022195.2_RagTag" & window_pos_1 >= 8.1e6 & window_pos_1 <=9.25e6)$count_diffs)/
  sum(filter(pi, pop=="WTX" & chromosome == "CM022195.2_RagTag"& window_pos_1 >= 8.1e6 & window_pos_1 <=9.25e6)$count_comparisons)
sum(filter(pi, pop=="AZ" & chromosome == "CM022195.2_RagTag"& window_pos_1 >= 8.1e6 & window_pos_1 <=9.25e6)$count_diffs)/
  sum(filter(pi, pop=="AZ" & chromosome == "CM022195.2_RagTag"& window_pos_1 >= 8.1e6 & window_pos_1 <=9.25e6)$count_comparisons)

sum(filter(pi, pop=="Sev" & chromosome == "CM022195.2_RagTag" & window_pos_1 >= 8.1e6 & window_pos_1 <=9.25e6)$count_diffs)/
  sum(filter(pi, pop=="Sev" & chromosome == "CM022195.2_RagTag"& window_pos_1 >= 8.1e6 & window_pos_1 <=9.25e6)$count_comparisons)
sum(filter(pi, pop=="Butte" & chromosome == "CM022195.2_RagTag"& window_pos_1 >= 8.1e6 & window_pos_1 <=9.25e6)$count_diffs)/
  sum(filter(pi, pop=="Butte" & chromosome == "CM022195.2_RagTag"& window_pos_1 >= 8.1e6 & window_pos_1 <=9.25e6)$count_comparisons)
