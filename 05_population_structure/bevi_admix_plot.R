##########
# By: Ethan Gyllenhaal
# Updated 1Feb2026
#
# R script used for generating admixture plots for Bell's Vireos.
# Makes Admixture bar plots
# Note that order and color of bars was changed in illustrator for clarity of demonstration.
########

library(viridis)

setwd('/path/to/workdir/')

# autosome only
noZ2=read.table("bevi_noZ_focal_75_012_sort.2.Q")
noZ3=read.table("bevi_noZ_focal_75_012_sort.3.Q")
noZ4=read.table("bevi_noZ_focal_75_012_sort.4.Q")

# Z only
Z2=read.table("bevi_Z_focal_75_012_sort.2.Q")
Z3=read.table("bevi_Z_focal_75_012_sort.3.Q")
Z4=read.table("bevi_Z_focal_75_012_sort.4.Q")

# plot autosome then Z in that order
par(mfrow=c(2,1), mar=c(1,2,0,0), oma=c(1,2,1,0))
barplot(t(as.matrix(noZ2)), col=c("#32dede", "#c83200"),ylab="Ancestry", border=NA)
barplot(t(as.matrix(Z2)), col=c("#32dede", "#c83200"),ylab="Ancestry", border=NA)
