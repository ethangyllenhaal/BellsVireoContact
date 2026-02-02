# By: Ethan Gyllenhaal
# Updated: 1Feb2026
#
# Script for making genomic PCAs and density plots
# First one (no Z) explains basic process in depth, others only have notes breaking from general patterns noted.

# Load packages, some may not be used here or automatically loaded as dependencies
library("adegenet")
library("ade4")
library("vcfR")
library("scales")
library("parallel")
library("viridis")
library(tidyverse)
library(gridExtra)
library(cowplot)

# set working directory
setwd('/path/to/workdir')

#formula for approximate p-values from "cell" in PCA
# p = 1 - exp( -[ c[cell]^2 /2 )
# 1.5 = 67%
# 2.5 = 95%

# set colors and labels for main PCA and regional PCAs
colors <- c("#c83200", "#c83264", "#c86496", "#646496", "#64b496", "#32deb4", "#32dede", "#6496fa", "#3232ff")
labels <- c("Arizona", "Gila", "Mimbres", "Elephant", "Sevilleta", "Pecos", "West Texas", "South Texas", "Missouri")
color_coarse <- c("#c83200", "#646496", "#32deb4")
abbrev_labels <- c("West", "Contact", "East")


### 75% focal noZ
bevi_focal_75_vcf <- read.vcfR("bevi_noZ_focal_75.vcf") # read VCF
bevi_focal_75_gl <- vcfR2genlight(bevi_focal_75_vcf) # convert to genlight
popmap <- read.csv('popmap.csv', sep = ',') # load relevant popmap
pop(bevi_focal_75_gl) <- popmap$pop # add population column to genlight object
pca_focal_75_bevi <- glPca(bevi_focal_75_gl, n.cores=4, nf=4) # run PCA on 4 cores with 4 retained PCs
scatter(pca_focal_75_bevi, cex=.25) # scatter plot for troubleshooting
# Make nice looking plot with s.class, note some options unique to this call
s.class(pca_focal_75_bevi$scores[,1:2], pop(bevi_focal_75_gl), col=colors, label=labels, 
        clab=1.5, cell=2.5, cpoint=1.5, axesell=0, pch=19)
# Make barplot of PC % variance explained
barplot(pca_focal_75_bevi$eig/sum(pca_focal_75_bevi$eig), main="eigenvalues", col=heat.colors(length(pca_focal_75_bevi$eig)))
# Print out scores
pca_focal_75_bevi$eig/sum(pca_focal_75_bevi$eig)

### 75% focal Z
bevi_Z_focal_75_vcf <- read.vcfR("bevi_Z_focal_75.vcf")
bevi_Z_focal_75_gl <- vcfR2genlight(bevi_Z_focal_75_vcf)
popmap_Z <- read.csv('popmap_Z.csv', sep = ',')
pop(bevi_Z_focal_75_gl) <- popmap_Z$pop
pca_Z_focal_75_bevi <- glPca(bevi_Z_focal_75_gl, n.cores=4, nf=4)
scatter(pca_Z_focal_75_bevi, cex=.1)
s.class(pca_Z_focal_75_bevi$scores[,1:2], pop(bevi_Z_focal_75_gl), col=colors, label=labels, 
        clab=1, cell=0, pch=19)
barplot(pca_Z_focal_75_bevi$eig/sum(pca_Z_focal_75_bevi$eig), main="eigenvalues", col=heat.colors(length(pca_Z_focal_75_bevi$eig)))

# specific regions

### 75% 13M Z

bevi_Z13M_focal_75_vcf <- read.vcfR("bevi_Z13M_focal_75.vcf")
bevi_Z13M_focal_75_gl <- vcfR2genlight(bevi_Z13M_focal_75_vcf)
popmap_Z <- read.csv('popmap_Z.csv', sep = ',')
pop(bevi_Z13M_focal_75_gl) <- popmap_Z$coarse
pca_Z13M_focal_75_bevi <- glPca(bevi_Z13M_focal_75_gl, n.cores=4, nf=4)

# make data frame of scores
Z13M_scores <- cbind(pca_Z13M_focal_75_bevi$scores[,1:4], bevi_Z13M_focal_75_gl$pop) %>% as.data.frame()

# PCA with ggplot from scores data frame
# geom point -> 90% ellipse -> remake points over top of ellipse for clarity
Z13M <- ggplot(data=filter(Z13M_scores, V5!=2), aes(x=as.numeric(PC1), y=as.numeric(PC2), color=as.factor(V5), fill=as.factor(V5))) +
  geom_point(stroke=NA, size=3) +
  scale_color_manual(values = color_coarse[c(1,3,2)], labels=c("West", "East", "Contact")) +
  scale_fill_manual(values = color_coarse[c(1,3,2)], labels=c("West", "East", "Contact")) +
  stat_ellipse(geom="polygon", alpha=0.1) +
  geom_point(data=filter(Z13M_scores, V5==2), size=3) + 
  ggtitle("ii") + theme_bw() + scale_x_reverse() +
  theme(axis.title=element_blank(), legend.position = "none", axis.text=element_blank(), axis.ticks = element_blank())
Z13M

# density plot for PC1 scores, note that X axis is flipped to match PCA and west-east axis
Z13M_dens <- ggplot(data=Z13M_scores, aes(x=-as.numeric(PC1))) +
  geom_density(fill="#03688C", adjust=0.3) + theme_bw() + scale_x_reverse() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text = element_text(size=10))
Z13M_dens

### 75% 16M Chr 1

bevi_Chr1_16M_focal_75_vcf <- read.vcfR("bevi_Chr1-16M_focal_75.vcf")
bevi_Chr1_16M_focal_75_gl <- vcfR2genlight(bevi_Chr1_16M_focal_75_vcf)
popmap <- read.csv('popmap.csv', sep = ',')
pop(bevi_Chr1_16M_focal_75_gl) <- popmap$coarse
pca_Chr1_16M_focal_75_bevi <- glPca(bevi_Chr1_16M_focal_75_gl, n.cores=4, nf=4)
Chr1_16M_scores <- cbind(pca_Chr1_16M_focal_75_bevi$scores[,1:4], bevi_Chr1_16M_focal_75_gl$pop) %>% as.data.frame()

Chr1_16M <- ggplot(data=filter(Chr1_16M_scores, V5!=2), aes(x=as.numeric(PC1), y=as.numeric(PC2), color=as.factor(V5), fill=as.factor(V5))) +
  geom_point(stroke=NA, size=3) +
  scale_color_manual(values = color_coarse[c(1,3,2)], labels=c("West", "East", "Contact")) +
  scale_fill_manual(values = color_coarse[c(1,3,2)], labels=c("West", "East", "Contact")) +
  stat_ellipse(geom="polygon", alpha=0.1) +
  geom_point(data=filter(Chr1_16M_scores, V5==2), size=3) + 
  ggtitle("i") + theme_bw() + scale_x_reverse() +
  theme(axis.title=element_blank(), legend.position = "none", axis.text=element_blank(), axis.ticks = element_blank())
Chr1_16M

Chr1_16M_dens <- ggplot(data=Chr1_16M_scores, aes(x=as.numeric(PC1))) +
  geom_density(fill="#CC9B2B", adjust=0.3) + theme_bw()+ scale_x_reverse() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text = element_text(size=10))
Chr1_16M_dens

### 75% 8M Chr 12

bevi_Chr12_8M_focal_75_vcf <- read.vcfR("bevi_Chr12-8M_focal_75.vcf")
bevi_Chr12_8M_focal_75_gl <- vcfR2genlight(bevi_Chr12_8M_focal_75_vcf)
popmap <- read.csv('popmap.csv', sep = ',')
pop(bevi_Chr12_8M_focal_75_gl) <- popmap$coarse
pca_Chr12_8M_focal_75_bevi <- glPca(bevi_Chr12_8M_focal_75_gl, n.cores=4, nf=4)
Chr12_8M_scores <- cbind(pca_Chr12_8M_focal_75_bevi$scores[,1:4], bevi_Chr12_8M_focal_75_gl$pop) %>% as.data.frame()

Chr12_8M <- ggplot(data=filter(Chr12_8M_scores, V5!=2), aes(x=as.numeric(PC1), y=as.numeric(PC2), color=as.factor(V5), fill=as.factor(V5))) +
  geom_point(stroke=NA, size=3) +
  scale_color_manual(values = color_coarse[c(1,3,2)], labels=c("West", "East", "Contact")) +
  scale_fill_manual(values = color_coarse[c(1,3,2)], labels=c("West", "East", "Contact")) +
  stat_ellipse(geom="polygon", alpha=0.1) +
  geom_point(data=filter(Z13M_scores, V5==2), size=3) + 
  ggtitle("iii") + theme_bw() + scale_x_reverse() +
  theme(axis.title=element_blank(), legend.position = "none", axis.text=element_blank(), axis.ticks = element_blank())
Chr12_8M

Chr12_8M_dens <- ggplot(data=Chr12_8M_scores, aes(x=as.numeric(PC1))) +
  geom_density(fill="#CC9B2B", adjust=0.3) + theme_bw() + scale_x_reverse() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text = element_text(size=10))
Chr12_8M_dens

# plot all regions everything together with different relative heights
plot_grid(Chr1_16M, Z13M, Chr12_8M, Chr1_16M_dens, Z13M_dens, Chr12_8M_dens, nrow=2, rel_heights = c(3/5,2/5))
