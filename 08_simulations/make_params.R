# simple script for making parameters for simulation
setwd("~/scratch/bevi")
# make data frame of grid
df <- expand.grid(c(0.0,0.01), c(0.0001), c(1.0, 0.5, 0.1), c(0.0, 0.5, 1), c(1:25))
# output grid
write.table(df, "params.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
