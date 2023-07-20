#!/usr/bin/Rscript

# Clear workspace
rm(list = ls())

# Load libraries
if(!require(igraph)) install.packages("igraph"); library(igraph)

# Read phenotypic data
info <- read.csv("https://raw.githubusercontent.com/zchuri/BilingualBrainNet/main/Data/Phenotypic/phenotypic.csv")
atlas <- read.csv("https://raw.githubusercontent.com/zchuri/BilingualBrainNet/main/Data/Atlas/atlas.csv")

# Create empty volume to store connectivity matrices
cmx <- array(data = as.numeric(NA),
             dim = c(nrow(atlas),nrow(atlas),nrow(info)))
# Read Brain data
cmx_path <- "https://raw.githubusercontent.com/zchuri/BilingualBrainNet/main/Data/Connectomes/"
for(ii in 1:nrow(info)){
  if(ii == 1) cat("Reading matrices:")
  cat(paste0("..",ii))
  # Generate filename
  num3char <- formatC(ii, width = 3, format = "d", flag = "0")
  filename <- paste0(cmx_path,"BrainNet_FCz_",num3char,".csv")
  # Read matrix
  cmx[,,ii] <- as.matrix(read.table(file = filename, sep = ","))
  if(ii == nrow(info)) cat("\n")
}
# Extract dimension
cmx_dim <- dim(cmx)

# Load graph theory formulas
source("https://raw.githubusercontent.com/zchuri/BilingualBrainNet/main/Code/AuxFunctions/inet_w.R")

# Remove negative values
cmx[cmx<0] <- 0
# Substract diagonal
for(ii in 1:cmx_dim[3]) diag(cmx[,,ii]) <- 0

# Inverse average distances (per node) - weighted
tic <- Sys.time()
w_E <- t(sapply(1:cmx_dim[3], function(x) iEff_w(cmx[,,x])))
toc <- Sys.time()
print("Efficiency:"); print(toc-tic)
info$E <- apply(w_E,1,mean)

# Modularity
tic <- Sys.time()
info$Q <- sapply(1:cmx_dim[3], function(x) iQ_w(ConnMx = cmx[,,x], modules = as.integer(as.factor(atlas$net))))
toc <- Sys.time()
print("Modularity:"); print(toc-tic)

# Save features
#write.csv(x = info,
#          file = "Data/Phenotypic/phenotypic.csv",
#          row.names = F)
