#!/usr/bin/Rscript

# Clear workspace
rm(list = ls())

# Load libraries
if(!require(NBR)) install.packages("NBR"); library(NBR)
source("https://raw.githubusercontent.com/zchuri/BilingualBrainNet/main/Code/AuxFunctions/valid_url.R")

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

# Remove negative values
cmx[cmx<0] <- 0
# Subtract diagonal
for(ii in 1:cmx_dim[3]) diag(cmx[,,ii]) <- 0

# Generate netwise object
net_n <- nlevels(as.factor(atlas$net))-1
net_mx <- array(as.numeric(NA), c(net_n,net_n,cmx_dim[3]))

# Extract functional connectivity
for(ii in 1:cmx_dim[3]){
  if(ii == 1) cat("Summarizing matrices:")
  cat(paste0("..",ii))
  # Summarize
  ex_cmx <- aggregate(x = cmx[,,ii],
                      by = list(atlas$net),
                      FUN = function(x) mean(x, na.rm=T))[-1]
  ex_cmx <- aggregate(x = t(ex_cmx),
                      by = list(atlas$net),
                      FUN = function(x) mean(x, na.rm=T))[-1]
  # Discard uncertain network
  ex_cmx <- as.matrix(ex_cmx)[,-13]
  net_mx[,,ii] <- ex_cmx[-13,]
  if(ii == cmx_dim[3]) cat("\n")
}

# Apply NBR
set.seed(18900217)
#outdir <-  file.path(getwd(),"Results","03-Modules")
#if(!dir.exists(outdir)) dir.create(path = outdir, recursive = T)
git_path <- "https://github.com/zchuri/BilingualBrainNet/raw/main/Results/03-Modules/"

# Simultaneous > Monolinguals
outfile <- paste0(git_path,"nbr_groups_mono-sim_gt01.rds")
# If file is not found, can be computed locally
if(!valid_url(outfile)){
  # Subsetting
  idx <- which(info$group=="Mono" | info$group=="Sim.B")
  phen <- info[idx,]
  phen$group <- factor(phen$group)
  sub_mx <- net_mx[,,idx]
  # Run NBR analysis
  tic <- Sys.time()
  nbr_result <- nbr_lm(net = sub_mx,
                         nnodes = net_n,
                         idata = phen,
                         mod = "~ group",
                         alternative = "greater",
                         diag = T,
                         thrP = 0.01/2,
                         nperm = 1000,
                         nudist = T)
  toc <- Sys.time()
  show(toc-tic)
  # Save object if desired
  #saveRDS(nbr_result, outfile)
  
}

# Early > Monolinguals
outfile <- paste0(git_path,"nbr_groups_mono-early_gt01.rds")
# If file is not found, can be computed locally
if(!valid_url(outfile)){
  # Subsetting
  idx <- which(info$group=="Mono" | info$group=="Early.B")
  phen <- info[idx,]
  phen$group <- factor(phen$group, levels = c("Mono","Early.B"))
  sub_mx <- net_mx[,,idx]
  # Run NBR analysis
  tic <- Sys.time()
  nbr_result <- nbr_lm(net = sub_mx,
                       nnodes = net_n,
                       idata = phen,
                       mod = "~ group",
                       alternative = "greater",
                       diag = T,
                       thrP = 0.01/2,
                       nperm = 1000,
                       nudist = T)
  toc <- Sys.time()
  show(toc-tic)
  # Save object if desired
  #saveRDS(nbr_result, outfile)
  
}

# Early-Late AoA
outfile <- paste0(git_path,"nbr_AoA_early-late_lt01.rds")
# If file is not found, can be computed locally
if(!valid_url(outfile)){
  # Subsetting
  idx <- which(info$group=="Late.B" | info$group=="Early.B")
  phen <- info[idx,]
  sub_mx <- net_mx[,,idx]
  # Run NBR analysis
  tic <- Sys.time()
  nbr_result <- nbr_lm(net = sub_mx,
                       nnodes = net_n,
                       idata = phen,
                       mod = "~ AOA",
                       alternative = "lower",
                       diag = T,
                       thrP = 0.01/2,
                       nperm = 1000,
                       nudist = T)
  toc <- Sys.time()
  show(toc-tic)
  # Save object if desired
  #saveRDS(nbr_result, outfile)
  
}

# Bilingual YoE
outfile <- paste0(git_path,"nbr_YoE_bili_gt01.rds")
# If file is not found, can be computed locally
if(!valid_url(outfile)){
  # Subsetting
  idx <- which(info$YOE>0)
  phen <- info[idx,]
  sub_mx <- net_mx[,,idx]
  
  # Run NBR analysis
  tic <- Sys.time()
  nbr_result <- nbr_lm(net = sub_mx,
                       nnodes = net_n,
                       idata = phen,
                       mod = "~ YOE",
                       alternative = "greater",
                       diag = T,
                       thrP = 0.01/2,
                       nperm = 1000,
                       nudist = T)
  toc <- Sys.time()
  show(toc-tic)
  # Save object if desired
  #saveRDS(nbr_result, outfile)
  
}

# Bilingual L2 proficiency
outfile <- paste0(git_path,"nbr_L2_bili_lt01.rds")
# If file is not found, can be computed locally
if(!valid_url(outfile)){
  # Subsetting
  idx <- which(!is.na(info$L2proficiency))
  phen <- info[idx,]
  sub_mx <- net_mx[,,idx]
  # Run NBR analysis
  tic <- Sys.time()
  nbr_result <- nbr_lm(net = sub_mx,
                       nnodes = net_n,
                       idata = phen,
                       mod = "~ rank(L2proficiency)",
                       alternative = "lower",
                       diag = T,
                       thrP = 0.01/2,
                       nperm = 1000,
                       nudist = T)
  toc <- Sys.time()
  show(toc-tic)
  # Save object if desired
  #saveRDS(nbr_result, outfile)
  
}

########################################################################
# Visualization
########################################################################

# Load 'corrplot', 'circlize' & 'scales' package
if(!require("corrplot")) install.packages("corrplot"); library(corrplot)
if(!require("circlize")) install.packages("circlize"); library(circlize)
if(!require("scales")) install.packages("scales"); library(scales)

# Input parameters
net_labs <- levels(as.factor(atlas$net))[-13]
colorcito <- colorRampPalette(c("blue4","blue","cyan","white","yellow","orange","red"))
# Chord diagram input parameters
net_fac <- as.factor(1:net_n)
levels(net_fac) <- net_labs

########################################################################
# Simultaneous > Monolinguals
idx <- which(info$group=="Mono" | info$group=="Sim.B")
phen <- info[idx,]
phen$group <- factor(phen$group)
sub_mx <- net_mx[,,idx]
# Read results
outfile <- paste0(git_path,"nbr_groups_mono-sim_gt01.rds")
nbr_result <- readRDS(url(outfile))

# Extract significant component first
fwe_tab <- nbr_result$fwe$groupsim
comp_idx <- fwe_tab[which(fwe_tab[,5]<0.05),1]
# Extract cluster connections
edge_tab <- nbr_result$components$groupsim
sig_tab <- edge_tab[which(edge_tab[,4]==comp_idx),]
# Create data.frame
sig_df <- data.frame(ROI1=vector("character",nrow(sig_tab)),
                     ROI2=vector("character",nrow(sig_tab)),
                     tval=vector("numeric",nrow(sig_tab)),
                     df=rep(nrow(phen)-2,nrow(sig_tab)),
                     pval=vector("numeric",nrow(sig_tab)),
                     stringsAsFactors = F)
# Fill edge info
for(ee in 1:nrow(sig_df)){
  # ROI1
  sig_df$ROI1[ee] <- net_labs[sig_tab[ee,2]]
  # ROI2
  sig_df$ROI2[ee] <- net_labs[sig_tab[ee,3]]
  # Fit LM
  fit_e <- summary(lm(sub_mx[sig_tab[ee,2],sig_tab[ee,3],]~phen$group))$coefficients
  sig_df$tval[ee] <- fit_e[2,3] # t-score
  sig_df$pval[ee] <- fit_e[2,4] # p-value
}# ee
# Print table
show(sig_df)
#outfile <- "nbr_mono-sim_gt01_edges.csv"
#write.csv(x = sig_df, file = outfile, row.names = F)

# Corrplot
mat_stat <- array(as.numeric(NA),c(net_n,net_n))
for(rr in 1:net_n) for(cc in 1:net_n) mat_stat[rr,cc] <- summary(lm(sub_mx[rr,cc,]~phen$group))$coefficients[2,3]
colim <- ceiling(max(abs(mat_stat),na.rm = T)/0.5)*0.5
colnames(mat_stat) <- row.names(mat_stat) <- net_labs
tthr <- qt(1-0.01/2,nrow(phen)-2) # All edges above this threshold are the ones in compoenent 2 of NBR object: nbr_result$fwe
# Plot
corrplot(mat_stat, is.corr = F, type = "lower",
         method = "square",
         p.mat = mat_stat>tthr,
         tl.col = "black", tl.cex = 0.7,
         insig = "pch", pch = "*", pch.cex = 2,
         cl.cex = 1, cl.length = 3, 
         col.lim = c(-colim,colim),
         diag = T, mar=c(3,1,1,1)+0.1,
         col = colorcito(100),
         title = "Sim.B > Mono")

# Chord diagram
# svg(outfile, 4, 4) # Export width/height to 4 for an optimal view
# Initialize the plot
circos.clear()
circos.initialize(factors = net_fac, xlim=c(-1,1))
# Build the regions of track #1
circos.trackPlotRegion(factors = net_fac, ylim = c(-1,1),
                       bg.col = "aliceblue",
                       bg.border = "black")
# Add labels
net_pos <- ux(16.9, "mm")*(1:nlevels(net_fac))
circos.text(net_pos,rep(0,nlevels(net_fac)),
            labels = net_labs,
            facing = "bending.inside", cex=1.4)
# Add a links between a point and another
# Find significant weigths
up_tri <- which(upper.tri(mat_stat), arr.ind = T)
sig_tri <- which(mat_stat[upper.tri(mat_stat)]>tthr)
# Draw links
for(ii in 1:length(sig_tri)){
  # Set central position
  rnum <- runif(1); p1 <- c(rnum,rnum-1)
  rnum <- runif(1); p2 <- c(rnum,rnum-1)
  # Set color
  tval <- mat_stat[up_tri[sig_tri[ii],1],up_tri[sig_tri[ii],2]]
  tper <- round(50*tval/colim)+50
  tcol <- colorcito(100)[tper]
  # Draw link
  circos.link(net_fac[up_tri[sig_tri[ii],1]], p1,
              net_fac[up_tri[sig_tri[ii],2]], p2,
              col = scales::alpha(tcol,.9),
              border = scales::alpha(tcol,.4),
              h.ratio=0.6)
}
# Add title
title("Sim.B > Mono")
# Clear circle parameters
circos.clear()
#dev.off()

########################################################################
# Early > Monolinguals
idx <- which(info$group=="Mono" | info$group=="Early.B")
phen <- info[idx,]
phen$group <- factor(phen$group, levels = c("Mono","Early.B"))
sub_mx <- net_mx[,,idx]
# Read results
outfile <- paste0(git_path,"nbr_groups_mono-early_gt01.rds")
nbr_result <- readRDS(url(outfile))

# Extract significant component first
fwe_tab <- nbr_result$fwe$groupearly
comp_idx <- fwe_tab[which(fwe_tab[,5]<0.05),1]
# Extract cluster connections
edge_tab <- nbr_result$components$groupearly
sig_tab <- edge_tab[which(edge_tab[,4]==comp_idx),]
# Create data.frame
sig_df <- data.frame(ROI1=vector("character",nrow(sig_tab)),
                     ROI2=vector("character",nrow(sig_tab)),
                     tval=vector("numeric",nrow(sig_tab)),
                     df=rep(nrow(phen)-2,nrow(sig_tab)),
                     pval=vector("numeric",nrow(sig_tab)),
                     stringsAsFactors = F)
# Fill edge info
for(ee in 1:nrow(sig_df)){
  # ROI1
  sig_df$ROI1[ee] <- net_labs[sig_tab[ee,2]]
  # ROI2
  sig_df$ROI2[ee] <- net_labs[sig_tab[ee,3]]
  # Fit LM
  fit_e <- summary(lm(sub_mx[sig_tab[ee,2],sig_tab[ee,3],]~phen$group))$coefficients
  sig_df$tval[ee] <- fit_e[2,3] # t-score
  sig_df$pval[ee] <- fit_e[2,4] # p-value
}# ee
# Print table
show(sig_df)
#outfile <- "nbr_mono-early_gt01_edges.csv"
#write.csv(x = sig_df, file = outfile, row.names = F)

# Corrplot
mat_stat <- array(as.numeric(NA),c(net_n,net_n))
for(rr in 1:net_n) for(cc in 1:net_n) mat_stat[rr,cc] <- summary(lm(sub_mx[rr,cc,]~phen$group))$coefficients[2,3]
colim <- ceiling(max(abs(mat_stat),na.rm = T)/0.5)*0.5
colnames(mat_stat) <- row.names(mat_stat) <- net_labs
tthr <- qt(1-0.01/2,nrow(phen)-2) # All edges above this threshold are the ones in compoenent 2 of NBR object: nbr_result$fwe
# Plot
corrplot(mat_stat, is.corr = F, type = "lower",
         method = "square",
         p.mat = mat_stat>tthr,
         tl.col = "black", tl.cex = 0.7,
         insig = "pch", pch = "*", pch.cex = 2,
         cl.cex = 1, cl.length = 3, 
         col.lim = c(-colim,colim),
         diag = T, mar=c(3,1,1,1)+0.1,
         col = colorcito(100),
         title = "Early.B > Mono")

# Chord diagram
# svg(outfile, 4, 4) # Export width/height to 4 for an optimal view
# Initialize the plot
circos.clear()
circos.initialize(factors = net_fac, xlim=c(-1,1))
# Build the regions of track #1
circos.trackPlotRegion(factors = net_fac, ylim = c(-1,1),
                       bg.col = "aliceblue",
                       bg.border = "black")
# Add labels
net_pos <- ux(16.9, "mm")*(1:nlevels(net_fac))
circos.text(net_pos,rep(0,nlevels(net_fac)),
            labels = net_labs,
            facing = "bending.inside", cex=1.4)
# Add a links between a point and another
# Find significant weigths
up_tri <- which(upper.tri(mat_stat), arr.ind = T)
sig_tri <- which(mat_stat[upper.tri(mat_stat)]>tthr)
# Draw links
for(ii in 1:length(sig_tri)){
  # Set central position
  rnum <- runif(1); p1 <- c(rnum,rnum-1)
  rnum <- runif(1); p2 <- c(rnum,rnum-1)
  # Set color
  tval <- mat_stat[up_tri[sig_tri[ii],1],up_tri[sig_tri[ii],2]]
  tper <- round(50*tval/colim)+50
  tcol <- colorcito(100)[tper]
  # Draw link
  circos.link(net_fac[up_tri[sig_tri[ii],1]], p1,
              net_fac[up_tri[sig_tri[ii],2]], p2,
              col = scales::alpha(tcol,.9),
              border = scales::alpha(tcol,.4),
              h.ratio=0.6)
}
# Add title
title("Early.B > Mono")
# Clear circle parameters
circos.clear()
#dev.off()
