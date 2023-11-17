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

########################################################################
# NBR analyses
########################################################################

# Apply NBR
set.seed(18900217)
#outdir <-  file.path(getwd(),"Results","03-Modules","Supplementary")
#if(!dir.exists(outdir)) dir.create(path = outdir, recursive = T)
git_path <- "https://github.com/zchuri/BilingualBrainNet/raw/main/Results/03-Modules/Supplementary/"

#----------------------------
# Modules
#----------------------------

# Simultaneous > Monolinguals
outfile <- paste0(git_path,"nbr_groups_mono-sim_gt05.rds")
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
                         thrP = 0.05/2,
                         nperm = 1000,
                         nudist = T)
  toc <- Sys.time()
  show(toc-tic)
  # Save object if desired
  #saveRDS(nbr_result, outfile)
  
}

# Early > Monolinguals
outfile <- paste0(git_path,"nbr_groups_mono-early_gt05.rds")
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
                       thrP = 0.05/2,
                       nperm = 1000,
                       nudist = T)
  toc <- Sys.time()
  show(toc-tic)
  # Save object if desired
  #saveRDS(nbr_result, outfile)
  
}

#----------------------------
# ROI
#----------------------------

#outdir <-  file.path(getwd(),"Results","04-ROI","Supplementary")
#if(!dir.exists(outdir)) dir.create(path = outdir, recursive = T)
git_path <- "https://github.com/zchuri/BilingualBrainNet/raw/main/Results/"

# Simultaneous > Monolinguals
outfile <- paste0(git_path,"04-ROI/Supplementary/nbr_groups_mono-sim_gt_net01_roi05.rds")

# Subsetting
idx <- which(info$group=="Mono" | info$group=="Sim.B")
phen <- info[idx,]
phen$group <- factor(phen$group)
# Read netwise NBR results
nbr_net <- readRDS(url(paste0(git_path,"03-Modules/nbr_groups_mono-sim_gt01.rds")))
# Extract FWE networks
net_labs <- levels(as.factor(atlas$net))[-13]
fwe_idx <- nbr_net$fwe$groupsim[which(nbr_net$fwe$groupsim[,5]<0.05),1]
fwe_edges <- which(nbr_net$components$groupsim[,4]==fwe_idx)
fwe_net <- sort(unique(c(nbr_net$components$groupsim[fwe_edges,2:3])))

# Create network restricted subnetwork
roi_idx <- which(!is.na(match(atlas$net,net_labs[fwe_net])))
sub_atlas <- atlas[roi_idx,]
sub_mx <- cmx[roi_idx,roi_idx,idx]

# ROI-to-ROI Network-Based R-Statistics
if(!valid_url(outfile)){
  set.seed(18900217)
  tic <- Sys.time()
  nbr_roi <- nbr_lm(net = sub_mx,
                    nnodes = length(roi_idx),
                    idata = phen,
                    mod = "~ group",
                    alternative = "greater",
                    thrP = 0.05/2,
                    nperm = 1000,
                    #cores = 32,
                    nudist = T)
  toc <- Sys.time()
  show(toc-tic)
  # Save object
  #saveRDS(nbr_roi, outfile)
  
} else nbr_roi <- readRDS(url(outfile))


########################################################################
# Visualization
########################################################################

#----------------------------
# Modules
#----------------------------

# Load 'corrplot', 'circlize' & 'scales' package
if(!require("corrplot")) install.packages("corrplot"); library(corrplot)
if(!require("circlize")) install.packages("circlize"); library(circlize)
if(!require("scales")) install.packages("scales"); library(scales)

# Input parameters
git_path <- "https://github.com/zchuri/BilingualBrainNet/raw/main/Results/03-Modules/Supplementary/"
net_labs <- levels(as.factor(atlas$net))[-13]
colorcito <- colorRampPalette(c("blue4","blue","cyan","white","yellow","orange","red"))
# Chord diagram input parameters
net_fac <- as.factor(1:net_n)
levels(net_fac) <- net_labs

# Simultaneous > Monolinguals
idx <- which(info$group=="Mono" | info$group=="Sim.B")
phen <- info[idx,]
phen$group <- factor(phen$group)
sub_mx <- net_mx[,,idx]
# Read results
outfile <- paste0(git_path,"nbr_groups_mono-sim_gt05.rds")
nbr_result <- readRDS(url(outfile))

# Extract significant component first
fwe_tab <- nbr_result$fwe$groupSim.B
comp_idx <- fwe_tab[which(fwe_tab[,5]<0.05),1]
# Extract cluster connections
edge_tab <- nbr_result$components$groupSim.B
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
#outfile <- "nbr_mono-sim_gt05_edges.csv"
#write.csv(x = sig_df, file = outfile, row.names = F)

# T-values matrix
mat_stat <- array(as.numeric(NA),c(net_n,net_n))
for(rr in 1:net_n) for(cc in 1:net_n){
  mat_stat[rr,cc] <- summary(lm(sub_mx[rr,cc,]~phen$group))$coefficients[2,3]  
} 
colim <- ceiling(max(abs(mat_stat),na.rm = T)/0.5)*0.5
colnames(mat_stat) <- row.names(mat_stat) <- net_labs
# FWE matrix
mat_fwe <- array(FALSE,c(net_n,net_n))
for(ss in 1:nrow(sig_tab)) mat_fwe[sig_tab[ss,2],sig_tab[ss,3]] <- mat_fwe[sig_tab[ss,3],sig_tab[ss,2]] <- TRUE
tthr <- qt(1-0.05/2,nrow(phen)-2)
# Plot
corrplot(mat_stat, is.corr = F, type = "lower",
         method = "square",
         p.mat = (mat_stat>tthr)*mat_fwe,
         tl.col = "black", tl.cex = 0.7,
         insig = "pch", pch = "*", pch.cex = 2,
         cl.cex = 1, cl.length = 3, 
         col.lim = c(-colim,colim),
         diag = T, mar=c(3,1,1,1)+0.1,
         col = colorcito(100),
         xlab = "t-value",
         title = "Sim.B > Mono (edge p < 0.05)")

# Early > Monolinguals
idx <- which(info$group=="Mono" | info$group=="Early.B")
phen <- info[idx,]
phen$group <- factor(phen$group, levels = c("Mono","Early.B"))
sub_mx <- net_mx[,,idx]
# Read results
outfile <- paste0(git_path,"nbr_groups_mono-early_gt05.rds")
nbr_result <- readRDS(url(outfile))

# Extract significant component first
fwe_tab <- nbr_result$fwe$groupEarly.B
comp_idx <- fwe_tab[which(fwe_tab[,5]<0.05),1]
# Extract cluster connections
edge_tab <- nbr_result$components$groupEarly.B
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
#outfile <- "nbr_mono-early_gt05_edges.csv"
#write.csv(x = sig_df, file = outfile, row.names = F)

# T-values matrix
mat_stat <- array(as.numeric(NA),c(net_n,net_n))
for(rr in 1:net_n) for(cc in 1:net_n) mat_stat[rr,cc] <- summary(lm(sub_mx[rr,cc,]~phen$group))$coefficients[2,3]
colim <- ceiling(max(abs(mat_stat),na.rm = T)/0.5)*0.5
colnames(mat_stat) <- row.names(mat_stat) <- net_labs
# FWE matrix
mat_fwe <- array(FALSE,c(net_n,net_n))
for(ss in 1:nrow(sig_tab)) mat_fwe[sig_tab[ss,2],sig_tab[ss,3]] <- mat_fwe[sig_tab[ss,3],sig_tab[ss,2]] <- TRUE

# Plot
corrplot(mat_stat, is.corr = F, type = "lower",
         method = "square",
         p.mat = (mat_stat>tthr)*mat_fwe,
         tl.col = "black", tl.cex = 0.7,
         insig = "pch", pch = "*", pch.cex = 2,
         cl.cex = 1, cl.length = 3, 
         col.lim = c(-colim,colim),
         diag = T, mar=c(3,1,1,1)+0.1,
         col = colorcito(100),
         title = "Early.B > Mono (edge p < 0.05)")

#----------------------------
# ROI
#----------------------------

# Read results
git_path <- "https://github.com/zchuri/BilingualBrainNet/raw/main/Results/"
outfile <- paste0(git_path,"04-ROI/Supplementary/nbr_groups_mono-sim_gt_net01_roi05.rds")
nbr_roi <- readRDS(url(outfile))

# Subsetting
idx <- which(info$group=="Mono" | info$group=="Sim.B")
phen <- info[idx,]
phen$group <- factor(phen$group)
# Read netwise NBR results
nbr_net <- readRDS(url(paste0(git_path,"03-Modules/nbr_groups_mono-sim_gt01.rds")))
# Extract FWE networks
net_labs <- levels(as.factor(atlas$net))[-13]
fwe_idx <- nbr_net$fwe$groupsim[which(nbr_net$fwe$groupsim[,5]<0.05),1]
fwe_edges <- which(nbr_net$components$groupsim[,4]==fwe_idx)
fwe_net <- sort(unique(c(nbr_net$components$groupsim[fwe_edges,2:3])))

# Create network restricted subnetwork
roi_idx <- which(!is.na(match(atlas$net,net_labs[fwe_net])))
sub_atlas <- atlas[roi_idx,]
sub_mx <- cmx[roi_idx,roi_idx,idx]

# Find significant component
fwe_comp <- which(nbr_roi$fwe[[1]][,5]<0.05)
if(length(fwe_comp)>0){
  
  # Describe results
  for(cc in fwe_comp){
    # Component label
    comp_idx <- nbr_roi$fwe[[1]][cc,1]
    # Extract edges
    res_mat <- nbr_roi$components[[1]]
    res_mat <- res_mat[which(res_mat[,4]==comp_idx),]
    res_df <- data.frame(ROI1=rep("",nrow(res_mat)),
                         ROI2=rep("",nrow(res_mat)),
                         net1=rep("",nrow(res_mat)),
                         net2=rep("",nrow(res_mat)),
                         idx1=rep(0,nrow(res_mat)),
                         idx2=rep(0,nrow(res_mat)),
                         tval=rep(0,nrow(res_mat)),
                         df=rep(0,nrow(res_mat)),
                         pval=rep(1,nrow(res_mat)),
                         stringsAsFactors = F)
    # Extract features
    for(ee in 1:nrow(res_mat)){
      # ROI names
      res_df$ROI1[ee] <- as.character(sub_atlas$name.full[res_mat[ee,2]])
      res_df$ROI2[ee] <- as.character(sub_atlas$name.full[res_mat[ee,3]])
      # Network names
      res_df$net1[ee] <- as.character(sub_atlas$net[res_mat[ee,2]])
      res_df$net2[ee] <- as.character(sub_atlas$net[res_mat[ee,3]])
      # Indices
      res_df$idx1[ee] <- roi_idx[res_mat[ee,2]]
      res_df$idx2[ee] <- roi_idx[res_mat[ee,3]]
      # Run model
      phen$y <- sub_mx[res_mat[ee,2],res_mat[ee,3],]
      res_mod <- lm(y~group, phen)
      # t-test
      res_stat <- summary(res_mod)
      res_df$tval[ee] <- res_stat$coefficients[2,3]
      res_df$df[ee] <- res_stat$df[2]
      res_df$pval[ee] <- res_stat$coefficients[2,4]
    }
    # Save result summary
    #tabfile <- "nbr_groups_mono-sim_gt_net01_roi01_table.csv"
    #write.csv(x = res_df,
    #          file = tabfile,
    #          row.names = F)
    
    # Create .edge file based on the results
    edge_mat <- matrix(data = 0,
                       nrow = nrow(atlas),
                       ncol = nrow(atlas))
    # Set t-value
    for(ee in 1:nrow(res_df)) edge_mat[res_df$idx1[ee],res_df$idx2[ee]] <- edge_mat[res_df$idx2[ee],res_df$idx1[ee]] <- res_df$tval[ee]
    #edgefile <- "nbr_groups_mono-sim_gt_net01_roi01_tval.edge"
    #write.table(x = edge_mat,
    #            file = edgefile,
    #            quote = F, sep = "\t",
    #            row.names = F, col.names = F)
    
    # Extract top ROIs
    roi_freq <- sort(table(c(res_df$idx1,res_df$idx2)), T)
    for(tt in 1:max(roi_freq)){
      # Extract top ROI index
      top_roi <- as.numeric(names(roi_freq)[tt])
      # Create .edge file for single ROI
      edge_roi <- matrix(0,nrow(atlas),nrow(atlas))
      # Add edges
      edge_roi[top_roi,] <- edge_mat[top_roi,]
      edge_roi[,top_roi] <- edge_mat[,top_roi]
      # Write
      #topfile <- paste0("nbr_groups_mono-sim_gt_net",thr_str[1],"_roi",thr_str[2],"_ROI-",top_roi,".edge")
      #write.table(x = edge_roi,
      #            file = topfile,
      #            quote = F, sep = "\t",
      #            row.names = F, col.names = F) 
    }# for tt in top ROIs
    
  }# for cc in fwe_comp
}# if length(fwe_comp)>0

# Load libraries
if(!require(ggraph)) install.packages("ggraph"); library(ggraph)
if(!require(igraph)) install.packages("igraph"); library(igraph)

# Plot chord diagram in an hemispheric fashion
# Mirrored order of ROIs by hemisphere and network
RH_idx <- which(atlas$hemi=="R")
LH_idx <- which(atlas$hemi=="L")
net_idx <- c(RH_idx[order(atlas$net[RH_idx])],
             LH_idx[rev(order(atlas$net[LH_idx]))])
# Remove uncertain network indices
net_idx <- setdiff(net_idx,which(atlas$net=="UNC"))
# Order by hemisphere/network
sort_atlas <- atlas[net_idx,]
sort_atlas$net <- factor(sort_atlas$net)
sort_edge <- edge_mat[net_idx,net_idx]
# Create graph from adjacency matrix
ig <- graph.adjacency(sort_edge, mode="undirected", weighted=T)
# Color limits and palette
colim <- ceiling(max(abs(sort_edge),na.rm = T)/0.5)*0.5
colorcito <- colorRampPalette(c("blue4","blue","cyan","white","yellow","orange","red"))
colnet <- c("blue3","brown3","chartreuse4",
                   "chocolate2","darkcyan","darkred",
                   "darkgoldenrod2","darkolivegreen1",
                   "darkorchid2","khaki1","cornflowerblue",
                   "lawngreen","yellow1","red1","black")
                   # A chord diagram
gg1 <- ggraph(ig, layout = 'linear', circular = TRUE) + 
  theme_void() +
  geom_edge_arc(aes(colour = E(ig)$weight),
                edge_width = 1) +
  scale_edge_colour_gradientn(colours = colorcito(100),
                              space = "Lab",
                              limits=c(-colim,colim),
                              na.value = "grey50",
                              guide = "edge_colourbar") +
  geom_node_point(aes(x = x*1.07, y=y*1.07, colour=sort_atlas$net), shape = 19, size=4) +
  scale_color_manual(values=colnet[1:nlevels(sort_atlas$net)]) +
  ggtitle("Sim.B > Mono (edge p < 0.05)")
show(gg1)
