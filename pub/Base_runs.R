rm(list=ls())
library(TMB)
library(devtools)
library(TMBhelper)
library(foreach)
library(doParallel)
library(RColorBrewer)
library(beanplot)


devtools::install_github("merrillrudd/LIME", build.vignettes=TRUE, dependencies=TRUE)
devtools::install_github("adrianhordyk/LBSPR", build.vignettes=TRUE, dependencies=TRUE)

library(LIME)
library(LBSPR)


### ----- directories and functions -----------###
# main_dir <- "C:\\Git_Projects\\LIME_sim"
main_dir <- "F:\\Merrill\\Git_Projects\\LIME_sim"

if(grepl("C:", main_dir)) ncores <- 2
if(grepl("F:", main_dir)) ncores <- 8

funs <- list.files(file.path(main_dir, "R"))
ignore <- sapply(1:length(funs), function(x) source(file.path(main_dir,"R", funs[x])))

fig_dir <- file.path(main_dir, "figs")
dir.create(fig_dir, showWarnings=FALSE)

res_dir <- file.path(main_dir, "calc_results")
dir.create(res_dir, showWarnings=FALSE)


############################################################
#### equilibrium - instantaneous LF
############################################################

### ----- models to run ----------- ###
lh_vec <- c("Short", "Medium", "Long")
Fdyn_vec <- "Constant"
Rdyn_vec <- "Constant"
data_vec <- c("Index_Catch_LC20", "Index_LC10", "Index_LC1", "Catch_LC10", "Catch_LC1", "LC10", "LC1", "LBSPR10", "LBSPR1")
SampleSize_vec <- 200
itervec <- 1:100

equil_modcombos <- expand.grid("Data_avail"=data_vec, "SampleSize"=paste0("SampleSize_", SampleSize_vec), "LH"=paste0("LH_",lh_vec), "Fdyn"=paste0("F_",Fdyn_vec), "Rdyn"=paste0("R_",Rdyn_vec), stringsAsFactors=FALSE)
equil_modcombos$C_opt <- rep(0, nrow(equil_modcombos))
	equil_modcombos$C_opt[which(grepl("Catch",equil_modcombos[,"Data_avail"]))] <- 2

### ----- equilibrium test ----------- ###

	equil_dir <- file.path(main_dir, "equil")
	# unlink(equil_dir, TRUE)
	dir.create(equil_dir, showWarnings=FALSE)	

	## setup equilibrium dirs
	equil_dir_vec <- model_paths(res_dir=equil_dir, modcombos=equil_modcombos[,-c(ncol(equil_modcombos))])
	ignore <- sapply(1:length(equil_dir_vec), function(m) sapply(itervec, function(x) dir.create(file.path(equil_dir_vec[m], x), showWarnings=FALSE)))	

	## setup scenario life history list
	## 1-parameter model for selectivity and maturity
	lh_list <- adj_variation(SigmaR=0.01, SigmaF=0.01, SigmaC=0.01, SigmaI=0.01, CVlen=0.1, rho=0, selex_param=2, mat_param=1, nseasons=1)
	lh_fig(lh_list)

	### ----- generate data -----------###
	### data rich cases only
	rich_dir <- equil_dir_vec[grepl("Index_Catch_LC20",equil_dir_vec)]
	rich_modcombos <- equil_modcombos[which(equil_modcombos$Data_avail=="Index_Catch_LC20"),]
	alt_dir <- equil_dir_vec[grepl("Index_Catch_LC20",equil_dir_vec)==FALSE]
	alt_modcombos <- equil_modcombos[which(equil_modcombos$Data_avail!="Index_Catch_LC20"),]

	### ----- use parallel cores -----------###
	registerDoParallel(cores=ncores)	

	start_gen <- Sys.time()
	foreach(loop=1:length(rich_dir), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=rich_dir[loop], data_avail=as.character(rich_modcombos[loop,"Data_avail"]), itervec=itervec, Fdynamics=as.character(strsplit(rich_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(rich_modcombos[loop,"Rdyn"],"_")[[1]][2]),  lh=lh_list[[as.character(strsplit(rich_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=20, Nyears_comp=20, comp_sample=as.numeric(strsplit(rich_modcombos[loop,"SampleSize"],"_")[[1]][2]), rewrite=FALSE, init_depl=c(0.05,0.95), pool=TRUE, mismatch=FALSE)
	end_gen <- Sys.time() - start_gen

	# copy data rich cases to data poor
	# when written using parallel cores, sometimes files not writing properly
	start_regen <- Sys.time()
	for(loop in 1:length(rich_dir)){
		copy_sim(fromdir=rich_dir[loop], fromcombos=rich_modcombos[loop,], todir=alt_dir, itervec=itervec, rewrite=FALSE, res_dir="equil")
	}
	end_regen <- Sys.time() - start_regen

	### ----- run LIME -----------###
	lime_dirs <- equil_dir_vec[which(grepl("LBSPR",equil_dir_vec)==FALSE)]
	lbspr_dirs <- equil_dir_vec[which(grepl("LBSPR",equil_dir_vec))]
	lime_combos <- equil_modcombos[which(grepl("LBSPR", equil_modcombos[,"Data_avail"])==FALSE),]
	lbspr_combos <- equil_modcombos[which(grepl("LBSPR", equil_modcombos[,"Data_avail"])),]	

	start_run <- Sys.time()
	foreach(loop=1:length(lime_dirs), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=lime_combos[loop,"C_opt"], LFdist=1, param_adjust=c("SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(0.737,0.2,0.2,0.2),  fix_param=FALSE, newtonsteps=3)
	end_run <- Sys.time() - start_run

	## pick up iterations not run in parallel
	for(loop in 1:length(lime_dirs)){
		run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=lime_combos[loop,"C_opt"], LFdist=1, param_adjust=c("SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(0.737,0.2,0.2,0.2),  fix_param=FALSE, newtonsteps=3)
	}

	### ----- run LBSPR -----------###	

	start_run <- Sys.time()
	foreach(loop=1:length(lbspr_dirs), .packages=c('LBSPR','LIME')) %dopar% run_LBSPR(modpath=lbspr_dirs[loop], lh=lh_list[[as.character(strsplit(lbspr_combos[loop,"LH"],"_")[[1]][2])]], itervec=itervec, species=NULL, rewrite=FALSE, simulation=TRUE)
	end_altrun <- Sys.time() - start_run


############################################################
#### base - with variation
############################################################
### ----- models to run ----------- ###
lh_vec <- c("Short", "Medium", "Long")
Fdyn_vec <- c("Ramp", "Increasing")#, "Decreasing")
Rdyn_vec <- c("AR")
data_vec <- c("Index_Catch_LC20", "Index_LC10", "Index_LC1", "Catch_LC10", "Catch_LC1", "LC10", "LC1", "LBSPR10", "LBSPR1")
SampleSize_vec <- 200
itervec <- 1:100

base_modcombos <- expand.grid("Data_avail"=data_vec, "SampleSize"=paste0("SampleSize_", SampleSize_vec), "LH"=paste0("LH_",lh_vec), "Fdyn"=paste0("F_",Fdyn_vec), "Rdyn"=paste0("R_",Rdyn_vec), stringsAsFactors=FALSE)
base_modcombos$C_opt <- rep(0, nrow(base_modcombos))
	base_modcombos$C_opt[which(grepl("Catch",base_modcombos[,"Data_avail"]))] <- 2

### ----- base runs ----------- ###
	base_dir <- file.path(main_dir, "base")
	# unlink(base_dir, TRUE)
	dir.create(base_dir, showWarnings=FALSE)	

	## setup equil dirs
	base_dir_vec <- model_paths(res_dir=base_dir, modcombos=base_modcombos[,-c(ncol(base_modcombos))])
	ignore <- sapply(1:length(base_dir_vec), function(m) sapply(itervec, function(x) dir.create(file.path(base_dir_vec[m], x), showWarnings=FALSE)))	

	## setup scenario life history list
	lh_list <- adj_variation(SigmaR=0.737, SigmaF=0.2, SigmaC=0.2, SigmaI=0.2, CVlen=0.1, rho=0.426, selex_param=2, mat_param=1,  nseasons=1)

	### ----- generate data -----------###
	### data rich cases only
	rich_dir <- base_dir_vec[grepl("Index_Catch_LC20",base_dir_vec)]
	rich_modcombos <- base_modcombos[which(base_modcombos$Data_avail=="Index_Catch_LC20"),]
	alt_dir <- base_dir_vec[grepl("Index_Catch_LC20",base_dir_vec)==FALSE]
	alt_modcombos <- base_modcombos[which(base_modcombos$Data_avail!="Index_Catch_LC20"),]

	### ----- use parallel cores -----------###
	# registerDoParallel(cores=ncores)	

	start_gen <- Sys.time()
	foreach(loop=1:length(rich_dir), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=rich_dir[loop], data_avail=as.character(rich_modcombos[loop,"Data_avail"]), itervec=itervec, Fdynamics=as.character(strsplit(rich_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(rich_modcombos[loop,"Rdyn"],"_")[[1]][2]),  lh=lh_list[[as.character(strsplit(rich_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=20, Nyears_comp=20, comp_sample=as.numeric(strsplit(rich_modcombos[loop,"SampleSize"],"_")[[1]][2]), rewrite=FALSE, init_depl=c(0.05,0.95))
	end_gen <- Sys.time() - start_gen

	# copy data rich cases to data poor
	# when written using parallel cores, sometimes files not writing properly
	start_regen <- Sys.time()
	for(loop in 1:length(rich_dir)){
		copy_sim(fromdir=rich_dir[loop], fromcombos=rich_modcombos[loop,], todir=alt_dir, itervec=itervec, rewrite=FALSE, res_dir="base")
	}
	end_regen <- Sys.time() - start_regen


	### ----- run LIME -----------###
	lime_dirs <- base_dir_vec[which(grepl("LBSPR",base_dir_vec)==FALSE)]
	lbspr_dirs <- base_dir_vec[which(grepl("LBSPR",base_dir_vec))]
	lime_combos <- base_modcombos[which(grepl("LBSPR", base_modcombos[,"Data_avail"])==FALSE),]
	lbspr_combos <- base_modcombos[which(grepl("LBSPR", base_modcombos[,"Data_avail"])),]	

	start_run <- Sys.time()
	foreach(loop=1:length(lime_dirs), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=lime_combos[loop,"C_opt"], LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.737,0.2,0.2),  fix_param=FALSE)
	end_run <- Sys.time() - start_run

	for(loop in 1:length(lime_dirs)){
		run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=lime_combos[loop,"C_opt"], LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.737,0.2,0.2),  fix_param=FALSE)
	}


	### ----- run LBSPR -----------###	

	start_run <- Sys.time()
	foreach(loop=1:length(lbspr_dirs), .packages=c('LBSPR','LIME')) %dopar% run_LBSPR(modpath=lbspr_dirs[loop], lh=lh_list[[as.character(strsplit(lbspr_combos[loop,"LH"],"_")[[1]][2])]], itervec=itervec, species=NULL, rewrite=FALSE, simulation=TRUE)
	end_altrun <- Sys.time() - start_run	

dirs1 <- c(equil_dir_vec, base_dir_vec)

# bp1 <- bias_precision(dirs1,itervec=1:100)
# saveRDS(bp1, file.path(res_dir, "equil_var_bias_precision.rds"))

bp1 <- readRDS(file.path(res_dir, "equil_var_bias_precision.rds"))

icol <- rev(brewer.pal(3,"Purples"))
ccol <- rev(brewer.pal(3, "Oranges"))
lcol <- rev(brewer.pal(3, "Blues"))
rcol <- rev(brewer.pal(3, "Greens"))
col_vec <- c(gray(0.3), icol[1:2], ccol[1:2], lcol[1:2], rcol[1:2])


png(file.path(fig_dir, "RelativeError_n200.png"), height=16, width=42, res=200, units="in")
par(mfrow=c(3,3), mar=c(0,4.5,0,0), omi=c(1.2,1.2,1,1), mgp=c(3,2,0))
## equilibrium
i2 <- which(grepl("Short", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("/equil/", dirs1))
# i2 <- i1[-c(which(grepl("LC5/",dirs1[i1]) | grepl("LC2/",dirs1[i1])))]
dirs1[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1,1.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(bp1$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(bp1$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=2.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=2.5)
abline(h=0, lty=2)
beanplot(as.data.frame(bp1$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
axis(2,at=seq(-0.5,1.5,by=0.5), las=2, cex.axis=3.5)
mtext(side=3, "Short-lived", line=1, cex=3)
print.letter("(a)", cex=3)

i2 <- which(grepl("Medium", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("/equil/", dirs1))
# i2 <- i1[-c(which(grepl("LC5/",dirs1[i1]) | grepl("LC2/",dirs1[i1])))]
dirs1[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1.5,2.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(bp1$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(bp1$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(bp1$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
mtext(side=3, "Medium-lived", line=1, cex=3)
axis(2,at=seq(-1,2.5,by=1), las=2, cex.axis=3.5)
print.letter("(b)", cex=3)

i2 <- which(grepl("Long", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("/equil/", dirs1))
# i2 <- i1[-c(which(grepl("LC5/",dirs1[i1]) | grepl("LC2/",dirs1[i1])))]
dirs1[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-2,3.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(bp1$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(bp1$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(bp1$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
mtext(side=3, "Longer-lived", line=1, cex=3)
mtext(side=4, "Equilibrium", line=2, cex=3)
axis(2,at=seq(-1,3,by=1), las=2, cex.axis=3.5)
print.letter("(c)", cex=3)

## base
i2 <- which(grepl("Short", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("base/", dirs1) & grepl("Ramp", dirs1))
# i2 <- i1[-c(which(grepl("LC5/",dirs1[i1]) | grepl("LC2/",dirs1[i1])))]
dirs1[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1,1.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(bp1$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(bp1$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(bp1$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
axis(2,at=seq(-0.5,1.5,by=0.5), las=2, cex.axis=3.5)
# axis(1, at=1:9, labels=c("Rich", "Index\nLC10", "Index\nLC1", "Catch\nLC10", "Catch\nLC1", "LC10", "LC1", "LBSPR\n10", "LBSPR\n1"), cex.axis=2.8)
print.letter("(d)", cex=3)

i2 <- which(grepl("Medium", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("base/", dirs1) & grepl("Ramp", dirs1))
# i2 <- i1[-c(which(grepl("LC5/",dirs1[i1]) | grepl("LC2/",dirs1[i1])))]
dirs1[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1.5,2.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(bp1$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(bp1$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(bp1$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
# axis(1, at=1:9, labels=c("Rich", "Index\nLC10", "Index\nLC1", "Catch\nLC10", "Catch\nLC1", "LC10", "LC1", "LBSPR\n10", "LBSPR\n1"), cex.axis=2.8)
axis(2,at=seq(-1,2.5,by=1), las=2, cex.axis=3.5)
print.letter("(e)", cex=3)

i2 <- which(grepl("Long", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("base/", dirs1) & grepl("Ramp", dirs1))
# i2 <- i1[-c(which(grepl("LC5/",dirs1[i1]) | grepl("LC2/",dirs1[i1])))]
dirs1[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-2,3.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(bp1$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(bp1$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
# axis(1, at=1:9, labels=c("Rich", "Index\nLC10", "Index\nLC1", "Catch\nLC10", "Catch\nLC1", "LC10", "LC1", "LBSPR\n10", "LBSPR\n1"), cex.axis=2.8)
axis(2,at=seq(-1,3,by=1), las=2, cex.axis=3.5)
beanplot(as.data.frame(bp1$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
mtext(side=4, "Two-way trip", line=2, cex=3)
print.letter("(f)", cex=3)


## base
i2 <- which(grepl("Short", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("base/", dirs1) & grepl("Increasing", dirs1))
# i2 <- i1[-c(which(grepl("LC5/",dirs1[i1]) | grepl("LC2/",dirs1[i1])))]
dirs1[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1,1.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(bp1$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(bp1$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(bp1$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
axis(2,at=seq(-0.5,1.5,by=0.5), las=2, cex.axis=3.5)
axis(1, at=1:length(i2), labels=c("Rich", "Index\nLC10", "Index\nLC1", "Catch\nLC10", "Catch\nLC1", "LC10", "LC1", "LBSPR\n10", "LBSPR\n1"), cex.axis=3)
print.letter("(g)", cex=3)

i2 <- which(grepl("Medium", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("base/", dirs1) & grepl("Increasing", dirs1))
# i2 <- i1[-c(which(grepl("LC5/",dirs1[i1]) | grepl("LC2/",dirs1[i1])))]
dirs1[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1.5,2.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(bp1$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(bp1$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(bp1$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
axis(1, at=1:length(i2), labels=c("Rich", "Index\nLC10", "Index\nLC1", "Catch\nLC10", "Catch\nLC1", "LC10", "LC1", "LBSPR\n10", "LBSPR\n1"), cex.axis=3)
axis(2,at=seq(-1,2.5,by=1), las=2, cex.axis=3.5)
print.letter("(h)", cex=3)

i2 <- which(grepl("Long", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("base/", dirs1) & grepl("Increasing", dirs1))
# i2 <- i1[-c(which(grepl("LC5/",dirs1[i1]) | grepl("LC2/",dirs1[i1])))]
dirs1[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-2,3.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(bp1$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(bp1$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
axis(1, at=1:length(i2), labels=c("Rich", "Index\nLC10", "Index\nLC1", "Catch\nLC10", "Catch\nLC1", "LC10", "LC1", "LBSPR\n10", "LBSPR\n1"), cex.axis=3)
axis(2,at=seq(-1,3,by=1), las=2, cex.axis=3.5)
beanplot(as.data.frame(bp1$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
mtext(side=4, "One-way trip", line=2, cex=3)
mtext(side=2, outer=TRUE, line=3, "Relative error", cex=3)
mtext(side=1, outer=TRUE, line=6, "Data availability scenario", cex=3)
print.letter("(i)", cex=3)
dev.off()




# ic1 <- interval_coverage(dirs1, itervec=1:100)
# saveRDS(ic1, file.path(res_dir, "equil_var_interval_coverage.rds"))

ic1 <- readRDS(file.path(res_dir, "equil_var_interval_coverage.rds"))



png(file.path(fig_dir, "IntervalCoverage.png"), height=15, width=28, res=200, units="in")
par(mfrow=c(3,3), mar=c(0,0,0,0), omi=c(1.2,1.8,1,1), mgp=c(3,2,0))
## equilibrium
i2 <- which(grepl("Short", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("/equil/", dirs1) & grepl("LBSPR", dirs1)==FALSE)
	dirs1[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,7.5),ylim=c(0,1.2))
print.letter("(a)", cex=3)
legend("bottomleft", pch=c(19,17), col=c("#11111175", "#AA000075"), legend=c("Interval coverage", "Convergence rate"), cex=3)
# cov <- sapply(1:length(i2), function(x) round(sum(all_ic$cover[which(all_ic$converge[,i2[x]]==1),i2[x]])/length(which(all_ic$converge[,i2[x]]==1)),2))
cov <- colSums(ic1$cover[,i2])/nrow(ic1$cover[,i2])
conv <- colSums(ic1$converge[,i2])/nrow(ic1$converge[,i2])
abline(lty=2,h=0.5, col=gray(0.3))
abline(lty=2,h=1, col="#AA0000")
points(x=1:7, y=conv, col="#AA000075", xpd=NA, cex=5, pch=17, lwd=3)
points(x=1:7, y=cov, col="#11111175", xpd=NA, cex=5, pch=19, lwd=3)
axis(2,at=seq(0.25,1,by=0.25), las=2,cex.axis=3.5)
mtext(side=3, "Short-lived", line=1, cex=3)
i2 <- which(grepl("Medium", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("/equil/", dirs1) & grepl("LBSPR", dirs1)==FALSE)
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,7.5),ylim=c(0,1.2))
print.letter("(b)", cex=3)
# cov <- sapply(1:length(i2), function(x) round(sum(ic1$cover[which(ic1$converge[,i2[x]]==1),i2[x]])/length(which(ic1$converge[,i2[x]]==1)),2))
cov <- colSums(ic1$cover[,i2])/nrow(ic1$cover[,i2])
conv <- colSums(ic1$converge[,i2])/nrow(ic1$converge[,i2])
abline(lty=2,h=0.5, col=gray(0.3))
abline(lty=2,h=1, col="#AA0000")
points(x=1:7, y=conv, col="#AA000075", xpd=NA, cex=5, pch=17, lwd=3)
points(x=1:7, y=cov, col="#11111175", xpd=NA, cex=5, pch=19, lwd=3)
mtext(side=3, "Medium-lived", line=1, cex=3)
i2 <- which(grepl("Long", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("/equil/", dirs1) & grepl("LBSPR", dirs1)==FALSE)
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,7.5),ylim=c(0,1.2))
print.letter("(c)", cex=3)
# cov <- sapply(1:length(i2), function(x) round(sum(ic1$cover[which(ic1$converge[,i2[x]]==1),i2[x]])/length(which(ic1$converge[,i2[x]]==1)),2))
cov <- colSums(ic1$cover[,i2])/nrow(ic1$cover[,i2])
conv <- colSums(ic1$converge[,i2])/nrow(ic1$converge[,i2])
abline(lty=2,h=0.5, col=gray(0.3))
abline(lty=2,h=1, col="#AA0000")
points(x=1:7, y=conv, col="#AA000075", xpd=NA, cex=5, pch=17, lwd=3)
points(x=1:7, y=cov, col="#11111175", xpd=NA, cex=5, pch=19, lwd=3)
mtext(side=3, "Longer-lived", line=1, cex=3)
mtext(side=4, "Equilibrium", line=2, cex=3)

## low var
i2 <- which(grepl("Short", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("/base/", dirs1) & grepl("Ramp", dirs1) & grepl("LBSPR", dirs1)==FALSE)
	dirs1[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,7.5),ylim=c(0,1.2))
print.letter("(d)", cex=3)
# cov <- sapply(1:length(i2), function(x) round(sum(ic1$cover[which(ic1$converge[,i2[x]]==1),i2[x]])/length(which(ic1$converge[,i2[x]]==1)),2))
cov <- colSums(ic1$cover[,i2])/nrow(ic1$cover[,i2])
conv <- colSums(ic1$converge[,i2])/nrow(ic1$converge[,i2])
abline(lty=2,h=0.5, col=gray(0.3))
abline(lty=2,h=1, col="#AA0000")
points(x=1:7, y=conv, col="#AA000075", xpd=NA, cex=5, pch=17, lwd=3)
points(x=1:7, y=cov, col="#11111175", xpd=NA, cex=5, pch=19, lwd=3)
axis(2,at=seq(0.25,1,by=0.25), las=2,cex.axis=3.5)
i2 <- which(grepl("Medium", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("/base/", dirs1) & grepl("Ramp", dirs1) & grepl("LBSPR", dirs1)==FALSE)
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,7.5),ylim=c(0,1.2))
print.letter("(e)", cex=3)
# cov <- sapply(1:length(i2), function(x) round(sum(ic1$cover[which(ic1$converge[,i2[x]]==1),i2[x]])/length(which(ic1$converge[,i2[x]]==1)),2))
cov <- colSums(ic1$cover[,i2])/nrow(ic1$cover[,i2])
conv <- colSums(ic1$converge[,i2])/nrow(ic1$converge[,i2])
abline(lty=2,h=0.5, col=gray(0.3))
abline(lty=2,h=1, col="#AA0000")
points(x=1:7, y=conv, col="#AA000075", xpd=NA, cex=5, pch=17, lwd=3)
points(x=1:7, y=cov, col="#11111175", xpd=NA, cex=5, pch=19, lwd=3)
i2 <- which(grepl("Long", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("/base/", dirs1) & grepl("Ramp", dirs1) & grepl("LBSPR", dirs1)==FALSE)
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,7.5),ylim=c(0,1.2))
print.letter("(f)", cex=3)
# cov <- sapply(1:length(i2), function(x) round(sum(ic1$cover[which(ic1$converge[,i2[x]]==1),i2[x]])/length(which(ic1$converge[,i2[x]]==1)),2))
cov <- colSums(ic1$cover[,i2])/nrow(ic1$cover[,i2])
conv <- colSums(ic1$converge[,i2])/nrow(ic1$converge[,i2])
abline(lty=2,h=0.5, col=gray(0.3))
abline(lty=2,h=1, col="#AA0000")
points(x=1:7, y=conv, col="#AA000075", xpd=NA, cex=5, pch=17, lwd=3)
points(x=1:7, y=cov, col="#11111175", xpd=NA, cex=5, pch=19, lwd=3)
mtext(side=4, "Two-way", line=2, cex=3)

## base
i2 <- which(grepl("Short", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("/base/", dirs1) & grepl("Increasing", dirs1) & grepl("LBSPR", dirs1)==FALSE)
	dirs1[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,7.5),ylim=c(0,1.2))
print.letter("(g)", cex=3)
# cov <- sapply(1:length(i2), function(x) round(sum(ic1$cover[which(ic1$converge[,i2[x]]==1),i2[x]])/length(which(ic1$converge[,i2[x]]==1)),2))
cov <- colSums(ic1$cover[,i2])/nrow(ic1$cover[,i2])
conv <- colSums(ic1$converge[,i2])/nrow(ic1$converge[,i2])
abline(lty=2,h=0.5, col=gray(0.3))
abline(lty=2,h=1, col="#AA0000")
points(x=1:7, y=conv, col="#AA000075", xpd=NA, cex=5, pch=17, lwd=3)
points(x=1:7, y=cov, col="#11111175", xpd=NA, cex=5, pch=19, lwd=3)
axis(2,at=seq(0.25,1,by=0.25), las=2,cex.axis=3.5)
axis(1, at=1:7, labels=c("Rich", "Index\nLC10", "Index\nLC1", "Catch\nLC10", "Catch\nLC1", "LC10", "LC1"), cex.axis=3)
i2 <- which(grepl("Medium", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("/base/", dirs1) & grepl("Increasing", dirs1) & grepl("LBSPR", dirs1)==FALSE)
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,7.5),ylim=c(0,1.2))
print.letter("(h)", cex=3)
# cov <- sapply(1:length(i2), function(x) round(sum(ic1$cover[which(ic1$converge[,i2[x]]==1),i2[x]])/length(which(ic1$converge[,i2[x]]==1)),2))
cov <- colSums(ic1$cover[,i2])/nrow(ic1$cover[,i2])
conv <- colSums(ic1$converge[,i2])/nrow(ic1$converge[,i2])
abline(lty=2,h=0.5, col=gray(0.3))
abline(lty=2,h=1, col="#AA0000")
points(x=1:7, y=conv, col="#AA000075", xpd=NA, cex=5, pch=17, lwd=3)
points(x=1:7, y=cov, col="#11111175", xpd=NA, cex=5, pch=19, lwd=3)
axis(1, at=1:7, labels=c("Rich", "Index\nLC10", "Index\nLC1", "Catch\nLC10", "Catch\nLC1", "LC10", "LC1"), cex.axis=3)
i2 <- which(grepl("Long", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("/base/", dirs1) & grepl("Increasing", dirs1) & grepl("LBSPR", dirs1)==FALSE)
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,7.5),ylim=c(0,1.2))
print.letter("(i)", cex=3)
# cov <- sapply(1:length(i2), function(x) round(sum(ic1$cover[which(ic1$converge[,i2[x]]==1),i2[x]])/length(which(ic1$converge[,i2[x]]==1)),2))
cov <- colSums(ic1$cover[,i2])/nrow(ic1$cover[,i2])
conv <- colSums(ic1$converge[,i2])/nrow(ic1$converge[,i2])
abline(lty=2,h=0.5, col=gray(0.3))
abline(lty=2,h=1, col="#AA0000")
points(x=1:7, y=conv, col="#AA000075", xpd=NA, cex=5, pch=17, lwd=3)
points(x=1:7, y=cov, col="#11111175", xpd=NA, cex=5, pch=19, lwd=3)
axis(1, at=1:7, labels=c("Rich", "Index\nLC10", "Index\nLC1", "Catch\nLC10", "Catch\nLC1", "LC10", "LC1"), cex.axis=3)
mtext(side=4, "One-way", line=2, cex=3)
mtext(side=2, outer=TRUE, line=8.5, "Proportion of iterations", cex=3)
mtext(side=1, outer=TRUE, line=6, "Data availability scenario", cex=3)
dev.off()









########### choose only length comp


lcol <- rev(brewer.pal(3, "Blues"))
rcol <- rev(brewer.pal(3, "Greens"))
col_vec <- c(lcol[1:2], rcol[1:2])


png(file.path(fig_dir, "LCequil_RelativeError_n200.png"), height=16, width=30, res=200, units="in")
par(mfrow=c(1,3), mar=c(0,4.5,0,0), omi=c(1.2,1.2,1,1), mgp=c(3,2,0))
## equilibrium
i2 <- which(grepl("Short", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("/equil/", dirs1) & grepl("Index", dirs1)==FALSE & grepl("Catch", dirs1)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",dirs1[i1]) | grepl("LC2/",dirs1[i1])))]
dirs1[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1,1.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(bp1$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(bp1$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=2.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=2.5)
abline(h=0, lty=2)
beanplot(as.data.frame(bp1$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
axis(2,at=seq(-0.5,1.5,by=0.5), las=2, cex.axis=4)
mtext(side=3, "Short-lived", line=1, cex=3)
axis(1, at=1:length(i2), labels=c("LC10", "LC1", "LBSPR\n10", "LBSPR\n1"), cex.axis=4)
# print.letter("(a)", cex=3)

i2 <- which(grepl("Medium", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("/equil/", dirs1) & grepl("Index", dirs1)==FALSE & grepl("Catch", dirs1)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",dirs1[i1]) | grepl("LC2/",dirs1[i1])))]
dirs1[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1.5,2.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(bp1$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(bp1$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(bp1$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
mtext(side=3, "Medium-lived", line=1, cex=3)
axis(2,at=seq(-1,2.5,by=1), las=2, cex.axis=4)
axis(1, at=1:length(i2), labels=c("LC10", "LC1", "LBSPR\n10", "LBSPR\n1"), cex.axis=4)
# print.letter("(b)", cex=3)

i2 <- which(grepl("Long", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("/equil/", dirs1) & grepl("Index", dirs1)==FALSE & grepl("Catch", dirs1)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",dirs1[i1]) | grepl("LC2/",dirs1[i1])))]
dirs1[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-2,3.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(bp1$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(bp1$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(bp1$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
mtext(side=3, "Longer-lived", line=1, cex=3)
# mtext(side=4, "Equilibrium", line=2, cex=3)
axis(2,at=seq(-1,3,by=1), las=2, cex.axis=4)
axis(1, at=1:length(i2), labels=c("LC10", "LC1", "LBSPR\n10", "LBSPR\n1"), cex.axis=4)
# print.letter("(c)", cex=3)
dev.off()
