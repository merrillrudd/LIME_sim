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
## LBSPR Operating model
############################################################
### ----- models to run ----------- ###
lh_vec <- c("Short", "Medium", "Long")
data_vec <- c("LC10", "LC1", "LBSPR10", "LBSPR1") #"Index_LC10", "Index_LC1", "Catch_LC10", "Catch_LC1", "LC10", "LC1", "LBSPR10", "LBSPR1")
itervec <- 1:100

LBSPR_modcombos <- expand.grid("Data_avail"=data_vec, "LH"=paste0("LH_",lh_vec), stringsAsFactors=FALSE)

	equil_LBSPR_dir <- file.path(main_dir, "equil_LBSPR")
	# unlink(equil_LBSPR_dir, TRUE)
	dir.create(equil_LBSPR_dir, showWarnings=FALSE)	

	## setup equilibrium dirs
	equil_LBSPR_dir_vec <- model_paths(res_dir=equil_LBSPR_dir, modcombos=LBSPR_modcombos)
	ignore <- sapply(1:length(equil_LBSPR_dir_vec), function(m) sapply(itervec, function(x) dir.create(file.path(equil_LBSPR_dir_vec[m], x), showWarnings=FALSE)))	

	## setup scenario life history list
	## 1-parameter model for selectivity and maturity
	lh_list <- adj_variation(SigmaR=0.01, SigmaF=0.01, SigmaC=0.01, SigmaI=0.01, CVlen=0.1, rho=0, selex_param=2, mat_param=1, nseasons=1)
	# lh_fig(lh_list, save=TRUE)

	### ----- use parallel cores -----------###
	for(loop in 1:length(equil_LBSPR_dir_vec)){
		left <- ifelse(grepl("LBSPR", LBSPR_modcombos[loop,"Data_avail"]), "LBSPR", "LC")
		ignore <- sim_LBSPR(dir=equil_LBSPR_dir_vec[loop], lh=lh_list[[as.character(strsplit(LBSPR_modcombos[loop,"LH"],"_")[[1]][2])]], itervec=itervec, Nyears=as.numeric(strsplit(LBSPR_modcombos[loop,"Data_avail"],left)[[1]][2]), sample_size=200)
	}

	### ----- run LIME -----------###
	lime_dirs <- equil_LBSPR_dir_vec[which(grepl("/LBSPR",equil_LBSPR_dir_vec)==FALSE)]
	lbspr_dirs <- equil_LBSPR_dir_vec[which(grepl("/LBSPR",equil_LBSPR_dir_vec))]
	lime_combos <- LBSPR_modcombos[which(grepl("LBSPR", LBSPR_modcombos[,"Data_avail"])==FALSE),]
	lbspr_combos <- LBSPR_modcombos[which(grepl("LBSPR", LBSPR_modcombos[,"Data_avail"])),]	

		registerDoParallel(cores=ncores)	


	start_run <- Sys.time()
	foreach(loop=1:length(lime_dirs), .packages=c('TMB','LIME')) %dopar%
	   run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(0.737,0.2,0.2,0.2),  fix_param=FALSE)
	end_run <- Sys.time() - start_run

	### ----- run LBSPR -----------###	

	start_run <- Sys.time()
	foreach(loop=1:length(lbspr_dirs), .packages=c('LBSPR','LIME')) %dopar% run_LBSPR(modpath=lbspr_dirs[loop], lh=lh_list[[as.character(strsplit(lbspr_combos[loop,"LH"],"_")[[1]][2])]], itervec=itervec, species=NULL, rewrite=FALSE, simulation=TRUE)
	end_altrun <- Sys.time() - start_run


############################################################
## monthly simulation, annual run
############################################################

	### ----- models to run ----------- ###
lh_vec <- c("Short", "Medium", "Long")
Fdyn_vec <- "Constant"
Rdyn_vec <- "Constant"
data_vec <- c("LC10", "LC1", "LBSPR10", "LBSPR1")
SampleSize_vec <- 200
itervec <- 1:100

equil_modcombos <- expand.grid("Data_avail"=data_vec, "SampleSize"=paste0("SampleSize_", SampleSize_vec), "LH"=paste0("LH_",lh_vec), "Fdyn"=paste0("F_",Fdyn_vec), "Rdyn"=paste0("R_",Rdyn_vec), stringsAsFactors=FALSE)
equil_modcombos$C_opt <- rep(0, nrow(equil_modcombos))
	equil_modcombos$C_opt[which(grepl("Catch",equil_modcombos[,"Data_avail"]))] <- 2

	equil_MoYr_dir <- file.path(main_dir, "equil_MoYr")
	# unlink(equil_MoYr_dir, TRUE)
	dir.create(equil_MoYr_dir, showWarnings=FALSE)	

	## setup equilibrium dirs
	equil_MoYr_dir_vec <- model_paths(res_dir=equil_MoYr_dir, modcombos=equil_modcombos[,-c(ncol(equil_modcombos))])
	ignore <- sapply(1:length(equil_MoYr_dir_vec), function(m) sapply(itervec, function(x) dir.create(file.path(equil_MoYr_dir_vec[m], x), showWarnings=FALSE)))	

	## setup scenario life history list
	## 1-parameter model for selectivity and maturity
	lh_list <- adj_variation(SigmaR=0.01, SigmaF=0.01, SigmaC=0.01, SigmaI=0.01, CVlen=0.1, rho=0, selex_param=2, mat_param=1, nseasons=12)
	# lh_fig(lh_list, save=TRUE)

	### ----- generate data -----------###
	### data rich cases only
	rich_dir <- equil_MoYr_dir_vec[grepl(data_vec[1],equil_MoYr_dir_vec)]
	rich_modcombos <- equil_modcombos[which(equil_modcombos$Data_avail==data_vec[1]),]
	alt_dir <- equil_MoYr_dir_vec[grepl(data_vec[1],equil_MoYr_dir_vec)==FALSE]
	alt_modcombos <- equil_modcombos[which(equil_modcombos$Data_avail!=data_vec[1]),]

	### ----- use parallel cores -----------###
	registerDoParallel(cores=ncores)	

	## monthly data collection but pool=TRUE
	start_gen <- Sys.time()
	foreach(loop=1:length(rich_dir), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=rich_dir[loop], data_avail=as.character(rich_modcombos[loop,"Data_avail"]), itervec=itervec, Fdynamics=as.character(strsplit(rich_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(rich_modcombos[loop,"Rdyn"],"_")[[1]][2]),  lh=lh_list[[as.character(strsplit(rich_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=20, Nyears_comp=20, comp_sample=as.numeric(strsplit(rich_modcombos[loop,"SampleSize"],"_")[[1]][2]), rewrite=FALSE, init_depl=c(0.05,0.95), pool=TRUE)
	end_gen <- Sys.time() - start_gen

	# copy data rich cases to data poor
	# when written using parallel cores, sometimes files not writing properly
	start_regen <- Sys.time()
	for(loop in 1:length(rich_dir)){
		copy_sim(fromdir=rich_dir[loop], fromcombos=rich_modcombos[loop,], todir=alt_dir, itervec=itervec, rewrite=FALSE, res_dir="equil_MoYr")
	}
	end_regen <- Sys.time() - start_regen

	### ----- run LIME -----------###
	lime_dirs <- equil_MoYr_dir_vec[which(grepl("LBSPR",equil_MoYr_dir_vec)==FALSE)]
	lbspr_dirs <- equil_MoYr_dir_vec[which(grepl("LBSPR",equil_MoYr_dir_vec))]
	lime_combos <- equil_modcombos[which(grepl("LBSPR", equil_modcombos[,"Data_avail"])==FALSE),]
	lbspr_combos <- equil_modcombos[which(grepl("LBSPR", equil_modcombos[,"Data_avail"])),]	

	lh_list <- adj_variation(SigmaR=0.01, SigmaF=0.01, SigmaC=0.01, SigmaI=0.01, CVlen=0.1, rho=0, selex_param=2, mat_param=1, nseasons=1)


	start_run <- Sys.time()
	foreach(loop=1:length(lime_dirs), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=lime_combos[loop,"C_opt"], LFdist=1, param_adjust=c("SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(0.737,0.2,0.2,0.2),  fix_param=FALSE)
	end_run <- Sys.time() - start_run


	### ----- run LBSPR -----------###	

	start_run <- Sys.time()
	foreach(loop=1:length(lbspr_dirs), .packages=c('LBSPR','LIME')) %dopar% run_LBSPR(modpath=lbspr_dirs[loop], lh=lh_list[[as.character(strsplit(lbspr_combos[loop,"LH"],"_")[[1]][2])]], itervec=itervec, species=NULL, rewrite=FALSE, simulation=TRUE)
	end_altrun <- Sys.time() - start_run



############################################################
## monthly simulation, monthly run
############################################################
### ----- equilibrium test ----------- ###
### ----- models to run ----------- ###
lh_vec <- c("Short")#, "Medium", "Long")
Fdyn_vec <- "Constant"
Rdyn_vec <- "Constant"
data_vec <- c("LC10", "LC1")#, "LBSPR10", "LBSPR1")
SampleSize_vec <- 200 #c(1000,200)#,20)
itervec <- 1:100

equil_modcombos <- expand.grid("Data_avail"=data_vec, "SampleSize"=paste0("SampleSize_", SampleSize_vec), "LH"=paste0("LH_",lh_vec), "Fdyn"=paste0("F_",Fdyn_vec), "Rdyn"=paste0("R_",Rdyn_vec), stringsAsFactors=FALSE)
equil_modcombos$C_opt <- rep(0, nrow(equil_modcombos))
	equil_modcombos$C_opt[which(grepl("Catch",equil_modcombos[,"Data_avail"]))] <- 2

	equil_runMonthly_dir <- file.path(main_dir, "equil_runMonthly")
	# unlink(equil_runMonthly_dir, TRUE)
	dir.create(equil_runMonthly_dir, showWarnings=FALSE)	

	## setup equilibrium dirs
	equil_runMonthly_dir_vec <- model_paths(res_dir=equil_runMonthly_dir, modcombos=equil_modcombos[,-c(ncol(equil_modcombos))])
	ignore <- sapply(1:length(equil_runMonthly_dir_vec), function(m) sapply(itervec, function(x) dir.create(file.path(equil_runMonthly_dir_vec[m], x), showWarnings=FALSE)))	

	## setup scenario life history list
	## 1-parameter model for selectivity and maturity
	lh_list <- adj_variation(SigmaR=0.01, SigmaF=0.01, SigmaC=0.01, SigmaI=0.01, CVlen=0.1, rho=0, selex_param=2, mat_param=1, nseasons=12)
	# lh_fig(lh_list, save=TRUE)

	### ----- generate data -----------###

	### ----- use parallel cores -----------###
	registerDoParallel(cores=ncores)	

	## key difference between other sims --- pool=FALSE (monthly length comp, not pooled annually)
	start_gen <- Sys.time()
	foreach(loop=1:length(equil_runMonthly_dir_vec), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=equil_runMonthly_dir_vec[loop], data_avail=as.character(equil_modcombos[loop,"Data_avail"]), itervec=itervec, Fdynamics=as.character(strsplit(equil_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(equil_modcombos[loop,"Rdyn"],"_")[[1]][2]),  lh=lh_list[[as.character(strsplit(equil_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=as.numeric(strsplit(equil_modcombos[loop,"Data_avail"],"LC")[[1]][2]), Nyears_comp=as.numeric(strsplit(equil_modcombos[loop,"Data_avail"],"LC")[[1]][2]), comp_sample=as.numeric(strsplit(equil_modcombos[loop,"SampleSize"],"_")[[1]][2]), rewrite=FALSE, init_depl=c(0.05,0.95), pool=FALSE)
	end_gen <- Sys.time() - start_gen


	### ----- run LIME -----------###
	lime_dirs <- equil_runMonthly_dir_vec[which(grepl("LBSPR",equil_runMonthly_dir_vec)==FALSE)]
	lime_combos <- equil_modcombos[which(grepl("LBSPR", equil_modcombos[,"Data_avail"])==FALSE),]

	start_run <- Sys.time()
	foreach(loop=1:length(lime_dirs), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=lime_combos[loop,"C_opt"], LFdist=1, param_adjust=c("SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(0.737,0.2,0.2,0.2),  fix_param=FALSE)
	end_run <- Sys.time() - start_run

all_dirs <- c(equil_LBSPR_dir_vec, equil_MoYr_dir_vec, equil_runMonthly_dir_vec)

# bp <- bias_precision(all_dirs,itervec=1:100)
# saveRDS(bp, file.path(res_dir, "explore_monthly_bias_precision.rds"))

bp <- readRDS(file.path(res_dir, "explore_monthly_bias_precision.rds"))

### equilibrium scenarios - LBSPR
dirnames <- sapply(1:length(all_dirs), function(x) strsplit(all_dirs[x], "equil_")[[1]][2])
inst_dirs <- which(grepl("LBSPR/", dirnames)==FALSE)
dirnames[inst_dirs] <- sapply(1:length(inst_dirs), function(x) strsplit(dirnames[inst_dirs[x]], "/F_Constant")[[1]][1])
bmat2 <- pmat2 <- rep(NA, length(dirnames))
names(bmat2) <- names(pmat2) <- dirnames

for(i in 1:length(dirnames)){
	b <- sapply(1:length(dirnames), function(x) median((bp$relerr[,i]), na.rm=TRUE))
	p <- sapply(1:length(dirnames), function(x) median(abs(bp$relerr[,i]), na.rm=TRUE))
	bmat2[i] <- b
	pmat2[i] <- p
	rm(b)
	rm(p)
}

write.csv(rbind(bmat2, pmat2), file.path(res_dir, "explore_monthly_bias_precision.csv"))

## calculations for paper
moyrd <- all_dirs[grepl("MoYr",all_dirs)]
mod <- all_dirs[grepl("runMonthly",all_dirs)]

idir <- which(grepl("runMonthly",all_dirs) & grepl("LC1/",all_dirs))
all_dirs[idir]
median(bp$relerr[,idir], na.rm=TRUE)

idir <- which(grepl("runMonthly",all_dirs) & grepl("LC10/",all_dirs))
all_dirs[idir]
median(bp$relerr[,idir], na.rm=TRUE)

idir <- which(grepl("MoYr",all_dirs) & grepl("LC1/",all_dirs) & grepl("Short",all_dirs))
all_dirs[idir]
median(bp$relerr[,idir], na.rm=TRUE)

idir <- which(grepl("MoYr",all_dirs) & grepl("LC10/",all_dirs) & grepl("Short",all_dirs))
all_dirs[idir]
median(bp$relerr[,idir], na.rm=TRUE)
## create figure for paper

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


dirs1 <- c(equil_dir_vec, base_dir_vec)

bp1 <- readRDS(file.path(res_dir, "equil_var_bias_precision.rds"))



lcol <- rev(brewer.pal(3, "Blues"))
rcol <- rev(brewer.pal(3, "Greens"))
col_vec <- c(lcol[1:2], rcol[1:2])




png(file.path(fig_dir, "EquilCompare_RelativeError_n200.png"), height=16, width=35, res=200, units="in")
par(mfrow=c(3,3), mar=c(0,7,0,0), omi=c(1.6,1.2,1,1.5), mgp=c(3,2,0))
## equilibrium
i2 <- which(grepl("Short", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("/equil/", dirs1) & grepl("Index", dirs1)==FALSE & grepl("Catch", dirs1)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",dirs1[i1]) | grepl("LC2/",dirs1[i1])))]
dirs1[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1,1))
b1 <- sapply(1:length(i2), function(x) round(median(abs(bp1$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(bp1$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=2.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=2.5)
abline(h=0, lty=2)
beanplot(as.data.frame(bp1$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
axis(2,at=seq(-0.5,1,by=0.5), las=2,cex.axis=3.5)
mtext(side=3, "Short-lived", line=1, cex=3)
print.letter(xy=c(0.05,0.95), "(a)", cex=3)

i2 <- which(grepl("Medium", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("/equil/", dirs1) & grepl("Index", dirs1)==FALSE & grepl("Catch", dirs1)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",dirs1[i1]) | grepl("LC2/",dirs1[i1])))]
dirs1[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1,1.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(bp1$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(bp1$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(bp1$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
mtext(side=3, "Medium-lived", line=1, cex=3)
axis(2,at=seq(-0.5,1.5,by=0.5), las=2,cex.axis=3.5)
print.letter(xy=c(0.05,0.95), "(b)", cex=3)

i2 <- which(grepl("Long", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("/equil/", dirs1) & grepl("Index", dirs1)==FALSE & grepl("Catch", dirs1)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",dirs1[i1]) | grepl("LC2/",dirs1[i1])))]
dirs1[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1.5,2.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(bp1$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(bp1$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(bp1$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
mtext(side=3, "Longer-lived", line=1, cex=3)
mtext(side=4, "Annual data,\nannual model", line=6.5, cex=3)
axis(2,at=seq(-1,2.5,by=0.5), las=2,cex.axis=3.5)
print.letter(xy=c(0.05,0.95), "(c)", cex=3)

## base
i2 <- which(grepl("Short", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("/equil_MoYr/", all_dirs) & grepl("Index", all_dirs)==FALSE & grepl("Catch", all_dirs)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1,4))
b1 <- sapply(1:length(i2), function(x) round(median(abs(bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
axis(2,at=seq(-1,4,by=1), las=2,cex.axis=3.5)
# axis(1, at=1:9, labels=c("Rich", "Index\nLC10", "Index\nLC1", "Catch\nLC10", "Catch\nLC1", "LC10", "LC1", "LBSPR\n10", "LBSPR\n1"), cex.axis=2.8)
print.letter(xy=c(0.05,0.95), "(d)", cex=3)

i2 <- which(grepl("Medium", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("/equil_MoYr/", all_dirs) & grepl("Index", all_dirs)==FALSE & grepl("Catch", all_dirs)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1.5,2))
b1 <- sapply(1:length(i2), function(x) round(median(abs(bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
# axis(1, at=1:9, labels=c("Rich", "Index\nLC10", "Index\nLC1", "Catch\nLC10", "Catch\nLC1", "LC10", "LC1", "LBSPR\n10", "LBSPR\n1"), cex.axis=2.8)
axis(2,at=seq(-1,2,by=0.5), las=2,cex.axis=3.5)
print.letter(xy=c(0.05,0.95), "(e)", cex=3)

i2 <- which(grepl("Long", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("/equil_MoYr/", all_dirs) & grepl("Index", all_dirs)==FALSE & grepl("Catch", all_dirs)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1.5,3.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
# axis(1, at=1:9, labels=c("Rich", "Index\nLC10", "Index\nLC1", "Catch\nLC10", "Catch\nLC1", "LC10", "LC1", "LBSPR\n10", "LBSPR\n1"), cex.axis=2.8)
axis(2,at=seq(-1,3.5,by=1), las=2,cex.axis=3.5)
beanplot(as.data.frame(bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
mtext(side=4, "Monthly data,\nannual model", line=6.5, cex=3)
print.letter(xy=c(0.05,0.95), "(f)", cex=3)


## base
i2 <- which(grepl("Short", all_dirs) & grepl("/equil_LBSPR/", all_dirs) & grepl("Index", all_dirs)==FALSE & grepl("Catch", all_dirs)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1,2))
b1 <- sapply(1:length(i2), function(x) round(median(abs(bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
axis(2,at=seq(-1,1.5,by=0.5), las=2,cex.axis=3.5)
axis(1, at=1:length(i2), labels=c("LC10","LC1", "LBSPR10", "LBSPR1"), cex.axis=3.5)
print.letter(xy=c(0.05,0.95), "(g)", cex=3)

i2 <- which(grepl("Medium", all_dirs) & grepl("/equil_LBSPR/", all_dirs) & grepl("Index", all_dirs)==FALSE & grepl("Catch", all_dirs)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1,1.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
axis(1, at=1:length(i2), labels=c("LC10","LC1", "LBSPR10", "LBSPR1"), cex.axis=3.5)
axis(2,at=seq(-1,1.5,by=0.5), las=2,cex.axis=3.5)
print.letter(xy=c(0.05,0.95), "(h)", cex=3)

i2 <- which(grepl("Long", all_dirs) & grepl("/equil_LBSPR/", all_dirs) & grepl("Index", all_dirs)==FALSE & grepl("Catch", all_dirs)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1.5,2))
b1 <- sapply(1:length(i2), function(x) round(median(abs(bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
axis(1, at=1:length(i2), labels=c("LC10","LC1", "LBSPR10", "LBSPR1"), cex.axis=3.5)
axis(2,at=seq(-1.5,2,by=0.5), las=2,cex.axis=3.5)
beanplot(as.data.frame(bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
mtext(side=4, "LB-SPR OM", line=2, cex=3)
mtext(side=2, outer=TRUE, line=3, "Relative error", cex=3)
mtext(side=1, outer=TRUE, line=8, "Model/Data availability scenario", cex=3)
print.letter(xy=c(0.05,0.95), "(i)", cex=3)
dev.off()



### appendix figure
## size and abundance at age for yearly vs. monthly time step
lh1 <- adj_variation(SigmaR=0.0001, SigmaF=0.0001, SigmaC=0.01, SigmaI=0.01, CVlen=0.1, rho=0, selex_param=2, mat_param=1, nseasons=1)[["Short"]]
lh2 <- adj_variation(SigmaR=0.0001, SigmaF=0.0001, SigmaC=0.01, SigmaI=0.01, CVlen=0.1, rho=0, selex_param=2, mat_param=1, nseasons=12)[["Short"]]

x1 <- generate_data(modpath=NULL, data_avail="LC", itervec=1, Fdynamics="Constant", Rdynamics="Constant",  lh=lh1, Nyears=20, Nyears_comp=20, comp_sample=200, init_depl=0.8, pool=TRUE)
x2 <- generate_data(modpath=NULL, data_avail="LC", itervec=1, Fdynamics="Constant", Rdynamics="Constant",  lh=lh2, Nyears=20, Nyears_comp=20, comp_sample=200, init_depl=0.8, pool=FALSE)

p1 <- t(x1$plba)
p2 <- t(x2$plba)
n1 <- x1$N_at
n2 <- x2$N_at
ages <- x1$ages
months <- x2$ages

orange_fun <- colorRampPalette(c("white", "orangered"))
oranges <- orange_fun(12)

png(file.path(fig_dir, "SizeAtAge.png"), height=10, width=8, res=200, units="in")
par(mfrow=c(5,1), mar=c(0,0,0,0), omi=c(1,1,1,1))
for(i in 1:length(ages)){
	if(i!=length(ages)) index <- which(months >= ages[i] & months < ages[i+1])
	if(i==length(ages)) index <- which(months >= ages[i])
	if(i==1){
		ylim=c(0,1)
		int=0.2
	}
	if(i!=1){
		ylim=c(0,0.2)
		int= 0.05
	}
	plot(p1[,i], col="blue", lwd=8, type="l", lty=1, ylim=ylim, xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")
	par(new=TRUE)
	matplot(p2[,index], col=oranges, lwd=2, type="l", lty=1, ylim=ylim, xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")
	axis(2, las=2, cex.axis=2, at=seq(0,(as.numeric(ylim[2])-int),by=int))
	print.letter(xy=c(0.9,0.9), paste0("Age = ", ages[i]), cex=2)
}
axis(1, cex.axis=2)
mtext("Length (cm)", cex=2, side=1, line=3.5)
mtext("Probability of being length given age", cex=2, side=2, line=5, outer=TRUE)
dev.off()

## old code -- need to rename
#### comparing equilibrium to monthly
png(file.path(fig_dir, "Compare_annual_monthly.png"), res=200, units="in", height=15, width=28)
par(mfrow=c(2,3), mar=c(0,4,2,2), omi=c(1,1,1,1))
i2 <- c(which(grepl("Short", dirs1) & grepl("SampleSize_200", dirs1) & grepl("Index", dirs1)==FALSE & grepl("Catch", dirs1)==FALSE & grepl("/equil/", dirs1)))
dirs1[i2]
b1 <- sapply(1:length(i2), function(x) median(bp1$relerr[,i2[x]], na.rm=TRUE))
p1 <- sapply(1:length(i2), function(x) median(abs(bp1$relerr[,i2[x]]), na.rm=TRUE))
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1,1))
abline(h=0, lty=2)
# axis(1, at=1:length(i2), labels=c("LC10", "LC1", "LBSPR10", "LBSPR1"), cex.axis=3.5)
axis(2,at=seq(-1,1,by=0.5), las=2,cex.axis=3.5)
beanplot(as.data.frame(bp1$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
mtext(side=3, "Annual data and model", line=1, cex=3)

i2 <- c(which(grepl("Short", all_dirs) & grepl("SampleSize_200", all_dirs) & grepl("Index", all_dirs)==FALSE & grepl("Catch", all_dirs)==FALSE & grepl("/equil_MoYr/", all_dirs)))
all_dirs[i2]
b1 <- sapply(1:length(i2), function(x) median(all_bp$relerr[,i2[x]], na.rm=TRUE))
p1 <- sapply(1:length(i2), function(x) median(abs(all_bp$relerr[,i2[x]]), na.rm=TRUE))
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1,4))
abline(h=0, lty=2)
# axis(1, at=1:length(i2), labels=c("LC10", "LC1", "LBSPR10", "LBSPR1"), cex.axis=3.5)
axis(2,at=seq(-1,4,by=1), las=2,cex.axis=3.5)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
mtext(side=3, "Monthly data, annual model", line=1, cex=3)


i2 <- c(which(grepl("Short", all_dirs) & grepl("SampleSize_200", all_dirs) & grepl("Index", all_dirs)==FALSE & grepl("Catch", all_dirs)==FALSE & grepl("equil", all_dirs) & grepl("runMonthly", all_dirs)))
all_dirs[i2]
b2 <- sapply(1:length(i2), function(x) median(all_bp$relerr[,i2[x]], na.rm=TRUE))
p2 <- sapply(1:length(i2), function(x) median(abs(all_bp$relerr[,i2[x]]), na.rm=TRUE))
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-2,2.5))
abline(h=0, lty=2)
# axis(1, at=1:length(i2), labels=c("LC10", "LC1", "LBSPR10", "LBSPR1"), cex.axis=3.5)
axis(2,at=seq(-2,2.5,by=1), las=2,cex.axis=3.5)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
mtext(side=3, "Monthly data and model", line=1, cex=3)
mtext(side=4, "Equilibrium", line=2, cex=3)


i2 <- c(which(grepl("Short", all_dirs) & grepl("SampleSize_200", all_dirs) & grepl("Index", all_dirs)==FALSE & grepl("Catch", all_dirs)==FALSE & grepl("Ramp", all_dirs) & grepl("/base/", all_dirs)))
all_dirs[i2]
b3 <- sapply(1:length(i2), function(x) median(all_bp$relerr[,i2[x]], na.rm=TRUE))
p3 <- sapply(1:length(i2), function(x) median(abs(all_bp$relerr[,i2[x]]), na.rm=TRUE))
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1,1))
abline(h=0, lty=2)
axis(1, at=1:length(i2), labels=c("LC10", "LC1", "LBSPR10", "LBSPR1"), cex.axis=3.5)
axis(2,at=seq(-1,1,by=0.5), las=2,cex.axis=3.5)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)

i2 <- c(which(grepl("Short", all_dirs) & grepl("SampleSize_200", all_dirs) & grepl("Index", all_dirs)==FALSE & grepl("Catch", all_dirs)==FALSE & grepl("/base_MoYr/", all_dirs)))
all_dirs[i2]
b1 <- sapply(1:length(i2), function(x) median(all_bp$relerr[,i2[x]], na.rm=TRUE))
p1 <- sapply(1:length(i2), function(x) median(abs(all_bp$relerr[,i2[x]]), na.rm=TRUE))
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1,4))
abline(h=0, lty=2)
axis(1, at=1:length(i2), labels=c("LC10", "LC1", "LBSPR10", "LBSPR1"), cex.axis=3.5)
axis(2,at=seq(-1,4,by=1), las=2,cex.axis=3.5)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)


i2 <- c(which(grepl("Short", all_dirs) & grepl("SampleSize_200", all_dirs) & grepl("Index", all_dirs)==FALSE & grepl("Catch", all_dirs)==FALSE & grepl("Ramp", all_dirs) & grepl("runMonthly", all_dirs)))
all_dirs[i2]
b4 <- sapply(1:length(i2), function(x) median(all_bp$relerr[,i2[x]], na.rm=TRUE))
p4 <- sapply(1:length(i2), function(x) median(abs(all_bp$relerr[,i2[x]]), na.rm=TRUE))
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-2,2.5))
abline(h=0, lty=2)
axis(1, at=1:length(i2), labels=c("LC10", "LC1", "LBSPR10", "LBSPR1"), cex.axis=3.5)
axis(2,at=seq(-2,2.5,by=1), las=2,cex.axis=3.5)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
mtext(side=4, "Variability", line=2, cex=3)
mtext("Relative error", side=2, outer=TRUE, line=3.5, cex=3)
mtext("Model type and data availability", side=1, outer=TRUE, line=5, cex=3)
dev.off()
