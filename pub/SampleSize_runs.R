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
#### SAMPLE SIZE equilibrium - instantaneous LF
############################################################
### ----- models to run ----------- ###
lh_vec <- c("Short", "Medium", "Long")
Fdyn_vec <- "Constant"
Rdyn_vec <- "Constant"
data_vec <- c("LC10", "LC1")
SampleSize_vec <- c(1000, 500, 200, 100, 50, 20)
itervec <- 1:100

equil_N_modcombos <- expand.grid("Data_avail"=data_vec, "SampleSize"=paste0("SampleSize_", SampleSize_vec), "LH"=paste0("LH_",lh_vec), "Fdyn"=paste0("F_",Fdyn_vec), "Rdyn"=paste0("R_",Rdyn_vec), stringsAsFactors=FALSE)
equil_N_modcombos$C_opt <- rep(0, nrow(equil_N_modcombos))
	equil_N_modcombos$C_opt[which(grepl("Catch",equil_N_modcombos[,"Data_avail"]))] <- 2

### ----- equilibrium test ----------- ###

	equil_N_dir <- file.path(main_dir, "equil_samplesize")
	# unlink(equil_N_dir, TRUE)
	dir.create(equil_N_dir, showWarnings=FALSE)	

	## setup equilibrium dirs
	equil_N_dir_vec <- model_paths(res_dir=equil_N_dir, modcombos=equil_N_modcombos[,-c(ncol(equil_N_modcombos))])
	ignore <- sapply(1:length(equil_N_dir_vec), function(m) sapply(itervec, function(x) dir.create(file.path(equil_N_dir_vec[m], x), showWarnings=FALSE)))	

	## setup scenario life history list
	## 1-parameter model for selectivity and maturity
	lh_list <- adj_variation(SigmaR=0.01, SigmaF=0.01, SigmaC=0.01, SigmaI=0.01, CVlen=0.1, rho=0, selex_param=2, mat_param=1, nseasons=1)
	lh_fig(lh_list)

	### ----- generate data -----------###
	### data rich cases only
	rich_dir <- equil_N_dir_vec[grepl(data_vec[1],equil_N_dir_vec)]
	rich_modcombos <- equil_N_modcombos[which(equil_N_modcombos$Data_avail==data_vec[1]),]
	alt_dir <- equil_N_dir_vec[grepl(data_vec[1],equil_N_dir_vec)==FALSE]
	alt_modcombos <- equil_N_modcombos[which(equil_N_modcombos$Data_avail!=data_vec[1]),]

	### ----- use parallel cores -----------###
	# registerDoParallel(cores=ncores)	

	start_gen <- Sys.time()
	foreach(loop=1:length(rich_dir), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=rich_dir[loop], data_avail=as.character(rich_modcombos[loop,"Data_avail"]), itervec=itervec, Fdynamics=as.character(strsplit(rich_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(rich_modcombos[loop,"Rdyn"],"_")[[1]][2]),  lh=lh_list[[as.character(strsplit(rich_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=20, Nyears_comp=20, comp_sample=as.numeric(strsplit(rich_modcombos[loop,"SampleSize"],"_")[[1]][2]), rewrite=FALSE, init_depl=c(0.05,0.95), pool=TRUE, mismatch=FALSE)
	end_gen <- Sys.time() - start_gen


	# copy data rich cases to data poor
	# when written using parallel cores, sometimes files not writing properly
	start_regen <- Sys.time()
	for(loop in 1:length(rich_dir)){
		copy_sim(fromdir=rich_dir[loop], fromcombos=rich_modcombos[loop,], todir=alt_dir, itervec=itervec, rewrite=FALSE, res_dir="equil_samplesize")
	}
	end_regen <- Sys.time() - start_regen

	### ----- run LIME -----------###
	lime_dirs <- equil_N_dir_vec[which(grepl("LBSPR",equil_N_dir_vec)==FALSE)]
	lbspr_dirs <- equil_N_dir_vec[which(grepl("LBSPR",equil_N_dir_vec))]
	lime_combos <- equil_N_modcombos[which(grepl("LBSPR", equil_N_modcombos[,"Data_avail"])==FALSE),]
	lbspr_combos <- equil_N_modcombos[which(grepl("LBSPR", equil_N_modcombos[,"Data_avail"])),]	

	start_run <- Sys.time()
	foreach(loop=1:length(lime_dirs), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=lime_combos[loop,"C_opt"], LFdist=1, param_adjust=c("SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(0.737,0.2,0.2,0.2),  fix_param=FALSE, newtonsteps=3)
	end_run <- Sys.time() - start_run

	## pick up iterations not run in parallel
	for(loop in 1:length(lime_dirs)){
		run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=lime_combos[loop,"C_opt"], LFdist=1, param_adjust=c("SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(0.737,0.2,0.2,0.2),  fix_param=FALSE, newtonsteps=3)
	}


############################################################
#### SAMPLE SIZE base - with variation
############################################################
### ----- models to run ----------- ###
lh_vec <- c("Short", "Medium", "Long")
Fdyn_vec <- c("Ramp", "Increasing")#, "Decreasing")
Rdyn_vec <- c("AR")
data_vec <- c("LC10", "LC1")
SampleSize_vec <- c(1000, 500, 100, 50, 20)
itervec <- 1:100

base_N_modcombos <- expand.grid("Data_avail"=data_vec, "SampleSize"=paste0("SampleSize_", SampleSize_vec), "LH"=paste0("LH_",lh_vec), "Fdyn"=paste0("F_",Fdyn_vec), "Rdyn"=paste0("R_",Rdyn_vec), stringsAsFactors=FALSE)
base_N_modcombos$C_opt <- rep(0, nrow(base_N_modcombos))
	base_N_modcombos$C_opt[which(grepl("Catch",base_N_modcombos[,"Data_avail"]))] <- 2

### ----- equilibrium test ----------- ###

	base_N_dir <- file.path(main_dir, "base_samplesize")
	# unlink(base_N_dir, TRUE)
	dir.create(base_N_dir, showWarnings=FALSE)	

	## setup equilibrium dirs
	base_N_dir_vec <- model_paths(res_dir=base_N_dir, modcombos=base_N_modcombos[,-c(ncol(base_N_modcombos))])
	ignore <- sapply(1:length(base_N_dir_vec), function(m) sapply(itervec, function(x) dir.create(file.path(base_N_dir_vec[m], x), showWarnings=FALSE)))	

	## setup scenario life history list
	lh_list <- adj_variation(SigmaR=0.737, SigmaF=0.2, SigmaC=0.2, SigmaI=0.2, CVlen=0.1, rho=0.426, selex_param=2, mat_param=1,  nseasons=1)
	# lh_fig(lh_list)

	### ----- generate data -----------###
	### data rich cases only
	rich_dir <- base_N_dir_vec[grepl(data_vec[1],base_N_dir_vec)]
	rich_modcombos <- base_N_modcombos[which(base_N_modcombos$Data_avail==data_vec[1]),]
	alt_dir <- base_N_dir_vec[grepl(data_vec[1],base_N_dir_vec)==FALSE]
	alt_modcombos <- base_N_modcombos[which(base_N_modcombos$Data_avail!=data_vec[1]),]

	### ----- use parallel cores -----------###
	registerDoParallel(cores=5)	

	start_gen <- Sys.time()
	foreach(loop=1:length(rich_dir), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=rich_dir[loop], data_avail=as.character(rich_modcombos[loop,"Data_avail"]), itervec=itervec, Fdynamics=as.character(strsplit(rich_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(rich_modcombos[loop,"Rdyn"],"_")[[1]][2]),  lh=lh_list[[as.character(strsplit(rich_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=20, Nyears_comp=20, comp_sample=as.numeric(strsplit(rich_modcombos[loop,"SampleSize"],"_")[[1]][2]), rewrite=FALSE, init_depl=c(0.05,0.95))
	end_gen <- Sys.time() - start_gen

	# copy data rich cases to data poor
	# when written using parallel cores, sometimes files not writing properly
	start_regen <- Sys.time()
	for(loop in 1:length(rich_dir)){
		copy_sim(fromdir=rich_dir[loop], fromcombos=rich_modcombos[loop,], todir=alt_dir, itervec=itervec, rewrite=TRUE, res_dir="base_samplesize")
	}
	end_regen <- Sys.time() - start_regen

	### ----- run LIME -----------###
	lime_dirs <- base_N_dir_vec[which(grepl("LBSPR",base_N_dir_vec)==FALSE)]
	lbspr_dirs <- base_N_dir_vec[which(grepl("LBSPR",base_N_dir_vec))]
	lime_combos <- base_N_modcombos[which(grepl("LBSPR", base_N_modcombos[,"Data_avail"])==FALSE),]
	lbspr_combos <- base_N_modcombos[which(grepl("LBSPR", base_N_modcombos[,"Data_avail"])),]	

	start_run <- Sys.time()
	foreach(loop=1:length(lime_dirs), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=lime_combos[loop,"C_opt"], LFdist=1, param_adjust=c("SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(0.737,0.2,0.2,0.2),  fix_param=FALSE)
	end_run <- Sys.time() - start_run
	
	## pick up iterations not run in parallel
	for(loop in 1:length(lime_dirs)){
		run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=lime_combos[loop,"C_opt"], LFdist=1, param_adjust=c("SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(0.737,0.2,0.2,0.2),  fix_param=FALSE)
	}

N_dirs <- c(equil_N_dir_vec, base_N_dir_vec)
# bp_N <- bias_precision(N_dirs, itervec=1:100)
# saveRDS(bp_N, file.path(res_dir, "sample_size_bias_precision.rds"))

bp_N <- readRDS(file.path(res_dir, "sample_size_bias_precision.rds"))


## sample size
data_vec <- paste0("SampleSize_", c(1000,500,100,50,20))
bmat <- pmat <- matrix(NA, nrow=length(data_vec), ncol=9)
rownames(bmat) <- rownames(pmat) <- data_vec
colnames(bmat) <- colnames(pmat) <- rep(c("equil", "two-way", "one-way"),3)

for(i in 1:length(data_vec)){
	idir <- c(which(grepl("Short", N_dirs) & grepl("/LC1/", N_dirs) & grepl(paste0(data_vec[i],"/"), N_dirs)), which(grepl("Medium", N_dirs) & grepl("/LC1/", N_dirs) & grepl(paste0(data_vec[i],"/"), N_dirs)), which(grepl("Long", N_dirs) & grepl("/LC1/", N_dirs) & grepl(paste0(data_vec[i],"/"), N_dirs)))
	b <- sapply(1:length(idir), function(x) median((bp_N$relerr[,idir[x]]), na.rm=TRUE))
	p <- sapply(1:length(idir), function(x) median(abs(bp_N$relerr[,idir[x]]), na.rm=TRUE))
	bmat[i,] <- b
	pmat[i,] <- p
	rm(b)
	rm(p)
}

write.csv(rbind(bmat,pmat), file.path(res_dir, "sample_size_biases_precision.csv"))

