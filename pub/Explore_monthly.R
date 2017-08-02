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
main_dir <- "C:\\Git_Projects\\LIME_sim"
# main_dir <- "F:\\Merrill\\Git_Projects\\LIME_sim"

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
		sim_LBSPR(dir=equil_LBSPR_dir_vec[loop], lh=lh_list[[as.character(strsplit(LBSPR_modcombos[loop,"LH"],"_")[[1]][2])]], itervec=itervec, Nyears=as.numeric(strsplit(LBSPR_modcombos[loop,"Data_avail"],"LC")[[1]][2]), sample_size=200)
	}

	### ----- run LIME -----------###
	lime_dirs <- equil_LBSPR_dir_vec[which(grepl("/LBSPR",equil_LBSPR_dir_vec)==FALSE)]
	lbspr_dirs <- equil_LBSPR_dir_vec[which(grepl("/LBSPR",equil_LBSPR_dir_vec))]
	lime_combos <- LBSPR_modcombos[which(grepl("LBSPR", LBSPR_modcombos[,"Data_avail"])==FALSE),]
	lbspr_combos <- LBSPR_modcombos[which(grepl("LBSPR", LBSPR_modcombos[,"Data_avail"])),]	

	start_run <- Sys.time()
	foreach(loop=1:length(lime_dirs), .packages=c('TMB','LIME')) %dopar
	   run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(0.737,0.2,0.2,0.2),  fix_param=FALSE)
	end_run <- Sys.time() - start_run

	### ----- run LBSPR -----------###	

	start_run <- Sys.time()
	foreach(loop=1:length(lbspr_dirs), .packages=c('LBSPR','LIME')) %dopar% run_LBSPR(modpath=lbspr_dirs[loop], lh=lh_list[[as.character(strsplit(lbspr_combos[loop,"LH"],"_")[[1]][2])]], itervec=itervec, species=NULL, rewrite=FALSE, simulation=TRUE)
	end_altrun <- Sys.time() - start_run


############################################################
## monthly simulation, annual run
############################################################
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
	# registerDoParallel(cores=ncores)	

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
lh_vec <- c("Short", "Medium", "Long")
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
	foreach(loop=1:length(equil_runMonthly_dir_vec), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=equil_runMonthly_dir_vec[loop], data_avail=as.character(equil_modcombos[loop,"Data_avail"]), itervec=itervec, Fdynamics=as.character(strsplit(equil_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(equil_modcombos[loop,"Rdyn"],"_")[[1]][2]),  lh=lh_list[[as.character(strsplit(equil_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=as.numeric(strsplit(equil_modcombos[loop,"Data_avail"],"LC")[[1]][2]), comp_sample=as.numeric(strsplit(equil_modcombos[loop,"SampleSize"],"_")[[1]][2]), rewrite=FALSE, init_depl=c(0.05,0.95), pool=FALSE)
	end_gen <- Sys.time() - start_gen


	### ----- run LIME -----------###
	lime_dirs <- equil_runMonthly_dir_vec[which(grepl("LBSPR",equil_runMonthly_dir_vec)==FALSE)]
	lbspr_dirs <- equil_runMonthly_dir_vec[which(grepl("LBSPR",equil_runMonthly_dir_vec))]
	lime_combos <- equil_modcombos[which(grepl("LBSPR", equil_modcombos[,"Data_avail"])==FALSE),]
	lbspr_combos <- equil_modcombos[which(grepl("LBSPR", equil_modcombos[,"Data_avail"])),]	

	start_run <- Sys.time()
	foreach(loop=1:length(lime_dirs), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=lime_combos[loop,"C_opt"], LFdist=1, param_adjust=c("SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(0.2,0.2,0.2,0.2),  fix_param=FALSE)
	end_run <- Sys.time() - start_run
