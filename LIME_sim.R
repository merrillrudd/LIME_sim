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
if(grepl("F:", main_dir)) ncores <- 10

funs <- list.files(file.path(main_dir, "R"))
ignore <- sapply(1:length(funs), function(x) source(file.path(main_dir,"R", funs[x])))

fig_dir <- file.path(main_dir, "figs")
dir.create(fig_dir, showWarnings=FALSE)


## look at possible stocks
# load(file=file.path(main_dir, "FAODataR"))
# full_data <- data.frame(Area,Country,ScientificName,Species,Catch,Measure,stringsAsFactors=FALSE)	#all measures are Quantity (tonnes)

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
	foreach(loop=1:length(rich_dir), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=rich_dir[loop], data_avail=as.character(rich_modcombos[loop,"Data_avail"]), itervec=itervec, Fdynamics=as.character(strsplit(rich_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(rich_modcombos[loop,"Rdyn"],"_")[[1]][2]),  lh=lh_list[[as.character(strsplit(rich_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=20, comp_sample=as.numeric(strsplit(rich_modcombos[loop,"SampleSize"],"_")[[1]][2]), rewrite=FALSE, init_depl=c(0.05,0.95))
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
	foreach(loop=1:length(lime_dirs), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=lime_combos[loop,"C_opt"], LFdist=1, param_adjust=c("SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(0.7,0.2,0.2,0.2),  fix_param=FALSE)
	end_run <- Sys.time() - start_run

	## pick up iterations not run in parallel
	for(loop in 1:length(lime_dirs)){
		run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=lime_combos[loop,"C_opt"], LFdist=1, param_adjust=c("SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(0.777,0.2,0.2,0.2),  fix_param=FALSE)
	}

	### ----- run LBSPR -----------###	

	start_run <- Sys.time()
	foreach(loop=1:length(lbspr_dirs), .packages=c('LBSPR','LIME')) %dopar% run_LBSPR(modpath=lbspr_dirs[loop], lh=lh_list[[as.character(strsplit(lbspr_combos[loop,"LH"],"_")[[1]][2])]], itervec=itervec, species=NULL, rewrite=FALSE, simulation=TRUE)
	end_altrun <- Sys.time() - start_run

######### compare sample size ##################
### ----- models to run ----------- ###
lh_vec <- c("Short", "Medium", "Long")
Fdyn_vec <- "Constant"
Rdyn_vec <- "Constant"
data_vec <- c("Index_Catch_LC20", "LC10", "LC1")
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
	rich_dir <- equil_N_dir_vec[grepl("Index_Catch_LC20",equil_N_dir_vec)]
	rich_modcombos <- equil_N_modcombos[which(equil_N_modcombos$Data_avail=="Index_Catch_LC20"),]
	alt_dir <- equil_N_dir_vec[grepl("Index_Catch_LC20",equil_N_dir_vec)==FALSE]
	alt_modcombos <- equil_N_modcombos[which(equil_N_modcombos$Data_avail!="Index_Catch_LC20"),]

	### ----- use parallel cores -----------###
	# registerDoParallel(cores=ncores)	

	start_gen <- Sys.time()
	foreach(loop=1:length(rich_dir), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=rich_dir[loop], data_avail=as.character(rich_modcombos[loop,"Data_avail"]), itervec=itervec, Fdynamics=as.character(strsplit(rich_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(rich_modcombos[loop,"Rdyn"],"_")[[1]][2]),  lh=lh_list[[as.character(strsplit(rich_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=20, comp_sample=as.numeric(strsplit(rich_modcombos[loop,"SampleSize"],"_")[[1]][2]), rewrite=FALSE, init_depl=c(0.05,0.95))
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
	foreach(loop=1:length(lime_dirs), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=lime_combos[loop,"C_opt"], LFdist=1, param_adjust=c("SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(0.7,0.2,0.2,0.2),  fix_param=FALSE)
	end_run <- Sys.time() - start_run

	## pick up iterations not run in parallel
	for(loop in 1:length(lime_dirs)){
		run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=lime_combos[loop,"C_opt"], LFdist=1, param_adjust=c("SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(0.777,0.2,0.2,0.2),  fix_param=FALSE)
	}

	# ### ----- run LBSPR -----------###	

	# start_run <- Sys.time()
	# foreach(loop=1:length(lbspr_dirs), .packages=c('LBSPR','LIME')) %dopar% run_LBSPR(modpath=lbspr_dirs[loop], lh=lh_list[[as.character(strsplit(lbspr_combos[loop,"LH"],"_")[[1]][2])]], itervec=itervec, species=NULL, rewrite=FALSE, simulation=TRUE)
	# end_altrun <- Sys.time() - start_run


### ----- equilibrium test ----------- ###

	# equil15_dir <- file.path(main_dir, "equil_15bins")
	# # unlink(equil15_dir, TRUE)
	# dir.create(equil15_dir, showWarnings=FALSE)	

	# ## setup equilibrium dirs
	# equil15_dir_vec <- model_paths(res_dir=equil15_dir, modcombos=equil_modcombos[,-c(ncol(equil_modcombos))])
	# ignore <- sapply(1:length(equil15_dir_vec), function(m) sapply(itervec, function(x) dir.create(file.path(equil15_dir_vec[m], x), showWarnings=FALSE)))	

	# ## setup scenario life history list
	# ## 1-parameter model for selectivity and maturity
	# lh_list <- adj_variation(SigmaR=0.01, SigmaF=0.01, SigmaC=0.01, SigmaI=0.01, CVlen=0.1, rho=0, selex_param=1, mat_param=1, nseasons=1, nbins=15)
	# # lh_fig(lh_list, save=TRUE)

	# ### ----- generate data -----------###
	# ### data rich cases only
	# rich_dir <- equil15_dir_vec[grepl("Index_Catch_LC20",equil15_dir_vec)]
	# rich_modcombos <- equil_modcombos[which(equil_modcombos$Data_avail=="Index_Catch_LC20"),]
	# alt_dir <- equil15_dir_vec[grepl("Index_Catch_LC20",equil15_dir_vec)==FALSE]
	# alt_modcombos <- equil_modcombos[which(equil_modcombos$Data_avail!="Index_Catch_LC20"),]

	# ### ----- use parallel cores -----------###
	# registerDoParallel(cores=ncores)	

	# start_gen <- Sys.time()
	# foreach(loop=1:length(rich_dir), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=rich_dir[loop], data_avail=as.character(rich_modcombos[loop,"Data_avail"]), itervec=itervec, Fdynamics=as.character(strsplit(rich_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(rich_modcombos[loop,"Rdyn"],"_")[[1]][2]),  lh=lh_list[[as.character(strsplit(rich_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=20, comp_sample=as.numeric(strsplit(rich_modcombos[loop,"SampleSize"],"_")[[1]][2]), rewrite=FALSE, init_depl=c(0.05,0.95))
	# end_gen <- Sys.time() - start_gen

	# # copy data rich cases to data poor
	# # when written using parallel cores, sometimes files not writing properly
	# start_regen <- Sys.time()
	# for(loop in 1:length(rich_dir)){
	# 	copy_sim(fromdir=rich_dir[loop], fromcombos=rich_modcombos[loop,], todir=alt_dir, itervec=itervec, rewrite=FALSE, res_dir="equil_15bins")
	# }
	# end_regen <- Sys.time() - start_regen

	# ### ----- run LIME -----------###
	# lime_dirs <- equil15_dir_vec[which(grepl("LBSPR",equil15_dir_vec)==FALSE)]
	# lbspr_dirs <- equil15_dir_vec[which(grepl("LBSPR",equil15_dir_vec))]
	# lime_combos <- equil_modcombos[which(grepl("LBSPR", equil_modcombos[,"Data_avail"])==FALSE),]
	# lbspr_combos <- equil_modcombos[which(grepl("LBSPR", equil_modcombos[,"Data_avail"])),]	

	# start_run <- Sys.time()
	# foreach(loop=1:length(lime_dirs), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=lime_combos[loop,"C_opt"], LFdist=1, param_adjust=c("SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(0.7,0.2,0.2,0.2),  fix_param=FALSE)
	# end_run <- Sys.time() - start_run

	# ### ----- run LBSPR -----------###	

	# start_run <- Sys.time()
	# foreach(loop=1:length(lbspr_dirs), .packages=c('LBSPR','LIME')) %dopar% run_LBSPR(modpath=lbspr_dirs[loop], lh=lh_list[[as.character(strsplit(lbspr_combos[loop,"LH"],"_")[[1]][2])]], itervec=itervec, species=NULL, rewrite=FALSE, simulation=TRUE)
	# end_altrun <- Sys.time() - start_run

### ----- equilibrium test ----------- ###

	equil_12seasons_dir <- file.path(main_dir, "equil_12seasons")
	# unlink(equil_12seasons_dir, TRUE)
	dir.create(equil_12seasons_dir, showWarnings=FALSE)	

	## setup equilibrium dirs
	equil_12seasons_dir_vec <- model_paths(res_dir=equil_12seasons_dir, modcombos=equil_modcombos[,-c(ncol(equil_modcombos))])
	ignore <- sapply(1:length(equil_12seasons_dir_vec), function(m) sapply(itervec, function(x) dir.create(file.path(equil_12seasons_dir_vec[m], x), showWarnings=FALSE)))	

	## setup scenario life history list
	## 1-parameter model for selectivity and maturity
	lh_list <- adj_variation(SigmaR=0.01, SigmaF=0.01, SigmaC=0.01, SigmaI=0.01, CVlen=0.1, rho=0, selex_param=2, mat_param=1, nseasons=12)
	# lh_fig(lh_list, save=TRUE)

	### ----- generate data -----------###
	### data rich cases only
	rich_dir <- equil_12seasons_dir_vec[grepl("Index_Catch_LC20",equil_12seasons_dir_vec)]
	rich_modcombos <- equil_modcombos[which(equil_modcombos$Data_avail=="Index_Catch_LC20"),]
	alt_dir <- equil_12seasons_dir_vec[grepl("Index_Catch_LC20",equil_12seasons_dir_vec)==FALSE]
	alt_modcombos <- equil_modcombos[which(equil_modcombos$Data_avail!="Index_Catch_LC20"),]

	### ----- use parallel cores -----------###
	# registerDoParallel(cores=ncores)	

	start_gen <- Sys.time()
	foreach(loop=1:length(rich_dir), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=rich_dir[loop], data_avail=as.character(rich_modcombos[loop,"Data_avail"]), itervec=itervec, Fdynamics=as.character(strsplit(rich_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(rich_modcombos[loop,"Rdyn"],"_")[[1]][2]),  lh=lh_list[[as.character(strsplit(rich_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=20, comp_sample=as.numeric(strsplit(rich_modcombos[loop,"SampleSize"],"_")[[1]][2]), rewrite=FALSE, init_depl=c(0.05,0.95))
	end_gen <- Sys.time() - start_gen

	# copy data rich cases to data poor
	# when written using parallel cores, sometimes files not writing properly
	start_regen <- Sys.time()
	for(loop in 1:length(rich_dir)){
		copy_sim(fromdir=rich_dir[loop], fromcombos=rich_modcombos[loop,], todir=alt_dir, itervec=itervec, rewrite=FALSE, res_dir="equil_12seasons")
	}
	end_regen <- Sys.time() - start_regen

	### ----- run LIME -----------###
	lime_dirs <- equil_12seasons_dir_vec[which(grepl("LBSPR",equil_12seasons_dir_vec)==FALSE)]
	lbspr_dirs <- equil_12seasons_dir_vec[which(grepl("LBSPR",equil_12seasons_dir_vec))]
	lime_combos <- equil_modcombos[which(grepl("LBSPR", equil_modcombos[,"Data_avail"])==FALSE),]
	lbspr_combos <- equil_modcombos[which(grepl("LBSPR", equil_modcombos[,"Data_avail"])),]	

	lh_list <- adj_variation(SigmaR=0.01, SigmaF=0.01, SigmaC=0.01, SigmaI=0.01, CVlen=0.1, rho=0, selex_param=2, mat_param=1, nseasons=1)


	start_run <- Sys.time()
	foreach(loop=1:length(lime_dirs), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=lime_combos[loop,"C_opt"], LFdist=1, param_adjust=c("SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(0.7,0.2,0.2,0.2),  fix_param=FALSE)
	end_run <- Sys.time() - start_run


	### ----- run LBSPR -----------###	

	start_run <- Sys.time()
	foreach(loop=1:length(lbspr_dirs), .packages=c('LBSPR','LIME')) %dopar% run_LBSPR(modpath=lbspr_dirs[loop], lh=lh_list[[as.character(strsplit(lbspr_combos[loop,"LH"],"_")[[1]][2])]], itervec=itervec, species=NULL, rewrite=FALSE, simulation=TRUE)
	end_altrun <- Sys.time() - start_run



### ----- equilibrium test ----------- ###
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
	registerDoParallel(cores=1)	

	# start_gen <- Sys.time()
	# foreach(loop=1:length(equil_LBSPR_dir_vec), .packages=c('TMB','LIME')) %dopar% sim_LBSPR(dir=equil_LBSPR_dir_vec[loop], lh=lh_list[[as.character(strsplit(LBSPR_modcombos[loop,"LH"],"_")[[1]][2])]], itervec=itervec, Nyears=1)
	# end_gen <- Sys.time() - start_gen

	for(loop in 1:length(equil_LBSPR_dir_vec)){
		sim_LBSPR(dir=equil_LBSPR_dir_vec[loop], lh=lh_list[[as.character(strsplit(LBSPR_modcombos[loop,"LH"],"_")[[1]][2])]], itervec=itervec, Nyears=1)
	}

	### ----- run LIME -----------###
	lime_dirs <- equil_LBSPR_dir_vec[which(grepl("/LBSPR",equil_LBSPR_dir_vec)==FALSE)]
	lbspr_dirs <- equil_LBSPR_dir_vec[which(grepl("/LBSPR",equil_LBSPR_dir_vec))]
	lime_combos <- LBSPR_modcombos[which(grepl("LBSPR", LBSPR_modcombos[,"Data_avail"])==FALSE),]
	lbspr_combos <- LBSPR_modcombos[which(grepl("LBSPR", LBSPR_modcombos[,"Data_avail"])),]	

	start_run <- Sys.time()
	foreach(loop=1:length(lime_dirs), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=0, LFdist=1, param_adjust=c("SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(0.7,0.2,0.2,0.2),  fix_param=FALSE, theta_type=0)
	end_run <- Sys.time() - start_run

	### ----- run LBSPR -----------###	

	start_run <- Sys.time()
	foreach(loop=1:length(lbspr_dirs), .packages=c('LBSPR','LIME')) %dopar% run_LBSPR(modpath=lbspr_dirs[loop], lh=lh_list[[as.character(strsplit(lbspr_combos[loop,"LH"],"_")[[1]][2])]], itervec=itervec, species=NULL, rewrite=FALSE, simulation=TRUE)
	end_altrun <- Sys.time() - start_run

# bp_check <- bias_precision(dirs=equil_LBSPR_dir_vec, itervec=itervec, param="SPR")


############################################################
#### equilibrium - selectivity
############################################################

### ----- models to run ----------- ###
lh_vec <- c("Short", "Medium", "Long")
Fdyn_vec <- "Constant"
Rdyn_vec <- "Constant"
data_vec <- c("Index_Catch_LC20", "LC10", "LC1", "LBSPR10", "LBSPR1")
SampleSize_vec <- 200
itervec <- 1:100

equil_pp_modcombos <- expand.grid("Data_avail"=data_vec, "SampleSize"=paste0("SampleSize_", SampleSize_vec), "LH"=paste0("LH_",lh_vec), "Fdyn"=paste0("F_",Fdyn_vec), "Rdyn"=paste0("R_",Rdyn_vec), stringsAsFactors=FALSE)
equil_pp_modcombos$C_opt <- rep(0, nrow(equil_pp_modcombos))
	equil_pp_modcombos$C_opt[which(grepl("Catch",equil_pp_modcombos[,"Data_avail"]))] <- 2

### ----- equil_1p1pibrium test ----------- ###

	equil_1p1p_dir <- file.path(main_dir, "equil_1p1p")
	# unlink(equil_1p1p_dir, TRUE)
	dir.create(equil_1p1p_dir, showWarnings=FALSE)	

	## setup equil_1p1pibrium dirs
	equil_1p1p_dir_vec <- model_paths(res_dir=equil_1p1p_dir, modcombos=equil_pp_modcombos[,-c(ncol(equil_pp_modcombos))])
	ignore <- sapply(1:length(equil_1p1p_dir_vec), function(m) sapply(itervec, function(x) dir.create(file.path(equil_1p1p_dir_vec[m], x), showWarnings=FALSE)))	

	## setup scenario life history list
	## 1-parameter model for selectivity and maturity
	lh_list <- adj_variation(SigmaR=0.01, SigmaF=0.01, SigmaC=0.01, SigmaI=0.01, CVlen=0.1, rho=0, selex_param=1, mat_param=1, nseasons=1)
	lh_fig(lh_list)

	### ----- generate data -----------###
	### data rich cases only
	rich_dir <- equil_1p1p_dir_vec[grepl("Index_Catch_LC20",equil_1p1p_dir_vec)]
	rich_modcombos <- equil_pp_modcombos[which(equil_pp_modcombos$Data_avail=="Index_Catch_LC20"),]
	alt_dir <- equil_1p1p_dir_vec[grepl("Index_Catch_LC20",equil_1p1p_dir_vec)==FALSE]
	alt_modcombos <- equil_pp_modcombos[which(equil_pp_modcombos$Data_avail!="Index_Catch_LC20"),]

	### ----- use parallel cores -----------###
	registerDoParallel(cores=ncores)	

	start_gen <- Sys.time()
	foreach(loop=1:length(rich_dir), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=rich_dir[loop], data_avail=as.character(rich_modcombos[loop,"Data_avail"]), itervec=itervec, Fdynamics=as.character(strsplit(rich_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(rich_modcombos[loop,"Rdyn"],"_")[[1]][2]),  lh=lh_list[[as.character(strsplit(rich_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=20, comp_sample=as.numeric(strsplit(rich_modcombos[loop,"SampleSize"],"_")[[1]][2]), rewrite=FALSE, init_depl=c(0.05,0.95))
	end_gen <- Sys.time() - start_gen

	# copy data rich cases to data poor
	# when written using parallel cores, sometimes files not writing properly
	start_regen <- Sys.time()
	for(loop in 1:length(rich_dir)){
		copy_sim(fromdir=rich_dir[loop], fromcombos=rich_modcombos[loop,], todir=alt_dir, itervec=itervec, rewrite=FALSE, res_dir="equil_1p1p")
	}
	end_regen <- Sys.time() - start_regen

	### ----- run LIME -----------###
	lime_dirs <- equil_1p1p_dir_vec[which(grepl("LBSPR",equil_1p1p_dir_vec)==FALSE)]
	lbspr_dirs <- equil_1p1p_dir_vec[which(grepl("LBSPR",equil_1p1p_dir_vec))]
	lime_combos <- equil_pp_modcombos[which(grepl("LBSPR", equil_pp_modcombos[,"Data_avail"])==FALSE),]
	lbspr_combos <- equil_pp_modcombos[which(grepl("LBSPR", equil_pp_modcombos[,"Data_avail"])),]	

	start_run <- Sys.time()
	foreach(loop=1:length(lime_dirs), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=lime_combos[loop,"C_opt"], LFdist=1, param_adjust=c("SigmaR","SigmaF","SigmaC","SigmaI", "Sdelta"), val_adjust=c(0.7,0.2,0.2,0.2,1.1),  fix_param="logSdelta")
	end_run <- Sys.time() - start_run

	for(loop in 1:length(lime_dirs)){
		run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=lime_combos[loop,"C_opt"], LFdist=1, param_adjust=c("SigmaR","SigmaF","SigmaC","SigmaI", "Sdelta"), val_adjust=c(0.7,0.2,0.2,0.2,1.1),  fix_param="logSdelta")
	}


	### ----- run LBSPR -----------###	

	start_run <- Sys.time()
	foreach(loop=1:length(lbspr_dirs), .packages=c('LBSPR','LIME')) %dopar% run_LBSPR(modpath=lbspr_dirs[loop], lh=lh_list[[as.character(strsplit(lbspr_combos[loop,"LH"],"_")[[1]][2])]], itervec=itervec, species=NULL, rewrite=FALSE, simulation=TRUE)
	end_altrun <- Sys.time() - start_run

### ----- equil_2p1pibrium test ----------- ###

	equil_2p1p_dir <- file.path(main_dir, "equil_2p1p")
	# unlink(equil_2p1p_dir, TRUE)
	dir.create(equil_2p1p_dir, showWarnings=FALSE)	

	## setup equil_2p1pibrium dirs
	equil_2p1p_dir_vec <- model_paths(res_dir=equil_2p1p_dir, modcombos=equil_pp_modcombos[,-c(ncol(equil_pp_modcombos))])
	ignore <- sapply(1:length(equil_2p1p_dir_vec), function(m) sapply(itervec, function(x) dir.create(file.path(equil_2p1p_dir_vec[m], x), showWarnings=FALSE)))	

	## setup scenario life history list
	## 1-parameter model for selectivity and maturity
	lh_list <- adj_variation(SigmaR=0.01, SigmaF=0.01, SigmaC=0.01, SigmaI=0.01, CVlen=0.1, rho=0, selex_param=2, mat_param=1, nseasons=1)
	lh_fig(lh_list)

	### ----- generate data -----------###
	### data rich cases only
	rich_dir <- equil_2p1p_dir_vec[grepl("Index_Catch_LC20",equil_2p1p_dir_vec)]
	rich_modcombos <- equil_pp_modcombos[which(equil_pp_modcombos$Data_avail=="Index_Catch_LC20"),]
	alt_dir <- equil_2p1p_dir_vec[grepl("Index_Catch_LC20",equil_2p1p_dir_vec)==FALSE]
	alt_modcombos <- equil_pp_modcombos[which(equil_pp_modcombos$Data_avail!="Index_Catch_LC20"),]

	### ----- use parallel cores -----------###
	# registerDoParallel(cores=ncores)	

	start_gen <- Sys.time()
	foreach(loop=1:length(rich_dir), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=rich_dir[loop], data_avail=as.character(rich_modcombos[loop,"Data_avail"]), itervec=itervec, Fdynamics=as.character(strsplit(rich_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(rich_modcombos[loop,"Rdyn"],"_")[[1]][2]),  lh=lh_list[[as.character(strsplit(rich_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=20, comp_sample=as.numeric(strsplit(rich_modcombos[loop,"SampleSize"],"_")[[1]][2]), rewrite=FALSE, init_depl=c(0.05,0.95))
	end_gen <- Sys.time() - start_gen

	# copy data rich cases to data poor
	# when written using parallel cores, sometimes files not writing properly
	start_regen <- Sys.time()
	for(loop in 1:length(rich_dir)){
		copy_sim(fromdir=rich_dir[loop], fromcombos=rich_modcombos[loop,], todir=alt_dir, itervec=itervec, rewrite=FALSE, res_dir="equil_2p1p")
	}
	end_regen <- Sys.time() - start_regen

	### ----- run LIME -----------###
	lime_dirs <- equil_2p1p_dir_vec[which(grepl("LBSPR",equil_2p1p_dir_vec)==FALSE)]
	lbspr_dirs <- equil_2p1p_dir_vec[which(grepl("LBSPR",equil_2p1p_dir_vec))]
	lime_combos <- equil_pp_modcombos[which(grepl("LBSPR", equil_pp_modcombos[,"Data_avail"])==FALSE),]
	lbspr_combos <- equil_pp_modcombos[which(grepl("LBSPR", equil_pp_modcombos[,"Data_avail"])),]	

	start_run <- Sys.time()
	foreach(loop=1:length(lime_dirs), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=lime_combos[loop,"C_opt"], LFdist=1, param_adjust=c("SigmaR","SigmaF","SigmaC","SigmaI", "Sdelta"), val_adjust=c(0.7,0.2,0.2,0.2,1.1),  fix_param="logSdelta")
	end_run <- Sys.time() - start_run


	### ----- run LBSPR -----------###	

	start_run <- Sys.time()
	foreach(loop=1:length(lbspr_dirs), .packages=c('LBSPR','LIME')) %dopar% run_LBSPR(modpath=lbspr_dirs[loop], lh=lh_list[[as.character(strsplit(lbspr_combos[loop,"LH"],"_")[[1]][2])]], itervec=itervec, species=NULL, rewrite=FALSE, simulation=TRUE)
	end_altrun <- Sys.time() - start_run


### ----- equil_1p2pibrium test ----------- ###

	equil_1p2p_dir <- file.path(main_dir, "equil_1p2p")
	# unlink(equil_1p2p_dir, TRUE)
	dir.create(equil_1p2p_dir, showWarnings=FALSE)	

	## setup equil_1p2pibrium dirs
	equil_1p2p_dir_vec <- model_paths(res_dir=equil_1p2p_dir, modcombos=equil_pp_modcombos[,-c(ncol(equil_pp_modcombos))])
	ignore <- sapply(1:length(equil_1p2p_dir_vec), function(m) sapply(itervec, function(x) dir.create(file.path(equil_1p2p_dir_vec[m], x), showWarnings=FALSE)))	

	## setup scenario life history list
	## 1-parameter model for selectivity and maturity
	lh_list <- adj_variation(SigmaR=0.01, SigmaF=0.01, SigmaC=0.01, SigmaI=0.01, CVlen=0.1, rho=0, selex_param=1, mat_param=1, nseasons=1)
	lh_fig(lh_list)

	### ----- generate data -----------###
	### data rich cases only
	rich_dir <- equil_1p2p_dir_vec[grepl("Index_Catch_LC20",equil_1p2p_dir_vec)]
	rich_modcombos <- equil_pp_modcombos[which(equil_pp_modcombos$Data_avail=="Index_Catch_LC20"),]
	alt_dir <- equil_1p2p_dir_vec[grepl("Index_Catch_LC20",equil_1p2p_dir_vec)==FALSE]
	alt_modcombos <- equil_pp_modcombos[which(equil_pp_modcombos$Data_avail!="Index_Catch_LC20"),]

	### ----- use parallel cores -----------###
	# registerDoParallel(cores=ncores)	

	start_gen <- Sys.time()
	foreach(loop=1:length(rich_dir), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=rich_dir[loop], data_avail=as.character(rich_modcombos[loop,"Data_avail"]), itervec=itervec, Fdynamics=as.character(strsplit(rich_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(rich_modcombos[loop,"Rdyn"],"_")[[1]][2]),  lh=lh_list[[as.character(strsplit(rich_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=20, comp_sample=as.numeric(strsplit(rich_modcombos[loop,"SampleSize"],"_")[[1]][2]), rewrite=FALSE, init_depl=c(0.05,0.95))
	end_gen <- Sys.time() - start_gen

	# copy data rich cases to data poor
	# when written using parallel cores, sometimes files not writing properly
	start_regen <- Sys.time()
	for(loop in 1:length(rich_dir)){
		copy_sim(fromdir=rich_dir[loop], fromcombos=rich_modcombos[loop,], todir=alt_dir, itervec=itervec, rewrite=FALSE, res_dir="equil_1p2p")
	}
	end_regen <- Sys.time() - start_regen

	### ----- run LIME -----------###
	lime_dirs <- equil_1p2p_dir_vec[which(grepl("LBSPR",equil_1p2p_dir_vec)==FALSE)]
	lbspr_dirs <- equil_1p2p_dir_vec[which(grepl("LBSPR",equil_1p2p_dir_vec))]
	lime_combos <- equil_pp_modcombos[which(grepl("LBSPR", equil_pp_modcombos[,"Data_avail"])==FALSE),]
	lbspr_combos <- equil_pp_modcombos[which(grepl("LBSPR", equil_pp_modcombos[,"Data_avail"])),]	

	start_run <- Sys.time()
	foreach(loop=1:length(lime_dirs), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=lime_combos[loop,"C_opt"], LFdist=1, param_adjust=c("SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(0.7,0.2,0.2,0.2),  fix_param=FALSE)
	end_run <- Sys.time() - start_run

	### ----- run LBSPR -----------###	

	start_run <- Sys.time()
	foreach(loop=1:length(lbspr_dirs), .packages=c('LBSPR','LIME')) %dopar% run_LBSPR(modpath=lbspr_dirs[loop], lh=lh_list[[as.character(strsplit(lbspr_combos[loop,"LH"],"_")[[1]][2])]], itervec=itervec, species=NULL, rewrite=FALSE, simulation=TRUE)
	end_altrun <- Sys.time() - start_run

### ----- equil_2p2pibrium test ----------- ###

	equil_2p2p_dir <- file.path(main_dir, "equil_2p2p")
	# unlink(equil_2p2p_dir, TRUE)
	dir.create(equil_2p2p_dir, showWarnings=FALSE)	

	## setup equil_2p2pibrium dirs
	equil_2p2p_dir_vec <- model_paths(res_dir=equil_2p2p_dir, modcombos=equil_pp_modcombos[,-c(ncol(equil_pp_modcombos))])
	ignore <- sapply(1:length(equil_2p2p_dir_vec), function(m) sapply(itervec, function(x) dir.create(file.path(equil_2p2p_dir_vec[m], x), showWarnings=FALSE)))	

	## setup scenario life history list
	## 1-parameter model for selectivity and maturity
	lh_list <- adj_variation(SigmaR=0.01, SigmaF=0.01, SigmaC=0.01, SigmaI=0.01, CVlen=0.1, rho=0, selex_param=2, mat_param=1, nseasons=1)
	lh_fig(lh_list)

	### ----- generate data -----------###
	### data rich cases only
	rich_dir <- equil_2p2p_dir_vec[grepl("Index_Catch_LC20",equil_2p2p_dir_vec)]
	rich_modcombos <- equil_pp_modcombos[which(equil_pp_modcombos$Data_avail=="Index_Catch_LC20"),]
	alt_dir <- equil_2p2p_dir_vec[grepl("Index_Catch_LC20",equil_2p2p_dir_vec)==FALSE]
	alt_modcombos <- equil_pp_modcombos[which(equil_pp_modcombos$Data_avail!="Index_Catch_LC20"),]

	### ----- use parallel cores -----------###
	# registerDoParallel(cores=ncores)	

	start_gen <- Sys.time()
	foreach(loop=1:length(rich_dir), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=rich_dir[loop], data_avail=as.character(rich_modcombos[loop,"Data_avail"]), itervec=itervec, Fdynamics=as.character(strsplit(rich_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(rich_modcombos[loop,"Rdyn"],"_")[[1]][2]),  lh=lh_list[[as.character(strsplit(rich_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=20, comp_sample=as.numeric(strsplit(rich_modcombos[loop,"SampleSize"],"_")[[1]][2]), rewrite=FALSE, init_depl=c(0.05,0.95))
	end_gen <- Sys.time() - start_gen

	# copy data rich cases to data poor
	# when written using parallel cores, sometimes files not writing properly
	start_regen <- Sys.time()
	for(loop in 1:length(rich_dir)){
		copy_sim(fromdir=rich_dir[loop], fromcombos=rich_modcombos[loop,], todir=alt_dir, itervec=itervec, rewrite=FALSE, res_dir="equil_2p2p")
	}
	end_regen <- Sys.time() - start_regen

	### ----- run LIME -----------###
	lime_dirs <- equil_2p2p_dir_vec[which(grepl("LBSPR",equil_2p2p_dir_vec)==FALSE)]
	lbspr_dirs <- equil_2p2p_dir_vec[which(grepl("LBSPR",equil_2p2p_dir_vec))]
	lime_combos <- equil_pp_modcombos[which(grepl("LBSPR", equil_pp_modcombos[,"Data_avail"])==FALSE),]
	lbspr_combos <- equil_pp_modcombos[which(grepl("LBSPR", equil_pp_modcombos[,"Data_avail"])),]	

	start_run <- Sys.time()
	foreach(loop=1:length(lime_dirs), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=lime_combos[loop,"C_opt"], LFdist=1, param_adjust=c("SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(0.7,0.2,0.2,0.2),  fix_param=FALSE)
	end_run <- Sys.time() - start_run

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
	foreach(loop=1:length(rich_dir), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=rich_dir[loop], data_avail=as.character(rich_modcombos[loop,"Data_avail"]), itervec=itervec, Fdynamics=as.character(strsplit(rich_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(rich_modcombos[loop,"Rdyn"],"_")[[1]][2]),  lh=lh_list[[as.character(strsplit(rich_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=20, comp_sample=as.numeric(strsplit(rich_modcombos[loop,"SampleSize"],"_")[[1]][2]), rewrite=FALSE, init_depl=c(0.05,0.95))
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
	foreach(loop=1:length(lime_dirs), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=lime_combos[loop,"C_opt"], LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.7,0.2,0.2),  fix_param=FALSE)
	end_run <- Sys.time() - start_run

	for(loop in 1:length(lime_dirs)){
		run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=lime_combos[loop,"C_opt"], LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.7,0.2,0.2),  fix_param=FALSE)
	}


	### ----- run LBSPR -----------###	

	start_run <- Sys.time()
	foreach(loop=1:length(lbspr_dirs), .packages=c('LBSPR','LIME')) %dopar% run_LBSPR(modpath=lbspr_dirs[loop], lh=lh_list[[as.character(strsplit(lbspr_combos[loop,"LH"],"_")[[1]][2])]], itervec=itervec, species=NULL, rewrite=FALSE, simulation=TRUE)
	end_altrun <- Sys.time() - start_run	


######### compare sample size ##################
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
	rich_dir <- base_N_dir_vec[grepl("LC10",base_N_dir_vec)]
	rich_modcombos <- base_N_modcombos[which(base_N_modcombos$Data_avail=="LC10"),]
	alt_dir <- base_N_dir_vec[grepl("LC10",base_N_dir_vec)==FALSE]
	alt_modcombos <- base_N_modcombos[which(base_N_modcombos$Data_avail!="LC10"),]

	### ----- use parallel cores -----------###
	registerDoParallel(cores=5)	

	start_gen <- Sys.time()
	foreach(loop=1:length(rich_dir), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=rich_dir[loop], data_avail=as.character(rich_modcombos[loop,"Data_avail"]), itervec=itervec, Fdynamics=as.character(strsplit(rich_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(rich_modcombos[loop,"Rdyn"],"_")[[1]][2]),  lh=lh_list[[as.character(strsplit(rich_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=20, comp_sample=as.numeric(strsplit(rich_modcombos[loop,"SampleSize"],"_")[[1]][2]), rewrite=FALSE, init_depl=c(0.05,0.95))
	end_gen <- Sys.time() - start_gen

	# copy data rich cases to data poor
	# when written using parallel cores, sometimes files not writing properly
	start_regen <- Sys.time()
	for(loop in 1:length(rich_dir)){
		copy_sim(fromdir=rich_dir[loop], fromcombos=rich_modcombos[loop,], todir=alt_dir, itervec=itervec, rewrite=FALSE, res_dir="base_samplesize")
	}
	end_regen <- Sys.time() - start_regen

	### ----- run LIME -----------###
	lime_dirs <- base_N_dir_vec[which(grepl("LBSPR",base_N_dir_vec)==FALSE)]
	lbspr_dirs <- base_N_dir_vec[which(grepl("LBSPR",base_N_dir_vec))]
	lime_combos <- base_N_modcombos[which(grepl("LBSPR", base_N_modcombos[,"Data_avail"])==FALSE),]
	lbspr_combos <- base_N_modcombos[which(grepl("LBSPR", base_N_modcombos[,"Data_avail"])),]	

	start_run <- Sys.time()
	foreach(loop=1:length(lime_dirs), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=lime_combos[loop,"C_opt"], LFdist=1, param_adjust=c("SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(0.7,0.2,0.2,0.2),  fix_param=FALSE)
	end_run <- Sys.time() - start_run
	
	## pick up iterations not run in parallel
	for(loop in 1:length(lime_dirs)){
		run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=lime_combos[loop,"C_opt"], LFdist=1, param_adjust=c("SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(0.777,0.2,0.2,0.2),  fix_param=FALSE)
	}
##########################################
### ----- base with 12 seasons assume 1 ----------- ###
	# base_12seasons_dir <- file.path(main_dir, "base_12seasons")
	# # unlink(base_12seasons_dir, TRUE)
	# dir.create(base_12seasons_dir, showWarnings=FALSE)	

	# ## setup  dirs
	# base_12seasons_dir_vec <- model_paths(res_dir=base_12seasons_dir, modcombos=base_modcombos[,-c(ncol(base_modcombos))])
	# ignore <- sapply(1:length(base_12seasons_dir_vec), function(m) sapply(itervec, function(x) dir.create(file.path(base_12seasons_dir_vec[m], x), showWarnings=FALSE)))	

	# ## setup scenario life history list
	# ## 1-parameter model for selectivity and maturity
	# lh_list <- adj_variation(SigmaR=0.737, SigmaF=0.2, SigmaC=0.2, SigmaI=0.2, CVlen=0.1, rho=0.426, selex_param=2, mat_param=1,  nseasons=12)
	# # lh_fig(lh_list, save=TRUE)

	# ### ----- generate data -----------###
	# ### data rich cases only
	# rich_dir <- base_12seasons_dir_vec[grepl("LC10",base_12seasons_dir_vec)]
	# rich_modcombos <- base_modcombos[which(base_modcombos$Data_avail=="LC10"),]
	# alt_dir <- base_12seasons_dir_vec[grepl("LC10",base_12seasons_dir_vec)==FALSE]
	# alt_modcombos <- base_modcombos[which(base_modcombos$Data_avail!="LC10"),]

	# ### ----- use parallel cores -----------###
	# # registerDoParallel(cores=ncores)	

	# start_gen <- Sys.time()
	# foreach(loop=1:length(rich_dir), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=rich_dir[loop], data_avail=as.character(rich_modcombos[loop,"Data_avail"]), itervec=itervec, Fdynamics=as.character(strsplit(rich_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(rich_modcombos[loop,"Rdyn"],"_")[[1]][2]),  lh=lh_list[[as.character(strsplit(rich_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=20, comp_sample=as.numeric(strsplit(rich_modcombos[loop,"SampleSize"],"_")[[1]][2]), rewrite=FALSE, init_depl=c(0.05,0.95), pool=TRUE)
	# end_gen <- Sys.time() - start_gen

	# for(loop in 1:length(rich_dir)){
	# 	generate_data(modpath=rich_dir[loop], data_avail=as.character(rich_modcombos[loop,"Data_avail"]), itervec=itervec, Fdynamics=as.character(strsplit(rich_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(rich_modcombos[loop,"Rdyn"],"_")[[1]][2]),  lh=lh_list[[as.character(strsplit(rich_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=20, comp_sample=as.numeric(strsplit(rich_modcombos[loop,"SampleSize"],"_")[[1]][2]), rewrite=FALSE, init_depl=c(0.05,0.95))
	# }

	# # copy data rich cases to data poor
	# # when written using parallel cores, sometimes files not writing properly
	# start_regen <- Sys.time()
	# for(loop in 1:length(rich_dir)){
	# 	copy_sim(fromdir=rich_dir[loop], fromcombos=rich_modcombos[loop,], todir=alt_dir, itervec=itervec, rewrite=FALSE, res_dir="base_12seasons")
	# }
	# end_regen <- Sys.time() - start_regen

	# ### ----- run LIME -----------###
	# lime_dirs <- base_12seasons_dir_vec[which(grepl("LBSPR",base_12seasons_dir_vec)==FALSE)]
	# lbspr_dirs <- base_12seasons_dir_vec[which(grepl("LBSPR",base_12seasons_dir_vec))]
	# lime_combos <- base_modcombos[which(grepl("LBSPR", base_modcombos[,"Data_avail"])==FALSE),]
	# lbspr_combos <- base_modcombos[which(grepl("LBSPR", base_modcombos[,"Data_avail"])),]	

	# lh_list <- adj_variation(SigmaR=0.737, SigmaF=0.2, SigmaC=0.2, SigmaI=0.2, CVlen=0.1, rho=0.426, selex_param=2, mat_param=1,  nseasons=1)

	# start_run <- Sys.time()
	# foreach(loop=1:length(lime_dirs), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=lime_combos[loop,"C_opt"], LFdist=1, param_adjust=c("SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(0.7,0.2,0.2,0.2),  fix_param=FALSE)
	# end_run <- Sys.time() - start_run


	# ### ----- run LBSPR -----------###	

	# start_run <- Sys.time()
	# foreach(loop=1:length(lbspr_dirs), .packages=c('LBSPR','LIME')) %dopar% run_LBSPR(modpath=lbspr_dirs[loop], lh=lh_list[[as.character(strsplit(lbspr_combos[loop,"LH"],"_")[[1]][2])]], itervec=itervec, species=NULL, rewrite=FALSE, simulation=TRUE)
	# end_altrun <- Sys.time() - start_run


############################################################
## monthly
############################################################
### ----- equilibrium test ----------- ###
### ----- models to run ----------- ###
lh_vec <- c("Short", "Medium", "Long")
Fdyn_vec <- "Constant"
Rdyn_vec <- "Constant"
data_vec <- c("LC10", "LC1")#, "LBSPR10", "LBSPR1")
SampleSize_vec <- 200 #c(1000,200)#,20)
itervec <- 1:25

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

	# for(loop in 1:length(lime_dirs)){
	# 	run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=lime_combos[loop,"C_opt"], LFdist=1, param_adjust=c("SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(0.7,0.2,0.2,0.2),  fix_param=FALSE)
	# }

	### ----- run LBSPR -----------###	

	# start_run <- Sys.time()
	# foreach(loop=1:length(lbspr_dirs), .packages=c('LBSPR','LIME')) %dopar% run_LBSPR(modpath=lbspr_dirs[loop], lh=lh_list[[as.character(strsplit(lbspr_combos[loop,"LH"],"_")[[1]][2])]], itervec=itervec, species=NULL, rewrite=FALSE, simulation=TRUE)
	# end_altrun <- Sys.time() - start_run
	# for(loop in 1:length(lbspr_dirs)){
	# 	run_LBSPR(modpath=lbspr_dirs[loop], lh=lh_list[[as.character(strsplit(lbspr_combos[loop,"LH"],"_")[[1]][2])]], itervec=itervec, species=NULL, rewrite=FALSE, simulation=TRUE)
	# }


################# base -- with 2-way F variability
### ----- models to run ----------- ###
lh_vec <- c("Short")#, "Medium", "Long")
Fdyn_vec <- "Ramp"
Rdyn_vec <- "AR"
data_vec <- c("LC10", "LC1", "LBSPR10", "LBSPR1")
SampleSize_vec <- 200 #c(1000,200)#,20)
itervec <- 1:100

base_modcombos <- expand.grid("Data_avail"=data_vec, "SampleSize"=paste0("SampleSize_", SampleSize_vec), "LH"=paste0("LH_",lh_vec), "Fdyn"=paste0("F_",Fdyn_vec), "Rdyn"=paste0("R_",Rdyn_vec), stringsAsFactors=FALSE)
base_modcombos$C_opt <- rep(0, nrow(base_modcombos))
	base_modcombos$C_opt[which(grepl("Catch",base_modcombos[,"Data_avail"]))] <- 2

	base_runMonthly_dir <- file.path(main_dir, "base_runMonthly")
	# unlink(base_runMonthly_dir, TRUE)
	dir.create(base_runMonthly_dir, showWarnings=FALSE)	

	## setup equilibrium dirs
	base_runMonthly_dir_vec <- model_paths(res_dir=base_runMonthly_dir, modcombos=base_modcombos[,-c(ncol(base_modcombos))])
	ignore <- sapply(1:length(base_runMonthly_dir_vec), function(m) sapply(itervec, function(x) dir.create(file.path(base_runMonthly_dir_vec[m], x), showWarnings=FALSE)))	

	## setup scenario life history list
	## 1-parameter model for selectivity and maturity
	lh_list <- adj_variation(SigmaR=0.737, SigmaF=0.2, SigmaC=0.2, SigmaI=0.2, CVlen=0.1, rho=0.426, selex_param=2, mat_param=1,  nseasons=12)
	# lh_fig(lh_list, save=TRUE)

	### ----- generate data -----------###
	### data rich cases only
	rich_dir <- base_runMonthly_dir_vec[grepl("LC10/",base_runMonthly_dir_vec)]
	rich_modcombos <- base_modcombos[which(base_modcombos$Data_avail=="LC10"),]
	alt_dir <- base_runMonthly_dir_vec[grepl("LC10/",base_runMonthly_dir_vec)==FALSE]
	alt_modcombos <- base_modcombos[which(base_modcombos$Data_avail!="LC10"),]

	### ----- use parallel cores -----------###
	# registerDoParallel(cores=ncores)	

	start_gen <- Sys.time()
	foreach(loop=1:length(rich_dir), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=rich_dir[loop], data_avail=as.character(rich_modcombos[loop,"Data_avail"]), itervec=itervec, Fdynamics=as.character(strsplit(rich_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(rich_modcombos[loop,"Rdyn"],"_")[[1]][2]),  lh=lh_list[[as.character(strsplit(rich_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=10, comp_sample=as.numeric(strsplit(rich_modcombos[loop,"SampleSize"],"_")[[1]][2]), rewrite=FALSE, init_depl=c(0.05,0.95), pool=FALSE)
	end_gen <- Sys.time() - start_gen

	# copy data rich cases to data poor
	# when written using parallel cores, sometimes files not writing properly
	start_regen <- Sys.time()
	for(loop in 1:length(rich_dir)){
		copy_sim(fromdir=rich_dir[loop], fromcombos=rich_modcombos[loop,], todir=alt_dir, itervec=itervec, rewrite=FALSE, res_dir="base_runMonthly")
	}
	end_regen <- Sys.time() - start_regen

	### ----- run LIME -----------###
	lime_dirs <- base_runMonthly_dir_vec[which(grepl("LBSPR",base_runMonthly_dir_vec)==FALSE)]
	lbspr_dirs <- base_runMonthly_dir_vec[which(grepl("LBSPR",base_runMonthly_dir_vec))]
	lime_combos <- base_modcombos[which(grepl("LBSPR", base_modcombos[,"Data_avail"])==FALSE),]
	lbspr_combos <- base_modcombos[which(grepl("LBSPR", base_modcombos[,"Data_avail"])),]	

	# start_run <- Sys.time()
	# foreach(loop=1:length(lime_dirs), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=lime_combos[loop,"C_opt"], LFdist=1, param_adjust=c("SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(0.7,0.2,0.2,0.2),  fix_param=FALSE)
	# end_run <- Sys.time() - start_run

	for(loop in 1:length(lime_dirs)){
		run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=lime_combos[loop,"C_opt"], LFdist=1, param_adjust=c("SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(0.7,0.2,0.2,0.2),  fix_param=FALSE)
	}


	### ----- run LBSPR -----------###	

	# start_run <- Sys.time()
	# foreach(loop=1:length(lbspr_dirs), .packages=c('LBSPR','LIME')) %dopar% run_LBSPR(modpath=lbspr_dirs[loop], lh=lh_list[[as.character(strsplit(lbspr_combos[loop,"LH"],"_")[[1]][2])]], itervec=itervec, species=NULL, rewrite=FALSE, simulation=TRUE)
	# end_altrun <- Sys.time() - start_run

	for(loop in 1:length(lbspr_dirs)){
		run_LBSPR(modpath=lbspr_dirs[loop], lh=lh_list[[as.character(strsplit(lbspr_combos[loop,"LH"],"_")[[1]][2])]], itervec=itervec, species=NULL, rewrite=FALSE, simulation=TRUE)
	}

############################################################
#### sensitivities - equilibrium
############################################################
### ----- models to run ----------- ###
lh_vec <- c("Short", "Medium", "Long")
Fdyn_vec <- "Constant"
Rdyn_vec <- "Constant" 
data_vec <- c("LC10")
SampleSize_vec <- 200
itervec <- 1:100
input_vec <- c("M", "linf", "vbk", "ML50", "CVlen")
val_vec <- c("high", "low")

sens_equil_modcombos <- expand.grid("Param"=input_vec, "Adjust"=val_vec, "Data_avail"=data_vec, "SampleSize"=paste0("SampleSize_", SampleSize_vec), "LH"=paste0("LH_",lh_vec), "Fdyn"=paste0("F_",Fdyn_vec), "Rdyn"=paste0("R_",Rdyn_vec), stringsAsFactors=FALSE)
sens_equil_modcombos$C_opt <- rep(0, nrow(sens_equil_modcombos))
	sens_equil_modcombos$C_opt[which(grepl("Catch",sens_equil_modcombos[,"Data_avail"]))] <- 2

### ----- sens runs ----------- ###
	sens_equil_dir <- file.path(main_dir, "sens_equil")
	# unlink(sens_equil_dir, TRUE)
	dir.create(sens_equil_dir, showWarnings=FALSE)	

	## setup dirs
	sens_equil_dir_vec <- model_paths(res_dir=sens_equil_dir, modcombos=sens_equil_modcombos[,-ncol(sens_equil_modcombos)])
	ignore <- sapply(1:length(sens_equil_dir_vec), function(m) sapply(itervec, function(x) dir.create(file.path(sens_equil_dir_vec[m], x), showWarnings=FALSE)))	

	## equilibrium
	lh_list <- adj_variation(SigmaR=0.01, SigmaF=0.01, SigmaC=0.01, SigmaI=0.01, CVlen=0.1, rho=0, selex_param=2, mat_param=1, nseasons=1)

	## matrix of sensitivity values
	sens_vals <- NULL
	for(ll in 1:length(lh_vec)){
		sens_vals[[lh_vec[ll]]] <- matrix(NA, nrow=length(val_vec), ncol=length(input_vec))
		rownames(sens_vals[[lh_vec[ll]]]) <- val_vec
		colnames(sens_vals[[lh_vec[ll]]]) <- input_vec
		for(i in 1:length(input_vec)){
			sens_vals[[lh_vec[ll]]]["low",input_vec[i]] <- lh_list[[lh_vec[ll]]][[input_vec[i]]] - 0.25*lh_list[[lh_vec[ll]]][[input_vec[i]]]
			sens_vals[[lh_vec[ll]]]["high",input_vec[i]] <- lh_list[[lh_vec[ll]]][[input_vec[i]]] + 0.25*lh_list[[lh_vec[ll]]][[input_vec[i]]]
		}
	}


	### ----- use parallel cores -----------###
	registerDoParallel(cores=5)	

	### ----- generate data -----------###
	## copy data files from similar previous directory
	init_dir <- equil_dir_vec[which(equil_modcombos$Data_avail %in% data_vec & equil_modcombos$SampleSize %in% paste0("SampleSize_", SampleSize_vec))]
	init_modcombos <- equil_modcombos[which(equil_modcombos$Data_avail %in% data_vec & equil_modcombos$SampleSize %in% paste0("SampleSize_", SampleSize_vec)),]
	alt_dir <- sens_equil_dir_vec
	alt_modcombos <- sens_equil_modcombos

	start_regen <- Sys.time()
	for(loop in 1:length(init_dir)){
		copy_sim(fromdir=init_dir[loop], fromcombos=init_modcombos[loop,], todir=alt_dir, itervec=itervec, rewrite=FALSE, res_dir="sens_equil", sensitivity=TRUE)
	}
	end_regen <- Sys.time() - start_regen

	## run LIME
	start_run <- Sys.time()
	foreach(loop=1:length(sens_equil_dir_vec), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=sens_equil_dir_vec[loop], lh=lh_list[[as.character(strsplit(sens_equil_modcombos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(sens_equil_modcombos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=sens_equil_modcombos[loop,"C_opt"], param_adjust=c(sens_equil_modcombos[loop,"Param"],"SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(sens_vals[[as.character(strsplit(sens_equil_modcombos[loop,"LH"],"_")[[1]][2])]][sens_equil_modcombos[loop,"Adjust"],sens_equil_modcombos[loop,"Param"]],0.7,0.2,0.2,0.2), LFdist=1)
	end_run <- Sys.time() - start_run

############################################################
#### sensitivities - base
############################################################
### ----- models to run ----------- ###
lh_vec <- c("Short", "Medium", "Long")
Fdyn_vec <- c("Ramp", "Increasing")#, "Decreasing")
Rdyn_vec <- "AR" 
data_vec <- c("LC10")
SampleSize_vec <- 200
itervec <- 1:100
input_vec <- c("M", "linf", "vbk", "ML50", "CVlen")
val_vec <- c("high", "low")

sens_base_modcombos <- expand.grid("Param"=input_vec, "Adjust"=val_vec, "Data_avail"=data_vec, "SampleSize"=paste0("SampleSize_", SampleSize_vec), "LH"=paste0("LH_",lh_vec), "Fdyn"=paste0("F_",Fdyn_vec), "Rdyn"=paste0("R_",Rdyn_vec), stringsAsFactors=FALSE)
sens_base_modcombos$C_opt <- rep(0, nrow(sens_base_modcombos))
	sens_base_modcombos$C_opt[which(grepl("Catch",sens_base_modcombos[,"Data_avail"]))] <- 2

### ----- sens runs ----------- ###
	sens_base_dir <- file.path(main_dir, "sens_base")
	# unlink(sens_base_dir, TRUE)
	dir.create(sens_base_dir, showWarnings=FALSE)	

	## setup dirs
	sens_base_dir_vec <- model_paths(res_dir=sens_base_dir, modcombos=sens_base_modcombos[,-ncol(sens_base_modcombos)])
	ignore <- sapply(1:length(sens_base_dir_vec), function(m) sapply(itervec, function(x) dir.create(file.path(sens_base_dir_vec[m], x), showWarnings=FALSE)))	

	## base
	lh_list <- adj_variation(SigmaR=0.737, SigmaF=0.2, SigmaC=0.2, SigmaI=0.2, CVlen=0.1, rho=0.426, selex_param=2, mat_param=1, nseasons=1)

	## matrix of sensitivity values
	sens_vals <- NULL
	for(ll in 1:length(lh_vec)){
		sens_vals[[lh_vec[ll]]] <- matrix(NA, nrow=length(val_vec), ncol=length(input_vec))
		rownames(sens_vals[[lh_vec[ll]]]) <- val_vec
		colnames(sens_vals[[lh_vec[ll]]]) <- input_vec
		for(i in 1:length(input_vec)){
			sens_vals[[lh_vec[ll]]]["low",input_vec[i]] <- lh_list[[lh_vec[ll]]][[input_vec[i]]] - 0.25*lh_list[[lh_vec[ll]]][[input_vec[i]]]
			sens_vals[[lh_vec[ll]]]["high",input_vec[i]] <- lh_list[[lh_vec[ll]]][[input_vec[i]]] + 0.25*lh_list[[lh_vec[ll]]][[input_vec[i]]]
		}
	}


	### ----- use parallel cores -----------###
	# registerDoParallel(cores=4)	

	### ----- generate data -----------###
	## copy data files from similar previous directory
	init_dir <- base_dir_vec[which(base_modcombos$Data_avail %in% data_vec & base_modcombos$SampleSize %in% paste0("SampleSize_", SampleSize_vec))]
	init_modcombos <- base_modcombos[which(base_modcombos$Data_avail %in% data_vec & base_modcombos$SampleSize %in% paste0("SampleSize_", SampleSize_vec)),]
	alt_dir <- sens_base_dir_vec
	alt_modcombos <- sens_base_modcombos

	start_regen <- Sys.time()
	for(loop in 1:length(init_dir)){
		copy_sim(fromdir=init_dir[loop], fromcombos=init_modcombos[loop,], todir=alt_dir, itervec=itervec, rewrite=FALSE, res_dir="sens_base", sensitivity=TRUE)
	}
	end_regen <- Sys.time() - start_regen

	## run LIME
	start_run <- Sys.time()
	foreach(loop=1:length(sens_base_dir_vec), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=sens_base_dir_vec[loop], lh=lh_list[[as.character(strsplit(sens_base_modcombos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(sens_base_modcombos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=sens_base_modcombos[loop,"C_opt"], param_adjust=c(sens_base_modcombos[loop,"Param"],"SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(sens_vals[[as.character(strsplit(sens_base_modcombos[loop,"LH"],"_")[[1]][2])]][sens_base_modcombos[loop,"Adjust"],sens_base_modcombos[loop,"Param"]],0.7,0.2,0.2,0.2), LFdist=1)
	end_run <- Sys.time() - start_run

############################################################
#### sensitivities - equil dome
############################################################
### ----- models to run ----------- ###
lh_vec <- c("Short", "Medium", "Long")
Fdyn_vec <- "Constant"
Rdyn_vec <- "Constant" 
data_vec <- c("LC10")
SampleSize_vec <- 200
itervec <- 1:100
val_vec <- c("high", "low")

dome_equil_modcombos <- expand.grid("Adjust"=val_vec, "Data_avail"=data_vec, "SampleSize"=paste0("SampleSize_", SampleSize_vec), "LH"=paste0("LH_",lh_vec), "Fdyn"=paste0("F_",Fdyn_vec), "Rdyn"=paste0("R_",Rdyn_vec), stringsAsFactors=FALSE)
dome_equil_modcombos$C_opt <- rep(0, nrow(dome_equil_modcombos))
	dome_equil_modcombos$C_opt[which(grepl("Catch",dome_equil_modcombos[,"Data_avail"]))] <- 2

### ----- sens runs ----------- ###
	dome_equil_dir <- file.path(main_dir, "dome_equil")
	# unlink(dome_equil_dir, TRUE)
	dir.create(dome_equil_dir, showWarnings=FALSE)	

	## setup dirs
	dome_equil_dir_vec <- model_paths(res_dir=dome_equil_dir, modcombos=dome_equil_modcombos[,-ncol(dome_equil_modcombos)])
	ignore <- sapply(1:length(dome_equil_dir_vec), function(m) sapply(itervec, function(x) dir.create(file.path(dome_equil_dir_vec[m], x), showWarnings=FALSE)))	

	## equilibrium
	lh_list <- adj_variation(SigmaR=0.01, SigmaF=0.01, SigmaC=0.01, SigmaI=0.01, CVlen=0.1, rho=0, selex_param=2, mat_param=1)

	## matrix of sensitivity values
	dome_vals <- NULL
	dome_vals[["Short"]] <- c(24, 50)
	names(dome_vals[["Short"]]) <- val_vec
	dome_vals[["Medium"]] <- c(37, 77)
	names(dome_vals[["Medium"]]) <- val_vec
	dome_vals[["Long"]] <- c(62, 140)
	names(dome_vals[["Long"]]) <- val_vec

	dome_low_dirs <- dome_equil_dir_vec[grepl("low", dome_equil_dir_vec)]
	dome_low_modcombos <- dome_equil_modcombos[which(dome_equil_modcombos$Adjust=="low"),]
	dome_high_dirs <- dome_equil_dir_vec[grepl("high", dome_equil_dir_vec)]
	dome_high_modcombos <- dome_equil_modcombos[which(dome_equil_modcombos$Adjust=="high"),]

	lh_dome_low <- NULL
	for(ll in 1:length(lh_list)){
		lh_dome_low[[lh_vec[ll]]] <- with(lh_list[[lh_vec[ll]]], create_lh_list(vbk=vbk, linf=linf, lwa=lwa, lwb=lwb, S50=SL50, M50=ML50, M=M, selex_input=selex_input, maturity_input=maturity_input, selex_type="dome", dome_sd=dome_vals[[lh_vec[ll]]]["low"], CVlen=CVlen, SigmaC=SigmaC, SigmaI=SigmaI, SigmaR=SigmaR, SigmaF=SigmaF, R0=R0, h=h, qcoef=qcoef, F1=F1, start_ages=0, rho=rho, theta=theta, AgeMax=AgeMax))
	}
	lh_dome_high <- NULL
	for(ll in 1:length(lh_list)){
		lh_dome_high[[lh_vec[ll]]] <- with(lh_list[[lh_vec[ll]]], create_lh_list(vbk=vbk, linf=linf, lwa=lwa, lwb=lwb, S50=SL50, M50=ML50, M=M, selex_input=selex_input, maturity_input=maturity_input, selex_type="dome", dome_sd=dome_vals[[lh_vec[ll]]]["high"], CVlen=CVlen, SigmaC=SigmaC, SigmaI=SigmaI, SigmaR=SigmaR, SigmaF=SigmaF, R0=R0, h=h, qcoef=qcoef, F1=F1, start_ages=0, rho=rho, theta=theta, AgeMax=AgeMax))
	}

	### ----- use parallel cores -----------###
	# registerDoParallel(cores=ncores)	

	start_gen <- Sys.time()
	foreach(loop=1:length(dome_low_dirs), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=dome_low_dirs[loop], data_avail=as.character(dome_low_modcombos[loop,"Data_avail"]), itervec=itervec, Fdynamics=as.character(strsplit(dome_low_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(dome_low_modcombos[loop,"Rdyn"],"_")[[1]][2]),  lh=lh_dome_low[[as.character(strsplit(dome_low_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=20, comp_sample=as.numeric(strsplit(dome_low_modcombos[loop,"SampleSize"],"_")[[1]][2]), rewrite=FALSE, init_depl=c(0.05,0.95))
	end_gen <- Sys.time() - start_gen

	start_gen <- Sys.time()
	foreach(loop=1:length(dome_high_dirs), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=dome_high_dirs[loop], data_avail=as.character(dome_high_modcombos[loop,"Data_avail"]), itervec=itervec, Fdynamics=as.character(strsplit(dome_high_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(dome_high_modcombos[loop,"Rdyn"],"_")[[1]][2]),  lh=lh_dome_high[[as.character(strsplit(dome_high_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=20, comp_sample=as.numeric(strsplit(dome_high_modcombos[loop,"SampleSize"],"_")[[1]][2]), rewrite=FALSE, init_depl=c(0.05,0.95))
	end_gen <- Sys.time() - start_gen

	## run LIME
	start_run <- Sys.time()
	foreach(loop=1:length(dome_equil_dir_vec), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=dome_equil_dir_vec[loop], lh=lh_list[[as.character(strsplit(dome_equil_modcombos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(dome_equil_modcombos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=dome_equil_modcombos[loop,"C_opt"], param_adjust=c("SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(0.7,0.2,0.2,0.2), LFdist=1)
	end_run <- Sys.time() - start_run

############################################################
#### sensitivities - variation dome
############################################################
### ----- models to run ----------- ###
lh_vec <- c("Short", "Medium", "Long")
Fdyn_vec <- c("Ramp", "Increasing")
Rdyn_vec <- "AR" 
data_vec <- c("LC10")
SampleSize_vec <- 200
itervec <- 1:100
val_vec <- c("high", "low")

dome_base_modcombos <- expand.grid("Adjust"=val_vec, "Data_avail"=data_vec, "SampleSize"=paste0("SampleSize_", SampleSize_vec), "LH"=paste0("LH_",lh_vec), "Fdyn"=paste0("F_",Fdyn_vec), "Rdyn"=paste0("R_",Rdyn_vec), stringsAsFactors=FALSE)
dome_base_modcombos$C_opt <- rep(0, nrow(dome_base_modcombos))
	dome_base_modcombos$C_opt[which(grepl("Catch",dome_base_modcombos[,"Data_avail"]))] <- 2

### ----- sens runs ----------- ###
	dome_base_dir <- file.path(main_dir, "dome_base")
	# unlink(dome_base_dir, TRUE)
	dir.create(dome_base_dir, showWarnings=FALSE)	

	## setup dirs
	dome_base_dir_vec <- model_paths(res_dir=dome_base_dir, modcombos=dome_base_modcombos[,-ncol(dome_base_modcombos)])
	ignore <- sapply(1:length(dome_base_dir_vec), function(m) sapply(itervec, function(x) dir.create(file.path(dome_base_dir_vec[m], x), showWarnings=FALSE)))	

	## variation
	lh_list <- adj_variation(SigmaR=0.737, SigmaF=0.2, SigmaC=0.2, SigmaI=0.2, CVlen=0.1, rho=0.426, selex_param=2, mat_param=1)

	## matrix of sensitivity values
	dome_vals <- NULL
	dome_vals[["Short"]] <- c(24, 50)
	names(dome_vals[["Short"]]) <- val_vec
	dome_vals[["Medium"]] <- c(37, 77)
	names(dome_vals[["Medium"]]) <- val_vec
	dome_vals[["Long"]] <- c(62, 140)
	names(dome_vals[["Long"]]) <- val_vec

	dome_low_dirs <- dome_base_dir_vec[grepl("low", dome_base_dir_vec)]
	dome_low_modcombos <- dome_base_modcombos[which(dome_base_modcombos$Adjust=="low"),]
	dome_high_dirs <- dome_base_dir_vec[grepl("high", dome_base_dir_vec)]
	dome_high_modcombos <- dome_base_modcombos[which(dome_base_modcombos$Adjust=="high"),]

	lh_dome_low <- NULL
	for(ll in 1:length(lh_list)){
		lh_dome_low[[lh_vec[ll]]] <- with(lh_list[[lh_vec[ll]]], create_lh_list(vbk=vbk, linf=linf, lwa=lwa, lwb=lwb, S50=SL50, M50=ML50, M=M, selex_input="length", maturity_input="length", selex_type="dome", dome_sd=dome_vals[[lh_vec[ll]]]["low"]), CVlen=CVlen, SigmaC=SigmaC, SigmaI=SigmaI, SigmaR=SigmaR, SigmaF=SigmaF, R0=R0, h=h, qcoef=qcoef, F1=F1, start_ages=start_ages, rho=rho, theta=theta)
	}
	lh_dome_high <- NULL
	for(ll in 1:length(lh_list)){
		lh_dome_high[[lh_vec[ll]]] <- with(lh_list[[lh_vec[ll]]], create_lh_list(vbk=vbk, linf=linf, lwa=lwa, lwb=lwb, S50=SL50, M50=ML50, M=M, selex_input="length", maturity_input="length", selex_type="dome", dome_sd=dome_vals[[lh_vec[ll]]]["high"]), CVlen=CVlen, SigmaC=SigmaC, SigmaI=SigmaI, SigmaR=SigmaR, SigmaF=SigmaF, R0=R0, h=h, qcoef=qcoef, F1=F1, start_ages=start_ages, rho=rho, theta=theta)
	}

	plot(lh_dome_low[[3]]$S_l, ylim=c(0,1))
	lh_dome_low[[3]]$S_l

	### ----- use parallel cores -----------###
	# registerDoParallel(cores=ncores)	

	start_gen <- Sys.time()
	foreach(loop=1:length(dome_low_dirs), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=dome_low_dirs[loop], data_avail=as.character(dome_low_modcombos[loop,"Data_avail"]), itervec=itervec, Fdynamics=as.character(strsplit(dome_low_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(dome_low_modcombos[loop,"Rdyn"],"_")[[1]][2]),  lh=lh_dome_low[[as.character(strsplit(dome_low_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=20, comp_sample=as.numeric(strsplit(dome_low_modcombos[loop,"SampleSize"],"_")[[1]][2]), rewrite=FALSE, init_depl=c(0.05,0.95))
	end_gen <- Sys.time() - start_gen

	start_gen <- Sys.time()
	foreach(loop=1:length(dome_high_dirs), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=dome_high_dirs[loop], data_avail=as.character(dome_high_modcombos[loop,"Data_avail"]), itervec=itervec, Fdynamics=as.character(strsplit(dome_high_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(dome_high_modcombos[loop,"Rdyn"],"_")[[1]][2]),  lh=lh_dome_high[[as.character(strsplit(dome_high_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=20, comp_sample=as.numeric(strsplit(dome_high_modcombos[loop,"SampleSize"],"_")[[1]][2]), rewrite=FALSE, init_depl=c(0.05,0.95))
	end_gen <- Sys.time() - start_gen

	## run LIME
	start_run <- Sys.time()
	foreach(loop=1:length(dome_base_dir_vec), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=dome_base_dir_vec[loop], lh=lh_list[[as.character(strsplit(dome_base_modcombos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(dome_base_modcombos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=dome_base_modcombos[loop,"C_opt"], param_adjust=c("SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(0.7,0.2,0.2,0.2), LFdist=1)
	end_run <- Sys.time() - start_run


