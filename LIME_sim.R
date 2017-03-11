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

############################################################
#### equilibrium
############################################################

### ----- models to run ----------- ###
lh_vec <- c("Short", "Medium", "Long")
Fdyn_vec <- "Constant"
Rdyn_vec <- "Constant"
data_vec <- c("Index_Catch_LC20", "Index_LC10", "Index_LC1", "Catch_LC10", "Catch_LC1", "LC10", "LC1", "LBSPR10", "LBSPR1")
SampleSize_vec <- c(1000, 500, 200, 50, 20)
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
	lh_list <- adj_variation(SigmaR=0.01, SigmaF=0.01, SigmaC=0.01, SigmaI=0.01, CVlen=0.1, rho=0, selex_param=1, mat_param=1)
	lh_fig(lh_list, save=TRUE)

	### ----- generate data -----------###
	### data rich cases only
	rich_dir <- equil_dir_vec[grepl("Index_Catch_LC20",equil_dir_vec)]
	rich_modcombos <- equil_modcombos[which(equil_modcombos$Data_avail=="Index_Catch_LC20"),]
	alt_dir <- equil_dir_vec[grepl("Index_Catch_LC20",equil_dir_vec)==FALSE]
	alt_modcombos <- equil_modcombos[which(equil_modcombos$Data_avail!="Index_Catch_LC20"),]

	### ----- use parallel cores -----------###
	registerDoParallel(cores=ncores)	

	start_gen <- Sys.time()
	foreach(loop=1:length(rich_dir), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=rich_dir[loop], data_avail=as.character(rich_modcombos[loop,"Data_avail"]), itervec=itervec, Fdynamics=as.character(strsplit(rich_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(rich_modcombos[loop,"Rdyn"],"_")[[1]][2]), write=TRUE, lh=lh_list[[as.character(strsplit(rich_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=20, comp_sample=as.numeric(strsplit(rich_modcombos[loop,"SampleSize"],"_")[[1]][2]), rewrite=FALSE, init_depl=c(0.05,0.95))
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
	foreach(loop=1:length(lime_dirs), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, REML=FALSE, f_true=FALSE, C_opt=lime_combos[loop,"C_opt"], LFdist=1, param_adjust=c("SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(0.7,0.2,0.2,0.2), write=TRUE, fix_param=FALSE)
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
Fdyn_vec <- c("Ramp", "Increasing", "Decreasing")
Rdyn_vec <- c("AR")
data_vec <- c("Index_Catch_LC20", "Index_LC10", "Index_LC1", "Catch_LC10", "Catch_LC1", "LC10", "LC1", "LBSPR10", "LBSPR1")
SampleSize_vec <- c(1000, 500, 200, 50, 20)
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
	lh_list <- adj_variation(SigmaR=0.737, SigmaF=0.2, SigmaC=0.2, SigmaI=0.2, CVlen=0.1, rho=0.426, selex_param=2, mat_param=2)
	# lh_fig(lh=lh_list, save=FALSE)

	### ----- generate data -----------###
	### data rich cases only
	rich_dir <- base_dir_vec[grepl("Index_Catch_LC20",base_dir_vec)]
	rich_modcombos <- base_modcombos[which(base_modcombos$Data_avail=="Index_Catch_LC20"),]
	alt_dir <- base_dir_vec[grepl("Index_Catch_LC20",base_dir_vec)==FALSE]
	alt_modcombos <- base_modcombos[which(base_modcombos$Data_avail!="Index_Catch_LC20"),]

	### ----- use parallel cores -----------###
	registerDoParallel(cores=ncores)	

	start_gen <- Sys.time()
	foreach(loop=1:length(rich_dir), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=rich_dir[loop], data_avail=as.character(rich_modcombos[loop,"Data_avail"]), itervec=itervec, Fdynamics=as.character(strsplit(rich_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(rich_modcombos[loop,"Rdyn"],"_")[[1]][2]), write=TRUE, lh=lh_list[[as.character(strsplit(rich_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=20, comp_sample=as.numeric(strsplit(rich_modcombos[loop,"SampleSize"],"_")[[1]][2]), rewrite=FALSE, init_depl=c(0.05,0.95))
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
	foreach(loop=1:length(lime_dirs), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, REML=FALSE, f_true=FALSE, C_opt=lime_combos[loop,"C_opt"], LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.7,0.2,0.2), write=TRUE, fix_param=FALSE)
	end_run <- Sys.time() - start_run


	### ----- run LBSPR -----------###	

	start_run <- Sys.time()
	foreach(loop=1:length(lbspr_dirs), .packages=c('LBSPR','LIME')) %dopar% run_LBSPR(modpath=lbspr_dirs[loop], lh=lh_list[[as.character(strsplit(lbspr_combos[loop,"LH"],"_")[[1]][2])]], itervec=itervec, species=NULL, rewrite=FALSE, simulation=TRUE)
	end_altrun <- Sys.time() - start_run	

	# ### ------- check output ------###	
	# dir_esslow <- base_dir_vec[grepl("SampleSize_100/", base_dir_vec)]
	# dir_esshigh <- base_dir_vec[grepl("SampleSize_1000", base_dir_vec)]
	
	# base_coverlow <- interval_coverage(dirs=dir_esslow[which(grepl("LBSPR",dir_esslow)==FALSE)], itervec=itervec)
	# base_coverhigh <- interval_coverage(dirs=dir_esshigh[which(grepl("LBSPR",dir_esshigh)==FALSE)], itervec=itervec)


	# png(file.path(fig_dir, "Base_SampleSizelow_RE.png"), height=10, width=25, res=200, units="in")
	# base_res <- plot_re(dirs=dir_esslow, modcombos=base_modcombos[which(base_modcombos$SampleSize=="SampleSize_100"),], itervec=itervec, compareToLH=paste0("F_", Fdyn_vec), lh_vec=lh_vec, cover=base_coverlow$cover, ylim=c(-1,1))
	# dev.off()

	# png(file.path(fig_dir, "Base_SampleSizehigh_RE.png"), height=10, width=25, res=200, units="in")
	# base_res <- plot_re(dirs=dir_esshigh, modcombos=base_modcombos[which(base_modcombos$SampleSize=="SampleSize_1000"),], itervec=itervec, compareToLH=paste0("F_", Fdyn_vec), lh_vec=lh_vec, cover=base_coverhigh$cover, ylim=c(-1,1))
	# dev.off()


# ### ----- base runs - sigR lower----------- ###
# 	base_lowsigR_dir <- file.path(main_dir, "base_lowsigR")
# 	# unlink(base_lowsigR_dir, TRUE)
# 	dir.create(base_lowsigR_dir, showWarnings=FALSE)	

# 	## setup base dirs
# 	base_lowsigR_dir_vec <- model_paths(res_dir=base_lowsigR_dir, modcombos=base_modcombos[,-c(ncol(base_modcombos))])
# 	ignore <- sapply(1:length(base_lowsigR_dir_vec), function(m) sapply(itervec, function(x) dir.create(file.path(base_lowsigR_dir_vec[m], x), showWarnings=FALSE)))	

# 	## setup scenario life history list
# 	lh_list <- adj_variation(SigmaR=0.384, SigmaF=0.2, SigmaC=0.2, SigmaI=0.2, CVlen=0.1, rho=0.426, selex_param=2, mat_param=2)
# 	# lh_fig(lh=lh_list, save=FALSE)

# 	### ----- generate data -----------###
# 	### data rich cases only
# 	rich_dir <- base_lowsigR_dir_vec[grepl("Index_Catch_LC20",base_lowsigR_dir_vec)]
# 	rich_modcombos <- base_modcombos[which(base_modcombos$Data_avail=="Index_Catch_LC20"),]
# 	alt_dir <- base_lowsigR_dir_vec[grepl("Index_Catch_LC20",base_lowsigR_dir_vec)==FALSE]
# 	alt_modcombos <- base_modcombos[which(base_modcombos$Data_avail!="Index_Catch_LC20"),]

# 	### ----- use parallel cores -----------###
# 	registerDoParallel(cores=ncores)	

# 	start_gen <- Sys.time()
# 	foreach(loop=1:length(rich_dir), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=rich_dir[loop], data_avail=as.character(rich_modcombos[loop,"Data_avail"]), itervec=itervec, Fdynamics=as.character(strsplit(rich_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(rich_modcombos[loop,"Rdyn"],"_")[[1]][2]), write=TRUE, lh=lh_list[[as.character(strsplit(rich_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=20, comp_sample=as.numeric(strsplit(rich_modcombos[loop,"SampleSize"],"_")[[1]][2]), rewrite=FALSE, init_depl=c(0.05,0.95))
# 	end_gen <- Sys.time() - start_gen

# 	# copy data rich cases to data poor
# 	# when written using parallel cores, sometimes files not writing properly
# 	start_regen <- Sys.time()
# 	for(loop in 1:length(rich_dir)){
# 		copy_sim(fromdir=rich_dir[loop], fromcombos=rich_modcombos[loop,], todir=alt_dir, itervec=itervec, rewrite=FALSE, res_dir="base_lowsigR")
# 	}
# 	end_regen <- Sys.time() - start_regen


# 	### ----- run LIME -----------###
# 	lime_dirs <- base_lowsigR_dir_vec[which(grepl("LBSPR",base_lowsigR_dir_vec)==FALSE)]
# 	lbspr_dirs <- base_lowsigR_dir_vec[which(grepl("LBSPR",base_lowsigR_dir_vec))]
# 	lime_combos <- base_modcombos[which(grepl("LBSPR", base_modcombos[,"Data_avail"])==FALSE),]
# 	lbspr_combos <- base_modcombos[which(grepl("LBSPR", base_modcombos[,"Data_avail"])),]	

# 	start_run <- Sys.time()
# 	foreach(loop=1:length(lime_dirs), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, REML=FALSE, f_true=FALSE, C_opt=lime_combos[loop,"C_opt"], LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.7,0.2,0.2), write=TRUE, fix_param=FALSE)
# 	end_run <- Sys.time() - start_run

# 	### ----- run LBSPR -----------###	

# 	start_run <- Sys.time()
# 	foreach(loop=1:length(lbspr_dirs), .packages=c('LBSPR','LIME')) %dopar% run_LBSPR(modpath=lbspr_dirs[loop], lh=lh_list[[as.character(strsplit(lbspr_combos[loop,"LH"],"_")[[1]][2])]], itervec=itervec, species=NULL, rewrite=FALSE, simulation=TRUE)
# 	end_altrun <- Sys.time() - start_run	

	### ------- check output ------###	
	# dir_esslow <- base_lowsigR_dir_vec[grepl("SampleSize_100/", base_lowsigR_dir_vec)]
	# dir_esshigh <- base_lowsigR_dir_vec[grepl("SampleSize_1000", base_lowsigR_dir_vec)]
	
	# base_coverlow <- interval_coverage(dirs=dir_esslow[which(grepl("LBSPR",dir_esslow)==FALSE)], itervec=itervec)
	# base_coverhigh <- interval_coverage(dirs=dir_esshigh[which(grepl("LBSPR",dir_esshigh)==FALSE)], itervec=itervec)


	# png(file.path(fig_dir, "Base_lowsigR_SampleSizelow_RE.png"), height=10, width=25, res=200, units="in")
	# base_res <- plot_re(dirs=dir_esslow, modcombos=base_modcombos[which(base_modcombos$SampleSize=="SampleSize_100"),], itervec=itervec, compareToLH=paste0("F_", Fdyn_vec), lh_vec=lh_vec, cover=base_coverlow$cover, ylim=c(-1,1))
	# dev.off()

	# png(file.path(fig_dir, "Base_lowsigR_SampleSizehigh_RE.png"), height=10, width=25, res=200, units="in")
	# base_res <- plot_re(dirs=dir_esshigh, modcombos=base_modcombos[which(base_modcombos$SampleSize=="SampleSize_1000"),], itervec=itervec, compareToLH=paste0("F_", Fdyn_vec), lh_vec=lh_vec, cover=base_coverhigh$cover, ylim=c(-1,1))
	# dev.off()



############################################################
#### sensitivities - equilibrium
############################################################
### ----- models to run ----------- ###
lh_vec <- c("Short", "Medium", "Long")
Fdyn_vec <- "Constant"
Rdyn_vec <- "Constant" 
data_vec <- c("LC10")
SampleSize_vec <- c(1000, 500, 200, 50, 20)
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
	lh_list <- adj_variation(SigmaR=0.01, SigmaF=0.01, SigmaC=0.01, SigmaI=0.01, CVlen=0.1, rho=0, selex_param=2, mat_param=2)

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
	registerDoParallel(cores=ncores)	

	### ----- generate data -----------###
	## copy data files from similar previous directory
	init_dir <- equil_dir_vec[grepl("/LC10/",equil_dir_vec)]
	init_modcombos <- equil_modcombos[which(equil_modcombos$Data_avail=="LC10"),]
	alt_dir <- sens_equil_dir_vec
	alt_modcombos <- sens_equil_modcombos

	start_regen <- Sys.time()
	for(loop in 1:length(init_dir)){
		copy_sim(fromdir=init_dir[loop], fromcombos=init_modcombos[loop,], todir=alt_dir, itervec=itervec, rewrite=FALSE, res_dir="sens_equil", sensitivity=TRUE)
	}
	end_regen <- Sys.time() - start_regen

	## run LIME
	start_run <- Sys.time()
	foreach(loop=1:length(sens_equil_dir_vec), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=sens_equil_dir_vec[loop], lh=lh_list[[as.character(strsplit(sens_equil_modcombos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(sens_equil_modcombos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, REML=FALSE, f_true=FALSE, C_opt=sens_equil_modcombos[loop,"C_opt"], param_adjust=c(sens_equil_modcombos[loop,"Param"],"SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(sens_vals[[as.character(strsplit(sens_equil_modcombos[loop,"LH"],"_")[[1]][2])]][sens_equil_modcombos[loop,"Adjust"],sens_equil_modcombos[loop,"Param"]],0.7,0.2,0.2,0.2), LFdist=1)
	end_run <- Sys.time() - start_run

############################################################
#### sensitivities - base
############################################################
### ----- models to run ----------- ###
lh_vec <- c("Short", "Medium", "Long")
Fdyn_vec <- c("Ramp", "Decreasing", "Increasing")
Rdyn_vec <- "AR" 
data_vec <- c("LC10")
SampleSize_vec <- c(1000, 500, 200, 50, 20)
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
	lh_list <- adj_variation(SigmaR=0.737, SigmaF=0.2, SigmaC=0.2, SigmaI=0.2, CVlen=0.1, rho=0.426, selex_param=2, mat_param=2)

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
	registerDoParallel(cores=ncores)	

	### ----- generate data -----------###
	## copy data files from similar previous directory
	init_dir <- base_dir_vec[grepl("/LC10/",base_dir_vec) & grepl("SampleSize_20",base_dir_vec)==FALSE]
	init_modcombos <- base_modcombos[which(base_modcombos$Data_avail=="LC10" & base_modcombos$SampleSize!="SampleSize_20"),]
	alt_dir <- sens_base_dir_vec
	alt_modcombos <- sens_base_modcombos

	start_regen <- Sys.time()
	for(loop in 1:length(init_dir)){
		copy_sim(fromdir=init_dir[loop], fromcombos=init_modcombos[loop,], todir=alt_dir, itervec=itervec, rewrite=FALSE, res_dir="sens_equil", sensitivity=TRUE)
	}
	end_regen <- Sys.time() - start_regen

	## run LIME
	start_run <- Sys.time()
	foreach(loop=1:length(sens_base_dir_vec), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=sens_base_dir_vec[loop], lh=lh_list[[as.character(strsplit(sens_base_modcombos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(sens_base_modcombos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, REML=FALSE, f_true=FALSE, C_opt=sens_base_modcombos[loop,"C_opt"], param_adjust=c(sens_base_modcombos[loop,"Param"],"SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(sens_vals[[as.character(strsplit(sens_base_modcombos[loop,"LH"],"_")[[1]][2])]][sens_base_modcombos[loop,"Adjust"],sens_base_modcombos[loop,"Param"]],0.7,0.2,0.2,0.2), LFdist=1)
	end_run <- Sys.time() - start_run


############################################################
#### sensitivities - equil dome
############################################################
### ----- models to run ----------- ###
lh_vec <- c("Short", "Medium", "Long")
Fdyn_vec <- "Constant"
Rdyn_vec <- "Constant" 
data_vec <- c("LC10")
SampleSize_vec <- c(1000, 500, 200, 50, 20)
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
	lh_list <- adj_variation(SigmaR=0.01, SigmaF=0.01, SigmaC=0.01, SigmaI=0.01, CVlen=0.1, rho=0, selex_param=2, mat_param=2)

	## matrix of sensitivity values
	dome_vals <- NULL
	dome_vals[["Short"]] <- c(20, 40)
	names(dome_vals[["Short"]]) <- val_vec
	dome_vals[["Medium"]] <- c(2, 5)
	names(dome_vals[["Medium"]]) <- val_vec
	dome_vals[["Long"]] <- c(15, 30)
	names(dome_vals[["Long"]]) <- val_vec

	dome_low_dirs <- dome_equil_dir_vec[grepl("low", dome_equil_dir_vec)]
	dome_low_modcombos <- dome_equil_modcombos[which(dome_equil_modcombos$Adjust=="low"),]
	dome_high_dirs <- dome_equil_dir_vec[grepl("high", dome_equil_dir_vec)]
	dome_high_modcombos <- dome_equil_modcombos[which(dome_equil_modcombos$Adjust=="high"),]

	lh_dome_low <- NULL
	for(ll in 1:length(lh_list)){
		lh_dome_low[[lh_vec[ll]]] <- with(lh_list[[lh_vec[ll]]], create_lh_list(vbk=vbk, linf=linf, lwa=lwa, lwb=lwb, S50=S50, M50=M50, M=M, selex_input=selex_input, maturity_input=maturity_input, selex_type="dome", dome_sd=dome_vals[[lh_vec[ll]]]["low"]), CVlen=CVlen, SigmaC=SigmaC, SigmaI=SigmaI, SigmaR=SigmaR, SigmaF=SigmaF, R0=R0, h=h, qcoef=qcoef, F1=F1, start_ages=start_ages, rho=rho, theta=theta)
	}
	lh_dome_high <- NULL
	for(ll in 1:length(lh_list)){
		lh_dome_high[[lh_vec[ll]]] <- with(lh_list[[lh_vec[ll]]], create_lh_list(vbk=vbk, linf=linf, lwa=lwa, lwb=lwb, S50=S50, M50=M50, M=M, selex_input=selex_input, maturity_input=maturity_input, selex_type="dome", dome_sd=dome_vals[[lh_vec[ll]]]["high"]), CVlen=CVlen, SigmaC=SigmaC, SigmaI=SigmaI, SigmaR=SigmaR, SigmaF=SigmaF, R0=R0, h=h, qcoef=qcoef, F1=F1, start_ages=start_ages, rho=rho, theta=theta)
	}


	### ----- use parallel cores -----------###
	registerDoParallel(cores=ncores)	

	start_gen <- Sys.time()
	foreach(loop=1:length(dome_low_dirs), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=dome_low_dirs[loop], data_avail=as.character(dome_low_modcombos[loop,"Data_avail"]), itervec=itervec, Fdynamics=as.character(strsplit(dome_low_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(dome_low_modcombos[loop,"Rdyn"],"_")[[1]][2]), write=TRUE, lh=lh_dome_low[[as.character(strsplit(dome_low_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=20, comp_sample=as.numeric(strsplit(dome_low_modcombos[loop,"SampleSize"],"_")[[1]][2]), rewrite=FALSE, init_depl=c(0.05,0.95))
	end_gen <- Sys.time() - start_gen

	start_gen <- Sys.time()
	foreach(loop=1:length(dome_high_dirs), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=dome_high_dirs[loop], data_avail=as.character(dome_high_modcombos[loop,"Data_avail"]), itervec=itervec, Fdynamics=as.character(strsplit(dome_high_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(dome_high_modcombos[loop,"Rdyn"],"_")[[1]][2]), write=TRUE, lh=lh_dome_high[[as.character(strsplit(dome_high_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=20, comp_sample=as.numeric(strsplit(dome_high_modcombos[loop,"SampleSize"],"_")[[1]][2]), rewrite=FALSE, init_depl=c(0.05,0.95))
	end_gen <- Sys.time() - start_gen

	## run LIME
	start_run <- Sys.time()
	foreach(loop=1:length(dome_equil_dir_vec), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=dome_equil_dir_vec[loop], lh=lh_list[[as.character(strsplit(dome_equil_modcombos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(dome_equil_modcombos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, REML=FALSE, f_true=FALSE, C_opt=dome_equil_modcombos[loop,"C_opt"], param_adjust=c("SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(0.7,0.2,0.2,0.2), LFdist=1)
	end_run <- Sys.time() - start_run

############################################################
#### sensitivities - variation dome
############################################################
### ----- models to run ----------- ###
lh_vec <- c("Short", "Medium", "Long")
Fdyn_vec <- c("Ramp", "Decreasing", "Increasing")
Rdyn_vec <- "AR" 
data_vec <- c("LC10")
SampleSize_vec <- c(1000, 500, 200, 50, 20)
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
	lh_list <- adj_variation(SigmaR=0.737, SigmaF=0.2, SigmaC=0.2, SigmaI=0.2, CVlen=0.1, rho=0.426, selex_param=2, mat_param=2)

	## matrix of sensitivity values
	dome_vals <- NULL
	dome_vals[["Short"]] <- c(20, 40)
	names(dome_vals[["Short"]]) <- val_vec
	dome_vals[["Medium"]] <- c(2, 5)
	names(dome_vals[["Medium"]]) <- val_vec
	dome_vals[["Long"]] <- c(15, 30)
	names(dome_vals[["Long"]]) <- val_vec

	dome_low_dirs <- dome_base_dir_vec[grepl("low", dome_base_dir_vec)]
	dome_low_modcombos <- dome_base_modcombos[which(dome_base_modcombos$Adjust=="low"),]
	dome_high_dirs <- dome_base_dir_vec[grepl("high", dome_base_dir_vec)]
	dome_high_modcombos <- dome_base_modcombos[which(dome_base_modcombos$Adjust=="high"),]

	lh_dome_low <- NULL
	for(ll in 1:length(lh_list)){
		lh_dome_low[[lh_vec[ll]]] <- with(lh_list[[lh_vec[ll]]], create_lh_list(vbk=vbk, linf=linf, lwa=lwa, lwb=lwb, S50=S50, M50=M50, M=M, selex_input=selex_input, maturity_input=maturity_input, selex_type="dome", dome_sd=dome_vals[[lh_vec[ll]]]["low"]), CVlen=CVlen, SigmaC=SigmaC, SigmaI=SigmaI, SigmaR=SigmaR, SigmaF=SigmaF, R0=R0, h=h, qcoef=qcoef, F1=F1, start_ages=start_ages, rho=rho, theta=theta)
	}
	lh_dome_high <- NULL
	for(ll in 1:length(lh_list)){
		lh_dome_high[[lh_vec[ll]]] <- with(lh_list[[lh_vec[ll]]], create_lh_list(vbk=vbk, linf=linf, lwa=lwa, lwb=lwb, S50=S50, M50=M50, M=M, selex_input=selex_input, maturity_input=maturity_input, selex_type="dome", dome_sd=dome_vals[[lh_vec[ll]]]["high"]), CVlen=CVlen, SigmaC=SigmaC, SigmaI=SigmaI, SigmaR=SigmaR, SigmaF=SigmaF, R0=R0, h=h, qcoef=qcoef, F1=F1, start_ages=start_ages, rho=rho, theta=theta)
	}


	### ----- use parallel cores -----------###
	registerDoParallel(cores=ncores)	

	start_gen <- Sys.time()
	foreach(loop=1:length(dome_low_dirs), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=dome_low_dirs[loop], data_avail=as.character(dome_low_modcombos[loop,"Data_avail"]), itervec=itervec, Fdynamics=as.character(strsplit(dome_low_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(dome_low_modcombos[loop,"Rdyn"],"_")[[1]][2]), write=TRUE, lh=lh_dome_low[[as.character(strsplit(dome_low_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=20, comp_sample=as.numeric(strsplit(dome_low_modcombos[loop,"SampleSize"],"_")[[1]][2]), rewrite=FALSE, init_depl=c(0.05,0.95))
	end_gen <- Sys.time() - start_gen

	start_gen <- Sys.time()
	foreach(loop=1:length(dome_high_dirs), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=dome_high_dirs[loop], data_avail=as.character(dome_high_modcombos[loop,"Data_avail"]), itervec=itervec, Fdynamics=as.character(strsplit(dome_high_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(dome_high_modcombos[loop,"Rdyn"],"_")[[1]][2]), write=TRUE, lh=lh_dome_high[[as.character(strsplit(dome_high_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=20, comp_sample=as.numeric(strsplit(dome_high_modcombos[loop,"SampleSize"],"_")[[1]][2]), rewrite=FALSE, init_depl=c(0.05,0.95))
	end_gen <- Sys.time() - start_gen

	## run LIME
	start_run <- Sys.time()
	foreach(loop=1:length(dome_base_dir_vec), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=dome_base_dir_vec[loop], lh=lh_list[[as.character(strsplit(dome_base_modcombos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(dome_base_modcombos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, REML=FALSE, f_true=FALSE, C_opt=dome_base_modcombos[loop,"C_opt"], param_adjust=c("SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(0.7,0.2,0.2,0.2), LFdist=1)
	end_run <- Sys.time() - start_run


