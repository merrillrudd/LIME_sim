rm(list=ls())
library(TMB)
library(devtools)
library(TMBhelper)
library(foreach)
library(doParallel)
library(RColorBrewer)

devtools::install_github("merrillrudd/LIME", build.vignettes=TRUE, dependencies=TRUE)
devtools::install_github("adrianhordyk/LBSPR", build.vignettes=TRUE, dependencies=TRUE)

library(LIME)
library(LBSPR)



### ----- directories and functions -----------###
# main_dir <- "C:\\Git_Projects\\LIME_sim"
main_dir <- "F:\\Merrill\\Git_Projects\\LIME_sim"
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
data_vec <- c("Index_Catch_LC20", "Index_LC10", "Index_LC1", "Catch_LC10", "Catch_LC1", "LC10", "LC5", "LC2", "LC1", "LBSPR10", "LBSPR1")
ESS_vec <- c(1000, 100)
itervec <- 1:100

equil_modcombos <- expand.grid("Data_avail"=data_vec, "ESS"=paste0("ESS_", ESS_vec), "LH"=paste0("LH_",lh_vec), "Fdyn"=paste0("F_",Fdyn_vec), "Rdyn"=paste0("R_",Rdyn_vec), stringsAsFactors=FALSE)
equil_modcombos$C_opt <- rep(0, nrow(equil_modcombos))
	equil_modcombos$C_opt[which(grepl("Catch",equil_modcombos[,"Data_avail"]))] <- 1

### ----- equilibrium test ----------- ###

	equil_dir <- file.path(main_dir, "equil")
	# unlink(equil_dir, TRUE)
	dir.create(equil_dir, showWarnings=FALSE)	

	## setup equilibrium dirs
	equil_dir_vec <- model_paths(res_dir=equil_dir, modcombos=equil_modcombos[,-c(ncol(equil_modcombos))])
	ignore <- sapply(1:length(equil_dir_vec), function(m) sapply(itervec, function(x) dir.create(file.path(equil_dir_vec[m], x), showWarnings=FALSE)))	

	## setup scenario life history list
	lh_list <- adj_variation(SigmaR=0.01, SigmaF=0.01, SigmaC=0.01, SigmaI=0.01, CVlen=0.01, rho=0)
	# lh_fig(lh=lh_list, save=FALSE)

	### ----- use parallel cores -----------###
	registerDoParallel(cores=8)	

	### ----- generate data -----------###
	start_gen <- Sys.time()
	foreach(loop=1:length(equil_dir_vec), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=equil_dir_vec[loop], data_avail=as.character(equil_modcombos[loop,"Data_avail"]), itervec=itervec, spatial=FALSE, Fdynamics=as.character(strsplit(equil_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(equil_modcombos[loop,"Rdyn"],"_")[[1]][2]), write=TRUE, lh=lh_list[[as.character(strsplit(equil_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=20, comp_sample=as.numeric(strsplit(equil_modcombos[loop,"ESS"],"_")[[1]][2]), rewrite=FALSE)
	end_gen <- Sys.time() - start_gen

	### ----- run LIME -----------###
	lime_dirs <- equil_dir_vec[which(grepl("LBSPR",equil_dir_vec)==FALSE)]
	lbspr_dirs <- equil_dir_vec[which(grepl("LBSPR",equil_dir_vec))]
	lime_combos <- equil_modcombos[which(grepl("LBSPR", equil_modcombos[,"Data_avail"])==FALSE),]
	lbspr_combos <- equil_modcombos[which(grepl("LBSPR", equil_modcombos[,"Data_avail"])),]	

	start_run <- Sys.time()
	foreach(loop=1:length(lime_dirs), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=TRUE, fix_f=0, simulation=TRUE, REML=FALSE, f_true=FALSE, C_opt=equil_modcombos[loop,"C_opt"])
	end_run <- Sys.time() - start_run


	### ----- run LBSPR -----------###	

	start_run <- Sys.time()
	foreach(loop=1:length(lbspr_dirs), .packages=c('LBSPR','LIME')) %dopar% run_LBSPR(modpath=lbspr_dirs[loop], lh=lh_list[[as.character(strsplit(lbspr_combos[loop,"LH"],"_")[[1]][2])]], itervec=itervec, species=NULL, rewrite=TRUE, simulation=TRUE)
	end_altrun <- Sys.time() - start_run


	### ------- check output ------###	
	dir_esslow <- equil_dir_vec[grepl("ESS_100/", equil_dir_vec)]
	dir_esshigh <- equil_dir_vec[grepl("ESS_1000", equil_dir_vec)]
	
	equil_coverlow <- interval_coverage(dirs=dir_esslow[which(grepl("LBSPR",dir_esslow)==FALSE)], itervec=itervec)
	equil_coverhigh <- interval_coverage(dirs=dir_esshigh[which(grepl("LBSPR",dir_esshigh)==FALSE)], itervec=itervec)

	png(file.path(fig_dir, "Equilibrium_ESSlow_RE.png"), height=10, width=25, res=200, units="in")
	equil_res <- plot_re(dirs=dir_esslow, modcombos=equil_modcombos[which(equil_modcombos$ESS=="ESS_100"),], itervec=itervec, compareToLH=paste0("F_", Fdyn_vec), lh_vec=lh_vec, cover=equil_coverlow$cover, ylim=c(-1,1))
	dev.off()

	png(file.path(fig_dir, "Equilibrium_ESShigh_RE.png"), height=10, width=25, res=200, units="in")
	equil_res <- plot_re(dirs=dir_esshigh, modcombos=equil_modcombos[which(equil_modcombos$ESS=="ESS_1000"),], itervec=itervec, compareToLH=paste0("F_", Fdyn_vec), lh_vec=lh_vec, cover=equil_coverhigh$cover, ylim=c(-1,1))
	dev.off()



############################################################
#### base - with variation
############################################################
### ----- models to run ----------- ###
lh_vec <- c("Short", "Medium", "Long")
Fdyn_vec <- c("Endogenous")
Rdyn_vec <- c("AR")
data_vec <- c("Index_Catch_LC20", "Index_LC10", "Index_LC1", "Catch_LC10", "Catch_LC1", "LC10", "LC5", "LC2", "LC1", "LBSPR10", "LBSPR1")
ESS_vec <- c(1000, 100)
itervec <- 1:100

base_modcombos <- expand.grid("Data_avail"=data_vec, "ESS"=paste0("ESS_", ESS_vec), "LH"=paste0("LH_",lh_vec), "Fdyn"=paste0("F_",Fdyn_vec), "Rdyn"=paste0("R_",Rdyn_vec), stringsAsFactors=FALSE)
base_modcombos$C_opt <- rep(0, nrow(base_modcombos))
	base_modcombos$C_opt[which(grepl("Catch",base_modcombos[,"Data_avail"]))] <- 1

### ----- base runs ----------- ###
	base_dir <- file.path(main_dir, "base")
	# unlink(base_dir, TRUE)
	dir.create(base_dir, showWarnings=FALSE)	

	## setup baseibrium dirs
	base_dir_vec <- model_paths(res_dir=base_dir, modcombos=base_modcombos[,-c(ncol(base_modcombos))])
	ignore <- sapply(1:length(base_dir_vec), function(m) sapply(itervec, function(x) dir.create(file.path(base_dir_vec[m], x), showWarnings=FALSE)))	

	## setup scenario life history list
	lh_list <- adj_variation(SigmaR=0.737, SigmaF=0.2, SigmaC=0.2, SigmaI=0.2, CVlen=0.1, rho=0.426)
	# lh_fig(lh=lh_list, save=FALSE)

	### ----- use parallel cores -----------###
	registerDoParallel(cores=8)	
	# registerDoParallel(cores=2)	

	### ----- generate data -----------###
	start_gen <- Sys.time()
	foreach(loop=1:length(base_dir_vec), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=base_dir_vec[loop], data_avail=as.character(base_modcombos[loop,"Data_avail"]), itervec=itervec, spatial=FALSE, Fdynamics=as.character(strsplit(base_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(base_modcombos[loop,"Rdyn"],"_")[[1]][2]), write=TRUE, lh=lh_list[[as.character(strsplit(base_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=20, comp_sample=as.numeric(strsplit(base_modcombos[loop,"ESS"],"_")[[1]][2]), rewrite=FALSE)
	end_gen <- Sys.time() - start_gen

	### ----- run LIME -----------###
	lime_dirs <- base_dir_vec[which(grepl("LBSPR",base_dir_vec)==FALSE)]
	lbspr_dirs <- base_dir_vec[which(grepl("LBSPR",base_dir_vec))]
	lime_combos <- base_modcombos[which(grepl("LBSPR", base_modcombos[,"Data_avail"])==FALSE),]
	lbspr_combos <- base_modcombos[which(grepl("LBSPR", base_modcombos[,"Data_avail"])),]	

	start_run <- Sys.time()
	foreach(loop=1:length(lime_dirs), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, fix_f=0, simulation=TRUE, REML=FALSE, f_true=FALSE, C_opt=base_modcombos[loop,"C_opt"])
	end_run <- Sys.time() - start_run

	### ----- run LBSPR -----------###	

	start_run <- Sys.time()
	foreach(loop=1:length(lbspr_dirs), .packages=c('LBSPR','LIME')) %dopar% run_LBSPR(modpath=lbspr_dirs[loop], lh=lh_list[[as.character(strsplit(lbspr_combos[loop,"LH"],"_")[[1]][2])]], itervec=itervec, species=NULL, rewrite=TRUE, simulation=TRUE)
	end_altrun <- Sys.time() - start_run	

	### ------- check output ------###	
	dir_esslow <- base_dir_vec[grepl("ESS_100/", base_dir_vec)]
	dir_esshigh <- base_dir_vec[grepl("ESS_1000", base_dir_vec)]
	
	base_coverlow <- interval_coverage(dirs=dir_esslow[which(grepl("LBSPR",dir_esslow)==FALSE)], itervec=itervec)
	base_coverhigh <- interval_coverage(dirs=dir_esshigh[which(grepl("LBSPR",dir_esshigh)==FALSE)], itervec=itervec)

	png(file.path(fig_dir, "Base_ESSlow_RE.png"), height=10, width=25, res=200, units="in")
	base_res <- plot_re(dirs=dir_esslow, modcombos=base_modcombos[which(base_modcombos$ESS=="ESS_100"),], itervec=itervec, compareToLH=paste0("F_", Fdyn_vec), lh_vec=lh_vec, cover=base_coverlow$cover, ylim=c(-1,1))
	dev.off()

	png(file.path(fig_dir, "Base_ESShigh_RE.png"), height=10, width=25, res=200, units="in")
	base_res <- plot_re(dirs=dir_esshigh, modcombos=base_modcombos[which(base_modcombos$ESS=="ESS_1000"),], itervec=itervec, compareToLH=paste0("F_", Fdyn_vec), lh_vec=lh_vec, cover=base_coverhigh$cover, ylim=c(-1,1))
	dev.off()


### ----- base runs - sigR lower----------- ###
	base_lowsigR_dir <- file.path(main_dir, "base_lowsigR")
	# unlink(base_lowsigR_dir, TRUE)
	dir.create(base_lowsigR_dir, showWarnings=FALSE)	

	## setup base dirs
	base_lowsigR_dir_vec <- model_paths(res_dir=base_lowsigR_dir, modcombos=base_modcombos[,-c(ncol(base_modcombos))])
	ignore <- sapply(1:length(base_lowsigR_dir_vec), function(m) sapply(itervec, function(x) dir.create(file.path(base_lowsigR_dir_vec[m], x), showWarnings=FALSE)))	

	## setup scenario life history list
	lh_list <- adj_variation(SigmaR=0.384, SigmaF=0.2, SigmaC=0.2, SigmaI=0.2, CVlen=0.1, rho=0.426)
	# lh_fig(lh=lh_list, save=FALSE)

	### ----- use parallel cores -----------###
	registerDoParallel(cores=8)	

	### ----- generate data -----------###
	start_gen <- Sys.time()
	foreach(loop=1:length(base_lowsigR_dir_vec), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=base_lowsigR_dir_vec[loop], data_avail=as.character(base_modcombos[loop,"Data_avail"]), itervec=itervec, spatial=FALSE, Fdynamics=as.character(strsplit(base_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(base_modcombos[loop,"Rdyn"],"_")[[1]][2]), write=TRUE, lh=lh_list[[as.character(strsplit(base_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=20, comp_sample=as.numeric(strsplit(base_modcombos[loop,"ESS"],"_")[[1]][2]), rewrite=FALSE)
	end_gen <- Sys.time() - start_gen

	### ----- run LIME -----------###
	lime_dirs <- base_lowsigR_dir_vec[which(grepl("LBSPR",base_lowsigR_dir_vec)==FALSE)]
	lbspr_dirs <- base_lowsigR_dir_vec[which(grepl("LBSPR",base_lowsigR_dir_vec))]
	lime_combos <- base_modcombos[which(grepl("LBSPR", base_modcombos[,"Data_avail"])==FALSE),]
	lbspr_combos <- base_modcombos[which(grepl("LBSPR", base_modcombos[,"Data_avail"])),]	

	start_run <- Sys.time()
	foreach(loop=1:length(lime_dirs), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, fix_f=0, simulation=TRUE, REML=FALSE, f_true=FALSE, C_opt=base_modcombos[loop,"C_opt"])
	end_run <- Sys.time() - start_run

	### ----- run LBSPR -----------###	

	start_run <- Sys.time()
	foreach(loop=1:length(lbspr_dirs), .packages=c('LBSPR','LIME')) %dopar% run_LBSPR(modpath=lbspr_dirs[loop], lh=lh_list[[as.character(strsplit(lbspr_combos[loop,"LH"],"_")[[1]][2])]], itervec=itervec, species=NULL, rewrite=TRUE, simulation=TRUE)
	end_altrun <- Sys.time() - start_run	

	### ------- check output ------###	

	dir_esslow <- base_lowsigR_dir_vec[grepl("ESS_100/", base_lowsigR_dir_vec)]
	dir_esshigh <- base_lowsigR_dir_vec[grepl("ESS_1000", base_lowsigR_dir_vec)]
	
	base_lowsigR_coverlow <- interval_coverage(dirs=dir_esslow[which(grepl("LBSPR",dir_esslow)==FALSE)], itervec=itervec)
	base_lowsigR_coverhigh <- interval_coverage(dirs=dir_esshigh[which(grepl("LBSPR",dir_esshigh)==FALSE)], itervec=itervec)

	png(file.path(fig_dir, "Base_LowSigR_ESSlow_RE.png"), height=10, width=25, res=200, units="in")
	base_lowsigR_res <- plot_re(dirs=dir_esslow, modcombos=base_modcombos[which(base_modcombos$ESS=="ESS_100"),], itervec=itervec, compareToLH=paste0("F_", Fdyn_vec), lh_vec=lh_vec, cover=base_lowsigR_coverlow$cover, ylim=c(-1,1))
	dev.off()

	png(file.path(fig_dir, "Base_LowSigR_ESShigh_RE.png"), height=10, width=25, res=200, units="in")
	base_lowsigR_res <- plot_re(dirs=dir_esshigh, modcombos=base_modcombos[which(base_modcombos$ESS=="ESS_1000"),], itervec=itervec, compareToLH=paste0("F_", Fdyn_vec), lh_vec=lh_vec, cover=base_lowsigR_coverhigh$cover, ylim=c(-1,1))
	dev.off()




############################################################
#### sensitivities - equilibrium
############################################################
### ----- models to run ----------- ###
lh_vec <- c("Short", "Medium", "Long")
Fdyn_vec <- "Constant"
Rdyn_vec <- "Constant" 
data_vec <- c("LC10")
ESS_vec <- c(1000, 100)
itervec <- 1:100
input_vec <- c("M", "linf", "vbk", "ML50", "CVlen")
val_vec <- c("high", "low")

sens_equil_modcombos <- expand.grid("Param"=input_vec, "Adjust"=val_vec, "Data_avail"=data_vec, "ESS"=paste0("ESS_", ESS_vec), "LH"=paste0("LH_",lh_vec), "Fdyn"=paste0("F_",Fdyn_vec), "Rdyn"=paste0("R_",Rdyn_vec), stringsAsFactors=FALSE)
sens_equil_modcombos$C_opt <- rep(0, nrow(sens_equil_modcombos))
	sens_equil_modcombos$C_opt[which(grepl("Catch",sens_equil_modcombos[,"Data_avail"]))] <- 1

### ----- sens runs ----------- ###
	sens_equil_dir <- file.path(main_dir, "sens_equil")
	# unlink(sens_equil_dir, TRUE)
	dir.create(sens_equil_dir, showWarnings=FALSE)	

	## setup dirs
	sens_equil_dir_vec <- model_paths(res_dir=sens_equil_dir, modcombos=sens_equil_modcombos[,-ncol(sens_equil_modcombos)])
	ignore <- sapply(1:length(sens_equil_dir_vec), function(m) sapply(itervec, function(x) dir.create(file.path(sens_equil_dir_vec[m], x), showWarnings=FALSE)))	

	## equilibrium
	lh_list <- adj_variation(SigmaR=0.01, SigmaF=0.01, SigmaC=0.01, SigmaI=0.01, CVlen=0.01, rho=0)

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
	registerDoParallel(cores=8)	

	### ----- generate data -----------###
	## copy data files from similar previous directory
	alt_dirs <- rep(NA, length(sens_equil_dir_vec))
	## find directory
	for(i in 1:length(alt_dirs)){
		xdir <- equil_dir_vec[grepl(paste0("/",sens_equil_modcombos[i,"Data_avail"]),equil_dir_vec)&grepl(paste0(sens_equil_modcombos[i,"ESS"],"/"),equil_dir_vec)&grepl(sens_equil_modcombos[i,"LH"],equil_dir_vec)&grepl(sens_equil_modcombos[i,"Fdyn"],equil_dir_vec)&grepl(sens_equil_modcombos[i,"Rdyn"],equil_dir_vec)]
		if(length(xdir)>1) break
		alt_dirs[i] <- xdir
	}
	## copy to directory
	for(i in 1:length(alt_dirs)){
		ignore <- sapply(1:length(itervec), function(x) file.copy(from=file.path(alt_dirs[i],x,"True.rds"), to=file.path(sens_equil_dir_vec[i],x)))
	}

	## run LIME
	start_run <- Sys.time()
	foreach(loop=1:length(sens_equil_dir_vec), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=sens_equil_dir_vec[loop], lh=lh_list[[as.character(strsplit(sens_equil_modcombos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(sens_equil_modcombos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, fix_f=0, simulation=TRUE, REML=FALSE, f_true=FALSE, C_opt=sens_equil_modcombos[loop,"C_opt"], param_adjust=sens_equil_modcombos[loop,"Param"], val_adjust=sens_vals[[as.character(strsplit(sens_equil_modcombos[loop,"LH"],"_")[[1]][2])]][sens_equil_modcombos[loop,"Adjust"],sens_equil_modcombos[loop,"Param"]])
	end_run <- Sys.time() - start_run

############################################################
#### sensitivities - base
############################################################
### ----- models to run ----------- ###
lh_vec <- c("Short", "Medium", "Long")
Fdyn_vec <- "Endogenous"
Rdyn_vec <- "AR" 
data_vec <- c("LC10")
ESS_vec <- c(1000, 100)
itervec <- 1:100
input_vec <- c("M", "linf", "vbk", "ML50", "CVlen")
val_vec <- c("high", "low")

sens_base_modcombos <- expand.grid("Param"=input_vec, "Adjust"=val_vec, "Data_avail"=data_vec, "ESS"=paste0("ESS_", ESS_vec), "LH"=paste0("LH_",lh_vec), "Fdyn"=paste0("F_",Fdyn_vec), "Rdyn"=paste0("R_",Rdyn_vec), stringsAsFactors=FALSE)
sens_base_modcombos$C_opt <- rep(0, nrow(sens_base_modcombos))
	sens_base_modcombos$C_opt[which(grepl("Catch",sens_base_modcombos[,"Data_avail"]))] <- 1

### ----- sens runs ----------- ###
	sens_base_dir <- file.path(main_dir, "sens_base")
	unlink(sens_base_dir, TRUE)
	dir.create(sens_base_dir, showWarnings=FALSE)	

	## setup dirs
	sens_base_dir_vec <- model_paths(res_dir=sens_base_dir, modcombos=sens_base_modcombos[,-ncol(sens_base_modcombos)])
	ignore <- sapply(1:length(sens_base_dir_vec), function(m) sapply(itervec, function(x) dir.create(file.path(sens_base_dir_vec[m], x), showWarnings=FALSE)))	

	## base
	lh_list <- adj_variation(SigmaR=0.737, SigmaF=0.2, SigmaC=0.2, SigmaI=0.2, CVlen=0.1, rho=0.426)

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
	registerDoParallel(cores=8)	

	### ----- generate data -----------###
	## copy data files from similar previous directory
	alt_dirs <- rep(NA, length(sens_base_dir_vec))
	## find directory
	for(i in 1:length(alt_dirs)){
		xdir <- base_dir_vec[grepl(paste0("/",sens_base_modcombos[i,"Data_avail"]),base_dir_vec)&grepl(paste0(sens_base_modcombos[i,"ESS"],"/"),base_dir_vec)&grepl(sens_base_modcombos[i,"LH"],base_dir_vec)&grepl(sens_base_modcombos[i,"Fdyn"],base_dir_vec)&grepl(sens_base_modcombos[i,"Rdyn"],base_dir_vec)]
		if(length(xdir)>1) break
		alt_dirs[i] <- xdir
	}
	## copy to directory
	for(i in 1:length(alt_dirs)){
		ignore <- sapply(1:length(itervec), function(x) file.copy(from=file.path(alt_dirs[i],x,"True.rds"), to=file.path(sens_base_dir_vec[i],x)))
	}

	## run LIME
	start_run <- Sys.time()
	foreach(loop=1:length(sens_base_dir_vec), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=sens_base_dir_vec[loop], lh=lh_list[[as.character(strsplit(sens_base_modcombos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(sens_base_modcombos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, fix_f=0, simulation=TRUE, REML=FALSE, f_true=FALSE, C_opt=sens_base_modcombos[loop,"C_opt"], param_adjust=sens_base_modcombos[loop,"Param"], val_adjust=sens_vals[[as.character(strsplit(sens_base_modcombos[loop,"LH"],"_")[[1]][2])]][sens_base_modcombos[loop,"Adjust"],sens_base_modcombos[loop,"Param"]])
	end_run <- Sys.time() - start_run







### ----- Piner plot - Catch_LC1 ----------- ###
	sens_dir1 <- file.path(main_dir, "sens1")
	# unlink(sens_dir1, TRUE)
	dir.create(sens_dir1, showWarnings=FALSE)	

	## setup base dirs
	sens_dir1_vec <- model_paths(res_dir=sens_dir1, modcombos=modcombos)
	dir.create(sens_dir1_vec, showWarnings=FALSE)
	# ignore <- sapply(1:length(sens_dir1_vec), function(m) sapply(itervec, function(x) dir.create(file.path(sens_dir1_vec[m], x), showWarnings=FALSE)))	

	## setup scenario life history list
	lh_list <- adj_variation(SigmaR=0.6, SigmaF=0.3, SigmaC=0.1, SigmaI=0.1)
	# lh_fig(lh=lh_list, save=FALSE)

	### ----- model info -----------###	

	info <- list("Nyears"=20, "comp_sample"=1000)	

	loop <- 1
	### ----- generate data -----------###

	R0_vec <- seq(0.5,2,by=0.1)
	for(i in 1:length(R0_vec)){
		dir <- file.path(sens_dir1_vec[loop], paste0("R0_", R0_vec[i]))
		dir.create(dir, showWarnings=FALSE)

		ignore <- generate_data(modpath=dir, data_avail=as.character(modcombos[loop,"Data_avail"]), itervec=itervec, spatial=FALSE, Fdynamics=as.character(modcombos[loop,"Fdyn"]), Rdynamics="Constant", write=TRUE, lh=lh_list[[as.character(modcombos[loop,"LH"])]], Nyears=info$Nyears, comp_sample=info$comp_sample, rewrite=FALSE)

	### ----- run LIME -----------###
		ignore <- run_LIME(modpath=dir, lh=lh_list[[as.character(modcombos[loop,"LH"])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(modcombos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, fix_f=0, simulation=TRUE, REML=FALSE, f_true=FALSE, write=TRUE, fix_param="beta", param_adjust="R0", val_adjust=R0_vec[i])
	}

	jnll_mat <- matrix(NA, nrow=length(R0_vec), ncol=3)
	for(i in 1:length(R0_vec)){
		rep <- readRDS(file.path(sens_dir1_vec, paste0("R0_", R0_vec[i]), 1, "Report.rds"))
		jnll_mat[i,1] <- rep$jnll
		jnll_mat[i,2] <- rep$jnll_comp[2]
		jnll_mat[i,3] <- rep$jnll_comp[5]
	}
	par(mfrow=c(1,3))
	plot(x=R0_vec, y=jnll_mat[,1], pch=19, type="o")
		abline(v=1, col="gray")
	plot(x=R0_vec, y=jnll_mat[,2], pch=19, col="steelblue", type="o")
		abline(v=1, col="gray")
	plot(x=R0_vec, y=jnll_mat[,3], pch=19, col="tomato", type="o")
		abline(v=1, col="gray")



	true <- readRDS(file.path(sens_dir1_vec[loop], 1, "True.rds"))
	inp <- readRDS(file.path(sens_dir1_vec[loop], 1, "Inputs.rds"))
	rep <- readRDS(file.path(sens_dir1_vec[loop], 1, "Report.rds"))
	sdrep <- readRDS(file.path(sens_dir1_vec[loop], 1, "Sdreport.rds"))

	FUN <- function(InputMat, log=TRUE, rel=FALSE){
          index <- which(is.na(InputMat[,2])==FALSE)
          if(log==TRUE) return(c( exp(InputMat[index,1]-1.96*InputMat[index,2]), rev(exp(InputMat[index,1]+1.96*InputMat[index,2]))))
          if(log==FALSE) return(c( InputMat[index,1]-1.96*InputMat[index,2], rev(InputMat[index,1]+1.96*InputMat[index,2])))
	} 
	
  plot(true$F_t, lwd=2, type="l", ylim=c(0, 1))
  lines(rep$F_t, lwd=2, col="blue", ylim=c(0, max(rep$F_t)*1.5), type="l")
  polygon( y=FUN(summary(sdrep)[which(rownames(summary(sdrep))=="lF_t"),], log=TRUE), x=c(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="lF_t"),2])==FALSE), rev(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="lF_t"),2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)

  exp(rep$beta)
  true$R0