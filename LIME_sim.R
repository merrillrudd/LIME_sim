rm(list=ls())
library(TMB)
library(devtools)
library(TMBhelper)
library(foreach)
library(doParallel)
library(RColorBrewer)

devtools::install_github("merrillrudd/LIME", build.vignettes=TRUE, dependencies=TRUE)
library(LIME)
devtools::install_github("adrianhordyk/LBSPR", build.vignettes=TRUE, dependencies=TRUE)
library(LBSPR)



### ----- directories and functions -----------###
main_dir <- "C:\\Git_Projects\\LIME_sim"
# main_dir <- "F:\\Merrill\\Git_Projects\\LIME_sim"
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
ESS_vec <- c(1000, 100)
itervec <- 1:100

equil_modcombos <- expand.grid("Data_avail"=data_vec, "ESS"=paste0("ESS_", ESS_vec), "LH"=paste0("LH_",lh_vec), "Fdyn"=paste0("F_",Fdyn_vec), "Rdyn"=paste0("R_",Rdyn_vec), stringsAsFactors=FALSE)

### ----- equilibrium test ----------- ###

	equil_dir <- file.path(main_dir, "equil")
	# unlink(equil_dir, TRUE)
	dir.create(equil_dir, showWarnings=FALSE)	

	## setup equilibrium dirs
	equil_dir_vec <- model_paths(res_dir=equil_dir, modcombos=equil_modcombos)
	ignore <- sapply(1:length(equil_dir_vec), function(m) sapply(itervec, function(x) dir.create(file.path(equil_dir_vec[m], x), showWarnings=FALSE)))	

	## setup scenario life history list
	lh_list <- adj_variation(SigmaR=0.01, SigmaF=0.01, SigmaC=0.01, SigmaI=0.01, CVlen=0.01, rho=0)
	# lh_fig(lh=lh_list, save=FALSE)

	### ----- use parallel cores -----------###
	# registerDoParallel(cores=8)	
	registerDoParallel(cores=2)	

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
	foreach(loop=1:length(lime_dirs), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(equil_modcombos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, fix_f=0, simulation=TRUE, REML=FALSE, f_true=FALSE)
	end_run <- Sys.time() - start_run

	### ----- run LBSPR -----------###	

	start_run <- Sys.time()
	foreach(loop=1:length(lbspr_dirs), .packages=c('LBSPR','LIME')) %dopar% run_LBSPR(modpath=lbspr_dirs[loop], lh=lh_list[[as.character(strsplit(equil_modcombos[loop,"LH"],"_")[[1]][2])]], itervec=itervec, species=NULL, rewrite=TRUE, simulation=TRUE)
	end_altrun <- Sys.time() - start_run	

	### ------- check output ------###	

	equil_cover <- interval_coverage(dirs=lime_dirs, itervec=itervec)

	png(file.path(fig_dir, "Equilibrium_RE.png"), height=10, width=25, res=200, units="in")
	equil_res <- plot_re(dirs=equil_dir_vec, modcombos=equil_modcombos, itervec=itervec, compareToLH=paste0("F_", Fdyn_vec), lh_vec=lh_vec, cover=equil_cover$cover, ylim=c(-1,1))
	dev.off()

	re_med <- sapply(1:ncol(equil_res$bias), function(x) median(equil_res$bias[,x], na.rm=TRUE))
	re_med_LIME <- re_med[which(grepl("LBSPR", equil_dir_vec)==FALSE)]	



############################################################
#### base - with variation
############################################################
### ----- models to run ----------- ###
lh_vec <- c("Short", "Medium", "Long")
Fdyn_vec <- c("Endogenous")
Rdyn_vec <- c("AR")
data_vec <- c("Index_Catch_LC20", "Index_LC10", "Index_LC1", "Catch_LC10", "Catch_LC1", "LC10", "LC1")#, "LBSPR10", "LBSPR1")
ESS_vec <- c(1000, 100)
itervec <- 1:100

base_modcombos <- expand.grid("Data_avail"=data_vec, "ESS"=paste0("ESS_", ESS_vec), "LH"=paste0("LH_",lh_vec), "Fdyn"=paste0("F_",Fdyn_vec), "Rdyn"=paste0("R_",Rdyn_vec), stringsAsFactors=FALSE)

### ----- base runs ----------- ###
	base_dir <- file.path(main_dir, "base")
	# unlink(base_dir, TRUE)
	dir.create(base_dir, showWarnings=FALSE)	

	## setup baseibrium dirs
	base_dir_vec <- model_paths(res_dir=base_dir, modcombos=base_modcombos)
	ignore <- sapply(1:length(base_dir_vec), function(m) sapply(itervec, function(x) dir.create(file.path(base_dir_vec[m], x), showWarnings=FALSE)))	

	## setup scenario life history list
	lh_list <- adj_variation(SigmaR=0.737, SigmaF=0.3, SigmaC=0.2, SigmaI=0.2, CVlen=0.1, rho=0.426)
	# lh_fig(lh=lh_list, save=FALSE)

	### ----- use parallel cores -----------###
	# registerDoParallel(cores=8)	
	registerDoParallel(cores=2)	

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
	foreach(loop=1:length(lime_dirs), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(base_modcombos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, fix_f=0, simulation=TRUE, REML=FALSE, f_true=FALSE)
	end_run <- Sys.time() - start_run

	### ----- run LBSPR -----------###	

	start_run <- Sys.time()
	foreach(loop=1:length(lbspr_dirs), .packages=c('LBSPR','LIME')) %dopar% run_LBSPR(modpath=lbspr_dirs[loop], lh=lh_list[[as.character(strsplit(base_modcombos[loop,"LH"],"_")[[1]][2])]], itervec=itervec, species=NULL, rewrite=TRUE, simulation=TRUE)
	end_altrun <- Sys.time() - start_run	

	### ------- check output ------###	

	base_cover <- interval_coverage(dirs=lime_dirs, itervec=itervec)

	png(file.path(fig_dir, "Equilibrium_RE.png"), height=10, width=25, res=200, units="in")
	base_res <- plot_re(dirs=base_dir_vec, modcombos=base_modcombos, itervec=itervec, compareToLH=paste0("F_", Fdyn_vec), lh_vec=lh_vec, cover=base_cover$cover, ylim=c(-1,1))
	dev.off()

	re_med <- sapply(1:ncol(base_res$bias), function(x) median(base_res$bias[,x], na.rm=TRUE))
	re_med_LIME <- re_med[which(grepl("LBSPR", base_dir_vec)==FALSE)]	

### ----- base runs - AR high----------- ###
	base_highAR_dir <- file.path(main_dir, "base_highAR")
	# unlink(base_highAR_dir, TRUE)
	dir.create(base_highAR_dir, showWarnings=FALSE)	

	## setup base dirs
	base_highAR_dir_vec <- model_paths(res_dir=base_highAR_dir, modcombos=modcombos)
	ignore <- sapply(1:length(base_highAR_dir_vec), function(m) sapply(itervec, function(x) dir.create(file.path(base_highAR_dir_vec[m], x), showWarnings=FALSE)))	

	## setup scenario life history list
	lh_list <- adj_variation(SigmaR=0.737, SigmaF=0.3, SigmaC=0.2, SigmaI=0.2, CVlen=0.1, rho=0.976)
	# lh_fig(lh=lh_list, save=FALSE)

	### ----- model info -----------###	

	info <- list("Nyears"=20, "comp_sample"=1000)	

	### ----- use parallel cores -----------###
	# registerDoParallel(cores=8)	
	registerDoParallel(cores=2)	

	### ----- generate data -----------###
	start_gen <- Sys.time()
	foreach(loop=1:length(base_highAR_dir_vec), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=base_highAR_dir_vec[loop], data_avail=as.character(modcombos[loop,"Data_avail"]), itervec=itervec, spatial=FALSE, Fdynamics=as.character(strsplit(modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(modcombos[loop,"Rdyn"],"_")[[1]][2]), write=TRUE, lh=lh_list[[as.character(strsplit(modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=info$Nyears, comp_sample=info$comp_sample, rewrite=FALSE)
	end_gen <- Sys.time() - start_gen

	### ----- run LIME -----------###
	lime_dirs <- base_highAR_dir_vec[which(grepl("LBSPR",base_highAR_dir_vec)==FALSE)]
	lbspr_dirs <- base_highAR_dir_vec[which(grepl("LBSPR",base_highAR_dir_vec))]
	lime_combos <- modcombos[which(grepl("LBSPR", modcombos[,"Data_avail"])==FALSE),]
	lbspr_combos <- modcombos[which(grepl("LBSPR", modcombos[,"Data_avail"])),]	

	start_run <- Sys.time()
	foreach(loop=1:length(lime_dirs), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(modcombos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, fix_f=0, simulation=TRUE, REML=FALSE, f_true=FALSE)
	end_run <- Sys.time() - start_run

	### ----- run LBSPR -----------###	

	start_run <- Sys.time()
	foreach(loop=1:length(lbspr_dirs), .packages=c('LBSPR','LIME')) %dopar% run_LBSPR(modpath=lbspr_dirs[loop], lh=lh_list[[as.character(strsplit(modcombos[loop,"LH"],"_")[[1]][2])]], itervec=itervec, species=NULL, rewrite=TRUE, simulation=TRUE)
	end_altrun <- Sys.time() - start_run	

	### ------- check output ------###	

	base_highAR_cover <- interval_coverage(dirs=lime_dirs, itervec=itervec)

	png(file.path(fig_dir, "Equilibrium_RE.png"), height=10, width=25, res=200, units="in")
	base_highAR_res <- plot_re(dirs=base_highAR_dir_vec, modcombos=modcombos, itervec=itervec, compareToLH=paste0("F_", Fdyn_vec), lh_vec=lh_vec, cover=base_highAR_cover$cover, ylim=c(-1,1))
	dev.off()

	re_med <- sapply(1:ncol(base_highAR_res$bias), function(x) median(base_highAR_res$bias[,x], na.rm=TRUE))
	re_med_LIME <- re_med[which(grepl("LBSPR", base_highAR_dir_vec)==FALSE)]	





############################################################
#### sensitivities
############################################################
### ----- models to run ----------- ###
lh_vec <- "Medium"
Fdyn_vec <- "Constant" 
data_vec <- "Catch_LC1"
itervec <- 1

modcombos <- expand.grid("Data_avail"=data_vec, "LH"=lh_vec, "Fdyn"=Fdyn_vec, stringsAsFactors=FALSE)


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