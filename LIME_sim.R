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
data_vec <- c("Index_Catch_LC20", "Index_LC10", "Index_LC1", "Catch_LC10", "Catch_LC1", "LC10", "LC5", "LC2", "LC1", "LBSPR10", "LBSPR1")
ESS_vec <- c(1000, 100, 20)
itervec <- 1:100

equil_modcombos <- expand.grid("Data_avail"=data_vec, "ESS"=paste0("ESS_", ESS_vec), "LH"=paste0("LH_",lh_vec), "Fdyn"=paste0("F_",Fdyn_vec), "Rdyn"=paste0("R_",Rdyn_vec), stringsAsFactors=FALSE)
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
	lh_list <- adj_variation(SigmaR=0.01, SigmaF=0.01, SigmaC=0.01, SigmaI=0.01, CVlen=0.1, rho=0)
	lh_fig(lh=lh_list, save=TRUE)


	### ----- generate data -----------###
	### data rich cases only
	rich_dir <- equil_dir_vec[grepl("Index_Catch_LC20",equil_dir_vec)]
	rich_modcombos <- equil_modcombos[which(equil_modcombos$Data_avail=="Index_Catch_LC20"),]
	alt_dir <- equil_dir_vec[grepl("Index_Catch_LC20",equil_dir_vec)==FALSE]
	alt_modcombos <- equil_modcombos[which(equil_modcombos$Data_avail!="Index_Catch_LC20"),]

	### ----- use parallel cores -----------###
	registerDoParallel(cores=ncores)	

	start_gen <- Sys.time()
	foreach(loop=1:length(rich_dir), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=rich_dir[loop], data_avail=as.character(rich_modcombos[loop,"Data_avail"]), itervec=itervec, Fdynamics=as.character(strsplit(rich_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(rich_modcombos[loop,"Rdyn"],"_")[[1]][2]), write=TRUE, lh=lh_list[[as.character(strsplit(rich_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=20, comp_sample=as.numeric(strsplit(rich_modcombos[loop,"ESS"],"_")[[1]][2]), rewrite=FALSE, init_depl=c(0.05,0.95))
	end_gen <- Sys.time() - start_gen

	# copy data rich cases to data poor
	# when written using parallel cores, sometimes files not writing properly
	start_regen <- Sys.time()
	for(loop in 1:length(rich_dir)){
		copy_sim(fromdir=rich_dir[loop], fromcombos=rich_modcombos[loop,], todir=alt_dir, itervec=itervec, rewrite=FALSE, res_dir="equil")
	}
	end_regen <- Sys.time() - start_regen


	# plot equilibrium scenarios
	png(file.path(fig_dir, "Equil_scenarios.png"), height=10, width=12, res=200, units="in")
	plot_scenarios(dirs=equil_dir_vec[grepl("ESS_1000/",equil_dir_vec) & grepl("Index_Catch_LC20",equil_dir_vec)], itervec=itervec)
	dev.off()

	### ----- run LIME -----------###
	lime_dirs <- equil_dir_vec[which(grepl("LBSPR",equil_dir_vec)==FALSE)]
	lbspr_dirs <- equil_dir_vec[which(grepl("LBSPR",equil_dir_vec))]
	lime_combos <- equil_modcombos[which(grepl("LBSPR", equil_modcombos[,"Data_avail"])==FALSE),]
	lbspr_combos <- equil_modcombos[which(grepl("LBSPR", equil_modcombos[,"Data_avail"])),]	

	start_run <- Sys.time()
	foreach(loop=1:length(lime_dirs), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, fix_f=0, simulation=TRUE, REML=FALSE, f_true=FALSE, C_opt=lime_combos[loop,"C_opt"], Sel0=1, LFdist=1, param_adjust=c("SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(0.7,0.2,0.2,0.2), write=TRUE, fix_param=FALSE)
	end_run <- Sys.time() - start_run


	### ----- run LBSPR -----------###	

	start_run <- Sys.time()
	foreach(loop=1:length(lbspr_dirs), .packages=c('LBSPR','LIME')) %dopar% run_LBSPR(modpath=lbspr_dirs[loop], lh=lh_list[[as.character(strsplit(lbspr_combos[loop,"LH"],"_")[[1]][2])]], itervec=itervec, species=NULL, rewrite=FALSE, simulation=TRUE)
	end_altrun <- Sys.time() - start_run


	### ------- check output ------###	
	dir_esslowest <- equil_dir_vec[grepl("ESS_20/", equil_dir_vec)]
	dir_esslow <- equil_dir_vec[grepl("ESS_100/", equil_dir_vec)]
	dir_esshigh <- equil_dir_vec[grepl("ESS_1000", equil_dir_vec)]
	
	equil_coverlowest <- interval_coverage(dirs=dir_esslowest[which(grepl("LBSPR",dir_esslowest)==FALSE)], itervec=itervec) ## 19
	equil_coverlow <- interval_coverage(dirs=dir_esslow[which(grepl("LBSPR",dir_esslow)==FALSE)], itervec=itervec) ## 13
	equil_coverhigh <- interval_coverage(dirs=dir_esshigh[which(grepl("LBSPR",dir_esshigh)==FALSE)], itervec=itervec) ## 22

	png(file.path(fig_dir, "Equilibrium_ESSlowest_RE.png"), height=10, width=25, res=200, units="in")
	equil_res <- plot_re(dirs=dir_esslowest, modcombos=equil_modcombos[which(equil_modcombos$ESS=="ESS_20"),], itervec=itervec, compareToLH=paste0("F_", Fdyn_vec), lh_vec=lh_vec, cover=equil_coverlowest$cover, ylim=c(-1,1))
	dev.off()

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
Fdyn_vec <- c("Ramp")#, "Increasing", "Decreasing")
Rdyn_vec <- c("AR")
data_vec <- c("Index_Catch_LC20", "Index_LC10", "Index_LC1", "Catch_LC10", "Catch_LC1", "LC10", "LC5", "LC2", "LC1", "LBSPR10", "LBSPR1")
ESS_vec <- c(1000, 100)
itervec <- 1:100

base_modcombos <- expand.grid("Data_avail"=data_vec, "ESS"=paste0("ESS_", ESS_vec), "LH"=paste0("LH_",lh_vec), "Fdyn"=paste0("F_",Fdyn_vec), "Rdyn"=paste0("R_",Rdyn_vec), stringsAsFactors=FALSE)
base_modcombos$C_opt <- rep(0, nrow(base_modcombos))
	base_modcombos$C_opt[which(grepl("Catch",base_modcombos[,"Data_avail"]))] <- 2

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

	### ----- generate data -----------###
	### data rich cases only
	rich_dir <- base_dir_vec[grepl("Index_Catch_LC20",base_dir_vec)]
	rich_modcombos <- base_modcombos[which(base_modcombos$Data_avail=="Index_Catch_LC20"),]
	alt_dir <- base_dir_vec[grepl("Index_Catch_LC20",base_dir_vec)==FALSE]
	alt_modcombos <- base_modcombos[which(base_modcombos$Data_avail!="Index_Catch_LC20"),]

	### ----- use parallel cores -----------###
	registerDoParallel(cores=ncores)	

	start_gen <- Sys.time()
	foreach(loop=1:length(rich_dir), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=rich_dir[loop], data_avail=as.character(rich_modcombos[loop,"Data_avail"]), itervec=itervec, Fdynamics=as.character(strsplit(rich_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(rich_modcombos[loop,"Rdyn"],"_")[[1]][2]), write=TRUE, lh=lh_list[[as.character(strsplit(rich_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=20, comp_sample=as.numeric(strsplit(rich_modcombos[loop,"ESS"],"_")[[1]][2]), rewrite=FALSE, init_depl=c(0.05,0.95))
	end_gen <- Sys.time() - start_gen

	# copy data rich cases to data poor
	# when written using parallel cores, sometimes files not writing properly
	start_regen <- Sys.time()
	for(loop in 1:length(rich_dir)){
		copy_sim(fromdir=rich_dir[loop], fromcombos=rich_modcombos[loop,], todir=alt_dir, itervec=itervec, rewrite=FALSE, res_dir="base")
	}
	end_regen <- Sys.time() - start_regen

	##plot base scenarios
	png(file.path(fig_dir, "Base_Ramp_scenarios.png"), height=10, width=12, res=200, units="in")
	plot_scenarios(dirs=base_dir_vec[grepl("ESS_1000/",base_dir_vec) & grepl("Index_Catch_LC20",base_dir_vec) & grepl("Ramp", base_dir_vec)], itervec=itervec)
	dev.off()

	# 	##plot base scenarios
	# png(file.path(fig_dir, "Base_Increasing_scenarios.png"), height=10, width=12, res=200, units="in")
	# plot_scenarios(dirs=base_dir_vec[grepl("ESS_1000/",base_dir_vec) & grepl("Index_Catch_LC20",base_dir_vec) & grepl("Increasing", base_dir_vec)], itervec=itervec)
	# dev.off()

	# 		##plot base scenarios
	# png(file.path(fig_dir, "Base_Decreasing_scenarios.png"), height=10, width=12, res=200, units="in")
	# plot_scenarios(dirs=base_dir_vec[grepl("ESS_1000/",base_dir_vec) & grepl("Index_Catch_LC20",base_dir_vec) & grepl("Decreasing", base_dir_vec)], itervec=itervec)
	# dev.off()


	### ----- run LIME -----------###
	lime_dirs <- base_dir_vec[which(grepl("LBSPR",base_dir_vec)==FALSE)]
	lbspr_dirs <- base_dir_vec[which(grepl("LBSPR",base_dir_vec))]
	lime_combos <- base_modcombos[which(grepl("LBSPR", base_modcombos[,"Data_avail"])==FALSE),]
	lbspr_combos <- base_modcombos[which(grepl("LBSPR", base_modcombos[,"Data_avail"])),]	

	start_run <- Sys.time()
	foreach(loop=1:length(lime_dirs), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, fix_f=0, simulation=TRUE, REML=FALSE, f_true=FALSE, C_opt=lime_combos[loop,"C_opt"], Sel0=1, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.7,0.2,0.2), write=TRUE, fix_param=FALSE)
	end_run <- Sys.time() - start_run


	### ----- run LBSPR -----------###	

	start_run <- Sys.time()
	foreach(loop=1:length(lbspr_dirs), .packages=c('LBSPR','LIME')) %dopar% run_LBSPR(modpath=lbspr_dirs[loop], lh=lh_list[[as.character(strsplit(lbspr_combos[loop,"LH"],"_")[[1]][2])]], itervec=itervec, species=NULL, rewrite=FALSE, simulation=TRUE)
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

	### plot base truths
	plot_vec <- base_dir_vec[grepl("ESS_1000/",base_dir_vec) & grepl("Index_Catch_LC20",base_dir_vec)]
	par(mfrow=c(3,3), mar=c(0,4,0,0), omi=c(1,1,0.5,0.5))
	##Fishing
	plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,20), xaxs="i", yaxs="i", ylab="", xaxt="n", yaxt="n")
	axis(2,las=2,cex.axis=1.5)
	mtext("Fishing\nmortality", side=2, line=3, cex=2)
	for(iter in itervec){
		true <- readRDS(file.path(plot_vec[1], iter, "True.rds"))
		lines(true$F_t, col="#AA000050")
	}
	plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,10), xaxs="i", yaxs="i", ylab="", xaxt="n", yaxt="n")
		axis(2,las=2,cex.axis=1.5)
	for(iter in itervec){
		true <- readRDS(file.path(plot_vec[2], iter, "True.rds"))
		lines(true$F_t, col="#AA000050")
	}
	plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,2), xaxs="i", yaxs="i", ylab="", xaxt="n", yaxt="n")
		axis(2,las=2,cex.axis=1.5)
	for(iter in itervec){
		true <- readRDS(file.path(plot_vec[3], iter, "True.rds"))
		lines(true$F_t, col="#AA000050")
	}
	## Recruitment
	plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,10), xaxs="i", yaxs="i", ylab="", xaxt="n", yaxt="n")
		axis(2,las=2,cex.axis=1.5)
		mtext("Recruitment", side=2, line=3, cex=2)
	for(iter in itervec){
		true <- readRDS(file.path(plot_vec[1], iter, "True.rds"))
		lines(true$R_t, col="#0000AA50")
	}
	plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,10), xaxs="i", yaxs="i", ylab="", xaxt="n", yaxt="n")
		axis(2,las=2,cex.axis=1.5)
	for(iter in itervec){
		true <- readRDS(file.path(plot_vec[2], iter, "True.rds"))
		lines(true$R_t, col="#0000AA50")
	}
	plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,10), xaxs="i", yaxs="i", ylab="", xaxt="n", yaxt="n")
		axis(2,las=2,cex.axis=1.5)
	for(iter in itervec){
		true <- readRDS(file.path(plot_vec[3], iter, "True.rds"))
		lines(true$R_t, col="#0000AA50")
	}
	## SPR
	plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,1), xaxs="i", yaxs="i", ylab="",  yaxt="n", cex.axis=1.5)
		axis(2,las=2,cex.axis=1.5)
		mtext("SPR", side=2, line=3, cex=2)
	for(iter in itervec){
		true <- readRDS(file.path(plot_vec[1], iter, "True.rds"))
		lines(true$SPR_t, col="#00AA0050")
	}
	plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,1), xaxs="i", yaxs="i", ylab="",  yaxt="n", cex.axis=1.5)
		axis(2,las=2,cex.axis=1.5)
	for(iter in itervec){
		true <- readRDS(file.path(plot_vec[2], iter, "True.rds"))
		lines(true$SPR_t, col="#00AA0050")
	}
	plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,1), xaxs="i", yaxs="i", ylab="",  yaxt="n", cex.axis=1.5)
		axis(2,las=2,cex.axis=1.5)
	for(iter in itervec){
		true <- readRDS(file.path(plot_vec[3], iter, "True.rds"))
		lines(true$SPR_t, col="#00AA0050")
	}
	mtext(side=1, "Year", outer=TRUE, line=3, cex=2)

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

	### ----- generate data -----------###
	### data rich cases only
	rich_dir <- base_lowsigR_dir_vec[grepl("Index_Catch_LC20",base_lowsigR_dir_vec)]
	rich_modcombos <- base_modcombos[which(base_modcombos$Data_avail=="Index_Catch_LC20"),]
	alt_dir <- base_lowsigR_dir_vec[grepl("Index_Catch_LC20",base_lowsigR_dir_vec)==FALSE]
	alt_modcombos <- base_modcombos[which(base_modcombos$Data_avail!="Index_Catch_LC20"),]

	### ----- use parallel cores -----------###
	registerDoParallel(cores=ncores)	

	start_gen <- Sys.time()
	foreach(loop=1:length(rich_dir), .packages=c('TMB','LIME')) %dopar% generate_data(modpath=rich_dir[loop], data_avail=as.character(rich_modcombos[loop,"Data_avail"]), itervec=itervec, Fdynamics=as.character(strsplit(rich_modcombos[loop,"Fdyn"],"_")[[1]][2]), Rdynamics=as.character(strsplit(rich_modcombos[loop,"Rdyn"],"_")[[1]][2]), write=TRUE, lh=lh_list[[as.character(strsplit(rich_modcombos[loop,"LH"],"_")[[1]][2])]], Nyears=20, comp_sample=as.numeric(strsplit(rich_modcombos[loop,"ESS"],"_")[[1]][2]), rewrite=FALSE, init_depl=c(0.05,0.95))
	end_gen <- Sys.time() - start_gen

	# copy data rich cases to data poor
	# when written using parallel cores, sometimes files not writing properly
	start_regen <- Sys.time()
	for(loop in 1:length(rich_dir)){
		copy_sim(fromdir=rich_dir[loop], fromcombos=rich_modcombos[loop,], todir=alt_dir, itervec=itervec, rewrite=FALSE, res_dir="base_lowsigR")
	}
	end_regen <- Sys.time() - start_regen


	### ----- run LIME -----------###
	lime_dirs <- base_lowsigR_dir_vec[which(grepl("LBSPR",base_lowsigR_dir_vec)==FALSE)]
	lbspr_dirs <- base_lowsigR_dir_vec[which(grepl("LBSPR",base_lowsigR_dir_vec))]
	lime_combos <- base_modcombos[which(grepl("LBSPR", base_modcombos[,"Data_avail"])==FALSE),]
	lbspr_combos <- base_modcombos[which(grepl("LBSPR", base_modcombos[,"Data_avail"])),]	

	start_run <- Sys.time()
	foreach(loop=1:length(lime_dirs), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=lime_dirs[loop], lh=lh_list[[as.character(strsplit(lime_combos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(lime_combos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, fix_f=0, simulation=TRUE, REML=FALSE, f_true=FALSE, C_opt=lime_combos[loop,"C_opt"], Sel0=1, LFdist=1, param_adjust=c("SigmaF","SigmaR","SigmaC","SigmaI"), val_adjust=c(0.2,0.7,0.2,0.2), write=TRUE, fix_param=FALSE)
	end_run <- Sys.time() - start_run

	### ----- run LBSPR -----------###	

	start_run <- Sys.time()
	foreach(loop=1:length(lbspr_dirs), .packages=c('LBSPR','LIME')) %dopar% run_LBSPR(modpath=lbspr_dirs[loop], lh=lh_list[[as.character(strsplit(lbspr_combos[loop,"LH"],"_")[[1]][2])]], itervec=itervec, species=NULL, rewrite=TRUE, simulation=TRUE)
	end_altrun <- Sys.time() - start_run	

	### ------- check output ------###	
	dir_esslow <- base_lowsigR_dir_vec[grepl("ESS_100/", base_lowsigR_dir_vec)]
	dir_esshigh <- base_lowsigR_dir_vec[grepl("ESS_1000", base_lowsigR_dir_vec)]
	
	base_coverlow <- interval_coverage(dirs=dir_esslow[which(grepl("LBSPR",dir_esslow)==FALSE)], itervec=itervec)
	base_coverhigh <- interval_coverage(dirs=dir_esshigh[which(grepl("LBSPR",dir_esshigh)==FALSE)], itervec=itervec)


	png(file.path(fig_dir, "Base_lowsigR_ESSlow_RE.png"), height=10, width=25, res=200, units="in")
	base_res <- plot_re(dirs=dir_esslow, modcombos=base_modcombos[which(base_modcombos$ESS=="ESS_100"),], itervec=itervec, compareToLH=paste0("F_", Fdyn_vec), lh_vec=lh_vec, cover=base_coverlow$cover, ylim=c(-1,1))
	dev.off()

	png(file.path(fig_dir, "Base_lowsigR_ESShigh_RE.png"), height=10, width=25, res=200, units="in")
	base_res <- plot_re(dirs=dir_esshigh, modcombos=base_modcombos[which(base_modcombos$ESS=="ESS_1000"),], itervec=itervec, compareToLH=paste0("F_", Fdyn_vec), lh_vec=lh_vec, cover=base_coverhigh$cover, ylim=c(-1,1))
	dev.off()



#************* NEED TO UPDATE WITH LATEST CODE
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
	sens_equil_modcombos$C_opt[which(grepl("Catch",sens_equil_modcombos[,"Data_avail"]))] <- 2

### ----- sens runs ----------- ###
	sens_equil_dir <- file.path(main_dir, "sens_equil")
	# unlink(sens_equil_dir, TRUE)
	dir.create(sens_equil_dir, showWarnings=FALSE)	

	## setup dirs
	sens_equil_dir_vec <- model_paths(res_dir=sens_equil_dir, modcombos=sens_equil_modcombos[,-ncol(sens_equil_modcombos)])
	ignore <- sapply(1:length(sens_equil_dir_vec), function(m) sapply(itervec, function(x) dir.create(file.path(sens_equil_dir_vec[m], x), showWarnings=FALSE)))	

	## equilibrium
	lh_list <- adj_variation(SigmaR=0.01, SigmaF=0.01, SigmaC=0.01, SigmaI=0.01, CVlen=0.1, rho=0)

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
	init_dir <- equil_dir_vec[grepl("/LC10/",equil_dir_vec) & grepl("ESS_20",equil_dir_vec)==FALSE]
	init_modcombos <- equil_modcombos[which(equil_modcombos$Data_avail=="LC10" & equil_modcombos$ESS!="ESS_20"),]
	alt_dir <- sens_equil_dir_vec
	alt_modcombos <- sens_equil_modcombos

	start_regen <- Sys.time()
	for(loop in 1:length(init_dir)){
		copy_sim(fromdir=init_dir[loop], fromcombos=init_modcombos[loop,], todir=alt_dir, itervec=itervec, rewrite=FALSE, res_dir="sens_equil", sensitivity=TRUE)
	}
	end_regen <- Sys.time() - start_regen

	## run LIME
	start_run <- Sys.time()
	foreach(loop=1:length(sens_equil_dir_vec), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=sens_equil_dir_vec[loop], lh=lh_list[[as.character(strsplit(sens_equil_modcombos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(sens_equil_modcombos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, fix_f=0, simulation=TRUE, REML=FALSE, f_true=FALSE, C_opt=sens_equil_modcombos[loop,"C_opt"], param_adjust=sens_equil_modcombos[loop,"Param"], val_adjust=sens_vals[[as.character(strsplit(sens_equil_modcombos[loop,"LH"],"_")[[1]][2])]][sens_equil_modcombos[loop,"Adjust"],sens_equil_modcombos[loop,"Param"]], LFdist=1, Sel0=1)
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
	registerDoParallel(cores=ncores)	

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



##### sensitivity plot
######## HIGH ESS
png(file.path(fig_dir, "SensM_ESShigh_LC10.png"), height=8, width=14, units="in", res=200)
par(mfrow=c(2,3), mar=c(0,0,0,0), omi=c(1,1,1,1))
########## Equilibrium
## short
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\M\\low\\LC10\\ESS_1000\\LH_Short\\F_Constant\\R_Constant"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LC10\\ESS_1000\\LH_Short\\F_Constant\\R_Constant"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\M\\high\\LC10\\ESS_1000\\LH_Short\\F_Constant\\R_Constant"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, xaxt="n", yaxt="n")
axis(2, cex.axis=2, las=2)
mtext("Short", cex=2, side=3, line=1)

## medium
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\M\\low\\LC10\\ESS_1000\\LH_Medium\\F_Constant\\R_Constant"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LC10\\ESS_1000\\LH_Medium\\F_Constant\\R_Constant"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\M\\high\\LC10\\ESS_1000\\LH_Medium\\F_Constant\\R_Constant"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, xaxt="n", yaxt="n")
mtext("Medium", cex=2, side=3, line=1)

## long
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\M\\low\\LC10\\ESS_1000\\LH_Long\\F_Constant\\R_Constant"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LC10\\ESS_1000\\LH_Long\\F_Constant\\R_Constant"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\M\\high\\LC10\\ESS_1000\\LH_Long\\F_Constant\\R_Constant"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, xaxt="n", yaxt="n")
mtext("Long", cex=2, side=3, line=1)
mtext("Equilibrium", cex=2, side=4, line=1)


########## Base
## short
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\M\\low\\LC10\\ESS_1000\\LH_Short\\F_Endogenous\\R_AR"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LC10\\ESS_1000\\LH_Short\\F_Endogenous\\R_AR"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\M\\high\\LC10\\ESS_1000\\LH_Short\\F_Endogenous\\R_AR"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, yaxt="n", xaxt="n")
axis(2, cex.axis=2, las=2)
axis(1, at=1:3, labels=c("-25%", "true", "+25%"), cex.axis=2)

## medium
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\M\\low\\LC10\\ESS_1000\\LH_Medium\\F_Endogenous\\R_AR"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LC10\\ESS_1000\\LH_Medium\\F_Endogenous\\R_AR"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\M\\high\\LC10\\ESS_1000\\LH_Medium\\F_Endogenous\\R_AR"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, yaxt="n", xaxt="n")
axis(1, at=1:3, labels=c("-25%", "true", "+25%"), cex.axis=2)

## long
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\M\\low\\LC10\\ESS_1000\\LH_Long\\F_Endogenous\\R_AR"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LC10\\ESS_1000\\LH_Long\\F_Endogenous\\R_AR"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\M\\high\\LC10\\ESS_1000\\LH_Long\\F_Endogenous\\R_AR"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, yaxt="n", xaxt="n")
axis(1, at=1:3, labels=c("-25%", "true", "+25%"), cex.axis=2)
mtext("Variability", cex=2, side=4, line=1)
mtext("Estimation error", outer=TRUE, side=2, line=3, cex=1.5)
mtext("Fixed value", outer=TRUE, side=1, line=3, cex=1.5)
dev.off()


png(file.path(fig_dir, "SensLinf_ESShigh_LC10.png"), height=8, width=14, units="in", res=200)
##### sensitivity plot
par(mfrow=c(2,3), mar=c(0,0,0,0), omi=c(1,1,1,1))
########## Equilibrium
## short
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\linf\\low\\LC10\\ESS_1000\\LH_Short\\F_Constant\\R_Constant"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LC10\\ESS_1000\\LH_Short\\F_Constant\\R_Constant"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\linf\\high\\LC10\\ESS_1000\\LH_Short\\F_Constant\\R_Constant"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, xaxt="n", yaxt="n")
axis(2, cex.axis=2, las=2)
mtext("Short", cex=2, side=3, line=1)

## medium
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\linf\\low\\LC10\\ESS_1000\\LH_Medium\\F_Constant\\R_Constant"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LC10\\ESS_1000\\LH_Medium\\F_Constant\\R_Constant"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\linf\\high\\LC10\\ESS_1000\\LH_Medium\\F_Constant\\R_Constant"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, xaxt="n", yaxt="n")
mtext("Medium", cex=2, side=3, line=1)

## long
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\linf\\low\\LC10\\ESS_1000\\LH_Long\\F_Constant\\R_Constant"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LC10\\ESS_1000\\LH_Long\\F_Constant\\R_Constant"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\linf\\high\\LC10\\ESS_1000\\LH_Long\\F_Constant\\R_Constant"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, xaxt="n", yaxt="n")
mtext("Long", cex=2, side=3, line=1)
mtext("Equilibrium", cex=2, side=4, line=1)


########## Base
## short
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\linf\\low\\LC10\\ESS_1000\\LH_Short\\F_Endogenous\\R_AR"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LC10\\ESS_1000\\LH_Short\\F_Endogenous\\R_AR"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\linf\\high\\LC10\\ESS_1000\\LH_Short\\F_Endogenous\\R_AR"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-10,4), itervec=itervec, yaxt="n", xaxt="n")
axis(2, cex.axis=2, las=2)
axis(1, at=1:3, labels=c("-25%", "true", "+25%"), cex.axis=2)

## medium
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\linf\\low\\LC10\\ESS_1000\\LH_Medium\\F_Endogenous\\R_AR"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LC10\\ESS_1000\\LH_Medium\\F_Endogenous\\R_AR"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\linf\\high\\LC10\\ESS_1000\\LH_Medium\\F_Endogenous\\R_AR"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-10,4), itervec=itervec, yaxt="n", xaxt="n")
axis(1, at=1:3, labels=c("-25%", "true", "+25%"), cex.axis=2)

## long
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\linf\\low\\LC10\\ESS_1000\\LH_Long\\F_Endogenous\\R_AR"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LC10\\ESS_1000\\LH_Long\\F_Endogenous\\R_AR"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\linf\\high\\LC10\\ESS_1000\\LH_Long\\F_Endogenous\\R_AR"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-10,4), itervec=itervec, yaxt="n", xaxt="n")
axis(1, at=1:3, labels=c("-25%", "true", "+25%"), cex.axis=2)
mtext("Variability", cex=2, side=4, line=1)
mtext("Estimation error", outer=TRUE, side=2, line=3, cex=1.5)
mtext("Fixed value", outer=TRUE, side=1, line=3, cex=1.5)
dev.off()


png(file.path(fig_dir, "SensVBK_ESShigh_LC10.png"), height=8, width=14, units="in", res=200)
##### sensitivity plot
par(mfrow=c(2,3), mar=c(0,0,0,0), omi=c(1,1,1,1))
########## Equilibrium
## short
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\vbk\\low\\LC10\\ESS_1000\\LH_Short\\F_Constant\\R_Constant"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LC10\\ESS_1000\\LH_Short\\F_Constant\\R_Constant"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\vbk\\high\\LC10\\ESS_1000\\LH_Short\\F_Constant\\R_Constant"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, xaxt="n", yaxt="n")
axis(2, cex.axis=2, las=2)
mtext("Short", cex=2, side=3, line=1)

## medium
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\vbk\\low\\LC10\\ESS_1000\\LH_Medium\\F_Constant\\R_Constant"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LC10\\ESS_1000\\LH_Medium\\F_Constant\\R_Constant"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\vbk\\high\\LC10\\ESS_1000\\LH_Medium\\F_Constant\\R_Constant"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, xaxt="n", yaxt="n")
mtext("Medium", cex=2, side=3, line=1)

## long
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\vbk\\low\\LC10\\ESS_1000\\LH_Long\\F_Constant\\R_Constant"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LC10\\ESS_1000\\LH_Long\\F_Constant\\R_Constant"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\vbk\\high\\LC10\\ESS_1000\\LH_Long\\F_Constant\\R_Constant"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, xaxt="n", yaxt="n")
mtext("Long", cex=2, side=3, line=1)
mtext("Equilibrium", cex=2, side=4, line=1)


########## Base
## short
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\vbk\\low\\LC10\\ESS_1000\\LH_Short\\F_Endogenous\\R_AR"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LC10\\ESS_1000\\LH_Short\\F_Endogenous\\R_AR"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\vbk\\high\\LC10\\ESS_1000\\LH_Short\\F_Endogenous\\R_AR"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, yaxt="n", xaxt="n")
axis(2, cex.axis=2, las=2)
axis(1, at=1:3, labels=c("-25%", "true", "+25%"), cex.axis=2)

## medium
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\vbk\\low\\LC10\\ESS_1000\\LH_Medium\\F_Endogenous\\R_AR"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LC10\\ESS_1000\\LH_Medium\\F_Endogenous\\R_AR"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\vbk\\high\\LC10\\ESS_1000\\LH_Medium\\F_Endogenous\\R_AR"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, yaxt="n", xaxt="n")
axis(1, at=1:3, labels=c("-25%", "true", "+25%"), cex.axis=2)

## long
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\vbk\\low\\LC10\\ESS_1000\\LH_Long\\F_Endogenous\\R_AR"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LC10\\ESS_1000\\LH_Long\\F_Endogenous\\R_AR"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\vbk\\high\\LC10\\ESS_1000\\LH_Long\\F_Endogenous\\R_AR"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, yaxt="n", xaxt="n")
axis(1, at=1:3, labels=c("-25%", "true", "+25%"), cex.axis=2)
mtext("Variability", cex=2, side=4, line=1)
mtext("Estimation error", outer=TRUE, side=2, line=3, cex=1.5)
mtext("Fixed value", outer=TRUE, side=1, line=3, cex=1.5)
dev.off()




png(file.path(fig_dir, "SensML50_ESShigh_LC10.png"), height=8, width=14, units="in", res=200)
##### sensitivity plot
par(mfrow=c(2,3), mar=c(0,0,0,0), omi=c(1,1,1,1))
########## Equilibrium
## short
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\ML50\\low\\LC10\\ESS_1000\\LH_Short\\F_Constant\\R_Constant"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LC10\\ESS_1000\\LH_Short\\F_Constant\\R_Constant"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\ML50\\high\\LC10\\ESS_1000\\LH_Short\\F_Constant\\R_Constant"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, xaxt="n", yaxt="n")
axis(2, cex.axis=2, las=2)
mtext("Short", cex=2, side=3, line=1)

## medium
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\ML50\\low\\LC10\\ESS_1000\\LH_Medium\\F_Constant\\R_Constant"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LC10\\ESS_1000\\LH_Medium\\F_Constant\\R_Constant"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\ML50\\high\\LC10\\ESS_1000\\LH_Medium\\F_Constant\\R_Constant"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, xaxt="n", yaxt="n")
mtext("Medium", cex=2, side=3, line=1)

## long
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\ML50\\low\\LC10\\ESS_1000\\LH_Long\\F_Constant\\R_Constant"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LC10\\ESS_1000\\LH_Long\\F_Constant\\R_Constant"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\ML50\\high\\LC10\\ESS_1000\\LH_Long\\F_Constant\\R_Constant"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, xaxt="n", yaxt="n")
mtext("Long", cex=2, side=3, line=1)
mtext("Equilibrium", cex=2, side=4, line=1)


########## Base
## short
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\ML50\\low\\LC10\\ESS_1000\\LH_Short\\F_Endogenous\\R_AR"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LC10\\ESS_1000\\LH_Short\\F_Endogenous\\R_AR"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\ML50\\high\\LC10\\ESS_1000\\LH_Short\\F_Endogenous\\R_AR"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, yaxt="n", xaxt="n")
axis(2, cex.axis=2, las=2)
axis(1, at=1:3, labels=c("-25%", "true", "+25%"), cex.axis=2)

## medium
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\ML50\\low\\LC10\\ESS_1000\\LH_Medium\\F_Endogenous\\R_AR"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LC10\\ESS_1000\\LH_Medium\\F_Endogenous\\R_AR"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\ML50\\high\\LC10\\ESS_1000\\LH_Medium\\F_Endogenous\\R_AR"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, yaxt="n", xaxt="n")
axis(1, at=1:3, labels=c("-25%", "true", "+25%"), cex.axis=2)

## long
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\ML50\\low\\LC10\\ESS_1000\\LH_Long\\F_Endogenous\\R_AR"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LC10\\ESS_1000\\LH_Long\\F_Endogenous\\R_AR"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\ML50\\high\\LC10\\ESS_1000\\LH_Long\\F_Endogenous\\R_AR"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, yaxt="n", xaxt="n")
axis(1, at=1:3, labels=c("-25%", "true", "+25%"), cex.axis=2)
mtext("Variability", cex=2, side=4, line=1)
mtext("Estimation error", outer=TRUE, side=2, line=3, cex=1.5)
mtext("Fixed value", outer=TRUE, side=1, line=3, cex=1.5)
dev.off()


png(file.path(fig_dir, "SensCVlen_ESShigh_LC10.png"), height=8, width=14, units="in", res=200)
##### sensitivity plot
par(mfrow=c(2,3), mar=c(0,0,0,0), omi=c(1,1,1,1))
########## Equilibrium
## short
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\CVlen\\low\\LC10\\ESS_1000\\LH_Short\\F_Constant\\R_Constant"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LC10\\ESS_1000\\LH_Short\\F_Constant\\R_Constant"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\CVlen\\high\\LC10\\ESS_1000\\LH_Short\\F_Constant\\R_Constant"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, xaxt="n", yaxt="n")
axis(2, cex.axis=2, las=2)
mtext("Short", cex=2, side=3, line=1)

## medium
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\CVlen\\low\\LC10\\ESS_1000\\LH_Medium\\F_Constant\\R_Constant"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LC10\\ESS_1000\\LH_Medium\\F_Constant\\R_Constant"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\CVlen\\high\\LC10\\ESS_1000\\LH_Medium\\F_Constant\\R_Constant"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, xaxt="n", yaxt="n")
mtext("Medium", cex=2, side=3, line=1)

## long
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\CVlen\\low\\LC10\\ESS_1000\\LH_Long\\F_Constant\\R_Constant"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LC10\\ESS_1000\\LH_Long\\F_Constant\\R_Constant"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\CVlen\\high\\LC10\\ESS_1000\\LH_Long\\F_Constant\\R_Constant"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, xaxt="n", yaxt="n")
mtext("Long", cex=2, side=3, line=1)
mtext("Equilibrium", cex=2, side=4, line=1)


########## Base
## short
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\CVlen\\low\\LC10\\ESS_1000\\LH_Short\\F_Endogenous\\R_AR"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LC10\\ESS_1000\\LH_Short\\F_Endogenous\\R_AR"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\CVlen\\high\\LC10\\ESS_1000\\LH_Short\\F_Endogenous\\R_AR"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, yaxt="n", xaxt="n")
axis(2, cex.axis=2, las=2)
axis(1, at=1:3, labels=c("-25%", "true", "+25%"), cex.axis=2)

## medium
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\CVlen\\low\\LC10\\ESS_1000\\LH_Medium\\F_Endogenous\\R_AR"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LC10\\ESS_1000\\LH_Medium\\F_Endogenous\\R_AR"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\CVlen\\high\\LC10\\ESS_1000\\LH_Medium\\F_Endogenous\\R_AR"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, yaxt="n", xaxt="n")
axis(1, at=1:3, labels=c("-25%", "true", "+25%"), cex.axis=2)

## long
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\CVlen\\low\\LC10\\ESS_1000\\LH_Long\\F_Endogenous\\R_AR"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LC10\\ESS_1000\\LH_Long\\F_Endogenous\\R_AR"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\CVlen\\high\\LC10\\ESS_1000\\LH_Long\\F_Endogenous\\R_AR"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, yaxt="n", xaxt="n")
axis(1, at=1:3, labels=c("-25%", "true", "+25%"), cex.axis=2)
mtext("Variability", cex=2, side=4, line=1)
mtext("Estimation error", outer=TRUE, side=2, line=3, cex=1.5)
mtext("Fixed value", outer=TRUE, side=1, line=3, cex=1.5)
dev.off()



#### LOW ESS
png(file.path(fig_dir, "SensM_ESSlow_LC10.png"), height=8, width=14, units="in", res=200)
##### sensitivity plot
par(mfrow=c(2,3), mar=c(0,0,0,0), omi=c(1,1,1,1))
########## Equilibrium
## short
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\M\\low\\LC10\\ESS_100\\LH_Short\\F_Constant\\R_Constant"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LC10\\ESS_100\\LH_Short\\F_Constant\\R_Constant"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\M\\high\\LC10\\ESS_100\\LH_Short\\F_Constant\\R_Constant"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, xaxt="n", yaxt="n")
axis(2, cex.axis=2, las=2)
mtext("Short", cex=2, side=3, line=1)

## medium
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\M\\low\\LC10\\ESS_100\\LH_Medium\\F_Constant\\R_Constant"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LC10\\ESS_100\\LH_Medium\\F_Constant\\R_Constant"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\M\\high\\LC10\\ESS_100\\LH_Medium\\F_Constant\\R_Constant"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, xaxt="n", yaxt="n")
mtext("Medium", cex=2, side=3, line=1)

## long
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\M\\low\\LC10\\ESS_100\\LH_Long\\F_Constant\\R_Constant"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LC10\\ESS_100\\LH_Long\\F_Constant\\R_Constant"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\M\\high\\LC10\\ESS_100\\LH_Long\\F_Constant\\R_Constant"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, xaxt="n", yaxt="n")
mtext("Long", cex=2, side=3, line=1)
mtext("Equilibrium", cex=2, side=4, line=1)


########## Base
## short
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\M\\low\\LC10\\ESS_100\\LH_Short\\F_Endogenous\\R_AR"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LC10\\ESS_100\\LH_Short\\F_Endogenous\\R_AR"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\M\\high\\LC10\\ESS_100\\LH_Short\\F_Endogenous\\R_AR"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, yaxt="n", xaxt="n")
axis(2, cex.axis=2, las=2)
axis(1, at=1:3, labels=c("-25%", "true", "+25%"), cex.axis=2)

## medium
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\M\\low\\LC10\\ESS_100\\LH_Medium\\F_Endogenous\\R_AR"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LC10\\ESS_100\\LH_Medium\\F_Endogenous\\R_AR"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\M\\high\\LC10\\ESS_100\\LH_Medium\\F_Endogenous\\R_AR"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, yaxt="n", xaxt="n")
axis(1, at=1:3, labels=c("-25%", "true", "+25%"), cex.axis=2)

## long
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\M\\low\\LC10\\ESS_100\\LH_Long\\F_Endogenous\\R_AR"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LC10\\ESS_100\\LH_Long\\F_Endogenous\\R_AR"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\M\\high\\LC10\\ESS_100\\LH_Long\\F_Endogenous\\R_AR"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, yaxt="n", xaxt="n")
axis(1, at=1:3, labels=c("-25%", "true", "+25%"), cex.axis=2)
mtext("Variability", cex=2, side=4, line=1)
mtext("Estimation error", outer=TRUE, side=2, line=3, cex=1.5)
mtext("Fixed value", outer=TRUE, side=1, line=3, cex=1.5)
dev.off()


png(file.path(fig_dir, "SensLinf_ESSlow_LC10.png"), height=8, width=14, units="in", res=200)
##### sensitivity plot
par(mfrow=c(2,3), mar=c(0,0,0,0), omi=c(1,1,1,1))
########## Equilibrium
## short
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\linf\\low\\LC10\\ESS_100\\LH_Short\\F_Constant\\R_Constant"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LC10\\ESS_100\\LH_Short\\F_Constant\\R_Constant"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\linf\\high\\LC10\\ESS_100\\LH_Short\\F_Constant\\R_Constant"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, xaxt="n", yaxt="n")
axis(2, cex.axis=2, las=2)
mtext("Short", cex=2, side=3, line=1)

## medium
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\linf\\low\\LC10\\ESS_100\\LH_Medium\\F_Constant\\R_Constant"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LC10\\ESS_100\\LH_Medium\\F_Constant\\R_Constant"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\linf\\high\\LC10\\ESS_100\\LH_Medium\\F_Constant\\R_Constant"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, xaxt="n", yaxt="n")
mtext("Medium", cex=2, side=3, line=1)

## long
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\linf\\low\\LC10\\ESS_100\\LH_Long\\F_Constant\\R_Constant"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LC10\\ESS_100\\LH_Long\\F_Constant\\R_Constant"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\linf\\high\\LC10\\ESS_100\\LH_Long\\F_Constant\\R_Constant"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, xaxt="n", yaxt="n")
mtext("Long", cex=2, side=3, line=1)
mtext("Equilibrium", cex=2, side=4, line=1)


########## Base
## short
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\linf\\low\\LC10\\ESS_100\\LH_Short\\F_Endogenous\\R_AR"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LC10\\ESS_100\\LH_Short\\F_Endogenous\\R_AR"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\linf\\high\\LC10\\ESS_100\\LH_Short\\F_Endogenous\\R_AR"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, yaxt="n", xaxt="n")
axis(2, cex.axis=2, las=2)
axis(1, at=1:3, labels=c("-25%", "true", "+25%"), cex.axis=2)

## medium
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\linf\\low\\LC10\\ESS_100\\LH_Medium\\F_Endogenous\\R_AR"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LC10\\ESS_100\\LH_Medium\\F_Endogenous\\R_AR"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\linf\\high\\LC10\\ESS_100\\LH_Medium\\F_Endogenous\\R_AR"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, yaxt="n", xaxt="n")
axis(1, at=1:3, labels=c("-25%", "true", "+25%"), cex.axis=2)

## long
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\linf\\low\\LC10\\ESS_100\\LH_Long\\F_Endogenous\\R_AR"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LC10\\ESS_100\\LH_Long\\F_Endogenous\\R_AR"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\linf\\high\\LC10\\ESS_100\\LH_Long\\F_Endogenous\\R_AR"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, yaxt="n", xaxt="n")
axis(1, at=1:3, labels=c("-25%", "true", "+25%"), cex.axis=2)
mtext("Variability", cex=2, side=4, line=1)
mtext("Estimation error", outer=TRUE, side=2, line=3, cex=1.5)
mtext("Fixed value", outer=TRUE, side=1, line=3, cex=1.5)
dev.off()


png(file.path(fig_dir, "SensVBK_ESSlow_LC10.png"), height=8, width=14, units="in", res=200)
##### sensitivity plot
par(mfrow=c(2,3), mar=c(0,0,0,0), omi=c(1,1,1,1))
########## Equilibrium
## short
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\vbk\\low\\LC10\\ESS_100\\LH_Short\\F_Constant\\R_Constant"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LC10\\ESS_100\\LH_Short\\F_Constant\\R_Constant"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\vbk\\high\\LC10\\ESS_100\\LH_Short\\F_Constant\\R_Constant"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, xaxt="n", yaxt="n")
axis(2, cex.axis=2, las=2)
mtext("Short", cex=2, side=3, line=1)

## medium
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\vbk\\low\\LC10\\ESS_100\\LH_Medium\\F_Constant\\R_Constant"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LC10\\ESS_100\\LH_Medium\\F_Constant\\R_Constant"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\vbk\\high\\LC10\\ESS_100\\LH_Medium\\F_Constant\\R_Constant"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, xaxt="n", yaxt="n")
mtext("Medium", cex=2, side=3, line=1)

## long
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\vbk\\low\\LC10\\ESS_100\\LH_Long\\F_Constant\\R_Constant"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LC10\\ESS_100\\LH_Long\\F_Constant\\R_Constant"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\vbk\\high\\LC10\\ESS_100\\LH_Long\\F_Constant\\R_Constant"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, xaxt="n", yaxt="n")
mtext("Long", cex=2, side=3, line=1)
mtext("Equilibrium", cex=2, side=4, line=1)


########## Base
## short
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\vbk\\low\\LC10\\ESS_100\\LH_Short\\F_Endogenous\\R_AR"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LC10\\ESS_100\\LH_Short\\F_Endogenous\\R_AR"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\vbk\\high\\LC10\\ESS_100\\LH_Short\\F_Endogenous\\R_AR"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, yaxt="n", xaxt="n")
axis(2, cex.axis=2, las=2)
axis(1, at=1:3, labels=c("-25%", "true", "+25%"), cex.axis=2)

## medium
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\vbk\\low\\LC10\\ESS_100\\LH_Medium\\F_Endogenous\\R_AR"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LC10\\ESS_100\\LH_Medium\\F_Endogenous\\R_AR"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\vbk\\high\\LC10\\ESS_100\\LH_Medium\\F_Endogenous\\R_AR"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, yaxt="n", xaxt="n")
axis(1, at=1:3, labels=c("-25%", "true", "+25%"), cex.axis=2)

## long
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\vbk\\low\\LC10\\ESS_100\\LH_Long\\F_Endogenous\\R_AR"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LC10\\ESS_100\\LH_Long\\F_Endogenous\\R_AR"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\vbk\\high\\LC10\\ESS_100\\LH_Long\\F_Endogenous\\R_AR"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, yaxt="n", xaxt="n")
axis(1, at=1:3, labels=c("-25%", "true", "+25%"), cex.axis=2)
mtext("Variability", cex=2, side=4, line=1)
mtext("Estimation error", outer=TRUE, side=2, line=3, cex=1.5)
mtext("Fixed value", outer=TRUE, side=1, line=3, cex=1.5)
dev.off()




png(file.path(fig_dir, "SensML50_ESSlow_LC10.png"), height=8, width=14, units="in", res=200)
##### sensitivity plot
par(mfrow=c(2,3), mar=c(0,0,0,0), omi=c(1,1,1,1))
########## Equilibrium
## short
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\ML50\\low\\LC10\\ESS_100\\LH_Short\\F_Constant\\R_Constant"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LC10\\ESS_100\\LH_Short\\F_Constant\\R_Constant"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\ML50\\high\\LC10\\ESS_100\\LH_Short\\F_Constant\\R_Constant"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, xaxt="n", yaxt="n")
axis(2, cex.axis=2, las=2)
mtext("Short", cex=2, side=3, line=1)

## medium
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\ML50\\low\\LC10\\ESS_100\\LH_Medium\\F_Constant\\R_Constant"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LC10\\ESS_100\\LH_Medium\\F_Constant\\R_Constant"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\ML50\\high\\LC10\\ESS_100\\LH_Medium\\F_Constant\\R_Constant"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, xaxt="n", yaxt="n")
mtext("Medium", cex=2, side=3, line=1)

## long
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\ML50\\low\\LC10\\ESS_100\\LH_Long\\F_Constant\\R_Constant"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LC10\\ESS_100\\LH_Long\\F_Constant\\R_Constant"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\ML50\\high\\LC10\\ESS_100\\LH_Long\\F_Constant\\R_Constant"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, xaxt="n", yaxt="n")
mtext("Long", cex=2, side=3, line=1)
mtext("Equilibrium", cex=2, side=4, line=1)


########## Base
## short
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\ML50\\low\\LC10\\ESS_100\\LH_Short\\F_Endogenous\\R_AR"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LC10\\ESS_100\\LH_Short\\F_Endogenous\\R_AR"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\ML50\\high\\LC10\\ESS_100\\LH_Short\\F_Endogenous\\R_AR"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, yaxt="n", xaxt="n")
axis(2, cex.axis=2, las=2)
axis(1, at=1:3, labels=c("-25%", "true", "+25%"), cex.axis=2)

## medium
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\ML50\\low\\LC10\\ESS_100\\LH_Medium\\F_Endogenous\\R_AR"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LC10\\ESS_100\\LH_Medium\\F_Endogenous\\R_AR"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\ML50\\high\\LC10\\ESS_100\\LH_Medium\\F_Endogenous\\R_AR"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, yaxt="n", xaxt="n")
axis(1, at=1:3, labels=c("-25%", "true", "+25%"), cex.axis=2)

## long
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\ML50\\low\\LC10\\ESS_100\\LH_Long\\F_Endogenous\\R_AR"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LC10\\ESS_100\\LH_Long\\F_Endogenous\\R_AR"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\ML50\\high\\LC10\\ESS_100\\LH_Long\\F_Endogenous\\R_AR"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, yaxt="n", xaxt="n")
axis(1, at=1:3, labels=c("-25%", "true", "+25%"), cex.axis=2)
mtext("Variability", cex=2, side=4, line=1)
mtext("Estimation error", outer=TRUE, side=2, line=3, cex=1.5)
mtext("Fixed value", outer=TRUE, side=1, line=3, cex=1.5)
dev.off()


png(file.path(fig_dir, "SensCVlen_ESSlow_LC10.png"), height=8, width=14, units="in", res=200)
##### sensitivity plot
par(mfrow=c(2,3), mar=c(0,0,0,0), omi=c(1,1,1,1))
########## Equilibrium
## short
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\CVlen\\low\\LC10\\ESS_100\\LH_Short\\F_Constant\\R_Constant"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LC10\\ESS_100\\LH_Short\\F_Constant\\R_Constant"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\CVlen\\high\\LC10\\ESS_100\\LH_Short\\F_Constant\\R_Constant"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, xaxt="n", yaxt="n")
axis(2, cex.axis=2, las=2)
mtext("Short", cex=2, side=3, line=1)

## medium
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\CVlen\\low\\LC10\\ESS_100\\LH_Medium\\F_Constant\\R_Constant"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LC10\\ESS_100\\LH_Medium\\F_Constant\\R_Constant"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\CVlen\\high\\LC10\\ESS_100\\LH_Medium\\F_Constant\\R_Constant"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, xaxt="n", yaxt="n")
mtext("Medium", cex=2, side=3, line=1)

## long
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\CVlen\\low\\LC10\\ESS_100\\LH_Long\\F_Constant\\R_Constant"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LC10\\ESS_100\\LH_Long\\F_Constant\\R_Constant"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_equil\\CVlen\\high\\LC10\\ESS_100\\LH_Long\\F_Constant\\R_Constant"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, xaxt="n", yaxt="n")
mtext("Long", cex=2, side=3, line=1)
mtext("Equilibrium", cex=2, side=4, line=1)


########## Base
## short
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\CVlen\\low\\LC10\\ESS_100\\LH_Short\\F_Endogenous\\R_AR"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LC10\\ESS_100\\LH_Short\\F_Endogenous\\R_AR"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\CVlen\\high\\LC10\\ESS_100\\LH_Short\\F_Endogenous\\R_AR"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, yaxt="n", xaxt="n")
axis(2, cex.axis=2, las=2)
axis(1, at=1:3, labels=c("-25%", "true", "+25%"), cex.axis=2)

## medium
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\CVlen\\low\\LC10\\ESS_100\\LH_Medium\\F_Endogenous\\R_AR"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LC10\\ESS_100\\LH_Medium\\F_Endogenous\\R_AR"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\CVlen\\high\\LC10\\ESS_100\\LH_Medium\\F_Endogenous\\R_AR"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, yaxt="n", xaxt="n")
axis(1, at=1:3, labels=c("-25%", "true", "+25%"), cex.axis=2)

## long
path1 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\CVlen\\low\\LC10\\ESS_100\\LH_Long\\F_Endogenous\\R_AR"
path2 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LC10\\ESS_100\\LH_Long\\F_Endogenous\\R_AR"
path3 <- "F:\\Merrill\\Git_Projects\\LIME_sim\\sens_base\\CVlen\\high\\LC10\\ESS_100\\LH_Long\\F_Endogenous\\R_AR"
plot_re2(dirs=c(path1,path2,path3), ylim=c(-2.5,2.5), itervec=itervec, yaxt="n", xaxt="n")
axis(1, at=1:3, labels=c("-25%", "true", "+25%"), cex.axis=2)
mtext("Variability", cex=2, side=4, line=1)
mtext("Estimation error", outer=TRUE, side=2, line=3, cex=1.5)
mtext("Fixed value", outer=TRUE, side=1, line=3, cex=1.5)
dev.off()

############################################################
#### check model fits
############################################################

	set.seed(143)
	iter <- sample(1:100,1)
########### high ESS
##### Medium-lived
## equilibrium
LIME_dir <- file.path("F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LC10\\ESS_1000\\LH_Medium\\F_Constant\\R_Constant")
LBSPR_dir <- file.path("F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LBSPR10\\ESS_1000\\LH_Medium\\F_Constant\\R_Constant")

res <- plot_output(LIME_dir=LIME_dir, LBSPR_dir=LBSPR_dir, sim=TRUE, iter=iter)

## base
LIME_dir <- file.path("F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LC10\\ESS_1000\\LH_Medium\\F_Endogenous\\R_AR")
LBSPR_dir <- file.path("F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LBSPR10\\ESS_1000\\LH_Medium\\F_Endogenous\\R_AR")

res <- plot_output(LIME_dir=LIME_dir, LBSPR_dir=LBSPR_dir, sim=TRUE, iter=iter)

## base - lower sigmaR
LIME_dir <- file.path("F:\\Merrill\\Git_Projects\\LIME_sim\\base_lowsigR\\LC10\\ESS_1000\\LH_Medium\\F_Endogenous\\R_AR")
LBSPR_dir <- file.path("F:\\Merrill\\Git_Projects\\LIME_sim\\base_lowsigR\\LBSPR10\\ESS_1000\\LH_Medium\\F_Endogenous\\R_AR")

res <- plot_output(LIME_dir=LIME_dir, LBSPR_dir=LBSPR_dir, sim=TRUE, iter=iter)


##### Short-lived
## equilibrium
LIME_dir <- file.path("F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LC10\\ESS_1000\\LH_Short\\F_Constant\\R_Constant")
LBSPR_dir <- file.path("F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LBSPR10\\ESS_1000\\LH_Short\\F_Constant\\R_Constant")

res <- plot_output(LIME_dir=LIME_dir, LBSPR_dir=LBSPR_dir, sim=TRUE, iter=3)

## base ---*** HITTING UPPER BOUND ON ESTIMATION OF SIGMAR, NON-CONVERGENCE
LIME_dir <- file.path("F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LC10\\ESS_1000\\LH_Short\\F_Endogenous\\R_AR")
LBSPR_dir <- file.path("F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LBSPR10\\ESS_1000\\LH_Short\\F_Endogenous\\R_AR")

res <- plot_output(LIME_dir=LIME_dir, LBSPR_dir=LBSPR_dir, sim=TRUE, iter=15)

## base - lower sigmaR
LIME_dir <- file.path("F:\\Merrill\\Git_Projects\\LIME_sim\\base_lowsigR\\LC10\\ESS_1000\\LH_Short\\F_Endogenous\\R_AR")
LBSPR_dir <- file.path("F:\\Merrill\\Git_Projects\\LIME_sim\\base_lowsigR\\LBSPR10\\ESS_1000\\LH_Short\\F_Endogenous\\R_AR")

res <- plot_output(LIME_dir=LIME_dir, LBSPR_dir=LBSPR_dir, sim=TRUE, iter=15)


##### Long-lived
## equilibrium
LIME_dir <- file.path("F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LC10\\ESS_1000\\LH_Long\\F_Constant\\R_Constant")
LBSPR_dir <- file.path("F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LBSPR10\\ESS_1000\\LH_Long\\F_Constant\\R_Constant")

res <- plot_output(LIME_dir=LIME_dir, LBSPR_dir=LBSPR_dir, sim=TRUE, iter=iter)

## base
LIME_dir <- file.path("F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LC10\\ESS_1000\\LH_Long\\F_Endogenous\\R_AR")
LBSPR_dir <- file.path("F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LBSPR10\\ESS_1000\\LH_Long\\F_Endogenous\\R_AR")

res <- plot_output(LIME_dir=LIME_dir, LBSPR_dir=LBSPR_dir, sim=TRUE, iter=iter)

## base - lower sigmaR
LIME_dir <- file.path("F:\\Merrill\\Git_Projects\\LIME_sim\\base_lowsigR\\LC10\\ESS_1000\\LH_Long\\F_Endogenous\\R_AR")
LBSPR_dir <- file.path("F:\\Merrill\\Git_Projects\\LIME_sim\\base_lowsigR\\LBSPR10\\ESS_1000\\LH_Long\\F_Endogenous\\R_AR")

res <- plot_output(LIME_dir=LIME_dir, LBSPR_dir=LBSPR_dir, sim=TRUE, iter=iter)


########### lower ESS
##### Medium-lived
## equilibrium
LIME_dir <- file.path("F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LC10\\ESS_100\\LH_Medium\\F_Constant\\R_Constant")
LBSPR_dir <- file.path("F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LBSPR10\\ESS_100\\LH_Medium\\F_Constant\\R_Constant")

res <- plot_output(LIME_dir=LIME_dir, LBSPR_dir=LBSPR_dir, sim=TRUE, iter=iter)

## base
LIME_dir <- file.path("F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LC10\\ESS_100\\LH_Medium\\F_Endogenous\\R_AR")
LBSPR_dir <- file.path("F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LBSPR10\\ESS_100\\LH_Medium\\F_Endogenous\\R_AR")

res <- plot_output(LIME_dir=LIME_dir, LBSPR_dir=LBSPR_dir, sim=TRUE, iter=iter)

## base - lower sigmaR
LIME_dir <- file.path("F:\\Merrill\\Git_Projects\\LIME_sim\\base_lowsigR\\LC10\\ESS_100\\LH_Medium\\F_Endogenous\\R_AR")
LBSPR_dir <- file.path("F:\\Merrill\\Git_Projects\\LIME_sim\\base_lowsigR\\LBSPR10\\ESS_100\\LH_Medium\\F_Endogenous\\R_AR")

res <- plot_output(LIME_dir=LIME_dir, LBSPR_dir=LBSPR_dir, sim=TRUE, iter=iter)


##### Short-lived
## equilibrium
LIME_dir <- file.path("F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LC10\\ESS_100\\LH_Short\\F_Constant\\R_Constant")
LBSPR_dir <- file.path("F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LBSPR10\\ESS_100\\LH_Short\\F_Constant\\R_Constant")

res <- plot_output(LIME_dir=LIME_dir, LBSPR_dir=LBSPR_dir, sim=TRUE, iter=iter)

## base ---*** HITTING UPPER BOUND ON ESTIMATION OF SIGMAR, ESTIMATING F VERY HIGH (although truth is very high)
LIME_dir <- file.path("F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LC10\\ESS_100\\LH_Short\\F_Endogenous\\R_AR")
LBSPR_dir <- file.path("F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LBSPR10\\ESS_100\\LH_Short\\F_Endogenous\\R_AR")

res <- plot_output(LIME_dir=LIME_dir, LBSPR_dir=LBSPR_dir, sim=TRUE, iter=iter)

## base - lower sigmaR
LIME_dir <- file.path("F:\\Merrill\\Git_Projects\\LIME_sim\\base_lowsigR\\LC10\\ESS_100\\LH_Short\\F_Endogenous\\R_AR")
LBSPR_dir <- file.path("F:\\Merrill\\Git_Projects\\LIME_sim\\base_lowsigR\\LBSPR10\\ESS_100\\LH_Short\\F_Endogenous\\R_AR")

res <- plot_output(LIME_dir=LIME_dir, LBSPR_dir=LBSPR_dir, sim=TRUE, iter=iter)


##### Long-lived
## equilibrium
LIME_dir <- file.path("F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LC10\\ESS_100\\LH_Long\\F_Constant\\R_Constant")
LBSPR_dir <- file.path("F:\\Merrill\\Git_Projects\\LIME_sim\\equil\\LBSPR10\\ESS_100\\LH_Long\\F_Constant\\R_Constant")

res <- plot_output(LIME_dir=LIME_dir, LBSPR_dir=LBSPR_dir, sim=TRUE, iter=iter)

## base
LIME_dir <- file.path("F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LC10\\ESS_100\\LH_Long\\F_Endogenous\\R_AR")
LBSPR_dir <- file.path("F:\\Merrill\\Git_Projects\\LIME_sim\\base\\LBSPR10\\ESS_100\\LH_Long\\F_Endogenous\\R_AR")

res <- plot_output(LIME_dir=LIME_dir, LBSPR_dir=LBSPR_dir, sim=TRUE, iter=iter)

## base - lower sigmaR
LIME_dir <- file.path("F:\\Merrill\\Git_Projects\\LIME_sim\\base_lowsigR\\LC10\\ESS_100\\LH_Long\\F_Endogenous\\R_AR")
LBSPR_dir <- file.path("F:\\Merrill\\Git_Projects\\LIME_sim\\base_lowsigR\\LBSPR10\\ESS_100\\LH_Long\\F_Endogenous\\R_AR")

res <- plot_output(LIME_dir=LIME_dir, LBSPR_dir=LBSPR_dir, sim=TRUE, iter=iter)








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

		ignore <- generate_data(modpath=dir, data_avail=as.character(modcombos[loop,"Data_avail"]), itervec=itervec, Fdynamics=as.character(modcombos[loop,"Fdyn"]), Rdynamics="Constant", write=TRUE, lh=lh_list[[as.character(modcombos[loop,"LH"])]], Nyears=info$Nyears, comp_sample=info$comp_sample, rewrite=FALSE)

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