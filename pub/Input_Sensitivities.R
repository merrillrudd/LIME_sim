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
if(grepl("F:", main_dir)) ncores <- 5

funs <- list.files(file.path(main_dir, "R"))
ignore <- sapply(1:length(funs), function(x) source(file.path(main_dir,"R", funs[x])))

fig_dir <- file.path(main_dir, "figs")
dir.create(fig_dir, showWarnings=FALSE)

res_dir <- file.path(main_dir, "calc_results")
dir.create(res_dir, showWarnings=FALSE)


############################################################
#### find base run directories
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
	registerDoParallel(cores=ncores)	

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
	foreach(loop=1:length(sens_equil_dir_vec), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=sens_equil_dir_vec[loop], lh=lh_list[[as.character(strsplit(sens_equil_modcombos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(sens_equil_modcombos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=sens_equil_modcombos[loop,"C_opt"], param_adjust=c(sens_equil_modcombos[loop,"Param"],"SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(sens_vals[[as.character(strsplit(sens_equil_modcombos[loop,"LH"],"_")[[1]][2])]][sens_equil_modcombos[loop,"Adjust"],sens_equil_modcombos[loop,"Param"]],0.7,0.2,0.2,0.2), LFdist=1, newtonsteps=3)
	end_run <- Sys.time() - start_run

############################################################
#### find base directories
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

############################################################
#### sensitivities with variation
############################################################
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
	foreach(loop=1:length(sens_base_dir_vec), .packages=c('TMB','LIME')) %dopar% run_LIME(modpath=sens_base_dir_vec[loop], lh=lh_list[[as.character(strsplit(sens_base_modcombos[loop,"LH"],"_")[[1]][2])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(sens_base_modcombos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, f_true=FALSE, C_opt=sens_base_modcombos[loop,"C_opt"], param_adjust=c(sens_base_modcombos[loop,"Param"],"SigmaR","SigmaF","SigmaC","SigmaI"), val_adjust=c(sens_vals[[as.character(strsplit(sens_base_modcombos[loop,"LH"],"_")[[1]][2])]][sens_base_modcombos[loop,"Adjust"],sens_base_modcombos[loop,"Param"]],0.7,0.2,0.2,0.2), LFdist=1, newtonsteps=3)
	end_run <- Sys.time() - start_run


all_sens_dir <- c(sens_equil_dir_vec, sens_base_dir_vec)

sens_bp <- bias_precision(dirs=all_sens_dir, itervec=itervec)
saveRDS(sens_bp, file.path(res_dir, "sens_bias_precision.rds"))
sens_bp <- readRDS(file.path(res_dir, "sens_bias_precision.rds"))

dirs1 <- c(equil_dir_vec, base_dir_vec)
bp1 <- readRDS(file.path(res_dir, "equil_var_bias_precision.rds"))


col_vec <- c(rep("steelblue",3), rep("tomato",3))

png(file.path(fig_dir, "ParamSensitivities.png"), height=16, width=38, res=200, units="in")
par(mfrow=c(3,5), mar=c(0,0,0,0), omi=c(1.2,1.2,1,1), mgp=c(3,2,0))
## short
i1 <- which(grepl("Short", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("equil/", all_sens_dir) & grepl("/M/", all_sens_dir))
	all_sens_dir[i1]
i2 <- which(grepl("Short", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("equil/", dirs1)& grepl("/LC10/", dirs1))
	dirs1[i2]

i3 <- which(grepl("Short", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("Ramp/", all_sens_dir) & grepl("/M/", all_sens_dir))
	all_sens_dir[i3]
i4 <- which(grepl("Short", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("Ramp/", dirs1)& grepl("/LC10/", dirs1))
	dirs1[i4]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,6.5),ylim=c(-1,1))
abline(h=0, lty=2)
pdf <- cbind.data.frame(sens_bp$relerr[,i1[2]], bp1$relerr[,i2], sens_bp$relerr[,i1[1]], sens_bp$relerr[,i3[2]], bp1$relerr[,i4], sens_bp$relerr[,i3[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
axis(2,at=seq(-0.5,1,by=0.5), las=2,cex.axis=3.5)
mtext(side=3, "Natural mortality rate", line=1, cex=3)
print.letter("(a)", cex=3)

i1 <- which(grepl("Short", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("equil/", all_sens_dir) & grepl("/linf/", all_sens_dir))
	all_sens_dir[i1]
i2 <- which(grepl("Short", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("equil/", dirs1)& grepl("/LC10/", dirs1))
	dirs1[i2]

i3 <- which(grepl("Short", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("Ramp/", all_sens_dir) & grepl("/linf/", all_sens_dir))
	all_sens_dir[i3]
i4 <- which(grepl("Short", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("Ramp/", dirs1)& grepl("/LC10/", dirs1))
	dirs1[i4]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,6.5),ylim=c(-1,1))
abline(h=0, lty=2)
pdf <- cbind.data.frame(sens_bp$relerr[,i1[2]], bp1$relerr[,i2], sens_bp$relerr[,i1[1]], sens_bp$relerr[,i3[2]], bp1$relerr[,i4], sens_bp$relerr[,i3[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
mtext(side=3, "Asymptotic length", line=1, cex=3)
print.letter("(b)", cex=3)

i1 <- which(grepl("Short", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("equil/", all_sens_dir) & grepl("/vbk/", all_sens_dir))
	all_sens_dir[i1]
i2 <- which(grepl("Short", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("equil/", dirs1)& grepl("/LC10/", dirs1))
	dirs1[i2]

i3 <- which(grepl("Short", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("Ramp/", all_sens_dir) & grepl("/vbk/", all_sens_dir))
	all_sens_dir[i3]
i4 <- which(grepl("Short", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("Ramp/", dirs1)& grepl("/LC10/", dirs1))
	dirs1[i4]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,6.5),ylim=c(-1,1))
abline(h=0, lty=2)
pdf <- cbind.data.frame(sens_bp$relerr[,i1[2]], bp1$relerr[,i2], sens_bp$relerr[,i1[1]], sens_bp$relerr[,i3[2]], bp1$relerr[,i4], sens_bp$relerr[,i3[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
mtext(side=3, "von Bertalanffy k", line=1, cex=3)
print.letter("(c)", cex=3)


i1 <- which(grepl("Short", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("equil/", all_sens_dir) & grepl("/ML50/", all_sens_dir))
	all_sens_dir[i1]
i2 <- which(grepl("Short", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("equil/", dirs1)& grepl("/LC10/", dirs1))
	dirs1[i2]

i3 <- which(grepl("Short", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("Ramp/", all_sens_dir) & grepl("/ML50/", all_sens_dir))
	all_sens_dir[i3]
i4 <- which(grepl("Short", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("Ramp/", dirs1)& grepl("/LC10/", dirs1))
	dirs1[i4]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,6.5),ylim=c(-1,1))
abline(h=0, lty=2)
pdf <- cbind.data.frame(sens_bp$relerr[,i1[2]], bp1$relerr[,i2], sens_bp$relerr[,i1[1]], sens_bp$relerr[,i3[2]], bp1$relerr[,i4], sens_bp$relerr[,i3[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
mtext(side=3, "Length at 50% maturity", line=1, cex=3)
print.letter("(d)", cex=3)

i1 <- which(grepl("Short", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("equil/", all_sens_dir) & grepl("/CVlen/", all_sens_dir))
	all_sens_dir[i1]
i2 <- which(grepl("Short", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("equil/", dirs1)& grepl("/LC10/", dirs1))
	dirs1[i2]

i3 <- which(grepl("Short", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("Ramp/", all_sens_dir) & grepl("/CVlen/", all_sens_dir))
	all_sens_dir[i3]
i4 <- which(grepl("Short", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("Ramp/", dirs1)& grepl("/LC10/", dirs1))
	dirs1[i4]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,6.5),ylim=c(-1,1))
abline(h=0, lty=2)
pdf <- cbind.data.frame(sens_bp$relerr[,i1[2]], bp1$relerr[,i2], sens_bp$relerr[,i1[1]], sens_bp$relerr[,i3[2]], bp1$relerr[,i4], sens_bp$relerr[,i3[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
legend("topright", legend=c("Equilibrium", "Two-way F with variation"), col=c("steelblue", "tomato"), pch=15, cex=3.5)
mtext(side=3, "Age-length CV", line=1, cex=3)
mtext(side=4, "Short-lived", line=2, cex=3)
print.letter("(e)", cex=3)

## medium
i1 <- which(grepl("Medium", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("equil/", all_sens_dir) & grepl("/M/", all_sens_dir))
	all_sens_dir[i1]
i2 <- which(grepl("Medium", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("equil/", dirs1)& grepl("/LC10/", dirs1))
	dirs1[i2]

i3 <- which(grepl("Medium", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("Ramp/", all_sens_dir) & grepl("/M/", all_sens_dir))
	all_sens_dir[i3]
i4 <- which(grepl("Medium", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("Ramp/", dirs1)& grepl("/LC10/", dirs1))
	dirs1[i4]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,6.5),ylim=c(-2,2))
abline(h=0, lty=2)
pdf <- cbind.data.frame(sens_bp$relerr[,i1[2]], bp1$relerr[,i2], sens_bp$relerr[,i1[1]], sens_bp$relerr[,i3[2]], bp1$relerr[,i4], sens_bp$relerr[,i3[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
axis(2,at=seq(-1,2,by=1), las=2,cex.axis=3.5)
print.letter("(f)", cex=3)

i1 <- which(grepl("Medium", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("equil/", all_sens_dir) & grepl("/linf/", all_sens_dir))
	all_sens_dir[i1]
i2 <- which(grepl("Medium", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("equil/", dirs1)& grepl("/LC10/", dirs1))
	dirs1[i2]

i3 <- which(grepl("Medium", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("Ramp/", all_sens_dir) & grepl("/linf/", all_sens_dir))
	all_sens_dir[i3]
i4 <- which(grepl("Medium", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("Ramp/", dirs1)& grepl("/LC10/", dirs1))
	dirs1[i4]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,6.5),ylim=c(-2,2))
abline(h=0, lty=2)
pdf <- cbind.data.frame(sens_bp$relerr[,i1[2]], bp1$relerr[,i2], sens_bp$relerr[,i1[1]], sens_bp$relerr[,i3[2]], bp1$relerr[,i4], sens_bp$relerr[,i3[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
print.letter("(g)", cex=3)

i1 <- which(grepl("Medium", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("equil/", all_sens_dir) & grepl("/vbk/", all_sens_dir))
	all_sens_dir[i1]
i2 <- which(grepl("Medium", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("equil/", dirs1)& grepl("/LC10/", dirs1))
	dirs1[i2]

i3 <- which(grepl("Medium", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("Ramp/", all_sens_dir) & grepl("/vbk/", all_sens_dir))
	all_sens_dir[i3]
i4 <- which(grepl("Medium", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("Ramp/", dirs1)& grepl("/LC10/", dirs1))
	dirs1[i4]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,6.5),ylim=c(-2,2))
abline(h=0, lty=2)
pdf <- cbind.data.frame(sens_bp$relerr[,i1[2]], bp1$relerr[,i2], sens_bp$relerr[,i1[1]], sens_bp$relerr[,i3[2]], bp1$relerr[,i4], sens_bp$relerr[,i3[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
print.letter("(h)", cex=3)


i1 <- which(grepl("Medium", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("equil/", all_sens_dir) & grepl("/ML50/", all_sens_dir))
	all_sens_dir[i1]
i2 <- which(grepl("Medium", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("equil/", dirs1)& grepl("/LC10/", dirs1))
	dirs1[i2]

i3 <- which(grepl("Medium", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("Ramp/", all_sens_dir) & grepl("/ML50/", all_sens_dir))
	all_sens_dir[i3]
i4 <- which(grepl("Medium", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("Ramp/", dirs1)& grepl("/LC10/", dirs1))
	dirs1[i4]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,6.5),ylim=c(-2,2))
abline(h=0, lty=2)
pdf <- cbind.data.frame(sens_bp$relerr[,i1[2]], bp1$relerr[,i2], sens_bp$relerr[,i1[1]], sens_bp$relerr[,i3[2]], bp1$relerr[,i4], sens_bp$relerr[,i3[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
print.letter("(i)", cex=3)

i1 <- which(grepl("Medium", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("equil/", all_sens_dir) & grepl("/CVlen/", all_sens_dir))
	all_sens_dir[i1]
i2 <- which(grepl("Medium", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("equil/", dirs1)& grepl("/LC10/", dirs1))
	dirs1[i2]

i3 <- which(grepl("Medium", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("Ramp/", all_sens_dir) & grepl("/CVlen/", all_sens_dir))
	all_sens_dir[i3]
i4 <- which(grepl("Medium", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("Ramp/", dirs1)& grepl("/LC10/", dirs1))
	dirs1[i4]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,6.5),ylim=c(-2,2))
abline(h=0, lty=2)
pdf <- cbind.data.frame(sens_bp$relerr[,i1[2]], bp1$relerr[,i2], sens_bp$relerr[,i1[1]], sens_bp$relerr[,i3[2]], bp1$relerr[,i4], sens_bp$relerr[,i3[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
mtext(side=4, "Medium-lived", line=2, cex=3)
print.letter("(j)", cex=3)

## long
i1 <- which(grepl("Long", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("equil/", all_sens_dir) & grepl("/M/", all_sens_dir))
	all_sens_dir[i1]
i2 <- which(grepl("Long", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("equil/", dirs1)& grepl("/LC10/", dirs1))
	dirs1[i2]

i3 <- which(grepl("Long", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("Ramp/", all_sens_dir) & grepl("/M/", all_sens_dir))
	all_sens_dir[i3]
i4 <- which(grepl("Long", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("Ramp/", dirs1)& grepl("/LC10/", dirs1))
	dirs1[i4]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,6.5),ylim=c(-3,3))
abline(h=0, lty=2)
pdf <- cbind.data.frame(sens_bp$relerr[,i1[2]], bp1$relerr[,i2], sens_bp$relerr[,i1[1]], sens_bp$relerr[,i3[2]], bp1$relerr[,i4], sens_bp$relerr[,i3[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
axis(2,at=seq(-2,3,by=1), las=2,cex.axis=3.5)
axis(1, at=1:6, labels=c("-25%","true","+25%","-25%","true","+25%"), cex.axis=3.2, col=c(rep("steelblue",3), rep("tomato",3)))
print.letter("(k)", cex=3)

i1 <- which(grepl("Long", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("equil/", all_sens_dir) & grepl("/linf/", all_sens_dir))
	all_sens_dir[i1]
i2 <- which(grepl("Long", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("equil/", dirs1)& grepl("/LC10/", dirs1))
	dirs1[i2]

i3 <- which(grepl("Long", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("Ramp/", all_sens_dir) & grepl("/linf/", all_sens_dir))
	all_sens_dir[i3]
i4 <- which(grepl("Long", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("Ramp/", dirs1)& grepl("/LC10/", dirs1))
	dirs1[i4]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,6.5),ylim=c(-3,3))
abline(h=0, lty=2)
pdf <- cbind.data.frame(sens_bp$relerr[,i1[2]], bp1$relerr[,i2], sens_bp$relerr[,i1[1]], sens_bp$relerr[,i3[2]], bp1$relerr[,i4], sens_bp$relerr[,i3[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
print.letter("(l)", cex=3)

axis(1, at=1:6, labels=c("-25%","true","+25%","-25%","true","+25%"), cex.axis=3.2, col=c(rep("steelblue",3), rep("tomato",3)))

i1 <- which(grepl("Long", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("equil/", all_sens_dir) & grepl("/vbk/", all_sens_dir))
	all_sens_dir[i1]
i2 <- which(grepl("Long", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("equil/", dirs1)& grepl("/LC10/", dirs1))
	dirs1[i2]

i3 <- which(grepl("Long", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("Ramp/", all_sens_dir) & grepl("/vbk/", all_sens_dir))
	all_sens_dir[i3]
i4 <- which(grepl("Long", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("Ramp/", dirs1)& grepl("/LC10/", dirs1))
	dirs1[i4]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,6.5),ylim=c(-3,3))
abline(h=0, lty=2)
pdf <- cbind.data.frame(sens_bp$relerr[,i1[2]], bp1$relerr[,i2], sens_bp$relerr[,i1[1]], sens_bp$relerr[,i3[2]], bp1$relerr[,i4], sens_bp$relerr[,i3[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
print.letter("(m)", cex=3)

axis(1, at=1:6, labels=c("-25%","true","+25%","-25%","true","+25%"), cex.axis=3.2, col=c(rep("steelblue",3), rep("tomato",3)))


i1 <- which(grepl("Long", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("equil/", all_sens_dir) & grepl("/ML50/", all_sens_dir))
	all_sens_dir[i1]
i2 <- which(grepl("Long", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("equil/", dirs1)& grepl("/LC10/", dirs1))
	dirs1[i2]

i3 <- which(grepl("Long", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("Ramp/", all_sens_dir) & grepl("/ML50/", all_sens_dir))
	all_sens_dir[i3]
i4 <- which(grepl("Long", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("Ramp/", dirs1)& grepl("/LC10/", dirs1))
	dirs1[i4]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,6.5),ylim=c(-3,3))
abline(h=0, lty=2)
pdf <- cbind.data.frame(sens_bp$relerr[,i1[2]], bp1$relerr[,i2], sens_bp$relerr[,i1[1]], sens_bp$relerr[,i3[2]], bp1$relerr[,i4], sens_bp$relerr[,i3[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
print.letter("(n)", cex=3)

axis(1, at=1:6, labels=c("-25%","true","+25%","-25%","true","+25%"), cex.axis=3.2, col=c(rep("steelblue",3), rep("tomato",3)))

i1 <- which(grepl("Long", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("equil/", all_sens_dir) & grepl("/CVlen/", all_sens_dir))
	all_sens_dir[i1]
i2 <- which(grepl("Long", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("equil/", dirs1)& grepl("/LC10/", dirs1))
	dirs1[i2]

i3 <- which(grepl("Long", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("Ramp/", all_sens_dir) & grepl("/CVlen/", all_sens_dir))
	all_sens_dir[i3]
i4 <- which(grepl("Long", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("Ramp/", dirs1)& grepl("/LC10/", dirs1))
	dirs1[i4]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,6.5),ylim=c(-3,3))
abline(h=0, lty=2)
pdf <- cbind.data.frame(sens_bp$relerr[,i1[2]], bp1$relerr[,i2], sens_bp$relerr[,i1[1]], sens_bp$relerr[,i3[2]], bp1$relerr[,i4], sens_bp$relerr[,i3[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
print.letter("(o)", cex=3)

axis(1, at=1:6, labels=c("-25%","true","+25%","-25%","true","+25%"), cex.axis=3.2, col=c(rep("steelblue",3), rep("tomato",3)))

mtext(side=4, "Longer-lived", line=2, cex=3)
mtext(side=2, outer=TRUE, line=6, "Relative error", cex=3)
mtext(side=1, outer=TRUE, line=6, "Fixed value", cex=3)
dev.off()











png(file.path(fig_dir, "ParamSensitivities_presentation.png"), height=16, width=38, res=200, units="in")
par(mfrow=c(1,5), mar=c(0,0,0,0), omi=c(1.2,1.2,1,1), mgp=c(3,2,0))
## medium
i1 <- which(grepl("Medium", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("equil/", all_sens_dir) & grepl("/M/", all_sens_dir))
	all_sens_dir[i1]
i2 <- which(grepl("Medium", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("equil/", dirs1)& grepl("/LC10/", dirs1))
	dirs1[i2]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,3.5),ylim=c(-1.2,2.5))
abline(h=0, lty=2)
pdf <- cbind.data.frame(sens_bp$relerr[,i1[2]], bp1$relerr[,i2], sens_bp$relerr[,i1[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
axis(2,at=seq(-1,2,by=1), las=2,cex.axis=3.5)
axis(1, at=1:3, labels=c("-25%","true","+25%"), cex.axis=3.2, col=c(rep("steelblue",3)))
mtext(side=3, "Natural mortality", line=1, cex=3)

# print.letter("(f)", cex=3)

i1 <- which(grepl("Medium", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("equil/", all_sens_dir) & grepl("/linf/", all_sens_dir))
	all_sens_dir[i1]
i2 <- which(grepl("Medium", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("equil/", dirs1)& grepl("/LC10/", dirs1))
	dirs1[i2]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,3.5),ylim=c(-1.2,2.5))
abline(h=0, lty=2)
pdf <- cbind.data.frame(sens_bp$relerr[,i1[2]], bp1$relerr[,i2], sens_bp$relerr[,i1[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
axis(1, at=1:3, labels=c("-25%","true","+25%"), cex.axis=3.2, col=c(rep("steelblue",3)))
mtext(side=3, "Asymptotic length", line=1, cex=3)
# print.letter("(g)", cex=3)

i1 <- which(grepl("Medium", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("equil/", all_sens_dir) & grepl("/vbk/", all_sens_dir))
	all_sens_dir[i1]
i2 <- which(grepl("Medium", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("equil/", dirs1)& grepl("/LC10/", dirs1))
	dirs1[i2]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,3.5),ylim=c(-1.2,2.5))
abline(h=0, lty=2)
pdf <- cbind.data.frame(sens_bp$relerr[,i1[2]], bp1$relerr[,i2], sens_bp$relerr[,i1[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
axis(1, at=1:3, labels=c("-25%","true","+25%"), cex.axis=3.2, col=c(rep("steelblue",3)))
mtext(side=3, "von Bertalanffy k", line=1, cex=3)
# print.letter("(h)", cex=3)


i1 <- which(grepl("Medium", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("equil/", all_sens_dir) & grepl("/ML50/", all_sens_dir))
	all_sens_dir[i1]
i2 <- which(grepl("Medium", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("equil/", dirs1)& grepl("/LC10/", dirs1))
	dirs1[i2]


plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,3.5),ylim=c(-1.2,2.5))
abline(h=0, lty=2)
pdf <- cbind.data.frame(sens_bp$relerr[,i1[2]], bp1$relerr[,i2], sens_bp$relerr[,i1[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
axis(1, at=1:3, labels=c("-25%","true","+25%"), cex.axis=3.2, col=c(rep("steelblue",3)))
mtext(side=3, "Length at 50% maturity", line=1, cex=3)

# print.letter("(i)", cex=3)

i1 <- which(grepl("Medium", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("equil/", all_sens_dir) & grepl("/CVlen/", all_sens_dir))
	all_sens_dir[i1]
i2 <- which(grepl("Medium", dirs1) & grepl("SampleSize_200/", dirs1) & grepl("equil/", dirs1)& grepl("/LC10/", dirs1))
	dirs1[i2]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,3.5),ylim=c(-1.2,2.5))
abline(h=0, lty=2)
pdf <- cbind.data.frame(sens_bp$relerr[,i1[2]], bp1$relerr[,i2], sens_bp$relerr[,i1[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
axis(1, at=1:3, labels=c("-25%","true","+25%"), cex.axis=3.2, col=c(rep("steelblue",3)))
mtext(side=4, "Medium-lived", line=2, cex=3)
mtext(side=3, "Age-length CV", line=1, cex=3)
# print.letter("(j)", cex=3)
dev.off()
