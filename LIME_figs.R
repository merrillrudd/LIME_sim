rm(list=ls())
library(RColorBrewer)
library(LIME)
library(beanplot)

### ----- directories and functions -----------###
# main_dir <- "C:\\Git_Projects\\LIME_sim"
main_dir <- "F:\\Merrill\\Git_Projects\\LIME_sim"


funs <- list.files(file.path(main_dir, "R"))
ignore <- sapply(1:length(funs), function(x) source(file.path(main_dir,"R", funs[x])))

main_dir <- "F:\\Merrill\\Git_Projects\\LIME_sim\\v6"

fig_dir <- file.path(main_dir, "figs")
dir.create(fig_dir, showWarnings=FALSE)

res_dir <- file.path(main_dir, "calc_results")
dir.create(res_dir, showWarnings=FALSE)
###########################################################
#### Figure 1. Confounding F and R
##########################################################

lh_med0 <- create_lh_list(vbk=0.21, linf=64.58, lwa=0.0245, lwb=2.790, S50=4, M50=34, selex_input="age", maturity_input="length", M=0.43, binwidth=1, t0=-0.01, CVlen=0.1, SigmaC=0.2, SigmaI=0.2, SigmaR=0.01, SigmaF=0.01, R0=1,  h=1, qcoef=1e-5, F1=0.25, Fequil=0.25, Frate=0.17, Fmax=0.7, start_ages=0, rho=0, nseasons=12)
lh_med1 <- create_lh_list(vbk=0.21, linf=64.58, lwa=0.0245, lwb=2.790, S50=4, M50=34, selex_input="age", maturity_input="length", M=0.43, binwidth=1, t0=-0.01, CVlen=0.1, SigmaC=0.2, SigmaI=0.2, SigmaR=0.7, SigmaF=0.01, R0=1,  h=1, qcoef=1e-5, F1=0.25, Fequil=0.25, Frate=0.17, Fmax=0.7, start_ages=0, rho=0.426, nseasons=12)

finc <- sim_pop(lh=lh_med0, Nyears=20, Fdynamics="Increasing", Rdynamics="Constant", Nyears_comp=20, comp_sample=100000, init_depl=0.5, nburn=1, seed=123, modname="ex")
rpulse <- sim_pop(lh=lh_med1, Nyears=20, Fdynamics="Constant", Rdynamics="AR", Nyears_comp=20, comp_sample=100000, init_depl=0.5, nburn=1, seed=123, modname="ex")
both <- sim_pop(lh=lh_med1, Nyears=20, Fdynamics="Increasing", Rdynamics="AR", Nyears_comp=20, comp_sample=100000, init_depl=0.5, nburn=1, seed=123, modname="ex")

png(file.path(fig_dir, "ConfoundingFR.png"), width=25, height=10, res=200, units="in")
# dev.new()
cols <- brewer.pal(5, "Set1")
par(mfrow=c(3,5), mar=c(0,5,0,0), omi=c(1,1,1,1))
plot(finc$F_t, col=cols[1], type="o", lwd=4, xlab="", ylab="", xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylim=c(0, 1))
axis(2, at=seq(0.2,1.2,by=0.2), las=2, cex.axis=2)
mtext(side=3, "Fishing mortality", cex=2, line=1)
print.letter("(a)", cex=3)
plot(finc$R_t, col=cols[2], type="o", lwd=4, xlab="", ylab="", xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylim=c(0,2))
axis(2, at=seq(0.5, 3.5, by=0.5), las=2, cex.axis=2)
mtext(side=3, "Recruitment", cex=2, line=1)
print.letter("(b)", cex=3)
plot(finc$ML_t, col=cols[3], type="o", lwd=4, xlab="", ylab="", xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylim=c(0,max(finc$ML_t)*1.5))
axis(2, at=seq(10,50,by=10), las=2, cex.axis=2)
mtext(side=3, "Mean length", cex=2, line=1)
print.letter("(c)", cex=3)
barplot(finc$LF0[nrow(finc$LF0),]/sum(finc$LF0[nrow(finc$LF0),]), col="gray", xlab="", ylab="", xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylim=c(0, max(finc$LF0[nrow(finc$LF0),]/sum(finc$LF0[nrow(finc$LF0),]))*2))
box()
par(new=TRUE)
barplot(finc$LF[nrow(finc$LF),]/sum(finc$LF[nrow(finc$LF),]), col=cols[4], xlab="", ylab="", xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylim=c(0, max(finc$LF0[nrow(finc$LF0),]/sum(finc$LF0[nrow(finc$LF0),]))*2))
axis(2, at=seq(0.01,0.06,by=0.01), las=2, cex.axis=2)
mtext(side=3, "Length composition", cex=2, line=1)
print.letter("(d)", cex=3)
plot(finc$D_t, col=cols[5], type="o", lwd=4, xlab="", ylab="", xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylim=c(0,max(rpulse$D_t)*1.5))
axis(2, at=seq(0.5,2.5,by=0.5), las=2, cex.axis=2)
mtext(side=3, "Relative biomass", cex=2, line=1)
print.letter("(e)", cex=3)

plot(rpulse$F_t, col=cols[1], type="o", lwd=4, xlab="", ylab="", xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylim=c(0, 1))
axis(2, at=seq(0.2,1.2,by=0.2), las=2, cex.axis=2)
print.letter("(f)", cex=3)
plot(rpulse$R_t, col=cols[2], type="o", lwd=4, xlab="", ylab="", xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylim=c(0,2))
axis(2, at=seq(0.5, 3.5, by=0.5), las=2, cex.axis=2)
print.letter("(g)", cex=3)
plot(rpulse$ML_t, col=cols[3], type="o", lwd=4, xlab="", ylab="", xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylim=c(0,max(finc$ML_t)*1.5))
axis(2, at=seq(10,50,by=10), las=2, cex.axis=2)
print.letter("(h)", cex=3)
barplot(rpulse$LF0[nrow(rpulse$LF0),]/sum(rpulse$LF0[nrow(rpulse$LF0),]), col="gray", xlab="", ylab="", xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylim=c(0, max(finc$LF0[nrow(finc$LF0),]/sum(finc$LF0[nrow(finc$LF0),]))*2))
box()
par(new=TRUE)
barplot(rpulse$LF[nrow(rpulse$LF),]/sum(rpulse$LF[nrow(rpulse$LF),]), col=cols[4], xlab="", ylab="", xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylim=c(0, max(finc$LF0[nrow(finc$LF0),]/sum(finc$LF0[nrow(finc$LF0),]))*2))
axis(2, at=seq(0.01,0.06,by=0.01), las=2, cex.axis=2)
print.letter("(i)", cex=3)
plot(rpulse$D_t, col=cols[5], type="o", lwd=4, xlab="", ylab="", xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylim=c(0,max(rpulse$D_t)*1.5))
axis(2, at=seq(0.5,2.5,by=0.5), las=2, cex.axis=2)
print.letter("(j)", cex=3)

plot(both$F_t, col=cols[1], type="o", lwd=4, xlab="", ylab="", xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylim=c(0, 1))
axis(2, at=seq(0,1.2,by=0.2), las=2, cex.axis=2)
axis(1, cex.axis=2)
mtext(side=1, "Year", line=3.5, cex=2)
print.letter("(k)", cex=3)
plot(both$R_t, col=cols[2], type="o", lwd=4, xlab="", ylab="", xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylim=c(0,2))
axis(2, at=seq(0, 3.5, by=0.5), las=2, cex.axis=2)
axis(1, cex.axis=2)
mtext(side=1, "Year", line=3.5, cex=2)
print.letter("(l)", cex=3)
plot(both$ML_t, col=cols[3], type="o", lwd=4, xlab="", ylab="", xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylim=c(0,max(finc$ML_t)*1.5))
axis(2, at=seq(0,50,by=10), las=2, cex.axis=2)
axis(1, cex.axis=2)
mtext(side=1, "Year", line=3.5, cex=2)
print.letter("(m)", cex=3)
barplot(both$LF0[nrow(both$LF0),]/sum(both$LF0[nrow(both$LF0),]), col="gray", xlab="", ylab="", xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylim=c(0, max(finc$LF0[nrow(finc$LF0),]/sum(finc$LF0[nrow(finc$LF0),]))*2))
box()
par(new=TRUE)
barplot(both$LF[nrow(both$LF),]/sum(both$LF[nrow(both$LF),]), col=cols[4], xlab="", ylab="", xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylim=c(0, max(finc$LF0[nrow(finc$LF0),]/sum(finc$LF0[nrow(finc$LF0),]))*2))
axis(2, at=seq(0,0.06,by=0.01), las=2, cex.axis=2)
axis(1, cex.axis=2)
mtext(side=1, "Length bin", line=3.5, cex=2)
print.letter("(n)", cex=3)
plot(both$D_t, col=cols[5], type="o", lwd=4, xlab="", ylab="", xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylim=c(0,max(rpulse$D_t)*1.5))
axis(2, at=seq(0,2.5,by=0.5), las=2, cex.axis=2)
axis(1, cex.axis=2)
mtext(side=1, "Year", line=3.5, cex=2)
print.letter("(o)", cex=3)
dev.off()

###########################################################
#### Figure 2. life history
############################################################
	## equilibrium
	lh_equil <- adj_variation(SigmaR=0.01, SigmaF=0.01, SigmaC=0.01, SigmaI=0.01, CVlen=0.1, rho=0, selex_param=1, mat_param=1, nseasons=1)

	## base
	lh_base <- adj_variation(SigmaR=0.737, SigmaF=0.2, SigmaC=0.2, SigmaI=0.2, CVlen=0.1, rho=0.426, selex_param=1, mat_param=1, nseasons=1)

	png(file.path(fig_dir, "life_history_comparison.png"), height=6, width=10, res=200, units="in")
	lh_fig(lh=lh_base)
	dev.off()


###########################################################
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
	## setup equilibrium dirs
	equil_dir_vec <- model_paths(res_dir=equil_dir, modcombos=equil_modcombos[,-c(ncol(equil_modcombos))])

	equil_12seasons_dir <- file.path(main_dir, "equil_12seasons")
	# unlink(equil_12seasons_dir, TRUE)
	dir.create(equil_12seasons_dir, showWarnings=FALSE)	

	## setup equilibrium dirs
	equil_12seasons_dir_vec <- model_paths(res_dir=equil_12seasons_dir, modcombos=equil_modcombos[,-c(ncol(equil_modcombos))])

lh_vec <- c("Short", "Medium", "Long")
data_vec <- c("LC10", "LC1", "LBSPR10", "LBSPR1") #"Index_LC10", "Index_LC1", "Catch_LC10", "Catch_LC1", "LC10", "LC1", "LBSPR10", "LBSPR1")
itervec <- 1:100

LBSPR_modcombos <- expand.grid("Data_avail"=data_vec, "LH"=paste0("LH_",lh_vec), stringsAsFactors=FALSE)

	equil_LBSPR_dir <- file.path(main_dir, "equil_LBSPR")
	# unlink(equil_LBSPR_dir, TRUE)
	dir.create(equil_LBSPR_dir, showWarnings=FALSE)	

	## setup equilibrium dirs
	equil_LBSPR_dir_vec <- model_paths(res_dir=equil_LBSPR_dir, modcombos=LBSPR_modcombos)





############################################################
#### base - with variation
############################################################
### ----- models to run ----------- ###
############################################################
#### base - with variation
############################################################
### ----- models to run ----------- ###
lh_vec <- c("Short", "Medium", "Long")
Fdyn_vec <- c("Ramp", "Increasing")#, "Decreasing")
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







############################################################
#### SUMMARY
############################################################

all_dirs <- c(equil_dir_vec, base_dir_vec, equil_12seasons_dir_vec, equil_LBSPR_dir_vec)


icol <- rev(brewer.pal(3,"Purples"))
ccol <- rev(brewer.pal(3, "Oranges"))
lcol <- rev(brewer.pal(3, "Blues"))
# rcol <- rev(brewer.pal(3, "Greens"))
col_vec <- c(gray(0.3), icol[1:2], ccol[1:2], lcol[1:2])#, rcol[1:2])
## example of model fits
	FUN <- function(InputMat, log=TRUE, rel=FALSE){
        index <- which(is.na(InputMat[,2])==FALSE)
        if(log==TRUE) return(c( exp(InputMat[index,1]-1.96*InputMat[index,2]), rev(exp(InputMat[index,1]+1.96*InputMat[index,2]))))
        if(log==FALSE) return(c( InputMat[index,1]-1.96*InputMat[index,2], rev(InputMat[index,1]+1.96*InputMat[index,2])))
	}
choose <- all_dirs[which(grepl("Medium", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("base/", all_dirs) & grepl("LBSPR", all_dirs)==FALSE & grepl("LC5/", all_dirs)==FALSE & grepl("LC2/", all_dirs)==FALSE & grepl("Ramp", all_dirs)) ]
names <- c("Rich", "Index+LC10","Index+LC1","Catch+LC10","Catch+LC1","LC10","LC1")
# set.seed(143)
iter <- 1
png(file.path(fig_dir, "Iter_example.png"), width=25, height=16, res=200, units="in")
par(mfrow=c(3,7), mar=c(0,0,0,0), omi=c(1,1,2,1))
for(i in 1:length(choose)){
	true <- readRDS(file.path(choose[i],iter,"True.rds"))
	rep <- readRDS(file.path(choose[i],iter,"Report.rds"))
	sdrep <- readRDS(file.path(choose[i],iter,"Sdreport.rds"))

	plot(x=1,y=1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",xaxs="i",yaxs="i",xlim=c(1,20),ylim=c(0,2.5))
    if(all(is.na(sdrep))==FALSE) polygon( y=FUN(summary(sdrep)[which(rownames(summary(sdrep))=="lF_t"),], log=TRUE), x=c(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="lF_t"),2])==FALSE), rev(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="lF_t"),2])==FALSE))), col=paste0(col_vec[i],50), border=NA)
    lines(rep$F_t, col=col_vec[i], lwd=2)
	lines(true$F_t, col="black", lwd=3, lty=2)
	if(i==1){
		axis(2, seq(0.5,2.5,by=0.5), las=2, cex.axis=2)
		mtext(side=2, "Fishing mortality", cex=3, line=4)
	}
	mtext(side=3, names[i], cex=3, line=2)
	print.letter(paste0("(", letters[i], ")"), cex=3)
}
for(i in 1:length(choose)){
	true <- readRDS(file.path(choose[i],iter,"True.rds"))
	rep <- readRDS(file.path(choose[i],iter,"Report.rds"))
	sdrep <- readRDS(file.path(choose[i],iter,"Sdreport.rds"))

	plot(x=1,y=1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",xaxs="i",yaxs="i",xlim=c(1,20),ylim=c(0,5))
    if(all(is.na(sdrep))==FALSE) polygon( y=FUN(summary(sdrep)[which(rownames(summary(sdrep))=="lR_t"),], log=TRUE), x=c(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="lR_t"),2])==FALSE), rev(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="lR_t"),2])==FALSE))), col=paste0(col_vec[i],50), border=NA)
    lines(rep$R_t, col=col_vec[i], lwd=2)
	lines(true$R_t, col="black", lwd=3, lty=2)
	if(i==1){
		axis(2, seq(1,4,by=1), las=2, cex.axis=2)
		mtext(side=2, "Recruitment", cex=3, line=4)
	}
	print.letter(paste0("(", letters[i+length(choose)], ")"), cex=3)
}
for(i in 1:length(choose)){
	true <- readRDS(file.path(choose[i],iter,"True.rds"))
	rep <- readRDS(file.path(choose[i],iter,"Report.rds"))
	sdrep <- readRDS(file.path(choose[i],iter,"Sdreport.rds"))

	plot(x=1,y=1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",xaxs="i",yaxs="i",xlim=c(1,20),ylim=c(0,1))
    if(all(is.na(sdrep))==FALSE) polygon( y=FUN(summary(sdrep)[which(rownames(summary(sdrep))=="SPR_t"),], log=FALSE), x=c(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="SPR_t"),2])==FALSE), rev(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="SPR_t"),2])==FALSE))), col=paste0(col_vec[i],50), border=NA)
    lines(rep$SPR_t, col=col_vec[i], lwd=2)
	lines(true$SPR_t, col="black", lwd=3, lty=2)
	if(i==1){
		axis(2, seq(0.2,1,by=0.2), las=2, cex.axis=2)
		mtext(side=2, "SPR", cex=3, line=4)
	}
	axis(1,at=pretty(c(1,20)), cex.axis=2)
	print.letter(paste0("(", letters[i+length(choose)*2], ")"), cex=3)
}
mtext("Year", outer=TRUE, side=1, cex=3, line=4)
mtext("Data availability scenario", outer=TRUE, side=3, line=7, cex=3)
dev.off()




# all_bp_check <- bias_precision(dirs=all_dirs, itervec=itervec)
# saveRDS(all_bp_check, file.path(res_dir, "all_bias_precision.rds"))
all_bp <- readRDS(file.path(res_dir, "all_bias_precision.rds"))
all_bp$esterr[which(is.infinite(all_bp$esterr))] <- NA

## interval coverage
# all_ic <- interval_coverage(all_dirs, itervec=itervec)
# saveRDS(all_ic, file.path(res_dir, "all_interval_coverage.rds"))
all_ic_raw <- readRDS(file.path(res_dir, "all_interval_coverage.rds"))

## check that I'm not including non-converged runs in the calculation of coverage
all_ic <- all_ic_raw
all_ic$cover[which(all_ic$converge==0)] <- 0


icol <- rev(brewer.pal(3,"Purples"))
ccol <- rev(brewer.pal(3, "Oranges"))
lcol <- rev(brewer.pal(3, "Blues"))
rcol <- rev(brewer.pal(3, "Greens"))
col_vec <- c(gray(0.3), icol[1:2], ccol[1:2], lcol[1:2], rcol[1:2])


png(file.path(fig_dir, "RelativeError_n200.png"), height=16, width=42, res=200, units="in")
par(mfrow=c(3,3), mar=c(0,4.5,0,0), omi=c(1.2,1.2,1,1), mgp=c(3,2,0))
## equilibrium
i2 <- which(grepl("Short", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("/equil/", all_dirs))
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1,1.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(all_bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(all_bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=2.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=2.5)
abline(h=0, lty=2)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
axis(2,at=seq(-0.5,1.5,by=0.5), las=2, cex.axis=3.5)
mtext(side=3, "Short-lived", line=1, cex=3)
print.letter("(a)", cex=3)

i2 <- which(grepl("Medium", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("/equil/", all_dirs))
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1.5,2.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(all_bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(all_bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
mtext(side=3, "Medium-lived", line=1, cex=3)
axis(2,at=seq(-1,2.5,by=1), las=2, cex.axis=3.5)
print.letter("(b)", cex=3)

i2 <- which(grepl("Long", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("/equil/", all_dirs))
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-2,3.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(all_bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(all_bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
mtext(side=3, "Longer-lived", line=1, cex=3)
mtext(side=4, "Equilibrium", line=2, cex=3)
axis(2,at=seq(-1,3,by=1), las=2, cex.axis=3.5)
print.letter("(c)", cex=3)

## base
i2 <- which(grepl("Short", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("base/", all_dirs) & grepl("Ramp", all_dirs))
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1,1.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(all_bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(all_bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
axis(2,at=seq(-0.5,1.5,by=0.5), las=2, cex.axis=3.5)
# axis(1, at=1:9, labels=c("Rich", "Index\nLC10", "Index\nLC1", "Catch\nLC10", "Catch\nLC1", "LC10", "LC1", "LBSPR\n10", "LBSPR\n1"), cex.axis=2.8)
print.letter("(d)", cex=3)

i2 <- which(grepl("Medium", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("base/", all_dirs) & grepl("Ramp", all_dirs))
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1.5,2.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(all_bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(all_bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
# axis(1, at=1:9, labels=c("Rich", "Index\nLC10", "Index\nLC1", "Catch\nLC10", "Catch\nLC1", "LC10", "LC1", "LBSPR\n10", "LBSPR\n1"), cex.axis=2.8)
axis(2,at=seq(-1,2.5,by=1), las=2, cex.axis=3.5)
print.letter("(e)", cex=3)

i2 <- which(grepl("Long", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("base/", all_dirs) & grepl("Ramp", all_dirs))
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-2,3.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(all_bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(all_bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
# axis(1, at=1:9, labels=c("Rich", "Index\nLC10", "Index\nLC1", "Catch\nLC10", "Catch\nLC1", "LC10", "LC1", "LBSPR\n10", "LBSPR\n1"), cex.axis=2.8)
axis(2,at=seq(-1,3,by=1), las=2, cex.axis=3.5)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
mtext(side=4, "Two-way trip", line=2, cex=3)
print.letter("(f)", cex=3)


## base
i2 <- which(grepl("Short", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("base/", all_dirs) & grepl("Increasing", all_dirs))
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1,1.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(all_bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(all_bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
axis(2,at=seq(-0.5,1.5,by=0.5), las=2, cex.axis=3.5)
axis(1, at=1:length(i2), labels=c("Rich", "Index\nLC10", "Index\nLC1", "Catch\nLC10", "Catch\nLC1", "LC10", "LC1", "LBSPR\n10", "LBSPR\n1"), cex.axis=3)
print.letter("(g)", cex=3)

i2 <- which(grepl("Medium", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("base/", all_dirs) & grepl("Increasing", all_dirs))
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1.5,2.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(all_bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(all_bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
axis(1, at=1:length(i2), labels=c("Rich", "Index\nLC10", "Index\nLC1", "Catch\nLC10", "Catch\nLC1", "LC10", "LC1", "LBSPR\n10", "LBSPR\n1"), cex.axis=3)
axis(2,at=seq(-1,2.5,by=1), las=2, cex.axis=3.5)
print.letter("(h)", cex=3)

i2 <- which(grepl("Long", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("base/", all_dirs) & grepl("Increasing", all_dirs))
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-2,3.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(all_bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(all_bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
axis(1, at=1:length(i2), labels=c("Rich", "Index\nLC10", "Index\nLC1", "Catch\nLC10", "Catch\nLC1", "LC10", "LC1", "LBSPR\n10", "LBSPR\n1"), cex.axis=3)
axis(2,at=seq(-1,3,by=1), las=2, cex.axis=3.5)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
mtext(side=4, "One-way trip", line=2, cex=3)
mtext(side=2, outer=TRUE, line=3, "Relative error", cex=3)
mtext(side=1, outer=TRUE, line=6, "Data availability scenario", cex=3)
print.letter("(i)", cex=3)
dev.off()



png(file.path(fig_dir, "RelativeError_n1000.png"), height=16, width=35, res=200, units="in")
par(mfrow=c(3,3), mar=c(0,4.5,0,0), omi=c(1.2,1.2,1,1), mgp=c(3,2,0))
## equilibrium
i2 <- which(grepl("Short", all_dirs) & grepl("SampleSize_1000/", all_dirs) & grepl("/equil/", all_dirs) & grepl("LBSPR",all_dirs)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1,1.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(all_bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(all_bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=2.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=2.5)
abline(h=0, lty=2)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
axis(2,at=seq(-0.5,1.5,by=0.5), las=2,cex.axis=2.5)
mtext(side=3, "Short-lived", line=1, cex=3)
print.letter("(a)", cex=3)

i2 <- which(grepl("Medium", all_dirs) & grepl("SampleSize_1000/", all_dirs) & grepl("/equil/", all_dirs) & grepl("LBSPR",all_dirs)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1.5,2.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(all_bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(all_bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
mtext(side=3, "Medium-lived", line=1, cex=3)
axis(2,at=seq(-1,2.5,by=1), las=2,cex.axis=2.5)
print.letter("(b)", cex=3)

i2 <- which(grepl("Long", all_dirs) & grepl("SampleSize_1000/", all_dirs) & grepl("/equil/", all_dirs) & grepl("LBSPR",all_dirs)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-2,3.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(all_bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(all_bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
mtext(side=3, "Longer-lived", line=1, cex=3)
mtext(side=4, "Equilibrium", line=2, cex=3)
axis(2,at=seq(-1,3,by=1), las=2,cex.axis=2.5)
print.letter("(c)", cex=3)

## base
i2 <- which(grepl("Short", all_dirs) & grepl("SampleSize_1000/", all_dirs) & grepl("base/", all_dirs) & grepl("Ramp", all_dirs) & grepl("LBSPR",all_dirs)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1,1.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(all_bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(all_bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
axis(2,at=seq(-0.5,1.5,by=0.5), las=2,cex.axis=2.5)
# axis(1, at=1:9, labels=c("Rich", "Index\nLC10", "Index\nLC1", "Catch\nLC10", "Catch\nLC1", "LC10", "LC1", "LBSPR\n10", "LBSPR\n1"), cex.axis=2.8)
print.letter("(d)", cex=3)

i2 <- which(grepl("Medium", all_dirs) & grepl("SampleSize_1000/", all_dirs) & grepl("base/", all_dirs) & grepl("Ramp", all_dirs) & grepl("LBSPR",all_dirs)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1.5,2.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(all_bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(all_bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
# axis(1, at=1:9, labels=c("Rich", "Index\nLC10", "Index\nLC1", "Catch\nLC10", "Catch\nLC1", "LC10", "LC1", "LBSPR\n10", "LBSPR\n1"), cex.axis=2.8)
axis(2,at=seq(-1,2.5,by=1), las=2,cex.axis=2.5)
print.letter("(e)", cex=3)

i2 <- which(grepl("Long", all_dirs) & grepl("SampleSize_1000/", all_dirs) & grepl("base/", all_dirs) & grepl("Ramp", all_dirs) & grepl("LBSPR",all_dirs)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-2,3.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(all_bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(all_bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
# axis(1, at=1:9, labels=c("Rich", "Index\nLC10", "Index\nLC1", "Catch\nLC10", "Catch\nLC1", "LC10", "LC1", "LBSPR\n10", "LBSPR\n1"), cex.axis=2.8)
axis(2,at=seq(-1,3,by=1), las=2,cex.axis=2.5)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
mtext(side=4, "Two-way trip", line=2, cex=3)
print.letter("(f)", cex=3)


## base
i2 <- which(grepl("Short", all_dirs) & grepl("SampleSize_1000/", all_dirs) & grepl("base/", all_dirs) & grepl("Increasing", all_dirs) & grepl("LBSPR",all_dirs)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1,1.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(all_bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(all_bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
axis(2,at=seq(-0.5,1.5,by=0.5), las=2,cex.axis=2.5)
axis(1, at=1:length(i2), labels=c("Rich", "Index\nLC10", "Index\nLC1", "Catch\nLC10", "Catch\nLC1", "LC10", "LC1"), cex.axis=2.8)
print.letter("(g)", cex=3)

i2 <- which(grepl("Medium", all_dirs) & grepl("SampleSize_1000/", all_dirs) & grepl("base/", all_dirs) & grepl("Increasing", all_dirs) & grepl("LBSPR",all_dirs)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1.5,2.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(all_bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(all_bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
axis(1, at=1:length(i2), labels=c("Rich", "Index\nLC10", "Index\nLC1", "Catch\nLC10", "Catch\nLC1", "LC10", "LC1"), cex.axis=2.8)
axis(2,at=seq(-1,2.5,by=1), las=2,cex.axis=2.5)
print.letter("(h)", cex=3)

i2 <- which(grepl("Long", all_dirs) & grepl("SampleSize_1000/", all_dirs) & grepl("base/", all_dirs) & grepl("Increasing", all_dirs) & grepl("LBSPR",all_dirs)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-2,3.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(all_bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(all_bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
axis(1, at=1:length(i2), labels=c("Rich", "Index\nLC10", "Index\nLC1", "Catch\nLC10", "Catch\nLC1", "LC10", "LC1"), cex.axis=2.8)
axis(2,at=seq(-1,3,by=1), las=2,cex.axis=2.5)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
mtext(side=4, "One-way trip", line=2, cex=3)
mtext(side=2, outer=TRUE, line=3, "Relative error", cex=3)
mtext(side=1, outer=TRUE, line=6, "Data availability scenario", cex=3)
print.letter("(i)", cex=3)
dev.off()



lcol <- rev(brewer.pal(5, "Blues"))
col_vec <- lcol

### compare SampleSize
png(file.path(fig_dir, "Compare_Sample_Size.png"), height=16, width=35, res=200, units="in")
par(mfrow=c(3,3), mar=c(0,4.5,0,0), omi=c(1.2,1.2,1,1), mgp=c(3,2,0))
## equilibrium
i2 <- which(grepl("Short", all_dirs) & grepl("/LC1/", all_dirs) & grepl("/equil/", all_dirs) & grepl("LBSPR",all_dirs)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1,1.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(all_bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(all_bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=2.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=2.5)
abline(h=0, lty=2)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
axis(2,at=seq(-0.5,1.5,by=0.5), las=2,cex.axis=2.5)
mtext(side=3, "Short-lived", line=1, cex=3)
print.letter("(a)", cex=3)

i2 <- which(grepl("Medium", all_dirs) & grepl("/LC1/", all_dirs) & grepl("/equil/", all_dirs) & grepl("LBSPR",all_dirs)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1.5,2.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(all_bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(all_bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
mtext(side=3, "Medium-lived", line=1, cex=3)
axis(2,at=seq(-1,2.5,by=1), las=2,cex.axis=2.5)
print.letter("(b)", cex=3)

i2 <- which(grepl("Long", all_dirs) & grepl("/LC1/", all_dirs) & grepl("/equil/", all_dirs) & grepl("LBSPR",all_dirs)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-2,3.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(all_bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(all_bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
mtext(side=3, "Longer-lived", line=1, cex=3)
mtext(side=4, "Equilibrium", line=2, cex=3)
axis(2,at=seq(-1,3,by=1), las=2,cex.axis=2.5)
print.letter("(c)", cex=3)

## base
i2 <- which(grepl("Short", all_dirs) & grepl("/LC1/", all_dirs) & grepl("base/", all_dirs) & grepl("Ramp", all_dirs) & grepl("LBSPR",all_dirs)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1,1.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(all_bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(all_bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
axis(2,at=seq(-0.5,1.5,by=0.5), las=2,cex.axis=2.5)
# axis(1, at=1:9, labels=c("Rich", "Index\nLC10", "Index\nLC1", "Catch\nLC10", "Catch\nLC1", "LC10", "LC1", "LBSPR\n10", "LBSPR\n1"), cex.axis=2.8)
print.letter("(d)", cex=3)

i2 <- which(grepl("Medium", all_dirs) & grepl("/LC1/", all_dirs) & grepl("base/", all_dirs) & grepl("Ramp", all_dirs) & grepl("LBSPR",all_dirs)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1.5,2.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(all_bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(all_bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
# axis(1, at=1:9, labels=c("Rich", "Index\nLC10", "Index\nLC1", "Catch\nLC10", "Catch\nLC1", "LC10", "LC1", "LBSPR\n10", "LBSPR\n1"), cex.axis=2.8)
axis(2,at=seq(-1,2.5,by=1), las=2,cex.axis=2.5)
print.letter("(e)", cex=3)

i2 <- which(grepl("Long", all_dirs) & grepl("/LC1/", all_dirs) & grepl("base/", all_dirs) & grepl("Ramp", all_dirs) & grepl("LBSPR",all_dirs)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-2,3.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(all_bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(all_bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
# axis(1, at=1:9, labels=c("Rich", "Index\nLC10", "Index\nLC1", "Catch\nLC10", "Catch\nLC1", "LC10", "LC1", "LBSPR\n10", "LBSPR\n1"), cex.axis=2.8)
axis(2,at=seq(-1,3,by=1), las=2,cex.axis=2.5)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
mtext(side=4, "Two-way trip", line=2, cex=3)
print.letter("(f)", cex=3)


## base
i2 <- which(grepl("Short", all_dirs) & grepl("/LC1/", all_dirs) & grepl("base/", all_dirs) & grepl("Increasing", all_dirs) & grepl("LBSPR",all_dirs)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1,1.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(all_bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(all_bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
axis(2,at=seq(-0.5,1.5,by=0.5), las=2,cex.axis=2.5)
axis(1, at=1:length(i2), labels=c(1000, 500, 200, 50, 20), cex.axis=2.8)
print.letter("(g)", cex=3)

i2 <- which(grepl("Medium", all_dirs) & grepl("/LC1/", all_dirs) & grepl("base/", all_dirs) & grepl("Increasing", all_dirs) & grepl("LBSPR",all_dirs)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1.5,2.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(all_bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(all_bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
axis(1, at=1:length(i2), labels=c(1000, 500, 200, 50, 20), cex.axis=2.8)
axis(2,at=seq(-1,2.5,by=1), las=2,cex.axis=2.5)
print.letter("(h)", cex=3)

i2 <- which(grepl("Long", all_dirs) & grepl("/LC1/", all_dirs) & grepl("base/", all_dirs) & grepl("Increasing", all_dirs) & grepl("LBSPR",all_dirs)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-2,3.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(all_bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(all_bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
axis(1, at=1:length(i2), labels=c(1000, 500, 200, 50, 20), cex.axis=2.8)
axis(2,at=seq(-1,3,by=1), las=2,cex.axis=2.5)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
mtext(side=4, "One-way trip", line=2, cex=3)
mtext(side=2, outer=TRUE, line=3, "Relative error", cex=3)
mtext(side=1, outer=TRUE, line=6, "Annual sample size of length measurements", cex=3)
print.letter("(i)", cex=3)
dev.off()




png(file.path(fig_dir, "IntervalCoverage.png"), height=15, width=28, res=200, units="in")
par(mfrow=c(3,3), mar=c(0,0,0,0), omi=c(1.2,1.8,1,1), mgp=c(3,2,0))
## equilibrium
i2 <- which(grepl("Short", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("/equil/", all_dirs) & grepl("LBSPR", all_dirs)==FALSE)
	all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,7.5),ylim=c(0,1.2))
print.letter("(a)", cex=3)
legend("bottomleft", pch=c(19,17), col=c("#11111175", "#AA000075"), legend=c("Interval coverage", "Convergence rate"), cex=3)
# cov <- sapply(1:length(i2), function(x) round(sum(all_ic$cover[which(all_ic$converge[,i2[x]]==1),i2[x]])/length(which(all_ic$converge[,i2[x]]==1)),2))
cov <- colSums(all_ic$cover[,i2])/nrow(all_ic$cover[,i2])
conv <- colSums(all_ic$converge[,i2])/nrow(all_ic$converge[,i2])
abline(lty=2,h=0.5, col=gray(0.3))
abline(lty=2,h=1, col="#AA0000")
points(x=1:7, y=conv, col="#AA000075", xpd=NA, cex=5, pch=17, lwd=3)
points(x=1:7, y=cov, col="#11111175", xpd=NA, cex=5, pch=19, lwd=3)
axis(2,at=seq(0.25,1,by=0.25), las=2,cex.axis=3.5)
mtext(side=3, "Short-lived", line=1, cex=3)
i2 <- which(grepl("Medium", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("/equil/", all_dirs) & grepl("LBSPR", all_dirs)==FALSE)
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,7.5),ylim=c(0,1.2))
print.letter("(b)", cex=3)
# cov <- sapply(1:length(i2), function(x) round(sum(all_ic$cover[which(all_ic$converge[,i2[x]]==1),i2[x]])/length(which(all_ic$converge[,i2[x]]==1)),2))
cov <- colSums(all_ic$cover[,i2])/nrow(all_ic$cover[,i2])
conv <- colSums(all_ic$converge[,i2])/nrow(all_ic$converge[,i2])
abline(lty=2,h=0.5, col=gray(0.3))
abline(lty=2,h=1, col="#AA0000")
points(x=1:7, y=conv, col="#AA000075", xpd=NA, cex=5, pch=17, lwd=3)
points(x=1:7, y=cov, col="#11111175", xpd=NA, cex=5, pch=19, lwd=3)
mtext(side=3, "Medium-lived", line=1, cex=3)
i2 <- which(grepl("Long", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("/equil/", all_dirs) & grepl("LBSPR", all_dirs)==FALSE)
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,7.5),ylim=c(0,1.2))
print.letter("(c)", cex=3)
# cov <- sapply(1:length(i2), function(x) round(sum(all_ic$cover[which(all_ic$converge[,i2[x]]==1),i2[x]])/length(which(all_ic$converge[,i2[x]]==1)),2))
cov <- colSums(all_ic$cover[,i2])/nrow(all_ic$cover[,i2])
conv <- colSums(all_ic$converge[,i2])/nrow(all_ic$converge[,i2])
abline(lty=2,h=0.5, col=gray(0.3))
abline(lty=2,h=1, col="#AA0000")
points(x=1:7, y=conv, col="#AA000075", xpd=NA, cex=5, pch=17, lwd=3)
points(x=1:7, y=cov, col="#11111175", xpd=NA, cex=5, pch=19, lwd=3)
mtext(side=3, "Longer-lived", line=1, cex=3)
mtext(side=4, "Equilibrium", line=2, cex=3)

## low var
i2 <- which(grepl("Short", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("/base/", all_dirs) & grepl("Ramp", all_dirs) & grepl("LBSPR", all_dirs)==FALSE)
	all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,7.5),ylim=c(0,1.2))
print.letter("(d)", cex=3)
# cov <- sapply(1:length(i2), function(x) round(sum(all_ic$cover[which(all_ic$converge[,i2[x]]==1),i2[x]])/length(which(all_ic$converge[,i2[x]]==1)),2))
cov <- colSums(all_ic$cover[,i2])/nrow(all_ic$cover[,i2])
conv <- colSums(all_ic$converge[,i2])/nrow(all_ic$converge[,i2])
abline(lty=2,h=0.5, col=gray(0.3))
abline(lty=2,h=1, col="#AA0000")
points(x=1:7, y=conv, col="#AA000075", xpd=NA, cex=5, pch=17, lwd=3)
points(x=1:7, y=cov, col="#11111175", xpd=NA, cex=5, pch=19, lwd=3)
axis(2,at=seq(0.25,1,by=0.25), las=2,cex.axis=3.5)
i2 <- which(grepl("Medium", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("/base/", all_dirs) & grepl("Ramp", all_dirs) & grepl("LBSPR", all_dirs)==FALSE)
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,7.5),ylim=c(0,1.2))
print.letter("(e)", cex=3)
# cov <- sapply(1:length(i2), function(x) round(sum(all_ic$cover[which(all_ic$converge[,i2[x]]==1),i2[x]])/length(which(all_ic$converge[,i2[x]]==1)),2))
cov <- colSums(all_ic$cover[,i2])/nrow(all_ic$cover[,i2])
conv <- colSums(all_ic$converge[,i2])/nrow(all_ic$converge[,i2])
abline(lty=2,h=0.5, col=gray(0.3))
abline(lty=2,h=1, col="#AA0000")
points(x=1:7, y=conv, col="#AA000075", xpd=NA, cex=5, pch=17, lwd=3)
points(x=1:7, y=cov, col="#11111175", xpd=NA, cex=5, pch=19, lwd=3)
i2 <- which(grepl("Long", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("/base/", all_dirs) & grepl("Ramp", all_dirs) & grepl("LBSPR", all_dirs)==FALSE)
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,7.5),ylim=c(0,1.2))
print.letter("(f)", cex=3)
# cov <- sapply(1:length(i2), function(x) round(sum(all_ic$cover[which(all_ic$converge[,i2[x]]==1),i2[x]])/length(which(all_ic$converge[,i2[x]]==1)),2))
cov <- colSums(all_ic$cover[,i2])/nrow(all_ic$cover[,i2])
conv <- colSums(all_ic$converge[,i2])/nrow(all_ic$converge[,i2])
abline(lty=2,h=0.5, col=gray(0.3))
abline(lty=2,h=1, col="#AA0000")
points(x=1:7, y=conv, col="#AA000075", xpd=NA, cex=5, pch=17, lwd=3)
points(x=1:7, y=cov, col="#11111175", xpd=NA, cex=5, pch=19, lwd=3)
mtext(side=4, "Two-way", line=2, cex=3)

## base
i2 <- which(grepl("Short", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("/base/", all_dirs) & grepl("Increasing", all_dirs) & grepl("LBSPR", all_dirs)==FALSE)
	all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,7.5),ylim=c(0,1.2))
print.letter("(g)", cex=3)
# cov <- sapply(1:length(i2), function(x) round(sum(all_ic$cover[which(all_ic$converge[,i2[x]]==1),i2[x]])/length(which(all_ic$converge[,i2[x]]==1)),2))
cov <- colSums(all_ic$cover[,i2])/nrow(all_ic$cover[,i2])
conv <- colSums(all_ic$converge[,i2])/nrow(all_ic$converge[,i2])
abline(lty=2,h=0.5, col=gray(0.3))
abline(lty=2,h=1, col="#AA0000")
points(x=1:7, y=conv, col="#AA000075", xpd=NA, cex=5, pch=17, lwd=3)
points(x=1:7, y=cov, col="#11111175", xpd=NA, cex=5, pch=19, lwd=3)
axis(2,at=seq(0.25,1,by=0.25), las=2,cex.axis=3.5)
axis(1, at=1:7, labels=c("Rich", "Index\nLC10", "Index\nLC1", "Catch\nLC10", "Catch\nLC1", "LC10", "LC1"), cex.axis=3)
i2 <- which(grepl("Medium", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("/base/", all_dirs) & grepl("Increasing", all_dirs) & grepl("LBSPR", all_dirs)==FALSE)
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,7.5),ylim=c(0,1.2))
print.letter("(h)", cex=3)
# cov <- sapply(1:length(i2), function(x) round(sum(all_ic$cover[which(all_ic$converge[,i2[x]]==1),i2[x]])/length(which(all_ic$converge[,i2[x]]==1)),2))
cov <- colSums(all_ic$cover[,i2])/nrow(all_ic$cover[,i2])
conv <- colSums(all_ic$converge[,i2])/nrow(all_ic$converge[,i2])
abline(lty=2,h=0.5, col=gray(0.3))
abline(lty=2,h=1, col="#AA0000")
points(x=1:7, y=conv, col="#AA000075", xpd=NA, cex=5, pch=17, lwd=3)
points(x=1:7, y=cov, col="#11111175", xpd=NA, cex=5, pch=19, lwd=3)
axis(1, at=1:7, labels=c("Rich", "Index\nLC10", "Index\nLC1", "Catch\nLC10", "Catch\nLC1", "LC10", "LC1"), cex.axis=3)
i2 <- which(grepl("Long", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("/base/", all_dirs) & grepl("Increasing", all_dirs) & grepl("LBSPR", all_dirs)==FALSE)
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,7.5),ylim=c(0,1.2))
print.letter("(i)", cex=3)
# cov <- sapply(1:length(i2), function(x) round(sum(all_ic$cover[which(all_ic$converge[,i2[x]]==1),i2[x]])/length(which(all_ic$converge[,i2[x]]==1)),2))
cov <- colSums(all_ic$cover[,i2])/nrow(all_ic$cover[,i2])
conv <- colSums(all_ic$converge[,i2])/nrow(all_ic$converge[,i2])
abline(lty=2,h=0.5, col=gray(0.3))
abline(lty=2,h=1, col="#AA0000")
points(x=1:7, y=conv, col="#AA000075", xpd=NA, cex=5, pch=17, lwd=3)
points(x=1:7, y=cov, col="#11111175", xpd=NA, cex=5, pch=19, lwd=3)
axis(1, at=1:7, labels=c("Rich", "Index\nLC10", "Index\nLC1", "Catch\nLC10", "Catch\nLC1", "LC10", "LC1"), cex.axis=3)
mtext(side=4, "One-way", line=2, cex=3)
mtext(side=2, outer=TRUE, line=8.5, "Proportion of iterations", cex=3)
mtext(side=1, outer=TRUE, line=6, "Data availability scenario", cex=3)
dev.off()



## instantaneous vs. continuous vs. LBSPR
# icol <- rev(brewer.pal(3,"Purples"))
# ccol <- rev(brewer.pal(3, "Oranges"))
lcol <- rev(brewer.pal(3, "Blues"))
rcol <- rev(brewer.pal(3, "Greens"))
col_vec <- c(lcol[1:2], rcol[1:2])


png(file.path(fig_dir, "EquilCompare_RelativeError_n200.png"), height=16, width=35, res=200, units="in")
par(mfrow=c(3,3), mar=c(0,7,0,0), omi=c(1.6,1.2,1,1.5), mgp=c(3,2,0))
## equilibrium
i2 <- which(grepl("Short", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("/equil/", all_dirs) & grepl("Index", all_dirs)==FALSE & grepl("Catch", all_dirs)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1,1))
b1 <- sapply(1:length(i2), function(x) round(median(abs(all_bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(all_bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=2.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=2.5)
abline(h=0, lty=2)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
axis(2,at=seq(-0.5,1,by=0.5), las=2,cex.axis=3.5)
mtext(side=3, "Short-lived", line=1, cex=3)
print.letter("(a)", cex=3)

i2 <- which(grepl("Medium", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("/equil/", all_dirs) & grepl("Index", all_dirs)==FALSE & grepl("Catch", all_dirs)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1,1.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(all_bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(all_bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
mtext(side=3, "Medium-lived", line=1, cex=3)
axis(2,at=seq(-0.5,1.5,by=0.5), las=2,cex.axis=3.5)
print.letter("(b)", cex=3)

i2 <- which(grepl("Long", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("/equil/", all_dirs) & grepl("Index", all_dirs)==FALSE & grepl("Catch", all_dirs)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-2,2))
b1 <- sapply(1:length(i2), function(x) round(median(abs(all_bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(all_bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
mtext(side=3, "Longer-lived", line=1, cex=3)
mtext(side=4, "Instantaneous\nsampling", line=6.5, cex=3)
axis(2,at=seq(-1,2,by=1), las=2,cex.axis=3.5)
print.letter("(c)", cex=3)

## base
i2 <- which(grepl("Short", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("/equil_12seasons/", all_dirs) & grepl("Index", all_dirs)==FALSE & grepl("Catch", all_dirs)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1.5,3))
b1 <- sapply(1:length(i2), function(x) round(median(abs(all_bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(all_bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
axis(2,at=seq(-1,3,by=1), las=2,cex.axis=3.5)
# axis(1, at=1:9, labels=c("Rich", "Index\nLC10", "Index\nLC1", "Catch\nLC10", "Catch\nLC1", "LC10", "LC1", "LBSPR\n10", "LBSPR\n1"), cex.axis=2.8)
print.letter("(d)", cex=3)

i2 <- which(grepl("Medium", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("/equil_12seasons/", all_dirs) & grepl("Index", all_dirs)==FALSE & grepl("Catch", all_dirs)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1,1.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(all_bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(all_bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
# axis(1, at=1:9, labels=c("Rich", "Index\nLC10", "Index\nLC1", "Catch\nLC10", "Catch\nLC1", "LC10", "LC1", "LBSPR\n10", "LBSPR\n1"), cex.axis=2.8)
axis(2,at=seq(-0.5,1,by=0.5), las=2,cex.axis=3.5)
print.letter("(e)", cex=3)

i2 <- which(grepl("Long", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("/equil_12seasons/", all_dirs) & grepl("Index", all_dirs)==FALSE & grepl("Catch", all_dirs)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1.5,3))
b1 <- sapply(1:length(i2), function(x) round(median(abs(all_bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(all_bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
# axis(1, at=1:9, labels=c("Rich", "Index\nLC10", "Index\nLC1", "Catch\nLC10", "Catch\nLC1", "LC10", "LC1", "LBSPR\n10", "LBSPR\n1"), cex.axis=2.8)
axis(2,at=seq(-1,3,by=1), las=2,cex.axis=3.5)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
mtext(side=4, "Continuous\nsampling", line=6.5, cex=3)
print.letter("(f)", cex=3)


## base
i2 <- which(grepl("Short", all_dirs) & grepl("/equil_LBSPR/", all_dirs) & grepl("Index", all_dirs)==FALSE & grepl("Catch", all_dirs)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-0.5,1.5))
b1 <- sapply(1:length(i2), function(x) round(median(abs(all_bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(all_bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
axis(2,at=seq(-0.5,1.5,by=0.5), las=2,cex.axis=3.5)
axis(1, at=1:length(i2), labels=c("LC10", "LC1", "LBSPR10", "LBSPR1"), cex.axis=3.5)
print.letter("(g)", cex=3)

i2 <- which(grepl("Medium", all_dirs) & grepl("/equil_LBSPR/", all_dirs) & grepl("Index", all_dirs)==FALSE & grepl("Catch", all_dirs)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-0.5,1))
b1 <- sapply(1:length(i2), function(x) round(median(abs(all_bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(all_bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
axis(1, at=1:length(i2), labels=c("LC10", "LC1", "LBSPR10", "LBSPR1"), cex.axis=3.5)
axis(2,at=seq(-0.5,1,by=0.5), las=2,cex.axis=3.5)
print.letter("(h)", cex=3)

i2 <- which(grepl("Long", all_dirs) & grepl("/equil_LBSPR/", all_dirs) & grepl("Index", all_dirs)==FALSE & grepl("Catch", all_dirs)==FALSE)
# i2 <- i1[-c(which(grepl("LC5/",all_dirs[i1]) | grepl("LC2/",all_dirs[i1])))]
all_dirs[i2]
plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,length(i2)+0.5),ylim=c(-1,1))
b1 <- sapply(1:length(i2), function(x) round(median(abs(all_bp$dev[,i2[x]]), na.rm=TRUE),3))
p1 <- round(all_bp$precision[i2],3)
# text(x=1:length(b1), y=0.95*3, b1, cex=3.5)
# text(x=1:length(p1), y=0.8*3, paste0("(",p1,")"), cex=3.5)
abline(h=0, lty=2)
axis(1, at=1:length(i2), labels=c("LC10", "LC1", "LBSPR10", "LBSPR1"), cex.axis=3.5)
axis(2,at=seq(-1,1,by=0.5), las=2,cex.axis=3.5)
beanplot(as.data.frame(all_bp$relerr[,i2]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
mtext(side=4, "LB-SPR OM", line=2, cex=3)
mtext(side=2, outer=TRUE, line=3, "Relative error", cex=3)
mtext(side=1, outer=TRUE, line=8, "Model/Data availability scenario", cex=3)
print.letter("(i)", cex=3)
dev.off()


############################################################
#### tables and values for paper
############################################################
## data availability scenarios
data_vec <- c("Index_Catch_LC20","Index_LC10","Index_LC1","Catch_LC10","Catch_LC1","LC10","LC1","LBSPR10","LBSPR1")
bmat <- pmat <- matrix(NA, nrow=9, ncol=9)
rownames(bmat) <- rownames(pmat) <- data_vec
colnames(bmat) <- colnames(pmat) <- rep(c("equil", "two-way", "one-way"),3)

for(i in 1:length(data_vec)){
	idir <- c(which(grepl("Short", all_dirs) & grepl(paste0("/", data_vec[i],"/"), all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("equil_", all_dirs)==FALSE), which(grepl("Medium", all_dirs) & grepl(paste0("/", data_vec[i],"/"), all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("equil_", all_dirs)==FALSE), which(grepl("Long", all_dirs) & grepl(paste0("/", data_vec[i],"/"), all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("equil_", all_dirs)==FALSE))
	b <- sapply(1:length(idir), function(x) median((all_bp$relerr[,idir[x]]), na.rm=TRUE))
	p <- sapply(1:length(idir), function(x) median(abs(all_bp$relerr[,idir[x]]), na.rm=TRUE))
	bmat[i,] <- b
	pmat[i,] <- p
	rm(b)
	rm(p)
}

write.csv(rbind(bmat, pmat), file.path(res_dir, "base_biases_precision.csv"))


### equilibrium scenarios - LBSPR
data_vec <- c("LC10","LC1","LBSPR10","LBSPR1")
bmat2 <- pmat2 <- matrix(NA, nrow=length(data_vec), ncol=9)
rownames(bmat2) <- rownames(pmat2) <- data_vec
colnames(bmat2) <- colnames(pmat2) <- rep(c("instantaneous", "continuous", "LB-SPR"),length(lh_vec))

for(i in 1:length(data_vec)){
	idir <- c(which(grepl("Short", all_dirs) & grepl(paste0("/", data_vec[i],"/"), all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("equil",all_dirs)), which(grepl("Short", all_dirs) & grepl(paste0("/", data_vec[i],"/"), all_dirs) & grepl("equil_LBSPR",all_dirs)), which(grepl("Medium", all_dirs) & grepl(paste0("/", data_vec[i],"/"), all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("equil",all_dirs)), which(grepl("Medium", all_dirs) & grepl(paste0("/", data_vec[i],"/"), all_dirs) & grepl("equil_LBSPR",all_dirs)), which(grepl("Long", all_dirs) & grepl(paste0("/", data_vec[i],"/"), all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("equil",all_dirs)), which(grepl("Long", all_dirs) & grepl(paste0("/", data_vec[i],"/"), all_dirs) & grepl("equil_LBSPR",all_dirs)))
	b <- sapply(1:length(idir), function(x) median((all_bp$relerr[,idir[x]]), na.rm=TRUE))
	p <- sapply(1:length(idir), function(x) median(abs(all_bp$relerr[,idir[x]]), na.rm=TRUE))
	bmat2[i,] <- b
	pmat2[i,] <- p
	rm(b)
	rm(p)
}

write.csv(rbind(bmat2, pmat2), file.path(res_dir, "equil_biases_precision.csv"))



## sample size
data_vec <- paste0("SampleSize_", c(1000,500,200,50,20))
bmat <- pmat <- matrix(NA, nrow=5, ncol=9)
rownames(bmat) <- rownames(pmat) <- data_vec
colnames(bmat) <- colnames(pmat) <- rep(c("equil", "two-way", "one-way"),3)

for(i in 1:length(data_vec)){
	idir <- c(which(grepl("Short", all_dirs) & grepl("/LC10/", all_dirs) & grepl(paste0(data_vec[i],"/"), all_dirs) & grepl("equil_", all_dirs)==FALSE), which(grepl("Medium", all_dirs) & grepl("/LC10/", all_dirs) & grepl(paste0(data_vec[i],"/"), all_dirs)  & grepl("equil_", all_dirs)==FALSE), which(grepl("Long", all_dirs) & grepl("/LC10/", all_dirs) & grepl(paste0(data_vec[i],"/"), all_dirs)  & grepl("equil_", all_dirs)==FALSE))
	b <- sapply(1:length(idir), function(x) median((all_bp$relerr[,idir[x]]), na.rm=TRUE))
	p <- sapply(1:length(idir), function(x) median(abs(all_bp$relerr[,idir[x]]), na.rm=TRUE))
	bmat[i,] <- b
	pmat[i,] <- p
	rm(b)
	rm(p)
}

write.csv(rbind(bmat,pmat), file.path(res_dir, "sample_size_biases_precision.csv"))


## abstract - improvement using LIME from LB-SPR with 1 year length comp for variability scenarios
idir <- c(which(grepl("SampleSize_200", all_dirs) & grepl("/LC1/", all_dirs)))
all_dirs[idir]

xdir <- which(grepl("SampleSize_200", all_dirs) & grepl("/LBSPR1/", all_dirs))
all_dirs[xdir]
(median(all_bp$relerr[,xdir], na.rm=TRUE) - (median(all_bp$relerr[,idir], na.rm=TRUE)))/(median(all_bp$relerr[,idir], na.rm=TRUE))
(median(abs(all_bp$relerr[,xdir]), na.rm=TRUE) - (median(abs(all_bp$relerr[,idir]), na.rm=TRUE)))/(median(abs(all_bp$relerr[,idir]), na.rm=TRUE))


### median bias for each life history type, 1 year length comp, 200 samples across scenarios
idir <- which(grepl("Short", all_dirs) & grepl("SampleSize_200", all_dirs) & grepl("/LC1/", all_dirs) & grepl("equil_", all_dirs)==FALSE)
all_dirs[idir]
median((all_bp$relerr[,idir]), na.rm=TRUE)

idir <- which(grepl("Medium", all_dirs) & grepl("SampleSize_200", all_dirs) & grepl("/LC1/", all_dirs) & grepl("equil_", all_dirs)==FALSE)
all_dirs[idir]
median((all_bp$relerr[,idir]), na.rm=TRUE)


idir <- which(grepl("Long", all_dirs) & grepl("SampleSize_200", all_dirs) & grepl("/LC1/", all_dirs) & grepl("equil_", all_dirs)==FALSE)
all_dirs[idir]
median((all_bp$relerr[,idir]), na.rm=TRUE)

## median percent improvement from 1 yr to 10 years
idir <- which(grepl("Short", all_dirs) & grepl("SampleSize_200", all_dirs) & grepl("/LC1/", all_dirs) & grepl("equil_", all_dirs)==FALSE)
all_dirs[idir]

xdir <- which(grepl("Short", all_dirs) & grepl("SampleSize_200", all_dirs) & grepl("/LC10/", all_dirs) & grepl("equil_", all_dirs)==FALSE)
(abs(median(all_bp$relerr[,xdir], na.rm=TRUE)) - abs(median(all_bp$relerr[,idir], na.rm=TRUE)))/abs(median(all_bp$relerr[,idir], na.rm=TRUE))

idir <- which(grepl("Medium", all_dirs) & grepl("SampleSize_200", all_dirs) & grepl("/LC1/", all_dirs) & grepl("equil_", all_dirs)==FALSE)
all_dirs[idir]

xdir <- which(grepl("Medium", all_dirs) & grepl("SampleSize_200", all_dirs) & grepl("/LC10/", all_dirs) & grepl("equil_", all_dirs)==FALSE)
(abs(median(all_bp$relerr[,xdir], na.rm=TRUE)) - abs(median(all_bp$relerr[,idir], na.rm=TRUE)))/abs(median(all_bp$relerr[,idir], na.rm=TRUE))

idir <- which(grepl("Long", all_dirs) & grepl("SampleSize_200", all_dirs) & grepl("/LC1/", all_dirs) & grepl("equil_", all_dirs)==FALSE)
all_dirs[idir]

xdir <- which(grepl("Long", all_dirs) & grepl("SampleSize_200", all_dirs) & grepl("/LC10/", all_dirs) & grepl("equil_", all_dirs)==FALSE)
(abs(median(all_bp$relerr[,xdir], na.rm=TRUE)) - abs(median(all_bp$relerr[,idir], na.rm=TRUE)))/abs(median(all_bp$relerr[,idir], na.rm=TRUE))


## median percent improvement integrating an abundance index with 1 year of length comp
idir <- which(grepl("SampleSize_200", all_dirs) & grepl("/LC1/", all_dirs) & grepl("equil_", all_dirs)==FALSE)
all_dirs[idir]

xdir <- which(grepl("SampleSize_200", all_dirs) & grepl("/Index_LC1/", all_dirs) & grepl("equil_", all_dirs)==FALSE)
(abs(median(all_bp$relerr[,xdir], na.rm=TRUE)) - abs(median(all_bp$relerr[,idir], na.rm=TRUE)))/abs(median(all_bp$relerr[,idir], na.rm=TRUE))
(abs(median(abs(all_bp$relerr[,xdir]), na.rm=TRUE)) - abs(median(abs(all_bp$relerr[,idir]), na.rm=TRUE)))/abs(median(abs(all_bp$relerr[,idir]), na.rm=TRUE))


## median percent change integrating catch with 1 year of length comp - equilibrium
idir <- which(grepl("SampleSize_200", all_dirs) & grepl("/LC1/", all_dirs) & grepl("equil_", all_dirs)==FALSE & grepl("/equil/", all_dirs))
all_dirs[idir]

xdir <- which(grepl("SampleSize_200", all_dirs) & grepl("/Catch_LC1/", all_dirs) & grepl("equil_", all_dirs)==FALSE & grepl("/equil/", all_dirs))
(abs(median(all_bp$relerr[,xdir], na.rm=TRUE)) - abs(median(all_bp$relerr[,idir], na.rm=TRUE)))/abs(median(all_bp$relerr[,idir], na.rm=TRUE))
(abs(median(abs(all_bp$relerr[,xdir]), na.rm=TRUE)) - abs(median(abs(all_bp$relerr[,idir]), na.rm=TRUE)))/abs(median(abs(all_bp$relerr[,idir]), na.rm=TRUE))

### median bias for each life history type, comparing 1 year LC and catch+1
idir <- which(grepl("Short", all_dirs) & grepl("SampleSize_200", all_dirs) & grepl("/LC1/", all_dirs) & grepl("Ramp", all_dirs))
all_dirs[idir]
median((all_bp$relerr[,idir]), na.rm=TRUE)

idir <- which(grepl("Medium", all_dirs) & grepl("SampleSize_200", all_dirs) & grepl("/LC1/", all_dirs) & grepl("equil_", all_dirs)==FALSE)
all_dirs[idir]
median((all_bp$relerr[,idir]), na.rm=TRUE)

idir <- which(grepl("Short", all_dirs) & grepl("SampleSize_200", all_dirs) & grepl("/Catch_LC1/", all_dirs) & grepl("Ramp", all_dirs))
all_dirs[idir]
median((all_bp$relerr[,idir]), na.rm=TRUE)

idir <- which(grepl("Medium", all_dirs) & grepl("SampleSize_200", all_dirs) & grepl("/Catch_LC1/", all_dirs) & grepl("equil_", all_dirs)==FALSE)
all_dirs[idir]
median((all_bp$relerr[,idir]), na.rm=TRUE)

## median precision for 1000, 500, 200, 50, 20 samples
idir <- which(grepl("SampleSize_1000/", all_dirs) & grepl("/LC1/", all_dirs) & grepl("equil_", all_dirs)==FALSE)
all_dirs[idir]
median(all_bp$relerr[,idir], na.rm=TRUE)

idir <- which(grepl("SampleSize_500/", all_dirs) & grepl("/LC1/", all_dirs) & grepl("equil_", all_dirs)==FALSE)
all_dirs[idir]
median(all_bp$relerr[,idir], na.rm=TRUE)

idir <- which(grepl("SampleSize_200/", all_dirs) & grepl("/LC1/", all_dirs) & grepl("equil_", all_dirs)==FALSE)
all_dirs[idir]
median(all_bp$relerr[,idir], na.rm=TRUE)

idir <- which(grepl("SampleSize_50/", all_dirs) & grepl("/LC1/", all_dirs) & grepl("equil_", all_dirs)==FALSE)
all_dirs[idir]
median(all_bp$relerr[,idir], na.rm=TRUE)

idir <- which(grepl("SampleSize_20/", all_dirs) & grepl("/LC1/", all_dirs) & grepl("equil_", all_dirs)==FALSE)
all_dirs[idir]
median(all_bp$relerr[,idir], na.rm=TRUE)

## median precision for 1000, 20 samples by life history type
idir <- which(grepl("Short", all_dirs) & grepl("SampleSize_1000/", all_dirs) & grepl("/LC1/", all_dirs) & grepl("equil_", all_dirs)==FALSE)
all_dirs[idir]
median(all_bp$relerr[,idir], na.rm=TRUE)

idir <- which(grepl("Short", all_dirs) & grepl("SampleSize_20/", all_dirs) & grepl("/LC1/", all_dirs) & grepl("equil_", all_dirs)==FALSE)
all_dirs[idir]
median(all_bp$relerr[,idir], na.rm=TRUE)

idir <- which(grepl("Long", all_dirs) & grepl("SampleSize_1000/", all_dirs) & grepl("/LC1/", all_dirs) & grepl("equil_", all_dirs)==FALSE)
all_dirs[idir]
median(all_bp$relerr[,idir], na.rm=TRUE)

idir <- which(grepl("Long", all_dirs) & grepl("SampleSize_20/", all_dirs) & grepl("/LC1/", all_dirs) & grepl("equil_", all_dirs)==FALSE)
all_dirs[idir]
median(all_bp$relerr[,idir], na.rm=TRUE)



## convergence rate
idir <- which(grepl("LBSPR", all_dirs)==FALSE)
all_dirs[idir]
sum(all_ic$converge[,idir])/length(all_ic$converge[,idir])

idir <- which(grepl("LBSPR", all_dirs)==FALSE & grepl("Index_Catch_LC20", all_dirs))
all_dirs[idir]
sum(all_ic$converge[,idir])/length(all_ic$converge[,idir])

idir <- which(grepl("LBSPR", all_dirs)==FALSE & grepl("/LC1", all_dirs))
all_dirs[idir]
sum(all_ic$converge[,idir])/length(all_ic$converge[,idir])

idir <- which(grepl("LBSPR", all_dirs)==FALSE & grepl("/Catch", all_dirs))
all_dirs[idir]
sum(all_ic$converge[,idir])/length(all_ic$converge[,idir])


idir <- which(grepl("LBSPR", all_dirs)==FALSE & grepl("/Index_LC", all_dirs))
all_dirs[idir]
sum(all_ic$converge[,idir])/length(all_ic$converge[,idir])



## coverage interval
idir <- which(grepl("LBSPR", all_dirs)==FALSE)
all_dirs[idir]
sum(all_ic$cover[,idir])/length(all_ic$cover[,idir])

idir <- which(grepl("LBSPR", all_dirs)==FALSE & grepl("Index_Catch_LC20", all_dirs))
all_dirs[idir]
sum(all_ic$cover[,idir])/length(all_ic$cover[,idir])

idir <- which(grepl("LBSPR", all_dirs)==FALSE & grepl("/LC1", all_dirs))
all_dirs[idir]
sum(all_ic$cover[,idir])/length(all_ic$cover[,idir])

idir <- which(grepl("LBSPR", all_dirs)==FALSE & grepl("/Catch", all_dirs))
all_dirs[idir]
sum(all_ic$cover[,idir])/length(all_ic$cover[,idir])


idir <- which(grepl("LBSPR", all_dirs)==FALSE & grepl("/Index_LC", all_dirs))
all_dirs[idir]
sum(all_ic$cover[,idir])/length(all_ic$cover[,idir])



## coverage interval -- by variation
idir <- which(grepl("/equil/", all_dirs) & grepl("LBSPR", all_dirs)==FALSE)
all_dirs[idir]
sum(all_ic$cover[,idir])/length(all_ic$cover[,idir])

idir <- which(grepl("Ramp/", all_dirs) & grepl("LBSPR", all_dirs)==FALSE)
all_dirs[idir]
sum(all_ic$cover[,idir])/length(all_ic$cover[,idir])

idir <- which(grepl("Increasing/", all_dirs) & grepl("LBSPR", all_dirs)==FALSE)
all_dirs[idir]
sum(all_ic$cover[,idir])/length(all_ic$cover[,idir])

## coverage interval -- by life history
idir <- which(grepl("Short", all_dirs) & grepl("LBSPR", all_dirs)==FALSE)
all_dirs[idir]
sum(all_ic$cover[,idir])/length(all_ic$cover[,idir])

idir <- which(grepl("Medium", all_dirs) & grepl("LBSPR", all_dirs)==FALSE)
all_dirs[idir]
sum(all_ic$cover[,idir])/length(all_ic$cover[,idir])

idir <- which(grepl("Long", all_dirs) & grepl("LBSPR", all_dirs)==FALSE)
all_dirs[idir]
sum(all_ic$cover[,idir])/length(all_ic$cover[,idir])






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

all_sens_dir <- c(sens_equil_dir_vec, sens_base_dir_vec)

# sens_bp <- bias_precision(dirs=all_sens_dir, itervec=itervec)
# saveRDS(sens_bp, file.path(res_dir, "sens_bias_precision.rds"))
sens_bp <- readRDS(file.path(res_dir, "sens_bias_precision.rds"))

col_vec <- c(rep("steelblue",3), rep("tomato",3))

png(file.path(fig_dir, "ParamSensitivities.png"), height=16, width=38, res=200, units="in")
par(mfrow=c(3,5), mar=c(0,0,0,0), omi=c(1.2,1.2,1,1), mgp=c(3,2,0))
## short
i1 <- which(grepl("Short", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("equil/", all_sens_dir) & grepl("/M/", all_sens_dir))
	all_sens_dir[i1]
i2 <- which(grepl("Short", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("equil/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i2]

i3 <- which(grepl("Short", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("Ramp/", all_sens_dir) & grepl("/M/", all_sens_dir))
	all_sens_dir[i3]
i4 <- which(grepl("Short", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("Ramp/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i4]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,6.5),ylim=c(-1,1))
abline(h=0, lty=2)
pdf <- cbind.data.frame(sens_bp$relerr[,i1[2]], all_bp$relerr[,i2], sens_bp$relerr[,i1[1]], sens_bp$relerr[,i3[2]], all_bp$relerr[,i4], sens_bp$relerr[,i3[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
axis(2,at=seq(-0.5,1,by=0.5), las=2,cex.axis=3.5)
mtext(side=3, "Natural mortality rate", line=1, cex=3)
print.letter("(a)", cex=3)

i1 <- which(grepl("Short", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("equil/", all_sens_dir) & grepl("/linf/", all_sens_dir))
	all_sens_dir[i1]
i2 <- which(grepl("Short", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("equil/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i2]

i3 <- which(grepl("Short", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("Ramp/", all_sens_dir) & grepl("/linf/", all_sens_dir))
	all_sens_dir[i3]
i4 <- which(grepl("Short", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("Ramp/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i4]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,6.5),ylim=c(-1,1))
abline(h=0, lty=2)
pdf <- cbind.data.frame(sens_bp$relerr[,i1[2]], all_bp$relerr[,i2], sens_bp$relerr[,i1[1]], sens_bp$relerr[,i3[2]], all_bp$relerr[,i4], sens_bp$relerr[,i3[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
mtext(side=3, "Asymptotic length", line=1, cex=3)
print.letter("(b)", cex=3)

i1 <- which(grepl("Short", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("equil/", all_sens_dir) & grepl("/vbk/", all_sens_dir))
	all_sens_dir[i1]
i2 <- which(grepl("Short", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("equil/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i2]

i3 <- which(grepl("Short", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("Ramp/", all_sens_dir) & grepl("/vbk/", all_sens_dir))
	all_sens_dir[i3]
i4 <- which(grepl("Short", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("Ramp/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i4]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,6.5),ylim=c(-1,1))
abline(h=0, lty=2)
pdf <- cbind.data.frame(sens_bp$relerr[,i1[2]], all_bp$relerr[,i2], sens_bp$relerr[,i1[1]], sens_bp$relerr[,i3[2]], all_bp$relerr[,i4], sens_bp$relerr[,i3[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
mtext(side=3, "von Bertalanffy k", line=1, cex=3)
print.letter("(c)", cex=3)


i1 <- which(grepl("Short", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("equil/", all_sens_dir) & grepl("/ML50/", all_sens_dir))
	all_sens_dir[i1]
i2 <- which(grepl("Short", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("equil/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i2]

i3 <- which(grepl("Short", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("Ramp/", all_sens_dir) & grepl("/ML50/", all_sens_dir))
	all_sens_dir[i3]
i4 <- which(grepl("Short", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("Ramp/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i4]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,6.5),ylim=c(-1,1))
abline(h=0, lty=2)
pdf <- cbind.data.frame(sens_bp$relerr[,i1[2]], all_bp$relerr[,i2], sens_bp$relerr[,i1[1]], sens_bp$relerr[,i3[2]], all_bp$relerr[,i4], sens_bp$relerr[,i3[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
mtext(side=3, "Length at 50% maturity", line=1, cex=3)
print.letter("(d)", cex=3)

i1 <- which(grepl("Short", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("equil/", all_sens_dir) & grepl("/CVlen/", all_sens_dir))
	all_sens_dir[i1]
i2 <- which(grepl("Short", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("equil/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i2]

i3 <- which(grepl("Short", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("Ramp/", all_sens_dir) & grepl("/CVlen/", all_sens_dir))
	all_sens_dir[i3]
i4 <- which(grepl("Short", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("Ramp/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i4]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,6.5),ylim=c(-1,1))
abline(h=0, lty=2)
pdf <- cbind.data.frame(sens_bp$relerr[,i1[2]], all_bp$relerr[,i2], sens_bp$relerr[,i1[1]], sens_bp$relerr[,i3[2]], all_bp$relerr[,i4], sens_bp$relerr[,i3[1]])
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
i2 <- which(grepl("Medium", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("equil/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i2]

i3 <- which(grepl("Medium", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("Ramp/", all_sens_dir) & grepl("/M/", all_sens_dir))
	all_sens_dir[i3]
i4 <- which(grepl("Medium", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("Ramp/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i4]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,6.5),ylim=c(-2,2))
abline(h=0, lty=2)
pdf <- cbind.data.frame(sens_bp$relerr[,i1[2]], all_bp$relerr[,i2], sens_bp$relerr[,i1[1]], sens_bp$relerr[,i3[2]], all_bp$relerr[,i4], sens_bp$relerr[,i3[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
axis(2,at=seq(-1,2,by=1), las=2,cex.axis=3.5)
print.letter("(f)", cex=3)

i1 <- which(grepl("Medium", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("equil/", all_sens_dir) & grepl("/linf/", all_sens_dir))
	all_sens_dir[i1]
i2 <- which(grepl("Medium", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("equil/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i2]

i3 <- which(grepl("Medium", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("Ramp/", all_sens_dir) & grepl("/linf/", all_sens_dir))
	all_sens_dir[i3]
i4 <- which(grepl("Medium", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("Ramp/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i4]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,6.5),ylim=c(-2,2))
abline(h=0, lty=2)
pdf <- cbind.data.frame(sens_bp$relerr[,i1[2]], all_bp$relerr[,i2], sens_bp$relerr[,i1[1]], sens_bp$relerr[,i3[2]], all_bp$relerr[,i4], sens_bp$relerr[,i3[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
print.letter("(g)", cex=3)

i1 <- which(grepl("Medium", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("equil/", all_sens_dir) & grepl("/vbk/", all_sens_dir))
	all_sens_dir[i1]
i2 <- which(grepl("Medium", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("equil/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i2]

i3 <- which(grepl("Medium", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("Ramp/", all_sens_dir) & grepl("/vbk/", all_sens_dir))
	all_sens_dir[i3]
i4 <- which(grepl("Medium", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("Ramp/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i4]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,6.5),ylim=c(-2,2))
abline(h=0, lty=2)
pdf <- cbind.data.frame(sens_bp$relerr[,i1[2]], all_bp$relerr[,i2], sens_bp$relerr[,i1[1]], sens_bp$relerr[,i3[2]], all_bp$relerr[,i4], sens_bp$relerr[,i3[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
print.letter("(h)", cex=3)


i1 <- which(grepl("Medium", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("equil/", all_sens_dir) & grepl("/ML50/", all_sens_dir))
	all_sens_dir[i1]
i2 <- which(grepl("Medium", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("equil/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i2]

i3 <- which(grepl("Medium", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("Ramp/", all_sens_dir) & grepl("/ML50/", all_sens_dir))
	all_sens_dir[i3]
i4 <- which(grepl("Medium", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("Ramp/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i4]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,6.5),ylim=c(-2,2))
abline(h=0, lty=2)
pdf <- cbind.data.frame(sens_bp$relerr[,i1[2]], all_bp$relerr[,i2], sens_bp$relerr[,i1[1]], sens_bp$relerr[,i3[2]], all_bp$relerr[,i4], sens_bp$relerr[,i3[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
print.letter("(i)", cex=3)

i1 <- which(grepl("Medium", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("equil/", all_sens_dir) & grepl("/CVlen/", all_sens_dir))
	all_sens_dir[i1]
i2 <- which(grepl("Medium", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("equil/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i2]

i3 <- which(grepl("Medium", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("Ramp/", all_sens_dir) & grepl("/CVlen/", all_sens_dir))
	all_sens_dir[i3]
i4 <- which(grepl("Medium", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("Ramp/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i4]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,6.5),ylim=c(-2,2))
abline(h=0, lty=2)
pdf <- cbind.data.frame(sens_bp$relerr[,i1[2]], all_bp$relerr[,i2], sens_bp$relerr[,i1[1]], sens_bp$relerr[,i3[2]], all_bp$relerr[,i4], sens_bp$relerr[,i3[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
mtext(side=4, "Medium-lived", line=2, cex=3)
print.letter("(j)", cex=3)


## long
i1 <- which(grepl("Long", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("equil/", all_sens_dir) & grepl("/M/", all_sens_dir))
	all_sens_dir[i1]
i2 <- which(grepl("Long", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("equil/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i2]

i3 <- which(grepl("Long", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("Ramp/", all_sens_dir) & grepl("/M/", all_sens_dir))
	all_sens_dir[i3]
i4 <- which(grepl("Long", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("Ramp/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i4]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,6.5),ylim=c(-3,3))
abline(h=0, lty=2)
pdf <- cbind.data.frame(sens_bp$relerr[,i1[2]], all_bp$relerr[,i2], sens_bp$relerr[,i1[1]], sens_bp$relerr[,i3[2]], all_bp$relerr[,i4], sens_bp$relerr[,i3[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
axis(2,at=seq(-2,3,by=1), las=2,cex.axis=3.5)
axis(1, at=1:6, labels=c("-25%","true","+25%","-25%","true","+25%"), cex.axis=3.2, col=c(rep("steelblue",3), rep("tomato",3)))
print.letter("(k)", cex=3)

i1 <- which(grepl("Long", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("equil/", all_sens_dir) & grepl("/linf/", all_sens_dir))
	all_sens_dir[i1]
i2 <- which(grepl("Long", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("equil/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i2]

i3 <- which(grepl("Long", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("Ramp/", all_sens_dir) & grepl("/linf/", all_sens_dir))
	all_sens_dir[i3]
i4 <- which(grepl("Long", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("Ramp/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i4]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,6.5),ylim=c(-3,3))
abline(h=0, lty=2)
pdf <- cbind.data.frame(sens_bp$relerr[,i1[2]], all_bp$relerr[,i2], sens_bp$relerr[,i1[1]], sens_bp$relerr[,i3[2]], all_bp$relerr[,i4], sens_bp$relerr[,i3[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
print.letter("(l)", cex=3)

axis(1, at=1:6, labels=c("-25%","true","+25%","-25%","true","+25%"), cex.axis=3.2, col=c(rep("steelblue",3), rep("tomato",3)))

i1 <- which(grepl("Long", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("equil/", all_sens_dir) & grepl("/vbk/", all_sens_dir))
	all_sens_dir[i1]
i2 <- which(grepl("Long", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("equil/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i2]

i3 <- which(grepl("Long", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("Ramp/", all_sens_dir) & grepl("/vbk/", all_sens_dir))
	all_sens_dir[i3]
i4 <- which(grepl("Long", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("Ramp/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i4]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,6.5),ylim=c(-3,3))
abline(h=0, lty=2)
pdf <- cbind.data.frame(sens_bp$relerr[,i1[2]], all_bp$relerr[,i2], sens_bp$relerr[,i1[1]], sens_bp$relerr[,i3[2]], all_bp$relerr[,i4], sens_bp$relerr[,i3[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
print.letter("(m)", cex=3)

axis(1, at=1:6, labels=c("-25%","true","+25%","-25%","true","+25%"), cex.axis=3.2, col=c(rep("steelblue",3), rep("tomato",3)))


i1 <- which(grepl("Long", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("equil/", all_sens_dir) & grepl("/ML50/", all_sens_dir))
	all_sens_dir[i1]
i2 <- which(grepl("Long", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("equil/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i2]

i3 <- which(grepl("Long", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("Ramp/", all_sens_dir) & grepl("/ML50/", all_sens_dir))
	all_sens_dir[i3]
i4 <- which(grepl("Long", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("Ramp/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i4]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,6.5),ylim=c(-3,3))
abline(h=0, lty=2)
pdf <- cbind.data.frame(sens_bp$relerr[,i1[2]], all_bp$relerr[,i2], sens_bp$relerr[,i1[1]], sens_bp$relerr[,i3[2]], all_bp$relerr[,i4], sens_bp$relerr[,i3[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
print.letter("(n)", cex=3)

axis(1, at=1:6, labels=c("-25%","true","+25%","-25%","true","+25%"), cex.axis=3.2, col=c(rep("steelblue",3), rep("tomato",3)))

i1 <- which(grepl("Long", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("equil/", all_sens_dir) & grepl("/CVlen/", all_sens_dir))
	all_sens_dir[i1]
i2 <- which(grepl("Long", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("equil/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i2]

i3 <- which(grepl("Long", all_sens_dir) & grepl("SampleSize_200/", all_sens_dir) & grepl("Ramp/", all_sens_dir) & grepl("/CVlen/", all_sens_dir))
	all_sens_dir[i3]
i4 <- which(grepl("Long", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("Ramp/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i4]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,6.5),ylim=c(-3,3))
abline(h=0, lty=2)
pdf <- cbind.data.frame(sens_bp$relerr[,i1[2]], all_bp$relerr[,i2], sens_bp$relerr[,i1[1]], sens_bp$relerr[,i3[2]], all_bp$relerr[,i4], sens_bp$relerr[,i3[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
print.letter("(o)", cex=3)

axis(1, at=1:6, labels=c("-25%","true","+25%","-25%","true","+25%"), cex.axis=3.2, col=c(rep("steelblue",3), rep("tomato",3)))

mtext(side=4, "Longer-lived", line=2, cex=3)
mtext(side=2, outer=TRUE, line=6, "Relative error", cex=3)
mtext(side=1, outer=TRUE, line=6, "Fixed value", cex=3)
dev.off()



############################################################
#### sensitivities -  dome
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
	## equilibrium
	lh_list <- adj_variation(SigmaR=0.01, SigmaF=0.01, SigmaC=0.01, SigmaI=0.01, CVlen=0.1, rho=0, selex_param=1, mat_param=1)

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



all_dome_dirs <- c(dome_equil_dir_vec, dome_base_dir_vec)
# dome_bp <- bias_precision(dirs=all_dome_dirs, itervec=itervec)
# saveRDS(dome_bp, file.path(res_dir, "dome_bias_precision.rds"))
dome_bp <- readRDS(file.path(res_dir, "dome_bias_precision.rds"))

png(file.path(fig_dir, "DomeSensitivity.png"), height=16, width=30, res=200, units="in")
mat <- matrix(c(1,2,3,
				1,2,3,
				1,2,3,
				4,5,6,
				7,8,9,
				7,8,9,
				7,8,9,
				10,11,12,
				10,11,12,
				10,11,12,
				13,14,15,
				13,14,15,
				13,14,15), nrow=13, ncol=3, byrow=TRUE)
nf <- layout(mat)
# layout.show(nf)
col_vec <- c("black", "steelblue", "tomato")
par(mar=c(0,6,0,0), omi=c(1.2,1.2,1,1), mgp=c(3,2,0))
plot(x=lh_list$Short$highs, y=lh_list$Short$S_l, lwd=3, ylim=c(0,1), cex.axis=2.5, yaxt="n", type="o", pch=19, cex=3, xaxs="i", yaxs="i", xpd=NA, xlab="", ylab="")
axis(2, cex.axis=2.5, las=2)
lines(x=lh_list$Short$highs, y=lh_dome_low$Short$S_l, lwd=3, col="steelblue", type="o", pch=19, cex=3, lty=2)
lines(x=lh_list$Short$highs, y=lh_dome_high$Short$S_l, lwd=3, col="tomato", type="o", pch=19, cex=3, lty=3)
mtext(side=2, "Selectivity", cex=3, line=6)
mtext(side=3, "Short-lived", cex=3, line=2)
print.letter("(a)", cex=3)

plot(x=lh_list$Medium$highs, y=lh_list$Medium$S_l, lwd=3, ylim=c(0,1), cex.axis=2.5, yaxt="n", type="o", pch=19, cex=3, xaxs="i", yaxs="i", xpd=NA, xlab="", ylab="")
axis(2, cex.axis=2.5, las=2)
lines(x=lh_list$Medium$highs, y=lh_dome_low$Medium$S_l, lwd=3, col="steelblue", type="o", pch=19, cex=3, lty=2)
lines(x=lh_list$Medium$highs, y=lh_dome_high$Medium$S_l, lwd=3, col="tomato", type="o", pch=19, cex=3, lty=3)
mtext(side=1, "Age", cex=3, line=5.5)
mtext(side=3, "Medium-lived", cex=3, line=2)
print.letter("(b)", cex=3)

plot(x=lh_list$Long$highs, y=lh_list$Long$S_l, lwd=3, ylim=c(0,1), cex.axis=2.5, yaxt="n", type="o", pch=19, cex=3, xaxs="i", yaxs="i", xpd=NA, xlab="", ylab="")
axis(2, cex.axis=2.5, las=2)
lines(x=lh_list$Long$highs, y=lh_dome_low$Long$S_l, lwd=3, col="steelblue", type="o", pch=19, cex=3, lty=2)
lines(x=lh_list$Long$highs, y=lh_dome_high$Long$S_l, lwd=3, col="tomato", type="o", pch=19, cex=3, lty=3)
mtext(side=3, "Longer-lived", cex=3, line=2)
print.letter("(c)", cex=3)

plot(x=1,y=1,type="n",axes=F,ann=F)
plot(x=1,y=1,type="n",axes=F,ann=F)
plot(x=1,y=1,type="n",axes=F,ann=F)

i1 <- which(grepl("Short", all_dome_dirs) & grepl("SampleSize_200/", all_dome_dirs) & grepl("equil/", all_dome_dirs))
	all_dome_dirs[i1]
i2 <- which(grepl("Short", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("/equil/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i2]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,3.5),ylim=c(-1,1))
abline(h=0, lty=2)
pdf <- cbind.data.frame(all_bp$relerr[,i2], dome_bp$relerr[,i1[2]], dome_bp$relerr[,i1[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
axis(2,at=seq(-0.5,1,by=0.5), las=2,cex.axis=2.5)
print.letter("(d)", cex=3)


i1 <- which(grepl("Medium", all_dome_dirs) & grepl("SampleSize_200/", all_dome_dirs) & grepl("equil/", all_dome_dirs))
	all_dome_dirs[i1]
i2 <- which(grepl("Medium", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("/equil/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i2]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,3.5),ylim=c(-2,2))
abline(h=0, lty=2)
pdf <- cbind.data.frame(all_bp$relerr[,i2], dome_bp$relerr[,i1[2]], dome_bp$relerr[,i1[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
axis(2,at=seq(-1,2,by=1), las=2,cex.axis=2.5)
print.letter("(e)", cex=3)


i1 <- which(grepl("Long", all_dome_dirs) & grepl("SampleSize_200/", all_dome_dirs) & grepl("equil", all_dome_dirs))
	all_dome_dirs[i1]
i2 <- which(grepl("Long", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("/equil/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i2]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,3.5),ylim=c(-3,3))
abline(h=0, lty=2)
pdf <- cbind.data.frame(all_bp$relerr[,i2], dome_bp$relerr[,i1[2]], dome_bp$relerr[,i1[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
axis(2,at=seq(-2,3,by=1), las=2,cex.axis=2.5)
mtext(side=4, "Equilibrium", line=2, cex=3)
print.letter("(f)", cex=3)

## base
i1 <- which(grepl("Short", all_dome_dirs) & grepl("SampleSize_200/", all_dome_dirs) & grepl("Ramp/", all_dome_dirs))
	all_dome_dirs[i1]
i2 <- which(grepl("Short", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("Ramp/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i2]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,3.5),ylim=c(-1,1))
abline(h=0, lty=2)
pdf <- cbind.data.frame(all_bp$relerr[,i2], dome_bp$relerr[,i1[2]], dome_bp$relerr[,i1[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
axis(2,at=seq(-0.5,1,by=0.5), las=2,cex.axis=2.5)
mtext(side=2, "Estimation error", cex=3, line=6)
print.letter("(g)", cex=3)

i1 <- which(grepl("Medium", all_dome_dirs) & grepl("SampleSize_200/", all_dome_dirs) & grepl("Ramp/", all_dome_dirs))
	all_dome_dirs[i1]
i2 <- which(grepl("Medium", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("Ramp/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i2]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,3.5),ylim=c(-2,2))
abline(h=0, lty=2)
pdf <- cbind.data.frame(all_bp$relerr[,i2], dome_bp$relerr[,i1[2]], dome_bp$relerr[,i1[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
axis(2,at=seq(-1,2,by=1), las=2,cex.axis=2.5)
print.letter("(h)", cex=3)


i1 <- which(grepl("Long", all_dome_dirs) & grepl("SampleSize_200/", all_dome_dirs) & grepl("Ramp/", all_dome_dirs))
	all_dome_dirs[i1]
i2 <- which(grepl("Long", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("Ramp/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i2]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,3.5),ylim=c(-3,3))
abline(h=0, lty=2)
pdf <- cbind.data.frame(all_bp$relerr[,i2], dome_bp$relerr[,i1[2]], dome_bp$relerr[,i1[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
axis(2,at=seq(-2,3,by=1), las=2,cex.axis=2.5)
print.letter("(i)", cex=3)
mtext(side=4, "Two-way", line=2, cex=3)

## base2
i1 <- which(grepl("Short", all_dome_dirs) & grepl("SampleSize_200/", all_dome_dirs) & grepl("Increasing/", all_dome_dirs))
	all_dome_dirs[i1]
i2 <- which(grepl("Short", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("Increasing/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i2]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,3.5),ylim=c(-1,1))
abline(h=0, lty=2)
pdf <- cbind.data.frame(all_bp$relerr[,i2], dome_bp$relerr[,i1[2]], dome_bp$relerr[,i1[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
axis(2,at=seq(-0.5,1,by=0.5), las=2,cex.axis=2.5)
axis(1,at=1:3, label=c("logistic", "low dome", "high dome"), cex.axis=3)
print.letter("(j)", cex=3)

i1 <- which(grepl("Medium", all_dome_dirs) & grepl("SampleSize_200/", all_dome_dirs) & grepl("Increasing/", all_dome_dirs))
	all_dome_dirs[i1]
i2 <- which(grepl("Medium", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("Increasing/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i2]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,3.5),ylim=c(-2,2))
abline(h=0, lty=2)
pdf <- cbind.data.frame(all_bp$relerr[,i2], dome_bp$relerr[,i1[2]], dome_bp$relerr[,i1[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
axis(2,at=seq(-1,2,by=1), las=2,cex.axis=2.5)
axis(1,at=1:3, label=c("logistic", "low dome", "high dome"), cex.axis=3)
print.letter("(k)", cex=3)


i1 <- which(grepl("Long", all_dome_dirs) & grepl("SampleSize_200/", all_dome_dirs) & grepl("Increasing/", all_dome_dirs))
	all_dome_dirs[i1]
i2 <- which(grepl("Long", all_dirs) & grepl("SampleSize_200/", all_dirs) & grepl("Increasing/", all_dirs)& grepl("/LC10/", all_dirs))
	all_dirs[i2]

plot(x=1,y=1,type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0.5,3.5),ylim=c(-3,3))
abline(h=0, lty=2)
pdf <- cbind.data.frame(all_bp$relerr[,i2], dome_bp$relerr[,i1[2]], dome_bp$relerr[,i1[1]])
beanplot(pdf, col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", xaxs="i", yaxs="i", lwd=3, na.rm=TRUE, what=c(0,1,1,0), beanlines="median", beanlinewd=3, add=TRUE)
allna <- sapply(1:ncol(pdf), function(x) all(is.na(pdf[,x])))
if(any(allna==TRUE)) points(x=which(allna==TRUE), y=rep(0,length(which(allna==TRUE))), pch=4, lwd=4, cex=5)
axis(2,at=seq(-2,3,by=1), las=2,cex.axis=2.5)
print.letter("(l)", cex=3)
axis(1,at=1:3, label=c("logistic", "low dome", "high dome"), cex.axis=3)
mtext(side=4, "One-way", line=2, cex=3)
mtext(side=1, "Selectivity type", line=7, cex=3, outer=TRUE)
dev.off()

############################################################
#### plot scenarios
############################################################
## all together
all_plot_vec <- c(equil_dir_vec[grepl("SampleSize_200/",equil_dir_vec) & grepl("Index_Catch_LC20",equil_dir_vec)], base_dir_vec[grepl("SampleSize_200/",base_dir_vec) & grepl("Index_Catch_LC20",base_dir_vec)])

cols <- brewer.pal(3, "Set1")

png(file.path(fig_dir, "FR_Scenarios.png"), height=11, width=20, res=200, units="in")
par(omi=c(1,1.5,1.5,1), mar=c(0,0,0,0))
mat <- matrix(c(1,1,1,2,2,2,3,3,3,4,5,5,5,6,6,6,7,7,7,8,9,9,9,10,10,10,11,11,11,
				12,12,12,13,13,13,14,14,14,15,16,16,16,17,17,17,18,18,18,19,20,20,20,21,21,21,22,22,22,
				23,23,23,24,24,24,25,25,25,26,27,27,27,28,28,28,29,29,29,30,31,31,31,32,32,32,33,33,33), ncol=29, nrow=3, byrow=TRUE)
tt <- layout(mat)
# layout.show(tt)
## fishing mortality
plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,4.5),xaxs="i",yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")
print.letter("(a)", cex=2)
axis(2, at=seq(1,4,by=1), las=2, cex.axis=2)
mtext("Fishing mortality", side=2, line=4, cex=2)
mtext("Equilibrium", side=3, line=1, cex=2)
pp <- sapply(1:length(itervec), function(x){
	true <- readRDS(file.path(all_plot_vec[1], itervec[x], "True.rds"))
	return(true$F_t)
})
up <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.90))
low <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.10))
polygon(x=c(1:length(up), length(up):1), y=c(low, rev(up)), col=paste0(cols[1],50), border=NA)
set.seed(208)
iters <- sample(itervec,3)
for(i in iters){
	true <- readRDS(file.path(all_plot_vec[1], i, "True.rds"))
	lines(true$F_t, col=cols[1], lwd=2)
}

plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,4.5),xaxs="i",yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")
print.letter("(b)", cex=2)
mtext("Two-way", side=3, line=1, cex=2)
mtext("Short-lived", side=3, line=4, cex=3.5, col=cols[1])
pp <- sapply(1:length(itervec), function(x){
	true <- readRDS(file.path(all_plot_vec[4], itervec[x], "True.rds"))
	return(true$F_t)
})
up <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.90))
low <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.10))
polygon(x=c(1:length(up), length(up):1), y=c(low, rev(up)), col=paste0(cols[1],50), border=NA)
for(i in iters){
	true <- readRDS(file.path(all_plot_vec[4], i, "True.rds"))
	lines(true$F_t, col=cols[1], lwd=2)
}

plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,4.5),xaxs="i",yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")
print.letter("(c)", cex=2)
mtext("One-way", side=3, line=1, cex=2)
pp <- sapply(1:length(itervec), function(x){
	true <- readRDS(file.path(all_plot_vec[7], itervec[x], "True.rds"))
	return(true$F_t)
})
up <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.90))
low <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.10))
polygon(x=c(1:length(up), length(up):1), y=c(low, rev(up)), col=paste0(cols[1],50), border=NA)
for(i in iters){
	true <- readRDS(file.path(all_plot_vec[7], i, "True.rds"))
	lines(true$F_t, col=cols[1], lwd=2)
}

plot(x=1,y=1,type="n",axes=F,ann=F)
plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,2),xaxs="i",yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")
print.letter("(d)", cex=2)
mtext("Equilibrium", side=3, line=1, cex=2)
axis(2, at=seq(0.5,2,by=0.5), las=2, cex.axis=2)
pp <- sapply(1:length(itervec), function(x){
	true <- readRDS(file.path(all_plot_vec[2], itervec[x], "True.rds"))
	return(true$F_t)
})
up <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.90))
low <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.10))
polygon(x=c(1:length(up), length(up):1), y=c(low, rev(up)), col=paste0(cols[2],50), border=NA)
for(i in iters){
	true <- readRDS(file.path(all_plot_vec[2], i, "True.rds"))
	lines(true$F_t, col=cols[2], lwd=2)
}

plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,2),xaxs="i",yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")
print.letter("(e)", cex=2)
mtext("Two-way", side=3, line=1, cex=2)
mtext("Medium-lived", side=3, line=4, cex=3.5, col=cols[2])
pp <- sapply(1:length(itervec), function(x){
	true <- readRDS(file.path(all_plot_vec[5], itervec[x], "True.rds"))
	return(true$F_t)
})
up <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.90))
low <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.10))
polygon(x=c(1:length(up), length(up):1), y=c(low, rev(up)), col=paste0(cols[2],50), border=NA)
for(i in iters){
	true <- readRDS(file.path(all_plot_vec[5], i, "True.rds"))
	lines(true$F_t, col=cols[2], lwd=2)
}

plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,2),xaxs="i",yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")
print.letter("(f)", cex=2)
mtext("One-way", side=3, line=1, cex=2)
pp <- sapply(1:length(itervec), function(x){
	true <- readRDS(file.path(all_plot_vec[8], itervec[x], "True.rds"))
	return(true$F_t)
})
up <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.90))
low <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.10))
polygon(x=c(1:length(up), length(up):1), y=c(low, rev(up)), col=paste0(cols[2],50), border=NA)
for(i in iters){
	true <- readRDS(file.path(all_plot_vec[8], i, "True.rds"))
	lines(true$F_t, col=cols[2], lwd=2)
}


plot(x=1,y=1,type="n",axes=F,ann=F)
plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,0.3),xaxs="i",yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")
print.letter("(g)", cex=2)
mtext("Equilibrium", side=3, line=1, cex=2)
axis(2, at=seq(0.1,0.3,by=0.1), las=2, cex.axis=2)
pp <- sapply(1:length(itervec), function(x){
	true <- readRDS(file.path(all_plot_vec[3], itervec[x], "True.rds"))
	return(true$F_t)
})
up <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.90))
low <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.10))
polygon(x=c(1:length(up), length(up):1), y=c(low, rev(up)), col=paste0(cols[3],50), border=NA)
for(i in iters){
	true <- readRDS(file.path(all_plot_vec[3], i, "True.rds"))
	lines(true$F_t, col=cols[3], lwd=2)
}


plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,0.5),xaxs="i",yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")
print.letter("(h)", cex=2)
mtext("Two-way", side=3, line=1, cex=2)
mtext("Longer-lived", side=3, line=4, cex=3.5, col=cols[3])
pp <- sapply(1:length(itervec), function(x){
	true <- readRDS(file.path(all_plot_vec[6], itervec[x], "True.rds"))
	return(true$F_t)
})
up <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.90))
low <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.10))
polygon(x=c(1:length(up), length(up):1), y=c(low, rev(up)), col=paste0(cols[3],50), border=NA)
for(i in iters){
	true <- readRDS(file.path(all_plot_vec[6], i, "True.rds"))
	lines(true$F_t, col=cols[3], lwd=2)
}

plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,0.5),xaxs="i",yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")
print.letter("(i)", cex=2)
mtext("One-way", side=3, line=1, cex=2)
pp <- sapply(1:length(itervec), function(x){
	true <- readRDS(file.path(all_plot_vec[9], itervec[x], "True.rds"))
	return(true$F_t)
})
up <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.90))
low <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.10))
polygon(x=c(1:length(up), length(up):1), y=c(low, rev(up)), col=paste0(cols[3],50), border=NA)
for(i in iters){
	true <- readRDS(file.path(all_plot_vec[9], i, "True.rds"))
	lines(true$F_t, col=cols[3], lwd=2)
}

## recruitment
plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,3.5),xaxs="i",yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")
print.letter("(j)", cex=2)
axis(2, at=seq(1,3.5,by=1), las=2, cex.axis=2)
mtext("Recruitment", side=2, line=4, cex=2)
pp <- sapply(1:length(itervec), function(x){
	true <- readRDS(file.path(all_plot_vec[1], itervec[x], "True.rds"))
	return(true$R_t)
})
up <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.90))
low <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.10))
polygon(x=c(1:length(up), length(up):1), y=c(low, rev(up)), col=paste0(cols[1],50), border=NA)
for(i in iters){
	true <- readRDS(file.path(all_plot_vec[1], i, "True.rds"))
	lines(true$R_t, col=cols[1], lwd=2)
}

plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,3.5),xaxs="i",yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")
print.letter("(k)", cex=2)
pp <- sapply(1:length(itervec), function(x){
	true <- readRDS(file.path(all_plot_vec[4], itervec[x], "True.rds"))
	return(true$R_t)
})
up <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.90))
low <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.10))
polygon(x=c(1:length(up), length(up):1), y=c(low, rev(up)), col=paste0(cols[1],50), border=NA)
for(i in iters){
	true <- readRDS(file.path(all_plot_vec[4], i, "True.rds"))
	lines(true$R_t, col=cols[1], lwd=2)
}

plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,3.5),xaxs="i",yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")
print.letter("(l)", cex=2)
pp <- sapply(1:length(itervec), function(x){
	true <- readRDS(file.path(all_plot_vec[4], itervec[x], "True.rds"))
	return(true$R_t)
})
up <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.90))
low <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.10))
polygon(x=c(1:length(up), length(up):1), y=c(low, rev(up)), col=paste0(cols[1],50), border=NA)
for(i in iters){
	true <- readRDS(file.path(all_plot_vec[4], i, "True.rds"))
	lines(true$R_t, col=cols[1], lwd=2)
}


plot(x=1,y=1,type="n",axes=F,ann=F)
plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,3.5),xaxs="i",yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")
print.letter("(m)", cex=2)
axis(2, at=seq(1,3.5,by=1), las=2, cex.axis=2)
pp <- sapply(1:length(itervec), function(x){
	true <- readRDS(file.path(all_plot_vec[2], itervec[x], "True.rds"))
	return(true$R_t)
})
up <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.90))
low <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.10))
polygon(x=c(1:length(up), length(up):1), y=c(low, rev(up)), col=paste0(cols[2],50), border=NA)
for(i in iters){
	true <- readRDS(file.path(all_plot_vec[2], i, "True.rds"))
	lines(true$R_t, col=cols[2], lwd=2)
}

plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,3.5),xaxs="i",yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")
print.letter("(n)", cex=2)
pp <- sapply(1:length(itervec), function(x){
	true <- readRDS(file.path(all_plot_vec[5], itervec[x], "True.rds"))
	return(true$R_t)
})
up <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.90))
low <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.10))
polygon(x=c(1:length(up), length(up):1), y=c(low, rev(up)), col=paste0(cols[2],50), border=NA)
for(i in iters){
	true <- readRDS(file.path(all_plot_vec[5], i, "True.rds"))
	lines(true$R_t, col=cols[2], lwd=2)
}

plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,3.5),xaxs="i",yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")
print.letter("(o)", cex=2)
pp <- sapply(1:length(itervec), function(x){
	true <- readRDS(file.path(all_plot_vec[8], itervec[x], "True.rds"))
	return(true$R_t)
})
up <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.90))
low <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.10))
polygon(x=c(1:length(up), length(up):1), y=c(low, rev(up)), col=paste0(cols[2],50), border=NA)
for(i in iters){
	true <- readRDS(file.path(all_plot_vec[8], i, "True.rds"))
	lines(true$R_t, col=cols[2], lwd=2)
}



plot(x=1,y=1,type="n",axes=F,ann=F)
plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,3.5),xaxs="i",yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")
print.letter("(p)", cex=2)
axis(2, at=seq(1,3.5,by=1), las=2, cex.axis=2)
pp <- sapply(1:length(itervec), function(x){
	true <- readRDS(file.path(all_plot_vec[3], itervec[x], "True.rds"))
	return(true$R_t)
})
up <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.90))
low <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.10))
polygon(x=c(1:length(up), length(up):1), y=c(low, rev(up)), col=paste0(cols[3],50), border=NA)
for(i in iters){
	true <- readRDS(file.path(all_plot_vec[3], i, "True.rds"))
	lines(true$R_t, col=cols[3], lwd=2)
}


plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,3.5),xaxs="i",yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")
print.letter("(q)", cex=2)
pp <- sapply(1:length(itervec), function(x){
	true <- readRDS(file.path(all_plot_vec[6], itervec[x], "True.rds"))
	return(true$R_t)
})
up <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.90))
low <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.10))
polygon(x=c(1:length(up), length(up):1), y=c(low, rev(up)), col=paste0(cols[3],50), border=NA)
for(i in iters){
	true <- readRDS(file.path(all_plot_vec[6], i, "True.rds"))
	lines(true$R_t, col=cols[3], lwd=2)
}

plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,3.5),xaxs="i",yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")
print.letter("(r)", cex=2)
pp <- sapply(1:length(itervec), function(x){
	true <- readRDS(file.path(all_plot_vec[9], itervec[x], "True.rds"))
	return(true$R_t)
})
up <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.90))
low <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.10))
polygon(x=c(1:length(up), length(up):1), y=c(low, rev(up)), col=paste0(cols[3],50), border=NA)
for(i in iters){
	true <- readRDS(file.path(all_plot_vec[9], i, "True.rds"))
	lines(true$R_t, col=cols[3], lwd=2)
}


## Relative biomass
plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,4),xaxs="i",yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")
print.letter("(s)", cex=2)
axis(2, las=2, cex.axis=2)
axis(1, cex.axis=2)
mtext("Relative biomass", side=2, line=4, cex=2)
pp <- sapply(1:length(itervec), function(x){
	true <- readRDS(file.path(all_plot_vec[1], itervec[x], "True.rds"))
	return(true$D_t)
})
up <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.90))
low <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.10))
polygon(x=c(1:length(up), length(up):1), y=c(low, rev(up)), col=paste0(cols[1],50), border=NA)
for(i in iters){
	true <- readRDS(file.path(all_plot_vec[1], i, "True.rds"))
	lines(true$D_t, col=cols[1], lwd=2)
}

plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,4),xaxs="i",yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")
print.letter("(t)", cex=2)
axis(1, cex.axis=2)
pp <- sapply(1:length(itervec), function(x){
	true <- readRDS(file.path(all_plot_vec[4], itervec[x], "True.rds"))
	return(true$D_t)
})
up <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.90))
low <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.10))
polygon(x=c(1:length(up), length(up):1), y=c(low, rev(up)), col=paste0(cols[1],50), border=NA)
for(i in iters){
	true <- readRDS(file.path(all_plot_vec[4], i, "True.rds"))
	lines(true$D_t, col=cols[1], lwd=2)
}


plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,4),xaxs="i",yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")
print.letter("(u)", cex=2)
axis(1, cex.axis=2)
pp <- sapply(1:length(itervec), function(x){
	true <- readRDS(file.path(all_plot_vec[7], itervec[x], "True.rds"))
	return(true$D_t)
})
up <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.90))
low <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.10))
polygon(x=c(1:length(up), length(up):1), y=c(low, rev(up)), col=paste0(cols[1],50), border=NA)
for(i in iters){
	true <- readRDS(file.path(all_plot_vec[7], i, "True.rds"))
	lines(true$D_t, col=cols[1], lwd=2)
}


plot(x=1,y=1,type="n",axes=F,ann=F)
plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,2.5),xaxs="i",yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")
print.letter("(v)", cex=2)
axis(1, cex.axis=2)
axis(2, las=2, cex.axis=2)
pp <- sapply(1:length(itervec), function(x){
	true <- readRDS(file.path(all_plot_vec[2], itervec[x], "True.rds"))
	return(true$D_t)
})
up <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.90))
low <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.10))
polygon(x=c(1:length(up), length(up):1), y=c(low, rev(up)), col=paste0(cols[2],50), border=NA)
for(i in iters){
	true <- readRDS(file.path(all_plot_vec[2], i, "True.rds"))
	lines(true$D_t, col=cols[2], lwd=2)
}

plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,2.5),xaxs="i",yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")
print.letter("(w)", cex=2)
axis(1, cex.axis=2)
pp <- sapply(1:length(itervec), function(x){
	true <- readRDS(file.path(all_plot_vec[5], itervec[x], "True.rds"))
	return(true$D_t)
})
up <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.90))
low <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.10))
polygon(x=c(1:length(up), length(up):1), y=c(low, rev(up)), col=paste0(cols[2],50), border=NA)
for(i in iters){
	true <- readRDS(file.path(all_plot_vec[5], i, "True.rds"))
	lines(true$D_t, col=cols[2], lwd=2)
}


plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,2.5),xaxs="i",yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")
print.letter("(x)", cex=2)
axis(1, cex.axis=2)
pp <- sapply(1:length(itervec), function(x){
	true <- readRDS(file.path(all_plot_vec[8], itervec[x], "True.rds"))
	return(true$D_t)
})
up <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.90))
low <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.10))
polygon(x=c(1:length(up), length(up):1), y=c(low, rev(up)), col=paste0(cols[2],50), border=NA)
for(i in iters){
	true <- readRDS(file.path(all_plot_vec[8], i, "True.rds"))
	lines(true$D_t, col=cols[2], lwd=2)
}

plot(x=1,y=1,type="n",axes=F,ann=F)
plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,2),xaxs="i",yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")
print.letter("(y)", cex=2)
axis(1, cex.axis=2)
axis(2, las=2, cex.axis=2)
pp <- sapply(1:length(itervec), function(x){
	true <- readRDS(file.path(all_plot_vec[3], itervec[x], "True.rds"))
	return(true$D_t)
})
up <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.90))
low <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.10))
polygon(x=c(1:length(up), length(up):1), y=c(low, rev(up)), col=paste0(cols[3],50), border=NA)
for(i in iters){
	true <- readRDS(file.path(all_plot_vec[3], i, "True.rds"))
	lines(true$D_t, col=cols[3], lwd=2)
}

plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,2),xaxs="i",yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")
print.letter("(z)", cex=2)
axis(1, cex.axis=2)
pp <- sapply(1:length(itervec), function(x){
	true <- readRDS(file.path(all_plot_vec[6], itervec[x], "True.rds"))
	return(true$D_t)
})
up <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.90))
low <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.10))
polygon(x=c(1:length(up), length(up):1), y=c(low, rev(up)), col=paste0(cols[3],50), border=NA)
for(i in iters){
	true <- readRDS(file.path(all_plot_vec[6], i, "True.rds"))
	lines(true$D_t, col=cols[3], lwd=2)
}

plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,2),xaxs="i",yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")
print.letter("(aa)", cex=2)
axis(1, cex.axis=2)
pp <- sapply(1:length(itervec), function(x){
	true <- readRDS(file.path(all_plot_vec[9], itervec[x], "True.rds"))
	return(true$D_t)
})
up <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.90))
low <- sapply(1:nrow(pp), function(x) quantile(pp[x,], prob=0.10))
polygon(x=c(1:length(up), length(up):1), y=c(low, rev(up)), col=paste0(cols[3],50), border=NA)
for(i in iters){
	true <- readRDS(file.path(all_plot_vec[9], i, "True.rds"))
	lines(true$D_t, col=cols[3], lwd=2)
}
mtext(side=1, outer=TRUE, "Year", cex=2, line=4)
dev.off()


	### plot equil truths
	png(file.path(fig_dir, "Equil_scenarios.png"), height=12, width=25, res=200, units="in")
	plot_vec <- equil_dir_vec[grepl("SampleSize_1000/",equil_dir_vec) & grepl("Index_Catch_LC20",equil_dir_vec)]
	par(mfrow=c(3,3), mar=c(0,4,0,0), omi=c(1,1,0.5,0.5))
	##Fishing
	plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,20), xaxs="i", yaxs="i", ylab="", xaxt="n", yaxt="n")
	axis(2,las=2,cex.axis=1.5)
	mtext("Fishing\nmortality", side=2, line=3, cex=2)
	for(iter in itervec){
		true <- readRDS(file.path(plot_vec[1], iter, "True.rds"))
		lines(true$F_t, col="#AA000075")
	}
	plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,10), xaxs="i", yaxs="i", ylab="", xaxt="n", yaxt="n")
		axis(2,las=2,cex.axis=1.5)
	for(iter in itervec){
		true <- readRDS(file.path(plot_vec[2], iter, "True.rds"))
		lines(true$F_t, col="#AA000075")
	}
	plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,2), xaxs="i", yaxs="i", ylab="", xaxt="n", yaxt="n")
		axis(2,las=2,cex.axis=1.5)
	for(iter in itervec){
		true <- readRDS(file.path(plot_vec[3], iter, "True.rds"))
		lines(true$F_t, col="#AA000075")
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
	dev.off()

	### plot base truths
	png(file.path(fig_dir, "Base_Ramp_scenarios.png"), height=12, width=25, res=200, units="in")
	plot_vec <- base_dir_vec[grepl("SampleSize_1000/",base_dir_vec) & grepl("Index_Catch_LC20",base_dir_vec)]
	par(mfrow=c(3,3), mar=c(0,4,0,0), omi=c(1,1,0.5,0.5))
	##Fishing
	plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,20), xaxs="i", yaxs="i", ylab="", xaxt="n", yaxt="n")
	axis(2,las=2,cex.axis=1.5)
	mtext("Fishing\nmortality", side=2, line=3, cex=2)
	for(iter in itervec){
		true <- readRDS(file.path(plot_vec[1], iter, "True.rds"))
		lines(true$F_t, col="#AA000075")
	}
	plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,10), xaxs="i", yaxs="i", ylab="", xaxt="n", yaxt="n")
		axis(2,las=2,cex.axis=1.5)
	for(iter in itervec){
		true <- readRDS(file.path(plot_vec[2], iter, "True.rds"))
		lines(true$F_t, col="#AA000075")
	}
	plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,2), xaxs="i", yaxs="i", ylab="", xaxt="n", yaxt="n")
		axis(2,las=2,cex.axis=1.5)
	for(iter in itervec){
		true <- readRDS(file.path(plot_vec[3], iter, "True.rds"))
		lines(true$F_t, col="#AA000075")
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
	dev.off()




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
		ignore <- run_LIME(modpath=dir, lh=lh_list[[as.character(modcombos[loop,"LH"])]], input_data=NULL, est_sigma=c("log_sigma_R"), data_avail=as.character(modcombos[loop,"Data_avail"]), itervec=itervec, rewrite=FALSE, simulation=TRUE, REML=FALSE, f_true=FALSE, write=TRUE, fix_param="beta", param_adjust="R0", val_adjust=R0_vec[i])
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
		abline(v=1, col=gray(0.3))
	plot(x=R0_vec, y=jnll_mat[,2], pch=19, col="steelblue", type="o")
		abline(v=1, col=gray(0.3))
	plot(x=R0_vec, y=jnll_mat[,3], pch=19, col="tomato", type="o")
		abline(v=1, col=gray(0.3))



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