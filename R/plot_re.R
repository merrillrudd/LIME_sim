plot_re <- function(dirs, modcombos, itervec, ylim=c(-1,1), compareToLH, lh_vec, recalc=TRUE, dev=NA, converge=NA, cover=NA, text=TRUE){


if(recalc==TRUE){
	### ------- check output ------###
	dev <- esterr <- matrix(NA, nrow=length(itervec), ncol=length(dirs))
	converge <- matrix(1, nrow=length(itervec), ncol=length(dirs))
	for(m in 1:length(dirs)){
		for(i in 1:length(itervec)){
			true <- readRDS(file.path(dirs[m], itervec[i], "True.rds"))	

			if(grepl("LBSPR", dirs[m])==FALSE){
				rep <- readRDS(file.path(dirs[m], itervec[i], "Report.rds"))
				sdrep <- readRDS(file.path(dirs[m], itervec[i], "Sdreport.rds"))
				inp <- readRDS(file.path(dirs[m], itervec[i], "Inputs.rds"))	
				der <- readRDS(file.path(dirs[m], itervec[i], "Derived_quants.rds"))	

				if(all(is.na(sdrep))){
					converge[i,m] <- 0
					next
				}	

				s <- summary(sdrep)[which(rownames(summary(sdrep))=="SPR_t")[length(which(rownames(summary(sdrep))=="SPR_t"))],]
				if(is.na(s[2])){
					converge[i,m] <- 0
					next
				}	

				dev[i,m] <- rep$SPR_t[length(rep$SPR_t)] - true$SPR_t[length(true$SPR_t)]	
				esterr[i,m] <- log(rep$SPR_t[length(rep$SPR_t)]) - log(true$SPR_t[length(true$SPR_t)])
			}	

			if(grepl("LBSPR", dirs[m])){
				if(file.exists(file.path(dirs[m], itervec[i], "non_convergence.txt"))){
					converge[i,m] <- 0
					next
				}
				rep <- readRDS(file.path(dirs[m], itervec[i], "LBSPR_results.rds"))
				dev[i,m] <- (rep$SPR[length(rep$SPR)] - true$SPR_t[length(true$SPR_t)])
				esterr[i,m] <- log(rep$SPR[length(rep$SPR)]) - log(true$SPR_t[length(true$SPR_t)])
			}	

		}
	}
}

nf <- length(compareToLH)
nl <- length(lh_vec)
# col_vec <- brewer.pal(3, "Set1")
names <- gsub("_", "\n", modcombos[,"Data_avail"])
bias <- sapply(1:length(dirs), function(x) median(abs(dev[,x]), na.rm=TRUE))
precision <- sapply(1:length(dirs), function(x) sqrt((1/nrow(dev))*sum(dev[,x]^2, na.rm=TRUE)))
lime_dirs <- dirs[which(grepl("LBSPR",dirs)==FALSE)]


par(mfrow=c(nf,nl), mar=c(0,0,0,0), omi=c(0.75,0.75,1,1))
for(ff in 1:nf){
	for(ll in 1:nl){
		index <- which(grepl(compareToLH[ff],dirs) & grepl(lh_vec[ll],dirs))
		lime_index <- which(grepl(compareToLH[ff],lime_dirs) & grepl(lh_vec[ll],lime_dirs))
		if(length(index)==1) xmax <- 1.5
		if(length(index)>1) xmax <- ncol(esterr[,index]) + 0.5

		# beanplot(as.data.frame(esterr[,index]), col=lapply(1:length(col_vec), function(x) c(col_vec[x],"black","black","black")), xaxt="n", yaxt="n", lwd=3, na.rm=TRUE, xaxs="i", yaxs="i", what=c(0,1,1,0), add=TRUE, beanlines="median", beanlinewd=3)
		plot(x=1, y=1, type="n", xaxs="i", yaxs="i", xlab="", ylab="", xlim=c(0.5,xmax), ylim=ylim, xaxt="n", yaxt="n")
		abline(h=0, lwd=5, col="red")
		par(new=TRUE)
		boxplot(esterr[,index], ylim=ylim, xaxt="n", yaxt="n", col="goldenrod", xlim=c(0.5, xmax), xaxs="i", yaxs="i")
		if(text==TRUE){
			text(x=1:length(index), y=0.8*ylim[2], round(precision[index],3), cex=2, col="red")
			text(x=1:length(index), y=0.7*ylim[2], round(bias[index],3), cex=2, col="blue")
			if(all(is.na(cover))==FALSE) text(x=1:length(lime_index), y=0.6*ylim[2], colSums(cover[,lime_index])/nrow(cover), cex=2, col="forestgreen")
			text(x=1:length(index), y=0.5*ylim[2], round(converge[index],2), cex=2, col="gray")
		}
		if(ff==nf) axis(1, at=1:length(index), labels=names[index], cex.axis=1.4)
		if(ll==1) axis(2, at=pretty(ylim), cex=1.2, las=2, cex.axis=1.4)
		if(ff==1) mtext(side=3, lh_vec[ll], font=2, cex=1.3, line=1)
		if(ll==nl) mtext(side=4, compareToLH[ff], font=2, cex=1.3, line=1)
	}
}
mtext("Data availability scenario", side=1, line=4, outer=TRUE, cex=1.5)
mtext("Estimation error", side=2, line=3, outer=TRUE, cex=1.5)
mtext("Life history type", side=3, line=4, outer=TRUE, cex=1.5)
mtext("Population dynamic scenario", side=4, line=4, outer=TRUE, cex=1.5)

# legend("topright", legend=c("Short-lived", "Medium-lived", "Long-lived"), title="Life history type", col=col_vec, pch=15)

Outs <- NULL
Outs$bias <- bias
Outs$precision <- precision
Outs$esterr <- esterr
Outs$converge <- converge
return(Outs)
}