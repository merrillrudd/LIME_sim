plot_scenarios <- function(dirs, itervec){
	par(mfrow=c(3,3), mar=c(0,0,0,0), omi=c(1,1,1,1))
	for(dd in 1:length(dirs)){
		plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,10), xaxt="n", yaxt="n", xaxs="i", yaxs="i", ann=F)
		for(ii in 1:length(itervec)){
			true <- readRDS(file.path(dirs[dd],itervec[ii],"True.rds"))
			lines(true$F_t, col="#AA000050", lwd=2)
		}
		if(dd==1){
			axis(2, las=2, cex=1.5)
			mtext(side=2, "Fishing mortality", cex=1.5, line=4)
		}
		mtext(side=3, lh_vec[dd], cex=1.5, line=1.5)
	}
	for(dd in 1:length(dirs)){
		plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,4), xaxt="n", yaxt="n", xaxs="i", yaxs="i", ann=F)
		for(ii in 1:length(itervec)){
			true <- readRDS(file.path(dirs[dd],itervec[ii],"True.rds"))
			lines(true$R_t, col="#0000AA50", lwd=2)
		}
		if(dd==1){
			axis(2, las=2, cex=1.5)
			mtext(side=2, "Recruitment", cex=1.5, line=4)
		}
	}
	for(dd in 1:length(dirs)){
		plot(x=1,y=1,type="n",xlim=c(1,20),ylim=c(0,1.5), xaxt="n", yaxt="n", xaxs="i", yaxs="i", ann=F)
		for(ii in 1:length(itervec)){
			true <- readRDS(file.path(dirs[dd],itervec[ii],"True.rds"))
			lines(true$D_t, col="#00AA0050", lwd=2)
		}
		if(dd==1){
			axis(2, las=2, cex=1.5)
			mtext(side=2, "Relative biomass", cex=1.5, line=4)
		}
		axis(1, cex=1.5)
	}
	mtext(side=1, "Year", outer=TRUE, cex=1.5, line=3.5)

}