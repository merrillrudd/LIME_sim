plot_re2 <- function(dirs, xlab, itervec, ...){
	esterr <- matrix(NA, nrow=length(itervec), ncol=length(dirs))
	for(i in 1:length(itervec)){
		rep1 <- readRDS(file.path(path1, itervec[i], "Report.rds"))
		rep2 <- readRDS(file.path(path2, itervec[i], "Report.rds"))
		rep3 <- readRDS(file.path(path3, itervec[i], "Report.rds"))
		true <- readRDS(file.path(path2, itervec[i], "True.rds"))
		esterr[itervec[i],1] <- (log(rep1$SPR_t[length(rep1$SPR_t)])-log(true$SPR_t[length(true$SPR_t)]))
		esterr[itervec[i],2] <- (log(rep2$SPR_t[length(rep2$SPR_t)])-log(true$SPR_t[length(true$SPR_t)]))
		esterr[itervec[i],3] <- (log(rep3$SPR_t[length(rep3$SPR_t)])-log(true$SPR_t[length(true$SPR_t)]))
	}
	plot(x=1,y=1,type="n", xlim=c(0,length(dirs)+1), xaxs="i", yaxs="i", ...)
	axis(2, las=2, cex.axis=2)
	if(xaxt!="n") axis(1, at=1:length(dirs), labels=xlab, cex.axis=2)
	abline(h=0, col=gray(0.4))
	par(new=TRUE)
	boxplot(esterr, xaxs="i", yaxs="i", xlim=c(0,length(dirs)+1), col="tomato", ...)
}