lh_fig <- function(lh){

	if(length(lh)==4) par(mfrow=c(2,2), mar=c(0,0,0,0), omi=c(1,1,1,1))
	if(length(lh)==3) par(mfrow=c(1,3), mar=c(0,0,0,0), omi=c(1,1,1,1))
	if(length(lh)==5) par(mfrow=c(3,2), mar=c(0,0,0,0), omi=c(1,1,1,1))

	for(ll in 1:length(lh)){
		xlim <- c(min(lh[[ll]]$ages), max(lh[[ll]]$ages))

		plot(x=lh[[ll]]$ages, y=lh[[ll]]$Mat_a, type="l", lwd=5, ylim=c(0,1), xaxs="i", yaxs="i", xlim=xlim, xpd=NA, xlab="", ylab="", xaxt="n", yaxt="n", xpd=NA)
		lines(x=lh[[ll]]$ages, y=lh[[ll]]$S_a, lty=3, lwd=5, col="gray", xpd=NA)
		if(ll==1) axis(2, at=seq(0,1,by=0.2), cex.axis=2, las=2)

		mtext(side=3, paste0("(", letters[ll],")"), font=2, cex=1.3, line=4)
		name <- names(lh)[ll]
		if(name=="Long") name <- "Longer"
		mtext(side=3, name, font=3, cex=1.4, line=1.5)


		par(new=TRUE)
		ylim <- c(0,max(unlist(sapply(1:length(lh), function(xx) lh[[xx]]$L_a))))
		plot(x=lh[[ll]]$ages, y=lh[[ll]]$L_a, lwd=5, type="l", col="blue", xpd=NA, xaxs="i", yaxs="i", ylim=ylim, xlim=xlim, xlab="", ylab="", xaxt="n", yaxt="n")
		if(ll==length(lh)) axis(4, at=pretty(ylim), las=2, cex.axis=2, col="blue", col.axis="blue")
		axis(1, at=pretty(xlim)[-length(pretty(xlim))], cex.axis=2)
	
		if(ll==length(lh)) legend("bottomright", legend=c("Maturity", "Selectivity", "Length"), col=c("black", "gray", "blue"), lty=c(1,3,1), lwd=4, cex=2)
	}
	mtext(side=2, "Proportion vulnerable/mature", cex=1.5, line=4, outer=TRUE)
	mtext(side=4, "Length (cm)", cex=1.5, line=4, outer=TRUE, col="blue")
	mtext(side=1, "Age", cex=1.5, line=4, outer=TRUE)

}