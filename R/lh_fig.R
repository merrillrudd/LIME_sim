lh_fig <- function(lh, save){

if(save==TRUE) png(file.path(fig_dir, "Life_history_comparison.png"), height=5, width=10, res=200, units="in")
	if(length(lh)==4) par(mfrow=c(2,2), mar=c(0,0,0,0), omi=c(1,1,1,1))
	if(length(lh)==3) par(mfrow=c(1,3), mar=c(0,0,0,0), omi=c(1,1,1,1))
	if(length(lh)==5) par(mfrow=c(3,2), mar=c(0,0,0,0), omi=c(1,1,1,1))

	for(ll in 1:length(lh)){
		plot(lh[[ll]]$Mat_a, type="l", lwd=4, ylim=c(0,1), xaxs="i", yaxs="i", xlim=c(1, max(unlist(sapply(1:length(lh), function(xx) lh[[xx]]$AgeMax + 1)))), xpd=NA, xlab="", ylab="", xaxt="n", yaxt="n")
		lines(lh[[ll]]$S_a, lty=3, lwd=4, col="gray", xpd=NA)
		if(length(lh)==4){
			if(ll==1) axis(2, at=seq(0.2,1,by=0.2), cex.axis=1.2, las=2)
			if(ll==3) axis(2, at=seq(0,1,by=0.2), cex.axis=1.2, las=2)
		}
		if(length(lh)==3 & ll==1) axis(2, at=seq(0,1,by=0.2), cex.axis=1.2, las=2)

		if(length(lh)==4){
			if(ll==1) text(x=7, y=0.1, "Type I: Long-lived, medium size", font=2, cex=1.3)
			if(ll==2) text(x=7, y=0.1, "Type II: Long-lived, smaller size", font=2, cex=1.3)
			if(ll==3) text(x=7, y=0.1, "Type III: Medium-lived, large size", font=2, cex=1.3)
			if(ll==4) text(x=7, y=0.1, "Type IV: Short-lived, medium size", font=2, cex=1.3)
		}
		if(length(lh)==3){
			if(names(lh[[ll]])=="Medium") print.letter("Snapper", xy=c(0.75,0.88), font=2, cex=1.3)
			if(names(lh[[ll]]=="Short") print.letter("Rabbitfish", xy=c(0.75,0.88), font=2, cex=1.3)
			if(names(lh[[ll]]=="Long")) print.letter("Hake", xy=c(0.75,0.88), font=2, cex=1.3)
		}


		par(new=TRUE)
		plot(lh[[ll]]$L_a, lwd=4, type="l", col="blue", xpd=NA, xaxs="i", yaxs="i", ylim=c(0,max(unlist(sapply(1:length(lh), function(xx) lh[[xx]]$L_a)))), xlim=c(1, max(unlist(sapply(1:length(lh), function(xx) lh[[xx]]$AgeMax + 1)))), xlab="", ylab="", xaxt="n", yaxt="n")
		if(length(lh)==4){
			if(ll==2) axis(4, at=seq(10,50,by=10), las=2, cex.axis=1.2, col="blue", col.axis="blue")
			if(ll==4) axis(4, at=seq(0,50,by=10), las=2, cex.axis=1.2, col="blue", col.axis="blue")
			if(ll==3) axis(1, at=seq(1,max(unlist(sapply(1:length(lh), function(xx) lh[[xx]]$AgeMax + 1)))), cex.axis=1.2)
			if(ll==4) axis(1, at=seq(1,max(unlist(sapply(1:length(lh), function(xx) lh[[xx]]$AgeMax + 1)))), cex.axis=1.2)
		}
		if(length(lh)==3){
			if(ll==3) axis(4, at=seq(0,100,by=10), las=2, cex.axis=1.2, col="blue", col.axis="blue")
			axis(1, at=seq(0,max(unlist(sapply(1:length(lh), function(xx) lh[[xx]]$AgeMax))), by=4), cex.axis=1.2)
		}
		if(ll==length(lh)) legend("bottomright", legend=c("Selectivity", "Maturity", "Length"), col=c("black", "gray", "blue"), lty=c(1,3,1), lwd=4, cex=1.2)
		print.letter(paste0("(", letters[ll], ")"), xy=c(0.75,0.93), cex=1.3, font=2)
	}
	mtext(side=2, "Proportion vulnerable/mature", cex=1.5, line=3, outer=TRUE)
	mtext(side=4, "Length (cm)", cex=1.5, line=3, outer=TRUE, col="blue")
	mtext(side=1, "Age", cex=1.5, line=3, outer=TRUE)
if(save==TRUE) dev.off()

}