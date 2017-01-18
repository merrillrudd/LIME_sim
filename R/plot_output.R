plot_output <- function(LIME_dir, LBSPR_dir=NULL, sim=TRUE, iter=1){
  
  if(sim==FALSE){
    inp <- readRDS(file.path(LIME_dir, "Inputs.rds"))
    rep <- readRDS(file.path(LIME_dir, "Report.rds"))
    sdrep <- readRDS(file.path(LIME_dir, "Sdreport.rds"))
    der <- readRDS(file.path(LIME_dir, "Derived_quants.rds"))
    if(is.null(LBSPR_dir)==FALSE) LBSPR_outs <- readRDS(file.path(LBSPR_dir, "LBSPR_results.rds"))
  }
  if(sim==TRUE){
    inp <- readRDS(file.path(LIME_dir, iter, "Inputs.rds"))
    rep <- readRDS(file.path(LIME_dir, iter, "Report.rds"))
    sdrep <- readRDS(file.path(LIME_dir, iter, "Sdreport.rds"))
    der <- readRDS(file.path(LIME_dir, iter, "Derived_quants.rds"))
    if(is.null(LBSPR_dir)==FALSE) LBSPR_outs <- readRDS(file.path(LBSPR_dir, iter, "LBSPR_results.rds"))
 
  }

if(sim==TRUE) true <- readRDS(file.path(LIME_dir, iter, "True.rds"))
if(sim==TRUE & is.null(LBSPR_dir)==FALSE){
	true2 <- readRDS(file.path(LBSPR_dir, iter, "True.rds"))
	if(identical(true$SPR,true2$SPR)==FALSE) break("True values for LIME & LBSPR are not the same")
}

Nyears_model <- length(rep$SPR_t)
Nyears_lf <- inp$Data$n_lc
Years_model <- 1:Nyears_model
Years_lf <- (Nyears_model - Nyears_lf + 1):Nyears_model

	FUN <- function(InputMat, log=TRUE, rel=FALSE){
        index <- which(is.na(InputMat[,2])==FALSE)
        if(log==TRUE) return(c( exp(InputMat[index,1]-1.96*InputMat[index,2]), rev(exp(InputMat[index,1]+1.96*InputMat[index,2]))))
        if(log==FALSE) return(c( InputMat[index,1]-1.96*InputMat[index,2], rev(InputMat[index,1]+1.96*InputMat[index,2])))
	}

      par(mfrow=c(2,2), mar=c(5,6,3,3))
      plot(rep$F_t, type="l", lwd=2, xlim=c(1,length(rep$F_t)), ylim=c(0, max(max(rep$F_t)*1.5,1)), col="blue", xaxs="i", yaxs="i", xlab="Year", ylab="", cex.axis=2, cex.lab=2, yaxt="n")
      axis(2, las=2, cex.axis=2)
      mtext(side=2, "Fishing mortality", cex=2, line=4)
      if(all(is.na(sdrep))==FALSE) polygon( y=FUN(summary(sdrep)[which(rownames(summary(sdrep))=="lF_t"),], log=TRUE), x=c(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="lF_t"),2])==FALSE), rev(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="lF_t"),2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
      if(is.null(LBSPR_dir)==FALSE) lines(x=which(Years_model %in% Years_lf), y=LBSPR_outs$FM*inp$Data$M, col="red", lwd=2)
  	  if(sim==TRUE) lines(true$F_t, lwd=2)
      
      plot(rep$R_t, type="l", lwd=2, xlim=c(1,length(rep$F_t)), ylim=c(0, max(max(rep$R_t)*1.5,2)), col="blue", xaxs="i", yaxs="i", xlab="Year", ylab="", cex.axis=2, cex.lab=2, yaxt="n")
      axis(2, las=2, cex.axis=2)
      mtext(side=2, "Recruitment", cex=2, line=4)
      if(all(is.na(sdrep))==FALSE) polygon( y=FUN(summary(sdrep)[which(rownames(summary(sdrep))=="lR_t"),], log=TRUE), x=c(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="lR_t"),2])==FALSE), rev(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="lR_t"),2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
      if(is.null(LBSPR_dir)==FALSE) lines(x=which(Years_model %in% Years_lf), y=rep(1,length(Years_lf)), col="red", lwd=2)
  	  if(sim==TRUE) lines(true$R_t, lwd=2)
      
      plot(x=1,y=1, type="n", xlim=c(1, length(rep$SPR_t)), ylim=c(0,1), xlab="Year", ylab="", xaxs="i", yaxs="i", cex.axis=2, cex.lab=2, yaxt="n")
      axis(2, las=2, cex.axis=2)
      mtext(side=2, "SPR", cex=2, line=4)
      polygon(x=c(0,Nyears_model, Nyears_model, 0), y=c(0.3,0.3,1,1), border=NA, col="#00AA0030")
      polygon(x=c(0,Nyears_model, Nyears_model, 0), y=c(0,0,0.3,0.3), border=NA, col="#AA000030")
      lines(x=rep$SPR_t, type="l", lwd=2, col="blue", xaxs="i", yaxs="i", xlab="Year", ylab="SPR")
      if(all(is.na(sdrep))==FALSE) polygon( y=FUN(summary(sdrep)[which(rownames(summary(sdrep))=="SPR_t"),], log=FALSE), x=c(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="SPR_t"),2])==FALSE), rev(which(is.na(summary(sdrep)[which(rownames(summary(sdrep))=="SPR_t"),2])==FALSE))), col=rgb(0,0,1,alpha=0.2), border=NA)
      if(is.null(LBSPR_dir)==FALSE) lines(x=which(Years_model %in% Years_lf), y=LBSPR_outs$SPR, col="red", lwd=2)
	  if(sim==TRUE) lines(true$SPR_t, lwd=2)


      plot(x=1, y=1, type="n", xlim=c(0,1), ylim=c(0,5), xaxs="i", yaxs="i", xlab="SPR", ylab="", cex.axis=2, cex.lab=2, yaxt="n")
      axis(2, las=2, cex.axis=2)
      mtext(side=2, "F/F30", cex=2, line=4)
      polygon(x=c(0.3,1,1,0.3), y=c(0,0,1,1), col="#00AA0030", border=NA)
      polygon(x=c(0,0.3,0.3,0), y=c(1,1,5,5), col="#AA000030", border=NA)
      abline(h=1, lty=2, lwd=3)
      abline(v=0.3, lty=2, lwd=3)
      points(x=rep$SPR_t[length(rep$SPR_t)], y=rep$F_t[length(rep$F_t)]/der$F30, col="blue", pch=19, cex=2)
      F30 <- tryCatch(with(rep, uniroot(calc_ref, lower=0, upper=50, Mat_a=Mat_a, W_a=W_a, M=M, S_a=S_a, ref=0.3)$root), error=function(e) NA)
      if(is.null(LBSPR_dir)==FALSE) if(is.na(F30)==FALSE) points(x=LBSPR_outs$SPR[length(LBSPR_outs$SPR)], y=(LBSPR_outs$FM[length(LBSPR_outs$FM)]*inp$Data$M)/F30, col="red", pch=19, cex=2)
      if(sim==TRUE) points(x=true$SPR_t[length(true$SPR_t)], y=true$FF30, col="black", pch=19, cex=2)

      if(iter==1){
        df <- read.csv(file.path(LIME_dir, "df.csv"), header=TRUE)
        return(df)
      }
      if(iter!=1){
        return(rep)
      }

}