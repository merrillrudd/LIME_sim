bias_precision <- function(dirs, itervec, param="SPR"){

	dev <- relerr <- esterr <- matrix(NA, nrow=length(itervec), ncol=length(dirs))
	converge <- matrix(1, nrow=length(itervec), ncol=length(dirs))
	bound <- matrix(0, nrow=length(itervec), ncol=length(dirs))

	for(m in 1:length(dirs)){
		for(i in 1:length(itervec)){
				if(file.exists(file.path(dirs[m], itervec[i], "high_final_gradient.txt"))){
					converge[i,m] <- 0
					next
				}	

			true <- readRDS(file.path(dirs[m], itervec[i], "True.rds"))

			if(grepl("/LBSPR", dirs[m])==FALSE){
				rep <- readRDS(file.path(dirs[m], itervec[i], "Report.rds"))
				sdrep <- readRDS(file.path(dirs[m], itervec[i], "Sdreport.rds"))

				if(all(is.na(sdrep))){
					converge[i,m] <- 0
					next
				}

				if(rep$sigma_R==2){
					bound[i,m] <- 1
					next
				}

				if(rep$S50==log(0.001)){
					bound[i,m] <- 1
					next
				}

				if(param=="SPR"){
					dev[i,m] <- rep$SPR_t[length(rep$SPR_t)] - true$SPR_t[length(true$SPR_t)]	
					esterr[i,m] <- log(rep$SPR_t[length(rep$SPR_t)]) - log(true$SPR_t[length(true$SPR_t)])
					relerr[i,m] <- (rep$SPR_t[length(rep$SPR_t)] - true$SPR_t[length(true$SPR_t)])/true$SPR_t[length(true$SPR_t)]
				}

				if(param=="SigmaR"){
					dev[i,m] <- rep$sigma_R - true$SigmaR
					esterr[i,m] <- log(rep$sigma_R) - log(true$SigmaR)
					relerr[i,m] <- (rep$sigma_R - true$SigmaR)/true$SigmaR
				}

				if(param=="F"){
					dev[i,m] <- rep$F_t[length(rep$F_t)] - true$F_t[length(true$F_t)]
					esterr[i,m] <- log(rep$F_t[length(rep$F_t)]) - log(true$F_t[length(true$F_t)])
					relerr[i,m] <- (rep$F_t[length(rep$F_t)] - true$F_t[length(true$F_t)])/true$F_t[length(true$F_t)]
				}
				if(param=="S50"){
					dev[i,m] <- rep$S50 - true$SL50
					esterr[i,m] <- log(rep$S50) - log(true$SL50)
					relerr[i,m] <- (rep$S50 - true$SL50)/true$SL50
				}
			}
			if(grepl("/LBSPR", dirs[m])){
				if(file.exists(file.path(dirs[m], itervec[i], "non_convergence.txt"))){
					converge[i,m] <- 0
					next
				}
				rep <- readRDS(file.path(dirs[m], itervec[i], "LBSPR_results.rds"))
				
				if(param=="SPR"){
					dev[i,m] <- (rep$SPR[length(rep$SPR)] - true$SPR_t[length(true$SPR_t)])
					esterr[i,m] <- log(rep$SPR[length(rep$SPR)]) - log(true$SPR_t[length(true$SPR_t)])
					relerr[i,m] <- (rep$SPR[length(rep$SPR)] - true$SPR_t[length(true$SPR_t)])/true$SPR_t[length(true$SPR_t)]
				}
				if(param=="S50"){
					dev[i,m] <- rep$SL50[length(rep$SL50)] - true$SL50
					esterr[i,m] <- log(rep$SL50[length(rep$SL50)]) - log(true$SL50)
					relerr[i,m] <- (rep$SL50[length(rep$SL50)] - true$SL50)/true$SL50
				}
				if(param=="F"){
					dev[i,m] <- rep$FM[length(rep$FM)]/true$M - true$F_t[length(true$F_t)]/true$M
					esterr[i,m] <- log(rep$FM[length(rep$FM)]/true$M) - log(true$F_t[length(true$F_t)]/true$M)
					relerr[i,m] <- ((rep$FM[length(rep$FM)]/true$M) - (true$F_t[length(true$F_t)]/true$M))/(true$F_t[length(true$F_t)]/true$M)
				}
			}	
		}
	}

	bias <- sapply(1:length(dirs), function(x) median(abs(dev[,x]), na.rm=TRUE))
	precision <- sapply(1:length(dirs), function(x) sqrt((1/length(which(is.na(dev[,x])==FALSE)))*sum(dev[,x]^2, na.rm=TRUE)))


	Outs <- NULL
	Outs$dev <- dev
	Outs$esterr <- esterr 
	Outs$bias <- bias
	Outs$precision <- precision
	Outs$bound <- bound
	Outs$converge <- converge
	Outs$relerr <- relerr
	return(Outs)
}