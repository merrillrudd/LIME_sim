interval_coverage <- function(dirs, itervec){

cover <- matrix(0, nrow=length(itervec), ncol=length(dirs))
converge <- matrix(1, nrow=length(itervec), ncol=length(dirs))
for(m in 1:length(dirs)){
	for(i in 1:length(itervec)){
			
			true <- readRDS(file.path(dirs[m], itervec[i], "True.rds"))
			if(file.exists(file.path(dirs[m], itervec[i], "Sdreport.rds"))) sdrep <- readRDS(file.path(dirs[m], itervec[i], "Sdreport.rds"))
			if(file.exists(file.path(dirs[m], itervec[i], "Sdreport.rds"))==FALSE){
				cover[i,m] <- NA
				converge[i,m] <- NA
				next
			} 


			if(all(is.na(sdrep))){
				converge[i,m] <- 0
				next
			}

			s <- summary(sdrep)[which(rownames(summary(sdrep))=="SPR_t")[length(which(rownames(summary(sdrep))=="SPR_t"))],]
			if(is.na(s[2])){
				converge[i,m] <- 0
				next
			}

			high50 <- s[1]+0.67449*s[2]
	        low50 <- s[1]-0.67449*s[2]
	        if(true$SPR_t[length(true$SPR_t)]< high50 & true$SPR_t[length(true$SPR_t)]>low50) cover[i,m] <- 1
	}
}

Outs <- NULL
Outs$cover <- cover
Outs$converge <- converge
return(Outs)

}