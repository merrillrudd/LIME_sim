sim_LBSPR <- function(dir, itervec, lh, Nyears, F=0.1, modtype="absel"){

require(LBSPR)
for(iter in 1:length(itervec)){
	set.seed(iter)
	F <- runif(1,0,1)
  	LB_pars <- new("LB_pars", verbose=FALSE)
	  LB_pars@MK <- lh$M/lh$vbk 
	  LB_pars@Linf <- lh$linf
	  LB_pars@CVLinf <- lh$CVlen
	  LB_pars@L50 <- lh$ML50 
	  LB_pars@L95 <- lh$ML95
	  LB_pars@Walpha <- lh$lwa
	  LB_pars@Wbeta <- lh$lwb
	  LB_pars@SL50 <- lh$SL50 
	  LB_pars@SL95 <- lh$SL95
	  LB_pars@FM <- F
	  LB_pars@BinWidth <- lh$binwidth
	  LB_pars@BinMax <- 1.5 * lh$linf
	  LB_pars@BinMin <- 0 
	  LB_pars@R0 <- lh$R0
	  sim <- LBSPRsim(LB_pars, Control=list(modtype=modtype))


	  true <- lh
	  true$DataScenario <- "LBSPR_test"
	  true$LF <- t(sim@pLCatch)
	  	colnames(true$LF) <- sim@LMids
	  	rownames(true$LF) <- Nyears
	  true$SPR_t <- rep(sim@SPR, Nyears)
	  true$Nyears <- Nyears
	  true$years <- 1:Nyears

	  saveRDS(true, file.path(dir, iter, "True.rds"))
}

}