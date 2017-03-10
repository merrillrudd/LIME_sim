#' run LBSPR within LIME framework
#'
#' \code{run_LBSPR} sets up LBSPR to run within LIME framework (adjusting for inputs, structure, etc.)
#'
#' @param modpath path to model directory
#' @param lh life history information
#' @param input_data list of input data including years of model and length frequency data
#' @param itervec vector of iterates to run for simulation
#' @param species name of species (will set as generic within LBSPR package if not specified)
#' @param  rewrite should the assessment be re-written in the directory if it already exists?
#' @param simulation is this a simulation or application to real data?
#' @param write write results to directory?

#' @return list of LBSPR output
#' @details Must load LBSPR package: devtools::install_github("AdrianHordyk/LBSPR")

run_LBSPR <- function(modpath, lh, input_data=NULL, itervec=NULL, species=NULL, rewrite=TRUE, simulation=TRUE, write=TRUE){

    if(simulation==FALSE) itervec <- 1

    for(iter in 1:length(itervec)){

        if(simulation==TRUE) iterpath <- file.path(modpath, iter)
        if(simulation==FALSE) iterpath <- modpath

        if(write==TRUE & rewrite==FALSE) if(any(grepl("LBSPR_results", list.files(path=iterpath)))) next
      if(write==TRUE & rewrite==TRUE){
        if(any(grepl("non_convergence", list.files(path=iterpath)))) unlink(file.path(iterpath, "non_convergence.txt"), TRUE)
      }

    if(write==TRUE & simulation==TRUE){
      sim <- readRDS(file.path(iterpath, "True.rds"))
      input_data <- list("years"=1:sim$Nyears, "LF"=sim$LF, "obs_per_year"=sim$obs_per_year)
    }

    inits <- create_inputs(lh=lh, input_data=input_data, param="h", val=0.99)
    Nyears <- inits$Nyears 
    Nyears_comp <- nrow(inits$LF)

        LB_pars <- new("LB_pars")

        if(simulation==TRUE) species <- "sim"
        if(simulation==FALSE) species <- species
        LB_pars@Species <- species

        LB_pars@MK <- inits$M/inits$vbk         

        LB_pars@Linf <- inits$linf

        LB_pars@CVLinf <- inits$CVlen

        LB_pars@L50 <- inits$ML50 
        LB_pars@L95 <- inits$ML95

        LB_pars@Walpha <- inits$lwa
        LB_pars@Wbeta <- inits$lwb

        LB_pars@SL50 <- inits$SL50 
        LB_pars@SL95 <- inits$SL95
        LB_pars@BinWidth <- inits$binwidth

        LB_pars@Steepness <- inits$h
        LB_pars@R0 <- inits$R0


        LB_lengths <- new("LB_lengths")
        LB_lengths@LMids <- inits$mids
        LB_lengths@LData <- t(inits$LF)
        LB_lengths@Years <- 1:Nyears_comp
        LB_lengths@NYears <- Nyears_comp
        
        lbspr_res <- tryCatch(LBSPRfit(LB_pars=LB_pars, LB_lengths=LB_lengths, yrs=NA, Control=list(modtype="absel")), error=function(e) NA)
        # lbspr_res <- tryCatch(LBSPRfit(LB_pars=LB_pars, LB_lengths=LB_lengths, yrs=NA, Control=list(modtype="GTG")), error=function(e) NA)


        if(isS4(lbspr_res)){
          LBSPR_outs <- list()
          LBSPR_outs$pLF <- lbspr_res@pLCatch
          LBSPR_outs$SL50 <- lbspr_res@Ests[,"SL50"]
          LBSPR_outs$SL95 <- lbspr_res@Ests[,"SL95"]
          LBSPR_outs$FM <- lbspr_res@Ests[,"FM"]
          LBSPR_outs$SPR <- lbspr_res@Ests[,"SPR"]
          if(write==TRUE) saveRDS(LBSPR_outs, file.path(iterpath, "LBSPR_results.rds"))
        }
        if(isS4(lbspr_res)==FALSE){
          if(write==TRUE) saveRDS("non_convergence", file.path(iterpath, "non_convergence.txt"))
          if(write==FALSE){
            LBSPR_outs <- NA
          }
        }


        rm(LB_lengths)
        rm(LB_pars)
        rm(lbspr_res)
    }


    if(write==TRUE) return(paste0(max(itervec), " iterates run in ", modpath))
    if(write==FALSE) return(LBSPR_outs)

}

