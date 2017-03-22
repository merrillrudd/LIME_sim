adj_variation <- function(SigmaR, SigmaF, SigmaC, SigmaI, CVlen=0.1, rho=0, selex_param=1, mat_param=1, nseasons=1, binwidth=1, nbins=NULL){

### ----- simulation settings -----------###
## F1 =0.04, Fequil=0.25, Frate=0.3

if(selex_param==1) S95 <- NULL
if(selex_param==2) S95 <- 22
if(mat_param==1) M95 <- NULL
if(mat_param==2) M95 <- 25
if(is.null(nbins)) bw <- binwidth
if(is.null(nbins)==FALSE) bw <- floor(36.2/nbins)
lh_short <- create_lh_list(vbk=0.87, linf=36.2, lwa=0.05970, lwb=2.754, S50=11.3, S95=S95, M50=20.2, M95=M95, selex_input="length", maturity_input="length", M=1.49, binwidth=bw, t0=-0.01, CVlen=CVlen, SigmaC=SigmaC, SigmaI=SigmaI, SigmaR=SigmaR, SigmaF=SigmaF, R0=1,  h=1, qcoef=1e-5, F1=0.25, Fequil=0.25, Frate=0.2, Fmax=0.7, start_ages=0, rho=rho, nseasons=nseasons)

## F1=0.009, Fequil=0.25, Frate=0.17
if(selex_param==1) S95 <- NULL
if(selex_param==2) S95 <- 40
if(mat_param==1) M95 <- NULL
if(mat_param==2) M95 <- 35
if(is.null(nbins)) bw <- binwidth
if(is.null(nbins)==FALSE) bw <- floor(64.58/nbins)
lh_med <- create_lh_list(vbk=0.21, linf=64.58, lwa=0.0245, lwb=2.790, S50=32, S95=S95, M50=34, M95=M95, selex_input="length", maturity_input="length", M=0.43, binwidth=bw, t0=-0.01, CVlen=CVlen, SigmaC=SigmaC, SigmaI=SigmaI, SigmaR=SigmaR, SigmaF=SigmaF, R0=1,  h=1, qcoef=1e-5, F1=0.25, Fequil=0.25, Frate=0.17, Fmax=0.7, start_ages=0, rho=rho, nseasons=nseasons, AgeMax=15)

## F10.002, Fequil=0.025, Frate=0.15
if(selex_param==1) S95 <- NULL
if(selex_param==2) S95 <- 40
if(mat_param==1) M95 <- NULL
if(mat_param==2) M95 <- 60
if(is.null(nbins)) bw <- binwidth
if(is.null(nbins)==FALSE) bw <- floor(100/nbins)
lh_long <- create_lh_list(vbk=0.13, linf=90, lwa=0.02635, lwb=2.9595, S50=25, S95=S95, M50=50, M95=M95, M=0.18, selex_input="length", maturity_input="length", binwidth=bw, t0=-0.01, CVlen=CVlen, SigmaC=SigmaC, SigmaI=SigmaI, SigmaR=SigmaR, SigmaF=SigmaF, R0=1,  h=1, qcoef=1e-5, F1=0.25, Fequil=0.25, Frate=0.15, Fmax=0.7, start_ages=0, rho=rho, nseasons=nseasons)


# par(mfrow=c(1,3))
# plot(lh_short$L_a)
# plot(lh_med$L_a)
# plot(lh_long$L_a)

lh_list <- list("Short"=lh_short, "Medium"=lh_med, "Long"=lh_long)
return(lh_list)

}