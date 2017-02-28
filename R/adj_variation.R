adj_variation <- function(SigmaR, SigmaF, SigmaC, SigmaI, CVlen=0.1, rho=0, selex_param=1, mat_param=1){

### ----- simulation settings -----------###
## F1 =0.04, Fequil=0.25, Frate=0.3

if(selex_param==1) S95 <- NULL
if(selex_param==2) S95 <- 22
if(mat_param==1) M95 <- NULL
if(mat_param==2) M95 <- 25
lh_short <- create_lh_list(vbk=0.87, linf=36.2, lwa=0.05970, lwb=2.754, S50=11.3, S95=S95, M50=20.2, M95=M95, selex_input="length", maturity_input="length", M=1.49, binwidth=1, t0=-0.01, CVlen=CVlen, SigmaC=SigmaC, SigmaI=SigmaI, SigmaR=SigmaR, SigmaF=SigmaF, R0=1,  h=1, qcoef=1e-5, F1=0.25, Fequil=0.25, Frate=0.2, Fmax=0.7, start_ages=0, rho=rho)

## F1=0.009, Fequil=0.25, Frate=0.17
if(selex_param==1) S95 <- NULL
if(selex_param==2) S95 <- 40
if(mat_param==1) M95 <- NULL
if(mat_param==2) M95 <- 35
lh_med <- create_lh_list(vbk=0.21, linf=64.58, lwa=0.0245, lwb=2.790, S50=37.0, S95=S95, M50=34, M95=M95, selex_input="length", maturity_input="length", M=0.43, binwidth=1, t0=-0.01, CVlen=CVlen, SigmaC=SigmaC, SigmaI=SigmaI, SigmaR=SigmaR, SigmaF=SigmaF, R0=1,  h=1, qcoef=1e-5, F1=0.25, Fequil=0.25, Frate=0.17, Fmax=0.7, start_ages=0, rho=rho)

## F10.002, Fequil=0.025, Frate=0.15
if(selex_param==1) S95 <- NULL
if(selex_param==2) S95 <- 40
if(mat_param==1) M95 <- NULL
if(mat_param==2) M95 <- 60
lh_long <- create_lh_list(vbk=0.14, linf=100, lwa=5e-6, lwb=3, S50=25.0, S95=S95, M50=47.4, M95=M95, M=0.15, selex_input="length", maturity_input="length", binwidth=1, t0=-0.01, CVlen=CVlen, SigmaC=SigmaC, SigmaI=SigmaI, SigmaR=SigmaR, SigmaF=SigmaF, R0=1,  h=1, qcoef=1e-5, F1=0.25, Fequil=0.25, Frate=0.15, Fmax=0.7, start_ages=0, rho=rho)


# par(mfrow=c(1,3))
# plot(lh_short$L_a)
# plot(lh_med$L_a)
# plot(lh_long$L_a)

lh_list <- list("Short"=lh_short, "Medium"=lh_med, "Long"=lh_long)
return(lh_list)

}