adj_variation <- function(SigmaR, SigmaF, SigmaC, SigmaI, CVlen=0.1, rho=0){

### ----- simulation settings -----------###
## F1 =0.04, Fequil=0.25, Frate=0.3
lh_short <- create_lh_list(vbk=0.87, linf=36.2, lwa=0.05970, lwb=2.754, S50=11.3, M50=20.2, selex_input="length", maturity_input="length", M=1.49, binwidth=1, t0=-0.01, CVlen=CVlen, SigmaC=SigmaC, SigmaI=SigmaI, SigmaR=SigmaR, SigmaF=SigmaF, R0=1,  h=1, qcoef=1e-5, F1=0.25, Fequil=0.25, Frate=0.2, Fmax=0.7, start_ages=0, rho=rho)

## F1=0.009, Fequil=0.25, Frate=0.17
lh_med <- create_lh_list(vbk=0.21, linf=64.58, lwa=0.0245, lwb=2.790, S50=4, M50=34, selex_input="age", maturity_input="length", M=0.43, binwidth=1, t0=-0.01, CVlen=CVlen, SigmaC=SigmaC, SigmaI=SigmaI, SigmaR=SigmaR, SigmaF=SigmaF, R0=1,  h=1, qcoef=1e-5, F1=0.25, Fequil=0.25, Frate=0.17, Fmax=0.7, start_ages=0, rho=rho)

## F10.002, Fequil=0.025, Frate=0.15
lh_long <- create_lh_list(vbk=0.14, linf=100, lwa=5e-6, lwb=3, S50=2, M50=47.4, selex_input="age", maturity_input="length", M=0.15, binwidth=1, t0=-0.01, CVlen=CVlen, SigmaC=SigmaC, SigmaI=SigmaI, SigmaR=SigmaR, SigmaF=SigmaF, R0=1,  h=1, qcoef=1e-5, F1=0.25, Fequil=0.25, Frate=0.15, Fmax=0.7, start_ages=0, rho=rho)

lh_list <- list("Medium"=lh_med, "Short"=lh_short, "Long"=lh_long)
return(lh_list)

}