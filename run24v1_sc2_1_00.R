#github trying
#Excellent references
#https://mc-stan.org/docs/reference-manual/logit-transform-jacobian.html
#Miller and Hyun

#
#sc1: scenario 1 = M constant by both length and year; 
#sc2_1: scenario 2.1 = M varies by length class alone;
#sc2_2: scenario 2.2 = M varies by year alone; 

#
#load("C:/Hyun/Proj_PICES/R4_lb_sc2/save231009.RData");

#
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd();

##################################################################
############ compile 
system("rm sbmack.o sbmack.cpp sbmack.dll");
source("cpp24v1_sc2_1.R");

##################################################################
############## updated data 
source("C:/Hyun/Proj5/SizeTMB2/data_import.R"); 

##################################################################
##################################################################
##################################################################
############ MakeADFun and nlminb
N1_code=0;     #how to treat N_{0,a}s; #0 = deviations; #1=equilibrium
select_code=0; #fishery selectivity;  #0=logistc; #1=domeshape; #2=double-logistic

# nages: number of age classes
nages=6;  #nages=7;      

#
CV_length_r=0.1;  #10%; 

#M_input=0.37;  #natural mortality from JH;

#likelihood weights 
lambda1=0.1; #0.25 #0.2; #0.1;  #for the multinomial
CV_N1styr_devs=0.15; #0.65; #0.6; #100%, #50%;  
  #CV_Ft=0.15;  #0.1; #0.15;  #30%;
CV_exp_F_devs=0.15; 
CV_yield=0.2; #0.25; #0.1;
CV_CPUE=0.1; #0.2; #0.1;
CV_mean_logLinf=c(0.07, log(51.67));  #0.07 good

#two parameters in the logistic selectivity
"log_L50"=3.33;        #this value is from JW's MS thesis
"log_gamma"=log(0.45); #this value is from JW's MS thesis

#four parameters in the dome-shape selectivity
bounds_domeinf1=c(5.0, 35.0);
bounds_domeinf2=c(35.0, 55.0);
bounds_domeg1=c(0.1, 10.0);
bounds_domeg2=c(0.1, 10.0);
logit_domeinf1_domain=c(0.2683);
logit_domeinf2_domain=c(-0.9694);
logit_domeg1_domain=c(-2.5759);
logit_domeg2_domain=c(-2.5759);

select_dome_estimates =c(1.392, 5.349, 0.355, 0.661);  #estimates of four parameters in the dome-selectivity;
                                                       #JW's study 
#age_sel=c(0.2486038, 0.8418758, 0.9618046, 0.8844492, 0.6289894, 0.2719256);

select_doublelogistic_pars_values=c(0.45, 27.5, 1.0, 0.4, 55, 0.0); 
                                    #slope alpha1, inflection beta1, omega 1, 
                                    #slope alpha2, inflection beta2, omega 2

#N0s=c(1.22E+06,	452476,	177484,	44762.6,	37686,	1208.15)*1000;
#N0s=c(452000, 177000, 44000, 37000, 1000)*1000;  #nages-1; #nages, 6
#N0s=c(452000, 177000, 44000, 37000, 1000, 1000)*1000; #nages-1; #nages, 7

N0s.f=function(nages=6) {
  if(nages==6)
    N0s=c(1.0e+06, 450000, 177000, 44000, 37000, 1000)*1000 #nages-1;
  else if(nages==7)
    N0s=c(1.0e+06,  450000, 177000, 44000, 37000, 1000, 500)*1000; #nages-1;
  return(N0s); 
};

#
N0s=N0s.f(nages=nages); 

#
N1.code.f=function(N1_code=0, nages=nages, 
                   val_logN1_pars_low=log(1.0e-5), val_logN1_pars_up=log(5.0),
                   logN01=log(500.0E+06), logF1=log(0.5), 
                   val_logN01_low =log(1.0e+05), val_logN01_up=log(1.0e+12), 
                   val_logF1_low=log(0.01),  val_logF1_up=log(2.0) ) {
  if(N1_code==0)  {
    logN1_pars=rep(log(1.0),nages);
    logN1_pars_low=rep(val_logN1_pars_low, nages);
    logN1_pars_up=rep(val_logN1_pars_up, nages);
  }  
  else if(N1_code==1) {
    logN1_pars=c(logN01, logF1);
    logN1_pars_low=c(val_logN01_low, val_logF1_low); 
    logN1_pars_up=c(val_logN01_up, val_logF1_up);
  };
  #
  return(list(N1_code=N1_code, logN1_pars=logN1_pars, logN1_pars_low=logN1_pars_low, logN1_pars_up=logN1_pars_up)); 
};

logN1_pars_f_output = N1.code.f(N1_code=N1_code, nages=nages, 
                      val_logF1_low=log(0.01),  val_logF1_up=log(2.5));
logN1_pars = logN1_pars_f_output$logN1_pars;
logN1_pars_low=logN1_pars_f_output$logN1_pars_low;
logN1_pars_up=logN1_pars_f_output$logN1_pars_up;
logN1_pars_N1_code=logN1_pars_f_output$N1_code;
#
# Nmat_init = read.csv(file="estimates_abundances.csv", header=F) #in thousands; 
#logN_re.df = log(Nmat_init[-1,-1]*1000); 
logNAA = matrix(rep(c(20.0, 19.0, 19.0, 18.0, 17.0), (nyrs-1)), nrow=nyrs-1, ncol=nages-1, byrow=T); 

#
F_input=seq(0.0, 2.0, by=0.1);

#
data=list(
          "N1_code"=N1_code,
          "select_code"=select_code,
          
          "h"=0.99, #0.7,  
          
          "mort_frac"=0.5,
          "nages"=nages,
          
          "F_input"=F_input, 
          
          "N0s"=N0s,
          "data_yieldCPUE"=data_yieldCPUE_mackerel,
          "CV_length_r"= CV_length_r,
          "x"=x,
          "data_length_freq"=data_length_freq_mackerel,
          "data_length_freq_first_half"=data_length_freq_first_half,
          "data_length_freq_second_half"=data_length_freq_second_half,
          "data_length_freq_eff_ss"=data_length_freq_eff_ss,
          "neff"=neff,
          
          "bounds_q"=bounds_q,
          
          "bounds_domeinf1"=bounds_domeinf1,
          "bounds_domeinf2"=bounds_domeinf2, 
          "bounds_domeg1"= bounds_domeg1, 
          "bounds_domeg2"=bounds_domeg2, 
          
          "data_LW"=data_LW,
          "log_alpha_LW"=log_alpha_LW,
          "log_beta_LW"=log_beta_LW,        
          "data_maturation"=data_maturation,
          "b0_mat"=b0_mat,
          "b1_mat"=b1_mat,
          "ratio_female"=ratio_female,
          
          #"M_input"=M_input, 
          
          "CV_N1styr_devs"=CV_N1styr_devs, 
          #"CV_Ft"=CV_Ft, 
          "CV_exp_F_devs"= CV_exp_F_devs, 
          "CV_yield"=CV_yield,
          "CV_CPUE"=CV_CPUE,
          "lambda1"=lambda1, 
          "CV_mean_logLinf" = CV_mean_logLinf, 
            
          "select_dome_estimates"=select_dome_estimates, 
          "select_doublelogistic_pars_values"= select_doublelogistic_pars_values
          
         );

#
initial_values=list("mu_lengths_r"=18.0,
                    
                    "logN1_pars" = logN1_pars,
                    
                    "logR0"=log(1e+9),
                    
                    # "log_SR_recruits"=rep(log(1e+9), nyrs),  #random effects;

                    #"log_sig_logRec" = log(0.15), 
                    
                    "logNAA"=logNAA,
                    
                    "log_sig_logNAA"= log(0.15),
                      
                    "logF1"=log(0.5),
                    "F_devs"=rep(0, nyrs-1),
    
                    #"logM"=log(0.35),
                    #"M_devs"=rep(0,nyrs-1), 
                    
                    "logit_q_domain"=c(-8.0),      
                 
                    "logLinf"=log(51.67),
                    
                    "log_kappa"=log(0.1053992),
                    
                    "log_sig_L"=log(0.99),

                    "log_L50"=log_L50,
                    "log_gamma"=log_gamma
                    
                    #"logit_domeinf1_domain"=logit_domeinf1_domain,
                    #"logit_domeinf2_domain"=logit_domeinf2_domain,
                    #"logit_domeg1_domain"=logit_domeg1_domain,
                    #"logit_domeg2_domain"=logit_domeg1_domain
                 );


#cal_output <- file("out_cal.txt")
#sink(cal_output, type = "output")

length_based=MakeADFun(data=data,
                       parameters=initial_values,
                       DLL="sbmack",
                       random=c(#"logN_1st_time",
                                #"log_SR_recruits", 
                         "logNAA"
                       ),
                       
                       map=list(
                          "mu_lengths_r"=factor(NA),
                          
                           #"log_kappa"=factor(NA),
                           
                           #"log_sig_logMx_re"=factor(NA), #"log_M_constant"=factor(NA), 
                          
                           "logLinf" = factor(NA),
                          
                           #"log_sig_L"=factor(NA)
                          
                           "log_L50"=factor(NA),
                           "log_gamma"=factor(NA)
                          
                           #"logit_domeinf1_domain"=factor(NA),
                           #"logit_domeinf2_domain"=factor(NA),
                           #"logit_domeg1_domain"=factor(NA),
                           #"logit_domeg2_domain"=factor(NA)
                           ),
                       
                       silent=T
                       );

## set bounds
Low_bound=list(
  #"mu_lengths_r"=5.0,
  
  "logN1_pars"=logN1_pars_low, 
  "logR0"=log(1000),
  #"log_SR_recruits" are random effects;
  #"log_sig_logRec" = log(0.01), 
                    
  #logNAA are treated as random effects;
  "log_sig_logNAA"= log(0.01),
  
  "logF1"=log(0.001),
  "F_devs"=rep(-2.0, nyrs-1),
  
  #"logM"=log(0.001),
  
  "logit_q_domain"=c(-20.0),   
  
  #"logLinf"=log(18.0),
  
  "log_kappa"=log(0.001),
 
  "log_sig_L"=log(0.001)
  
  #"log_L50"=log(18.0)
  #"log_gamma"=log(0.001)
  
  #"logit_domeinf1_domain"=c(-1.609),
  ##"logit_domeinf2_domain"=c(-1.735),
  #"logit_domeg1_domain"=c(-3.168)
  ##"logit_domeg2_domain"=c(-3.168)
  );

Up_bound=list(
  # "mu_length_r"  =30.0,
  
  "logN1_pars"=logN1_pars_up, 
  "logR0"=log(1.0e+16),
  #"log_SR_recruits" are random effects;
  #"log_sig_logRec" = log(2.0), 
               
  #logN_re are treated as random effects;
  "log_sig_logNAA"= log(2.0),
  
  "logF1"=log(5.0),
  "F_devs"=rep(2.0, nyrs-1),
  
  #"logM"=log(3.0),    

  "logit_q_domain"=c(30.0),     
   
  #"logLinf"=log(70.0),
  
  "log_kappa"=log(3.0),
  
  "log_sig_L"=log(5.0)
  
   #"log_L50"=log(51.67)
   #"log_gamma"=log(5.0)
  
   #"logit_domeinf1_domain"=c(2.197),
   ##"logit_domeinf2_domain"=c(1.099),
   #"logit_domeg1_domain"=c(1.374)
   ## "logit_domeg2_domain"=c(1.374)
  );

#
fit=nlminb(length_based$par, length_based$fn, length_based$gr, 
           upper=unlist(Up_bound), lower=unlist(Low_bound), 
           control=list(eval.max=1000, iter.max=1000)
           );
fit;

#sink(type = "output")
#close(cal_output)

sdrep<-sdreport(obj=length_based);
sdrep;

# 
summary_results_raw_neff <- summary( sdrep);    #summary( sdreport(length_based));
summary_results_raw_neff;

TrueOrFalse_Hessian= sdrep$pdHess;
TrueOrFalse_Hessian;

max.abs.grad=max(abs( sdrep$gradient.fixed)); #maximum of the *absolute values* of all gradients;
max.abs.grad;

reports=length_based$report(); 

names(reports); 

reports$select_length;

reports$NAA;
reports$NAA/1.0e+6;

reports$Catch;

reports$BAA;
reports$B; 
reports$SpawnerBiomass;
reports$Spawners; 

reports$M_tx;
reports$Ft; 

reports$AIC;

reports$SSB_y;
reports$logRec;
reports$SR_alpha;
reports$SR_beta;

length_based$report()$Mu;
length_based$report()$pp;      #3d array
length_based$report()$f_total; #3d array

#
par(mfcol=c(2,2))
plot(x, reports$M_tx[1,], type='l', ylab="M", ylim=c(0.0, 1.0));
for(t in 2:nyrs)  
  lines(x, reports$M_tx[t,], lty=t); 

#plot(year_t, reports$Mx[,t], type='l', ylab="Mt", ylim=c(0.0, 1.0));

plot(x, reports$select_length, ylab="fishery selectivity", ylim=c(0.0, 1.0));

plot(year_t, reports$Ft, type="l", ylab="Ft"); 

#
tem=length_based$report()$f_total;
par(mfrow=c(3,3))
for(m in 1:9) {
  tem2=tem[,,m]; #3rd index:  year 
  plot(x, tem2[1,], type='l');   
  for(a in 2:6)
     lines(x, tem2[a,], lty=a, col=a);
};

for(m in 10:18) {
  tem2=tem[,,m]; #3rd index:  year 
  plot(x, tem2[1,], type='l');   
  for(a in 2:6)
    lines(x, tem2[a,], lty=a, col=a);
};

for(m in 19:22) {
  tem2=tem[,,m]; #3rd index:  year 
  plot(x, tem2[1,], type='l');   
  for(a in 2:6)
    lines(x, tem2[a,], lty=a, col=a);
};

#####################################################################
# Goodness of fit as figures
#==========================
# yield and CPUE
Yieldhat=length_based$report()$Yieldhat; 
All_Yield=data_yieldCPUE_mackerel[,"All_yield"];  

predCPUE=length_based$report()$predcpue;
LPS_CPUE=data_yieldCPUE_mackerel[,"LPS_cpue"];

par(mfrow=c(2,1), mar=c(5,5,5,3), oma=c(0,0,0,0))
plot(year_t, All_Yield/10^6,  ylim=c(0.05, 0.25), 
     xlab="Year", ylab=expression(paste("MT (x ", 10^6, ")")), main="Yield", 
     cex.lab=1.5, cex.axis=1.5, cex.main=2.0)
lines(year_t, Yieldhat/10^3/10^6, col="red")

#MSY
abline(h=reports$MSY/10^3/10^6, col="red", lty=2, lwd=3); 

#
plot(year_t, LPS_CPUE, xlab="Year", ylab="MT/haul", main="CPUE", ylim=c(10, 40), cex.lab=1.5,
     cex.axis=1.5, cex.main=2.0)
lines(year_t, predCPUE, col="red")

#
#=========================================== 
#length frequency;
Catch_hat=length_based$report()$Catch

length_composition_hat<-matrix(data=NA, nrow=nyrs, ncol=nlengths)
for(i in 1:nyrs)
   length_composition_hat[i,]<-(Catch_hat[i,]/sum(Catch_hat[i,]))*neff;
#
cexx=1.5;
par(mfrow=c(3,3), oma=c(3,3,0.5,0.5), mar=c(2, 2, 1.5, 1), cex=cexx)

for(i in 1:9) {
  plot(x, data_length_freq_eff_ss[i,], type="h", ylim=c(0,20)); 
  lines(x, length_composition_hat[i,], col='red');
};

for(i in 10:18) {
  plot(x, data_length_freq_eff_ss[i,], type="h", ylim=c(0,20)); 
  lines(x, length_composition_hat[i,], col='red');
};

for(i in 19:22) {
  plot(x, data_length_freq_eff_ss[i,], type="h", ylim=c(0,20)); 
  lines(x, length_composition_hat[i,], col='red');
};

#
#===========================================
#length frequency by age class
tem3d=length_based$report()$f_total; #3d array
tem3d[,,1];  #year 2001
tem3d[,,22]; #year 2022;

rowSums(tem3d[,,1]);
rowSums(tem3d[,,3]);
rowSums(tem3d[,,22]);

par(mfrow=c(5,1))
for(y in c(1,5,10,15,20)) {
 tem=tem3d[,,y]; 
 
 plot(x, tem[1,], type='l');
 lines(x, tem[2,], col="red");
 lines(x, tem[3,], col="blue");
 lines(x, tem[4,], col="green");
 lines(x, tem[5,], col="pink");
 lines(x, tem[6,], col='orange');
};

#
par(mfrow=c(1,1))
plot(x, tem3d[1,,2014-1999], type='l'); 
lines(x, tem3d[2,,2015-1999], col="red"); 
lines(x, tem3d[3,,2016-1999], col="blue"); 
lines(x, tem3d[4,,2017-1999], col="green"); 
lines(x, tem3d[5,,2018-1999], col="pink"); 
lines(x, tem3d[6,,2019-1999], col="orange"); 

#
#===============================================================
#YPR & SPR
par(mfrow=c(2,1));
plot(F_input, reports$YPR);
plot(F_input, reports$spawnpotenratio);
abline(h=0.3,lty=2);

#
#As of 22 Jan 2024
#AS of 22 Jan 2024
#
