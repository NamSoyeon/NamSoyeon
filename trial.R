#NAA not random effects
#M by lengths (Charnov et al.)
#N1
#Estimate sigO
#Schunute It
#no f_total
#no growth in plus group downward calculation
#trying

#v2: log(R_t) = log(R_{t-1})+ error_t = log(f(alpha, beta, SSB_{t-2})) + error_t;
#v3: sig_logRec = sig_logNAA; 
#v4: penalized likelihood or prior with the other uncertainty being treated as free parameters
##DATA_SCALAR(CV_N1styr_devs); //less important because this assumption is applied to the first year alone
##DATA_SCALAR(CV_exp_F_devs);  //DATA_SCALAR(CV_Ft);
##//DATA_SCALAR(CV_yield);
##//DATA_SCALAR(CV_cpue);
#
#v5: (1) the uncertainty in both yield and cpue data are externally given;
#    (2) using fishing effort data from the large purse seine fishery; 
#        Fully selected F = q * Effort; 
#
#the Dirichlet-Multinomial likelihood for the length composition data; #as of 8 Feb 2024;
#xx0: the traditional multinomial for the length composition data,
#     where the uncertainty of cpue is NOT estimated;
#xx1: the dirichlet-multinomial for the length composition data, 
#     where an additional parameter, say theta is estimated;
#
#sc1: scenario 1 = estimate M, which is constant by both length and year; 
#sc2_1: scenario 2.1 = M varies by length class alone;
#sc2_2: scenario 2.2 = M varies by year alone; 

#
library(TMB);

#
sbmack <- '
 //A size-based (sb) model for the Korean mackerel (mack) stock assessment;
 //Adding R_T, SPR, YPR, and MSY, etc.
 //Authors: Gim and Hyun as of 10 Feb 2024; 
 
 #include <TMB.hpp>
 
 // pass missing values
 template <class Type>
 bool isNA(Type x){
    return R_IsNA(asDouble(x));
 }
  
 // square
 template <class Type>
 Type square(Type i) {
    return i*i; 
 }
 
 //dome shape selectivity by age, where parameter estimates are from Gim study.
 template <class Type> 
 vector<Type> dome_shape_sel_by_age_f(int nages, vector<Type> select_dome_estimates) {
 Type one = Type(1.0); 
 Type inf1=select_dome_estimates(0); //inf1
 Type inf2=select_dome_estimates(1); //inf2
 Type g1=select_dome_estimates(2); //g1
 Type g2=select_dome_estimates(3); //g2
 //
 vector<Type> age(nages);
 vector<Type> selectivity(nages); 
 for(int a_index=0;a_index<nages;a_index++) {
    age(a_index)=a_index+1;
    selectivity(a_index)=(one/(one+exp(( Type(-1.0)*age(a_index)+inf1 )/g1 )) )*(one/( one+exp( (age(a_index)-inf2 )/g2 )) );
 };
 return selectivity;
 };

 //dome shape selectivity by length;
 template <class Type> 
 Type dome_shape_sel_by_length_f(Type length, vector<Type> select_dome_pars) {
 Type one = Type(1.0);
 Type inf1=select_dome_pars(0); //inf1
 Type inf2=select_dome_pars(1); //inf2
 Type g1=select_dome_pars(2); //g1
 Type g2=select_dome_pars(3); //g2
 //
 Type selectivity=(one/(one+exp(( Type(-1.0)*length+inf1 )/g1 )) )*(one/( one+exp( (length-inf2 )/g2 )) );
 return selectivity;
 };
 
 //double logistic selectivity by length;
 template <class Type> 
 Type doublelogistic_sel_by_length_f(Type length, vector<Type> select_doublelogistic_pars) {
 Type one=Type(1.0); 
 Type al1=select_doublelogistic_pars(0); //slope alpha1
 Type be1=select_doublelogistic_pars(1); //inflection beta1
 Type om1=select_doublelogistic_pars(2); //omega 1
 Type al2=select_doublelogistic_pars(3); //slope alpha2
 Type be2=select_doublelogistic_pars(4); //inflection beta2
 Type om2=select_doublelogistic_pars(5); //omega 2
 //
 Type selectivity=om1/(one+exp(-al1*(length-be1 )))-(1-om2)/(one+exp(-al2*(length-be2)) );
 return(selectivity); 
 };
 
 //loglikelihood for Dirichlet-Multinomial likelihood for the length compositions in a year
 template <class Type>
 Type logDiriMultinom_f(int year, int nlengths, int SampleSize, Type theta, vector<Type> prop_obs, 
                      vector<Type> prop_pred) {
 Type zero=Type(0.0);
 Type one=Type(1.0);
 vector<Type> freq_obs(nlengths);
 vector<Type> freq_pred(nlengths); 
 freq_obs=SampleSize*prop_obs;
 freq_pred=SampleSize*prop_pred;
 Type ll=zero;
 ll+=lgamma(SampleSize+one)+lgamma(SampleSize*theta)-lgamma(SampleSize+SampleSize*theta);
 for(int il=0; il<nlengths; il++)
    ll+= -lgamma(freq_obs(il)+one) + lgamma(freq_obs(il)+freq_pred(il)*theta +Type(1.0e-10))-lgamma(freq_pred(il)*theta+Type(1.0e-10));
                                                              //Type(1.0e-10) to avoid Inf
 return ll; //return loglikelihood;
 };
 
 //objective function
 template <class Type>
 Type objective_function<Type>::operator() () {
  
 //Data section
 //DATA_INTEGER(N1_code);      //the option for the assumption of age-specific abundances in the 1st year, N_{0,a}s;
                             //0 = deviations; 1=equilibrium
 DATA_INTEGER(select_code);  //the option for the fishery selectivity;
                             //0= logistic; 1 = domeshape; 2 = double logistic 
 DATA_INTEGER(lengthcomp_code);  //0= multinomial; 1 = dirichlet-multinomial
 
 DATA_SCALAR(h);             //steepness 
 DATA_SCALAR(mort_frac);
                                                               //1= around the mean and with a prior 
 DATA_INTEGER(nages);     //number of imaginary age classes; 
 DATA_VECTOR(N0s);        //plausible abundances at the 1st year; //N01, N02,  ... ,N0A; //length=nages-1;
 
 DATA_MATRIX(data_yieldcpue); //yield and cpue data; 
                              //23 rows (23 years of 2000 - 2023) and 5 columns (year, LPS_Yield, TheOthers, All_yield, LPS_cpue);
 //length frequency data
 DATA_SCALAR(CV_length_r);        //CV of lengths at age 1
 DATA_VECTOR(x);                  //a vector of the respective midpoints of length classes (43 classes);
 DATA_MATRIX(data_length_freq);   //length frequency data (23 years x 43 classes); //Do NOT use IMATRIX;  //Doyul found it; 
 //DATA_MATRIX(data_length_freq_first_half);
 //DATA_MATRIX(data_length_freq_second_half);
 DATA_IVECTOR(SS_length);       //raw sample size for the length composition data; 

 //QAA = q * selectivity for the predicted cpue
 DATA_VECTOR(bounds_q);      //logit(q)

 //Estimates of the four parameters in the dome shape selectivity for the LPS (large purse seine) fishery
 //These values are from the previous study (J. Gim)
 DATA_VECTOR(bounds_domeinf1);
 DATA_VECTOR(bounds_domeinf2);
 DATA_VECTOR(bounds_domeg1);
 DATA_VECTOR(bounds_domeg2);
 
 //length-weigth relationship
 DATA_MATRIX(data_LW);        //the matrix of length-weight data (related to parameter log_alpha_LW, and log_beta_LW )
 DATA_SCALAR(log_alpha_LW);
 DATA_SCALAR(log_beta_LW);

 //length-maturation relationship
 DATA_MATRIX(data_maturation);       //the matrix of length-matruation ratio data (related to parameter b0_mat and b1_mat)
 DATA_SCALAR(b0_mat);
 DATA_SCALAR(b1_mat);
 DATA_SCALAR(ratio_female);
 
 //DATA_SCALAR(M_input); //natural mortality

 //CV 
 //DATA_SCALAR(CV_N1styr_devs); //less important because this assumption is applied to the first year alone
 //DATA_SCALAR(CV_exp_F_devs);  //DATA_SCALAR(CV_Ft);
 //DATA_SCALAR(CV_yield);
 //DATA_SCALAR(CV_cpue);

 //DATA_VECTOR(CV_mean_logLinf); 
 DATA_VECTOR(Shaperate_sig_L); //beta penalty

 
 //Estimates of parameters in the dome shape selectivity 
 //for the LPS (large purse seine) fishery; //These estimate are from the previous study, J. Gim.
 DATA_VECTOR(select_dome_estimates);  //dome_inf1, dome_inf2, dome_g1, dome_g2;

 //six parameters assumed for the double-logistic selectivity 
 DATA_VECTOR(select_doublelogistic_pars_values); //slope alpha1, inflection beta1, omega 1, slope alpha2, inflection beta2, omega 2

 //Fishing effort; 
 //DATA_VECTOR(effort);
 
 //for the upper bound of censored catch data; // observedCatch < trueCatch < p_upper_catch*observedCatch;
 //DATA_VECTOR(upper_catch); 
 
 //YPR, SPR
 DATA_VECTOR(F_input);
 int indnum_F =  F_input.size();
 
 Type zero = Type(0.0);
 Type one = Type(1.0); 
 Type two = Type(2.0);
 
 //std::cout<<"effort: "<<effort.transpose()<<std::endl;
 //PARAMETER(b0_mat); 
 //PARAMETER(b1_mat);  
 
 //==================================================================================
 //==================================================================================
 //Parameter section
 PARAMETER(mu_lengths_r);                     //the population mean of lengths at age 1;
 //PARAMETER(sig2_lengths_r);                   
 Type sig2_lengths_r=square(CV_length_r*mu_lengths_r);//the variance of lengths at age 1;
 
 //abundances in initial time
 //Age-specific abundances in the 1st year, following the NOAA ASAP practice where 
 //N_{0,a} = plausiable values x deviation, and log(deviations) become free parameters.
 //if(N1_code==0)
 //PARAMETER_VECTOR(logN1_pars)  //lengths = nages-1 or 2;
                                  //option 0: log(dev_N01), ..., log(dev_N0A); //log(dev)~Normal(0,sig2);
 //else if(N1_code==1)
 PARAMETER(logN1_par);   //option 1: one element, log(NAA(0,0)); 
                        
 PARAMETER(logR0);   //unexploited recruitment
 
 PARAMETER_VECTOR(logRec);  //log(recruits) as random effects;
 PARAMETER(log_sig_logRec);  
 //PARAMETER_MATRIX(logNAA);  //log(abundances)(nyrs-1,nages-1) as random effects, excluding cells at the first row and column 
 //PARAMETER(log_sig_logNAA); //log(sd of the process errors, e_{y,a}); //log(N_{y,a})=log(N_{y-1,a-1})-Z_{y-1,a-1}+e_{y,a}; 
 PARAMETER(logF1); 
 PARAMETER_VECTOR(logFt); 
 PARAMETER(log_sig_logF); 
 
 PARAMETER(sig_yield);
 PARAMETER(sig_cpue); 
 
 //PARAMETER(logq_F); //logFt = logq_F + log(effort)
 
 //natural mortality   
 //PARAMETER(logM);          //common M      
 //Type M=exp(logM);

 //QL(xind) = q*select_length(xind);   //cpue data are also from the LPS fishery;
 PARAMETER(logit_q_domain);   //logit_q_domain=logit((q-lower)/(upper-lower));
 
 //growth
 PARAMETER(logLinf);
 PARAMETER(log_kappa);  //kappa in von-Bertalanffy growth equation;
                        //Brody coefficient, Rho=exp(Type(-1.0)*kappa);
 PARAMETER(log_sig_L);  //sig_L, where L is the L_{a+1} equation; 
 
 int nlengths=x.size();                     //the number of length classes; // median length of maxFL = 52.5 cm; 
 vector<Type> select_length(nlengths);      //fishery selectivity 
 select_length.setZero();
 
 if(select_code==0) { //logistic fishery selectivity;
   PARAMETER(log_L50);                    //gear selectivity parameter;
   PARAMETER(log_gamma);                  //gear selectivity parameter;
  
   Type L50=exp(log_L50);              //logistic selectivity;
   Type gamma=exp(log_gamma);          //logistic selectivity;
   
   for(int xind=0; xind<nlengths; xind++)   //gear selectivity by length
       select_length(xind)= one/(one+exp(Type(-1.0)*gamma*(x(xind)-L50))); 
   }
 else if(select_code==1) { //dome fishery selectivity; //four parameters: inf1, inf2, g1, and g2
   PARAMETER(logit_domeinf1_domain);  
   PARAMETER(logit_domeinf2_domain);  
   PARAMETER(logit_domeg1_domain);  
   PARAMETER(logit_domeg2_domain);  
    
   Type domeinf1_lower=bounds_domeinf1(0);
   Type domeinf1_upper=bounds_domeinf1(1); 
   Type domeinf2_lower=bounds_domeinf2(0);
   Type domeinf2_upper=bounds_domeinf2(1); 
   Type domeg1_lower=bounds_domeg1(0);
   Type domeg1_upper=bounds_domeg1(1); 
   Type domeg2_lower=bounds_domeg2(0);
   Type domeg2_upper=bounds_domeg2(1); 
   
   Type domeinf1 = domeinf1_lower+(domeinf1_upper - domeinf1_lower)/(1+exp(-logit_domeinf1_domain));
   Type domeinf2 = domeinf2_lower+(domeinf2_upper - domeinf2_lower)/(1+exp(-logit_domeinf2_domain));
   Type domeg1 = domeg1_lower+ (domeg1_upper - domeg1_lower)/(1+exp(-logit_domeg1_domain));
   Type domeg2 = domeg2_lower+(domeg2_upper - domeg2_lower)/(1+exp(-logit_domeg2_domain));
   vector<Type> select_dome_pars(4);
   select_dome_pars(0)=domeinf1;
   select_dome_pars(1)=domeinf2;
   select_dome_pars(2)=domeg1;
   select_dome_pars(3)=domeg2;
   
   for(int xind=0; xind<nlengths; xind++)   //gear selectivity by length
      select_length(xind)=dome_shape_sel_by_length_f(x(xind), select_dome_pars);                   
 }
 else if(select_code==2) { //double logistic selectivity
                           //slope alpha1, inflection beta1, omega 1, slope alpha2, inflection beta2, omega 2
   for(int xind=0; xind<nlengths; xind++)   //gear selectivity by length
      select_length(xind)=doublelogistic_sel_by_length_f(x(xind), select_doublelogistic_pars_values);
      
   select_length/=max(select_length); 
 };
 
 //PARAMETER(logsig_yield); 
 //PARAMETER(logsig_cpue);   //sd of cpue, which are a lognormal RV
 
 //==================================================================================
 //==================================================================================
 //Derived quantities of data
 Type SSB0;   //spawning biomass under no fishing;
 
 //vector<Type> yield=data_yieldcpue.col(1);    //LPS fishery yield data in data_yieldcpue
                                                //Assuming that LPS fishery catch = the total fishery catch;  
 vector<Type> yield=data_yieldcpue.col(3);      //total yield data in data_yieldcpue
 vector<Type> cpue=data_yieldcpue.col(4);       //cpue data in data_yieldcpue
 vector<Type> n_obs=data_maturation.col(1);       //maturation level in maturation data
 vector<Type> Matdat=data_maturation.col(2);       //maturation level in maturation data

 int nyrs=yield.size();                           //the number of years
 int r=1;                                         //recruit is defined as the pop size at age 1; 
 
 vector<Type> L=x;                                //length after growth; 
 
 vector<Type> data_length_LW=data_LW.col(0);      //length data in data_LW
 vector<Type> data_weight_LW=data_LW.col(1);      //weight data in data_LW
 Type alpha_LW=exp(log_alpha_LW);         
 Type beta_LW=exp(log_beta_LW);           
 //int size_data_LW=data_length_LW.size();
 vector<Type> data_length_maturation=data_maturation.col(0);      //length data in data_maturation
 vector<Type> data_rate_maturation=data_maturation.col(1);        //maturation rate data in data_maturation
 //int size_data_maturation=data_length_maturation.size();

 Type q_lower=bounds_q(0);
 Type q_upper=bounds_q(1);

 Type kappa=exp(log_kappa);              //von-Bertalanffy;
 Type Rho=exp(Type(-1.0)*kappa);         //Brody coefficient;
 
 Type sig_L=exp(log_sig_L);
 
 //QL for predicted cpue = q * selectivity 
 vector<Type> QL(nlengths); 
 QL.setZero();
 
  //Derived quantities of parameters
 matrix<Type> logM_tx(nyrs,nlengths); //instantaneous natural mortality
 logM_tx.setZero();
 matrix<Type> M_tx(nyrs,nlengths);
 M_tx.setZero(); 
 
 vector<Type> pred_Wt(nlengths);                 //body weight by length classes;
 pred_Wt.setZero();
 vector<Type> pred_Maturation(nlengths);         //maturation rate by length classes; 
 pred_Maturation.setZero();
 
 vector<Type> Ft=exp(logFt);
 //Ft.setZero();
 matrix<Type> F_tx(nyrs,nlengths);               //instantaneous fishing mortality; considering selectivity
 F_tx.setZero();
 matrix<Type> Z_tx(nyrs,nlengths);               //instantaneous fishing mortality + instantaneous natural mortality
 Z_tx.setZero();
 matrix<Type> ExpZ_tx(nyrs,nlengths);            //survival rate
 ExpZ_tx.setZero();

 vector<Type> Mu(nlengths);                      //diffe by length
 Mu.setZero();
 
 vector<Type> SS(nages);                         //differ by age
 SS.setZero();
 vector<Type> p(nlengths);
 p.setZero();
 vector<Type> p_plus(nlengths);                  //the last age class; 
 p_plus.setZero(); 
 
 array<Type> f_total(nages, nlengths, nyrs);
 f_total.setZero();
 array<Type> pp(nlengths,nlengths,nages+1);        //pp(L,x,a); //where nlengths (row) x nlengths (columns) at each level;
 pp.setZero();  

 //Type logLinf=CV_mean_logLinf(1);
 Type Linf=exp(logLinf); 
 
 array<Type> NL(nyrs,nlengths,nages);  //at the level of each age
 NL.setZero();
   //array<Type> NL_re(nyrs,nlengths,nages);     //at the level of each age
   //NL_re.setZero();
 
 matrix<Type> NAA(nyrs,nages);      //prevously NAA(nages,nyrs);
 NAA.setZero();
 matrix<Type> pred_NAA(nyrs,nages); //pred_NAA = N_{y-1,a-1}*exp(-Z_{y-1,a-1}); 
 pred_NAA.setZero(); 
   //matrix<Type> N_save(nyrs, nages);//N_save(nyrs-1, nages-1);  //N_save = N_{y-1,a-1}*exp(-Z_{y-1,a-1}); 
   //N_save.setZero(); 
 
 matrix<Type> Spawners(nyrs,nages);   
 Spawners.setZero();
 array<Type> SpawnersL(nyrs,nlengths,nages); //at the level of each age
 SpawnersL.setZero();
 matrix<Type> SpawnerBiomass(nyrs,nages);   
 SpawnerBiomass.setZero();
 array<Type> SpawnerBiomassL(nyrs,nlengths,nages); //at the level of each age
 SpawnerBiomassL.setZero();
 vector<Type> SSB_y(nyrs);
 SSB_y.setZero();
 vector<Type> recruits_pred(nyrs);
 recruits_pred.setZero();
 
 array<Type> predcpue_3d(nyrs,nlengths,nages); 
 predcpue_3d.setZero(); 
 matrix<Type> predcpue_ya(nyrs,nages);
 predcpue_ya.setZero();
 vector<Type> predcpue(nyrs);
 predcpue.setZero();
 array<Type> Yield3D(nyrs,nlengths,nages);
 Yield3D.setZero();
                       
 Type CNum=zero;
 Type CWt=zero;
 vector<Type> TCatch(nyrs);
 TCatch.setZero();
 matrix<Type> Catch(nyrs,nlengths);
 Catch.setZero();
 vector<Type> Yieldhat(nyrs);
 Yieldhat.setZero();
 
 vector<Type> Pop(nyrs); 
 Pop.setZero();
 array<Type> ENx(nyrs,nlengths,nages); 
 ENx.setZero();
 vector<Type> EN(nyrs);
 EN.setZero();
 
 array<Type> BL(nyrs,nlengths,nages); //see NL, which is a 3d array
 BL.setZero(); 
 matrix<Type> BAA(nyrs,nages); 
 BAA.setZero();
 vector<Type> B(nyrs); 
 B.setZero();
 vector<Type> EB(nyrs); 
 EB.setZero();
 array<Type> BDY(nyrs,nlengths,nages); //Biomass for calculation of CPUE during the year
 BDY.setZero(); 

 vector<Type> Neff_lengthcomp(nyrs); //effective sample size for the length composition 
 Neff_lengthcomp.setZero(); 
 
 vector<Type> dispersion_lengthcomp(nyrs);
 dispersion_lengthcomp.setZero();
 
 //Weight and maturation rate by length
 for(int xind=0; xind<nlengths; xind++) {                     //dont forget zero index in TMB
    pred_Wt(xind)=alpha_LW*pow(x(xind),beta_LW)/Type(1000);   //the division of 1000 is to convert gram to kg;
    pred_Maturation(xind)=one/(one+exp(b0_mat-b1_mat*x(xind))); //maturation rate 
 }; 
 
 for(int m=0; m<nyrs; m++) {
   for(int xind=0; xind<nlengths; xind++) { 
     logM_tx(m,xind)=-1.5*(log(x(xind)) - log(Linf) )+log(kappa);  // Charnov et al. (2013)
     M_tx(m,xind)=exp(logM_tx(m,xind));
     //logM_tx(m,xind)=0.55-1.61*log(x(xind))+1.44*log(Linf)+log(kappa); //Gislason et al. (2010)
  }
 };
 
// std::cout<<"M_tx: "<<M_tx<<std::endl;

/*
 //Fishing mortality
 //When fishing effort is available;
 for(int t=0; t<nyrs; t++) {
    logFt(t)=logq_F+log(effort(t)); 
    Ft(t)=exp(logFt(t));
 };
*/

/* 
 //When fishing effort is missing
 for(int t=1; t<nyrs; t++)  {
    logFt(t)=logFt(t-1) + F_devs(t-1);
    Ft(t)=exp(logFt(t)); 
 };
*/
 
 for(int t=0; t<nyrs; t++) { //t is year; //it is m, month in Quinns code; 
    for(int xind=0; xind<nlengths; xind++){
       F_tx(t,xind)=select_length(xind)*Ft(t); 
       Z_tx(t,xind)=M_tx(t,xind)+F_tx(t,xind);
          //Z_tx(t,xind)=exp(logM_tx_re(xind))+F_tx(t,xind);
       ExpZ_tx(t,xind)=exp(Type(-1.0)*Z_tx(t,xind));  //survival; 
    };
 };
 
 vector<Type> ExpZ_tx_colsums(nlengths);
 ExpZ_tx_colsums=ExpZ_tx.colwise().sum();
 
 vector<Type> ExpZ_tx_mean(nlengths);
 ExpZ_tx_mean = ExpZ_tx_colsums/nyrs;
 
 //std::cout<<"ExpZ_tx_colsums: "<<ExpZ_tx_colsums<<std::endl;
 //std::cout<<"ExpZ_tx_mean: "<<ExpZ_tx_mean<<std::endl;
 
 //LVB body growth; 
 int a=0;
 SS(a)=sig2_lengths_r;             //SS(0): Var{lengths at age 1}
 for(int m=0;m<nyrs;m++)
   for(int xind=0;xind<nlengths;xind++)
     f_total(a,xind,m)=dnorm(x(xind), mu_lengths_r, sqrt(SS(a)))/sum(dnorm(x,mu_lengths_r,sqrt(SS(a))));  
 
 for(int xind=0; xind<nlengths; xind++) {
    Mu(xind)=Linf-(Linf-x(xind))*Rho; 
 };
 
 //std::cout<<"Mu(xind) at age 1 (i.e., a=0): "<<Mu.transpose()<<std::endl;
 
 for(int a=1; a<nages; a++) { //be careful the index in TMB; //i.e., a=1 means the 2nd element; 
    //this SS is from Cohen and Fishman (1980); //it was used for the shrimp in the Quinns paper;
    SS(a)=square(sig_L)*(one-pow(Rho,(two*(a+1)-two*r)))/(one-square(Rho))+pow(Rho,(two*(a+1)-two*r))*sig2_lengths_r; 
 };
 
 Type kkk;
 for(int a=0; a<nages; a++) {
    for(int xind=0; xind<nlengths; xind++) {
       kkk=Type(0.0); 
       for(int Lind=0; Lind<nlengths; Lind++) {
          pp(Lind,xind,a)=Type(0.0);
            
          if(Lind>=xind)  {
             pp(Lind,xind,a)=dnorm(L(Lind),Mu(xind),sqrt(SS(a)));  // f(L|x) in Quinn et al. (1998); 
             kkk=kkk+pp(Lind,xind,a);
             //std::cout<<"Mu(xind): "<<Mu(xind)<<std::endl; 
          };
       };
       for(int Lind=0; Lind<nlengths; Lind++) {
          pp(Lind,xind,a)=pp(Lind,xind,a)/kkk;  //normalize f(L|x);
       };   
    };
 };
 
 for(int Lind = 0; Lind < nlengths; Lind++) {
   for(int xind = 0; xind < nlengths; xind++) {
     if(Lind == xind) {
       pp(Lind, xind, nages) = 1.0; 
     } else {
       pp(Lind, xind, nages) = 0.0;
     };
   };
 };

 /*
 std::cout<<"pp, as a 3d array. Are the matrix (L,x) different by age?"<<std::endl;
 for(int Lind=0;Lind<nlengths;Lind++) {
    std::cout<<"pp(Lind,5,1)"<<std::endl;
    std::cout<<pp(Lind,5,1)<<std::endl;
    std::cout<<"pp(Lind,5,3)"<<std::endl; 
    std::cout<<pp(Lind,5,3)<<std::endl; 
 };
 */

 // R0, SSB0, h for calculation of the two parameters, say alpha and beta in the Beverton-Holt model; SPR
 //SPR_0 ???
 /// p_ref0 ///
 matrix<Type> p_ref0(nages,nlengths);
 vector<Type> p_refA0(nlengths);
 Type R0;
 Type SR_alpha;
 Type SR_beta;
 
 for(int Lind=0; Lind<nlengths; Lind++) {
      p_ref0(0,Lind) = f_total(0,Lind,nyrs-1); //0 = age 1; Lind = length index; nyrs-1 = terminal year
 };

 for(int iage=1; iage<nages; iage++) {
    for(int Lind=0; Lind<nlengths; Lind++) {
        p_ref0(iage,Lind)=zero;
        for(int xind=0; xind<nlengths; xind++) {
             p_ref0(iage,Lind)+=p_ref0(iage-1,xind)*exp(Type(-1.0)*M_tx(nyrs-1,xind))*pp(Lind,xind,iage);
        };  
    }; 
 };

 for(int Lind=0; Lind<nlengths; Lind++) {
      p_refA0(Lind)=zero;
      for(int xind=0; xind<nlengths; xind++) {  
         p_refA0(Lind)=p_ref0(nages-1,xind)*exp(Type(-1.0)*M_tx(nyrs-1,xind))*pp(Lind,xind,nages-1);
     };
 };

 for(int Lind=0; Lind<nlengths; Lind++) {
     p_ref0(nages-1,Lind)+=p_refA0(Lind);
 };

 vector<Type> SPR_0_vec(nlengths);
 SPR_0_vec.setZero();
 Type SPR_0; // scalar

 for(int Lind=0;Lind<nlengths;Lind++){
      SPR_0_vec(Lind)=p_ref0(0,Lind)*pred_Wt(Lind)*pred_Maturation(Lind)*ratio_female*exp(-mort_frac*M_tx(nyrs-1,Lind))+ // new 231222
                      p_ref0(1,Lind)*pred_Wt(Lind)*pred_Maturation(Lind)*ratio_female*exp(-mort_frac*M_tx(nyrs-1,Lind))+
                      p_ref0(2,Lind)*pred_Wt(Lind)*pred_Maturation(Lind)*ratio_female*exp(-mort_frac*M_tx(nyrs-1,Lind))+
                      p_ref0(3,Lind)*pred_Wt(Lind)*pred_Maturation(Lind)*ratio_female*exp(-mort_frac*M_tx(nyrs-1,Lind))+
                      p_ref0(4,Lind)*pred_Wt(Lind)*pred_Maturation(Lind)*ratio_female*exp(-mort_frac*M_tx(nyrs-1,Lind))+
                      p_ref0(5,Lind)*pred_Wt(Lind)*pred_Maturation(Lind)*ratio_female*exp(-mort_frac*M_tx(nyrs-1,Lind));
  };
  
  SPR_0=SPR_0_vec.sum();
 
  R0=exp(logR0);
  SSB0=SPR_0*R0; // SPR_0 required previously??; // I set the SPR_0 calculation above
  SR_alpha=Type(4.0)*h*R0/((5.0*h)-1.0);
  SR_beta=SSB0*(1.0-h)/(5*h-1.0);
  
  // ==========================================================================================
  // ====== S: length frequency at different age classes in the 1st time (e.g., 1st year)
  int m=0;
  for(int a=1; a<nages; a++) {  //the case of a=0 was made above; 
     for(int Lind=0; Lind<nlengths; Lind++) {
        for(int xind=0; xind<nlengths; xind++) {                 //note ExpZ_tx_mean below instead of ExpZ_tx; 
            f_total(a,Lind,m)+=f_total(a-1,xind,m)*ExpZ_tx_mean(xind)*pp(Lind,xind,a);   
        };
     };
  };
  
  //normalize f_total(a,Lind,0); 
  vector<Type> temsumbyage(nages-1); //the case of a=0 was made above; 
  for(int newa=0; newa<nages-1; newa++) { 
    Type temSum=0.0;
    for(int Lind=0;Lind<nlengths;Lind++) {
        temSum+=f_total(newa+1,Lind,0);
        temsumbyage(newa)=temSum;
    };
  };
  for(int newa=0; newa<nages-1; newa++)
     for(int Lind=0;Lind<nlengths;Lind++)
        f_total(newa+1,Lind,0)/=temsumbyage(newa);
 
  //std::cout<<"f_total(0,L,0): "<<std::endl; 
  //for(int Lind=0;Lind<nlengths;Lind++)
  //    std::cout<<f_total(0,Lind,0)<<std::endl;
     
  // ====== E: length frequency at different age classes in the 1st time (e.g., 1st year)
  // ==========================================================================================
 
 
  // ==========================================================================================
  // ====== S: f_total(0,xind,0), f_total(0,xind,1), ..., f_total(0,xind,m), ..., f_total(0,xind,nyrs-1)
  // N_{0,0}, N_{0,1}, N_{0,2}, ... , N_{0,A};  //vector<Type> N_init(nages-1);
  // N_{0,0}

  Type M_tx_sum=M_tx.row(0).sum(); 
  Type M_tx_mean=M_tx_sum/nlengths;
  
  Type F_tx_sum=F_tx.row(0).sum();
  Type F_tx_mean=F_tx_sum/nlengths;
 
  //std::cout<<"M_0x_mean: "<<M_0x_mean<<std::endl;
 
 /*
  for(int a=0;a<nages;a++) {  //N1styr_devs(1), N1styr_devs(2), ..., N1styr_devs(A);
    //if(N1_code==0)
    //   NAA(0,a)=N0s(a)*exp(logN1_pars(a));   //N0s(a)*N1styr_devs(a); 
        
    //else if(N1_code==1) {
    if(a==0) 
          NAA(0,a) = exp(logN1_par);     //NAA(0,0)
    else    {
          if(a == nages-1) 
              NAA(0,a) = NAA(0,a-1)/(one + exp( -M_tx(0,nlengths-1) - Ft(0) )   );
              //NAA(0,a) = NAA(0,a-1)/(one + exp(-1.0*M_tx_mean - exp(logN1_pars(1))  ) );
              //WHAM: NAA(0,a) = NAA(0,a-1)/(one + exp(-Mta(0,a) - exp(logN1_pars(1)) * Fta(0,a)/Fta(0,a-1)) );
          else 
              NAA(0,a) = NAA(0,a-1)*exp(-M_tx_mean - F_tx_mean);
              //NAA(0,a) = NAA(0,a-1)*exp(-1.0*M_tx_mean - exp(logN1_pars(1)) );
              //WHAM: NAA(0,a) = NAA(0,a-1)*exp(-Mta(0,a) -  exp(logN1_pars(1)) * Fta(0,a)/Fta(0,a-1));
    };
    //};
      
    pred_NAA(0,a)=NAA(0,a);
  };
*/


  //error("Stop hereeeeeeeeeeeee"); //ok upto this line;

  //logRec(0)  
  //need something
  a=0;
  NAA(0,a) = exp(logN1_par);     //NAA(0,0)

  for(int Lind=0;Lind<nlengths;Lind++){
    NL(0,Lind,a)=NAA(0,a)*f_total(a,Lind,0);
  };

  for(int a=1; a<nages; a++) {
     for(int Lind=0;Lind<nlengths;Lind++){
       for(int xind=0; xind<nlengths; xind++) {
       NL(0,Lind,a)+=NL(0,xind,a-1)*ExpZ_tx_mean(xind)*pp(Lind,xind,a);
       };
       pred_NAA(m,a) += NL(m,Lind,a);
     };
     NAA(m,a)=pred_NAA(m,a);
  };
  
  for(int a=0; a<nages; a++) {
     for(int Lind=0;Lind<nlengths;Lind++){
        SpawnersL(0,Lind,a)=NL(0,Lind,a)*pred_Maturation(Lind)*ratio_female*exp(-mort_frac*Z_tx(0,Lind));  
        SpawnerBiomassL(0,Lind,a)=SpawnersL(0,Lind,a)*pred_Wt(Lind);
        Spawners(0,a)+=SpawnersL(0,Lind,a);
        SpawnerBiomass(0,a)+=SpawnerBiomassL(0,Lind,a);
        };
        SSB_y(0)+=SpawnerBiomass(0,a);
  };
  
  //std::cout<<"SSB_y(0): "<<SSB_y(0)<<std::endl;
  
  //R_{t+1} = alpha*SSB_t/(beta+SSB_t)
  recruits_pred(0)=(SR_alpha*SSB_y(0))/(SR_beta+SSB_y(0)); 
     //time index 0 in SSB_y = year 0
     //time index 0 in recruits_pred = year 1
  //logRec(0)=log(recruits_pred(0));  //see the objective functions below
 
  // ====== E: f_total(0,xind,0), f_total(0,xind,1), ..., f_total(0,xind,m), ..., f_total(0,xind,nyrs-1)
  // ==========================================================================================
 
 //std::cout<<"before the loop over years"<<std::endl; 

 //loop over years instead of cohorts
 for(int m=1;m<nyrs;m++) {
   
   //S: a=0 means age 1
   a=0;
   for(int xind=0; xind<nlengths; xind++) {
       NL(m,xind,a)=exp(logRec(m))*f_total(a,xind,m);  // length of logRec is nyrs-1, initial recruitment is already estimated
       NAA(m,a) += NL(m,xind,a);
       
       SpawnersL(m,xind,a)=NL(m,xind,a)*pred_Maturation(xind)*ratio_female;
       SpawnerBiomassL(m,xind,a)=SpawnersL(m,xind,a)*pred_Wt(xind);
       Spawners(m,a)+=SpawnersL(m,xind,a);
       SpawnerBiomass(m,a)+=SpawnerBiomassL(m,xind,a);
   };
   SSB_y(m)+=SpawnerBiomass(m,a);
   
   pred_NAA(m,a)=NAA(m,a); 
   //E: a=0 means age 1;
 
   for(int a=1; a<nages; a++) {
    
      if(a!=(nages-1)) {
         Type SumP=zero;
         for(int Lind=0; Lind<nlengths; Lind++) {
            p(Lind)=zero;
            for(int xind=0; xind<nlengths; xind++) {
               p(Lind)+=f_total(a-1,xind,m-1)*ExpZ_tx(m-1,xind)*pp(Lind,xind,a); //double check with JW and Soyun 
                                                                                 //about a-1, m-1, and  a,  
            };  
            SumP+=p(Lind); 
         };      
     
         for(int Lind=0; Lind<nlengths; Lind++) {
            f_total(a,Lind,m)=p(Lind)/SumP;
            for(int xind=0; xind<nlengths; xind++) {
              NL(m,Lind,a)+=NL(m-1,xind,a-1)*ExpZ_tx(m-1,xind)*pp(Lind,xind,a); 
            };
            pred_NAA(m,a) += NL(m,Lind,a);
            //NAA(m,a)+=NL(m,Lind,a); 
            //N_save(m,a)+=NL(m,Lind,a);
         
            SpawnersL(m,Lind,a)=NL(m,Lind,a)*pred_Maturation(Lind)*ratio_female;
            SpawnerBiomassL(m,Lind,a)=SpawnersL(m,Lind,a)*pred_Wt(Lind);
            Spawners(m,a)+=SpawnersL(m,Lind,a);
            SpawnerBiomass(m,a)+=SpawnerBiomassL(m,Lind,a);
         };
          
         NAA(m,a)=pred_NAA(m,a);
         //pred_N_re(m-1,a-1)=NAA(m,a); 
      }
      else if(a==(nages-1))  {  //for the last age plus group, A
         //S: Part 1 for making f_total(A,xind,m);
         Type SumP=zero;
         for(int Lind=0; Lind<nlengths; Lind++) {
            p(Lind)=zero;
            for(int xind=0; xind<nlengths; xind++) {
               p(Lind)+=f_total(a-1,xind,m-1)*ExpZ_tx(m-1,xind)*pp(Lind,xind,a); //double check with JW and Soyun 
                                                                                 //about a-1, m-1, and  a,  
            };  
            SumP+=p(Lind); 
         };      
     
         for(int Lind=0; Lind<nlengths; Lind++) 
            f_total(a,Lind,m)=p(Lind)/SumP;  
         //E: Part 1 for making f_total(A,xind,m);
        
         for(int Lind=0; Lind<nlengths; Lind++) {
            p_plus(Lind)=zero;
            for(int xind=0; xind<nlengths; xind++) {
               p_plus(Lind)+=f_total(a,xind,m-1)*ExpZ_tx(m-1,xind)*pp(Lind,xind,a); //a means A
            };  
         };      
    
         for(int Lind=0; Lind<nlengths; Lind++) {
           for(int xind=0; xind<nlengths; xind++) {
             NL(m,Lind,a)+=NL(m-1,xind,a-1)*ExpZ_tx(m-1,xind)*pp(Lind,xind,a) + NL(m-1,xind,a)*ExpZ_tx(m-1,xind)*pp(Lind,xind,a+1);
           }; 
            pred_NAA(m,a) += NL(m,Lind,a);
            //NAA(m,a)+=NL(m,Lind,a); 
            //N_save(m,a)+=NL(m,Lind,a);
            
            SpawnersL(m,Lind,a)=NL(m,Lind,a)*pred_Maturation(Lind)*ratio_female;
            SpawnerBiomassL(m,Lind,a)=SpawnersL(m,Lind,a)*pred_Wt(Lind);
            Spawners(m,a)+=SpawnersL(m,Lind,a);
            SpawnerBiomass(m,a)+=SpawnerBiomassL(m,Lind,a);
         };
         
         NAA(m,a)=pred_NAA(m,a);
         //pred_N_re(m-1,a-1)=NAA(m,a); 
      };
      SSB_y(m)+=SpawnerBiomass(m,a);
   }; //age a
   
   recruits_pred(m)=(SR_alpha*SSB_y(m))/(SR_beta+SSB_y(m));
       //Recall that time index m is different between SSB_y and recruits_pred;
       
   //logRec(m)=log(recruits_pred(m));  //until nyrs + 1;  //see the objective functions below
 }; //m 
 


 /*
 std::cout<<"NL(0,0,0), NL(0,22,0): "<<std::endl; 
 std::cout<<NL(0,0,0)<<" "<<NL(0,22,0)<<std::endl; 
 std::cout<<"NL(1,0,0), NL(1,22,0): "<<std::endl; 
 std::cout<<NL(1,0,0)<<" "<<NL(1,22,0)<<std::endl; 
 
 std::cout<<"NL(1,1,1): "<<std::endl; 
 std::cout<<NL(1,1,1)<<std::endl; 
 
 std::cout<<"NAA(0,0): "<<std::endl; 
 std::cout<<NAA(0,0)<<std::endl; 
   
 std::cout<<"NAA(1,1): "<<std::endl; 
 std::cout<<NAA(1,1)<<std::endl; 
 */
 
 //
 for(int m=0;m<nyrs;m++) {
	  for(int a=0;a<nages;a++) {
	     for(int xind=0;xind<nlengths;xind++) {
		       CNum=NL(m,xind,a)*(F_tx(m,xind)/Z_tx(m,xind))*(one - ExpZ_tx(m,xind));
		       CWt=CNum*pred_Wt(xind);              //in kg
		       Catch(m,xind)+=CNum;                 //Catch(m,xind)+=CNum; 
		       TCatch(m)+=CNum; 
		       Yieldhat(m)+=CWt;
		       Yield3D(m,xind,a) = CWt;
		       BL(m,xind,a)=NL(m,xind,a)*pred_Wt(xind); 
		       BDY(m,xind,a)=0.5*(BL(m,xind,a)+(BL(m,xind,a)-Yield3D(m,xind,a)));
		       BAA(m,a)+=NL(m,xind,a)*pred_Wt(xind); 
		       B(m)+=NL(m,xind,a)*pred_Wt(xind);
	     };
    };  
 };  
 
 //std::cout<<"CNum: "<<CNum<<std::endl; 
 //std::cout<<"Catch: "<<Catch<<std::endl; 

 //For predicted cpue; //The cpue data for the chub mackerel fish were from the LPS fishery;
 Type q = q_lower+(q_upper-q_lower)/(1+exp(-logit_q_domain));
 for(int xind=0; xind<nlengths; xind++)
    QL(xind) = q*select_length(xind);   //cpue data are also from the LPS fishery;

 for(int m=0;m<nyrs;m++) 
	  for(int a=0;a<nages;a++) 
	     for(int xind=0;xind<nlengths;xind++) {
           predcpue_3d(m,xind,a)=BDY(m,xind,a)*QL(xind);
           predcpue_ya(m,a)+=predcpue_3d(m,xind,a);
	     };         
 predcpue=predcpue_ya.rowwise().sum(); 
 
 //====================================================================================
 //====================================================================================
 //The objective functions
 vector<Type> nll(8);     //elements of the objective funtion, which is the negative loglikelihood;
 nll.setZero(); 
 
 Type nll0 = zero,  nll1=zero, nll2=zero, nll3=zero, nll4=zero, nll5=zero, nll6=zero;

 //length-composition data
 matrix<Type> freq_out(nyrs, nlengths);
 matrix<Type> freq_obs(nyrs, nlengths);
 matrix<Type> freq_pred(nyrs, nlengths);
 
 if(lengthcomp_code==0)  {  //multinomial
   for(int m=0;m<nyrs;m++) {
     vector<Type> prop_obs = data_length_freq.row(m)/data_length_freq.row(m).sum();
     vector<Type> prop_pred = Catch.row(m)/Catch.row(m).sum();
     vector<Type> freq=prop_obs*SS_length(m);
     
     freq_out.row(m) = freq;
     freq_obs.row(m) = prop_obs;
     freq_pred.row(m) = prop_pred;
     
     
     nll(0)-= dmultinom(freq, prop_pred, true);
    
     Type temsum1=zero;
     Type temsum2=zero;
     for(int il=0; il<nlengths; il++) {
        temsum1 += prop_pred(il)*(1-prop_pred(il));
        temsum2 += square(prop_obs(il) - prop_pred (il) );
     };  
     Neff_lengthcomp(m)=temsum1/temsum2;
     //Neff_lengthcomp(m) = SS_length(m);
   };
 }
 else {  //Dirichlet-multinomial
   PARAMETER(logtheta);  //Dirichlet-Multinomial;
   Type theta=exp(logtheta); 
   for(int m=0;m<nyrs;m++)   {
     vector<Type> prop_obs = data_length_freq.row(m)/data_length_freq.row(m).sum();
     vector<Type> prop_pred = Catch.row(m)/Catch.row(m).sum();
   
     nll(0)-=logDiriMultinom_f(m, nlengths, SS_length(m), theta, prop_obs, prop_pred); 
    
     Neff_lengthcomp(m)=(one+theta*SS_length(m))/(one+theta); 
     dispersion_lengthcomp(m)=SS_length(m)*theta;
   };
 };
 nll0=nll(0);
 std::cout<<"nll(0): "<<nll(0)<<std::endl; 
 //std::cout<<"prop_obs: "<<prop_obs<<std::endl; 

 /*
 //likelihood for censored catch data
 DATA_SCALAR(CV_yield_cen);
 Type sig2_yield_cen =log(square(CV_yield_cen)+one  );  //yield~lognormal
 Type sig_yield_cen=sqrt(sig2_yield_cen); //Type sig_yield=exp(logsig_yield); 
 Type sig2_yield =log(square(CV_yield)+one  );  //yield~lognormal
 Type sig_yield=sqrt(sig2_yield); //Type sig_yield=exp(logsig_yield); 
 //note that pnorm below is the cumulatie probability  
 for(int m=0;m<8;m++)
    nll(1)-= log( pnorm(log(yield(m)), log(Yieldhat(m)/Type(1000)), sig_yield_cen)
                    -pnorm(log(yield(m)*0.9), log(Yieldhat(m)/Type(1000)), sig_yield_cen) + 1.0e-30); 
 for(int m=8;m<nyrs;m++)
    nll(1)-= dnorm(log(yield(m)), log(Yieldhat(m)/Type(1000)), sig_yield, true); 
 nll1=nll(1);
 std::cout<<"nll(1): "<<nll(1)<<std::endl; 
*/

 //lognormal for yield data (traditional way)
 //Type sig2_yield =log(square(CV_yield)+one  );  //yield~lognormal
 //Type sig_yield=sqrt(sig2_yield); //Type sig_yield=exp(logsig_yield); 
 for(int m=0;m<nyrs;m++) 
    nll(1)-= dnorm(log(yield(m)), log(Yieldhat(m)/Type(1000)), sig_yield, true);
 nll1=nll(1);
 std::cout<<"nll(1): "<<nll(1)<<std::endl; 

 //lognormal for cpue data
 //log(cpue) ~ normal(log(q*Bt), sig2_logcpue);
 //Type sig2_cpue =log(square(CV_cpue)+one  );  //cpue~lognormal
 //Type sig_cpue=sqrt(sig2_cpue);    //Type sig_cpue=exp(logsig_cpue);
 for(int m=0;m<nyrs;m++) 
    nll(2)-= dnorm(log(cpue(m)), log(predcpue(m)), sig_cpue, true);
 nll2=nll(2);
 std::cout<<"nll(2): "<<nll(2)<<std::endl; 

/*
  //binomial for maturation data
 for(int i=0;i<nlengths;i++){
   nll(3)-=dbinom(Type(Matdat(i)), Type(n_obs(i)), pred_Maturation(i), true);
 };
 */
 nll(3)=0.0;
 nll3=nll(3);
 std::cout<<"nll(3): "<<nll(3)<<std::endl; 

 //process of recruits as random effects
 //nll(4)-=dnorm(logRec(0),log(NAA(0,0)), exp(log_sig_logNAA),true); 
 nll(4)-=dnorm(logRec(0),log(NAA(0,0)), exp(log_sig_logRec),true); 
 for(int m=1;m<nyrs;m++)
     //nll(4)-=dnorm(logRec(m), log(recruits_pred(m-1)), exp(log_sig_logNAA),true);
     nll(4)-=dnorm(logRec(m), log(recruits_pred(m-1)), exp(log_sig_logRec),true);
 nll4=nll(4);
 std::cout<<"nll(4): "<<nll(4)<<std::endl; 
  
/*  
 //process of N_{y,a} as random effects; 
 for(int m=0; m<nyrs-1; m++)  {
   for(int a=0; a<nages-1; a++)
      nll(5)-=dnorm(logNAA(m,a), log(pred_NAA(m+1,a+1)+1.0e-15), log_sig_logNAA, true);
 };    
 nll5=nll(5);
 std::cout<<"nll(5): "<<nll(5)<<std::endl; 
*/

 //When fishing effort is missing;
 //deviataions around the fishing mortality;
 nll(6)-=dnorm(logFt(0), logF1, exp(log_sig_logF), true);
 for(int t=1; t<nyrs; t++) 
    nll(6)-=dnorm(logFt(t), logFt(t-1), exp(log_sig_logF), true);  
   
 nll6=nll(6);
 std::cout<<"nll(6): "<<nll(6)<<std::endl; 
 
 //Informative prior for sig_L; //sigo_Y ~  gamma 
 Type alpha_sig_L=Shaperate_sig_L(0);  
 Type beta_sig_L=Shaperate_sig_L(1); 
 nll(7)-= dgamma(exp(log_sig_L), alpha_sig_L, (1/beta_sig_L), true );
 
 Type jnll=nll.sum();

 

 
 //REPORT
 REPORT(nll0);
 REPORT(nll1);
 REPORT(nll2);
 REPORT(nll3);
 REPORT(nll4);
 REPORT(nll5);
 REPORT(nll6);

 REPORT(jnll);
 //REPORT(no_parameters); 
 REPORT(SR_alpha); 
 REPORT(SR_beta);
 REPORT(SSB_y);
 REPORT(recruits_pred);
 REPORT(logRec);
 REPORT(predcpue);
 REPORT(Catch); 
 REPORT(Yieldhat); 
 REPORT(NAA);
 REPORT(BAA);
 REPORT(B);
 REPORT(Spawners); 
 REPORT(SpawnerBiomass);
 REPORT(NL); 
 REPORT(select_length);
 REPORT(QL);
 REPORT(Mu);
 REPORT(f_total);
 REPORT(pp);
 REPORT(M_tx);
 REPORT(Ft); 
 REPORT(dispersion_lengthcomp);
 REPORT(Neff_lengthcomp);
 REPORT(freq_out);
 REPORT(freq_obs);
 REPORT(freq_pred);
 REPORT(Rho);
 REPORT(sig2_lengths_r);
 REPORT(SS);
 REPORT(Mu);
 //REPORT(log_sig_logNAA);
 REPORT(log_sig_logRec);
 REPORT(q);
 REPORT(logit_q_domain);
 REPORT(SSB0);
 REPORT(Catch);
 REPORT(pred_Wt);
 REPORT(pred_Maturation);
 REPORT(SPR_0);
 REPORT(kappa);
 REPORT(log_sig_L);
 REPORT(logR0);
 REPORT(sig_yield);
 REPORT(sig_cpue);
 REPORT(sig2_lengths_r);
 //REPORT(b0_mat);
 //REPORT(b1_mat);
 
 ADREPORT(q);
 ADREPORT(kappa);
 ADREPORT(R0);

 // ========================================================================================
 // ========================================================================================
 /////////////////////////////////////// YPR & SPR /////////////////////////////////////////
  
 array<Type> p_ref(nages,nlengths,indnum_F);
 p_ref.setZero();
 matrix<Type> p_refA(nlengths,indnum_F); // plus age group
 p_refA.setZero();
 matrix<Type> Zref(indnum_F,nlengths);
 Zref.setZero();
 matrix<Type> ExpZ_ref(indnum_F,nlengths);
 ExpZ_ref.setZero();

 // p_ref: Survival and growth from the cell at age 1 to the last age class along 
 //        one cohort in the terminal year under the equilibrium assumption
 //        This is comparable to part in red in JWs summary of SPR 
 for(int ind_F=0;ind_F<indnum_F;ind_F++){
    for(int xind=0;xind<nlengths;xind++){
      Zref(ind_F,xind)=F_input(ind_F)*select_length(xind)+M_tx(nyrs-1,xind);
      ExpZ_ref(ind_F,xind)=exp(Type(-1.0)*Zref(ind_F,xind));
    };
 };

 for(int ind_F=0;ind_F<indnum_F;ind_F++){
   for(int Lind=0; Lind<nlengths; Lind++) {
    p_ref(0,Lind,ind_F) = f_total(0,Lind,nyrs-1); // terminal year;
   };
 };

 for(int ind_F=0;ind_F<indnum_F;ind_F++){
   for(int iage=1; iage<nages; iage++) {
     for(int Lind=0; Lind<nlengths; Lind++) {
       p_ref(iage,Lind,ind_F)=zero;
       for(int xind=0; xind<nlengths; xind++) {
          p_ref(iage,Lind,ind_F)+=p_ref(iage-1,xind,ind_F)*ExpZ_ref(ind_F,xind)*pp(Lind,xind,iage); //double check with JW and Soyun 
       };
     }; 
   };
 };

 for(int ind_F=0;ind_F<indnum_F;ind_F++){
   for(int Lind=0; Lind<nlengths; Lind++) {
     p_refA(Lind,ind_F)=zero;
     for(int xind=0; xind<nlengths; xind++) {  
        p_refA(Lind,ind_F)=p_ref(nages-1,xind,ind_F)*ExpZ_ref(ind_F,xind)*pp(Lind,xind,nages-1); // ???
     };
   };
 };

 for(int ind_F=0;ind_F<indnum_F;ind_F++){
   for(int Lind=0; Lind<nlengths; Lind++) {
     p_ref(nages-1,Lind,ind_F)+=p_refA(Lind,ind_F);
   };
 };

 ///
 matrix<Type> YPR_mat(indnum_F,nlengths);
 vector<Type> YPR(indnum_F);

 for(int ind_F=0;ind_F<indnum_F;ind_F++){
   for(int Lind=0;Lind<nlengths;Lind++){
      YPR_mat(ind_F,Lind)=p_ref(0,Lind,ind_F)*pred_Wt(Lind)*(F_input(ind_F)*select_length(Lind)/Zref(ind_F,Lind))*(one-exp(-Zref(ind_F,Lind)))+
                          p_ref(1,Lind,ind_F)*pred_Wt(Lind)*(F_input(ind_F)*select_length(Lind)/Zref(ind_F,Lind))*(one-exp(-Zref(ind_F,Lind)))+
                          p_ref(2,Lind,ind_F)*pred_Wt(Lind)*(F_input(ind_F)*select_length(Lind)/Zref(ind_F,Lind))*(one-exp(-Zref(ind_F,Lind)))+
                          p_ref(3,Lind,ind_F)*pred_Wt(Lind)*(F_input(ind_F)*select_length(Lind)/Zref(ind_F,Lind))*(one-exp(-Zref(ind_F,Lind)))+  
                          p_ref(4,Lind,ind_F)*pred_Wt(Lind)*(F_input(ind_F)*select_length(Lind)/Zref(ind_F,Lind))*(one-exp(-Zref(ind_F,Lind)))+ 
                          p_ref(5,Lind,ind_F)*pred_Wt(Lind)*(F_input(ind_F)*select_length(Lind)/Zref(ind_F,Lind))*(one-exp(-Zref(ind_F,Lind)));
   };
 };  

 YPR=YPR_mat.rowwise().sum();  

 /////  
 matrix<Type> SPR_T_matrix(indnum_F,nlengths);
 SPR_T_matrix.setZero();
 vector<Type> SPR_T(indnum_F);
 SPR_T.setZero();

 vector<Type> spawnpotenratio(indnum_F);  //spawning potential ratio = SPR_T/SPR0   
 spawnpotenratio.setZero();

 for(int ind_F=0;ind_F<indnum_F;ind_F++){
   for(int Lind=0;Lind<nlengths;Lind++){
      SPR_T_matrix(ind_F,Lind)=p_ref(0,Lind,ind_F)*pred_Wt(Lind)*pred_Maturation(Lind)*ratio_female*exp(-mort_frac*Zref(ind_F,Lind))+ // add 231222
                               p_ref(1,Lind,ind_F)*pred_Wt(Lind)*pred_Maturation(Lind)*ratio_female*exp(-mort_frac*Zref(ind_F,Lind))+
                               p_ref(2,Lind,ind_F)*pred_Wt(Lind)*pred_Maturation(Lind)*ratio_female*exp(-mort_frac*Zref(ind_F,Lind))+
                               p_ref(3,Lind,ind_F)*pred_Wt(Lind)*pred_Maturation(Lind)*ratio_female*exp(-mort_frac*Zref(ind_F,Lind))+
                               p_ref(4,Lind,ind_F)*pred_Wt(Lind)*pred_Maturation(Lind)*ratio_female*exp(-mort_frac*Zref(ind_F,Lind))+
                               p_ref(5,Lind,ind_F)*pred_Wt(Lind)*pred_Maturation(Lind)*ratio_female*exp(-mort_frac*Zref(ind_F,Lind));
   };
 };

 SPR_T=SPR_T_matrix.rowwise().sum(); // sum over length??? // length of ind_F; function of F
 
 spawnpotenratio=SPR_T/SPR_0; // 
 
 //std::cout<<"YPR: "<<YPR<<std::endl;
 //std::cout<<"logRec(23): "<<logRec(23)<<std::endl;
 //error("Stop hereeeeeeeeeeeee");
 
 //Type temsumRec=exp(logRec(17))+ exp(logRec(18)) + exp(logRec(19)) + exp(logRec(20)) + exp(logRec(21)) + exp(logRec(22));
 //Type temmeanRec=temsumRec/6;
 
 //Type MSY=max(YPR)*temmeanRec;  
 
 vector<Type> SSB_F(indnum_F);
 SSB_F.setZero();
 vector<Type> Rec_F(indnum_F);
 Rec_F.setZero();
 vector<Type> Yield_F(indnum_F);
 Yield_F.setZero();
 Type MSY;
 
 //SR_alpha=Type(4.0)*h*R0/((5.0*h)-1.0);
 //SR_beta=SSB0*(1.0-h)/(5*h-1.0);
 
 SSB_F=SR_alpha*SPR_T-SR_beta;
 Rec_F=SSB_F/SPR_T;
 
 Yield_F=YPR*Rec_F; //yield;
 MSY=max(Yield_F);  //msy
 
 //REPORT
 REPORT(YPR);
 REPORT(spawnpotenratio);
 REPORT(Yield_F)
 REPORT(MSY);
 REPORT(Rec_F);
 REPORT(SSB_F);
 //  
 return jnll;
 }  
'

#
write(sbmack, file="sbmack.cpp"); 
compile("sbmack.cpp"); 
#compile("sbmack.cpp", "&> out.txt"); #note this 
dyn.load(dynlib("sbmack"));