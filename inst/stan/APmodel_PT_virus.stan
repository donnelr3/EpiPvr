// Description: Probability model of PT access period assay used in Bayesian analysis, see Donnelly Tankam Gilligan 2024 appendix S1 
// data : access period data NOTE that we assume there are acquisition access latent access and inoculation access components for PT virus
// when all 3 subassays are not available the dataset can simply pass 0 e.g. if there is no LAP then: 
// inputdata$D_LensL=0,inputdata$D_RepsL=0,inputdata$D_InfsL=0 can be passed as part of inputdata
//
//  stan file is used to estimate viral parameters: alpha beta gamma mu from inputdata 
// the probability model is formulated conditional on discretiseed insect lab survival 
// i.e. insect mortality is additionally estimated as a nuisance parameter
// user-inputted guess for upper bound on lab insect survival is required (D_lsPars) to structure the discretisation 
// if none is provided estimation will proceed based on upper bound of 50d


functions {  // 2xseries expansions and a linspace
real special_taylor(real pivotMnode, real t_in, int precis_in, int orderExpand){
  real tmp;
  real output;
  if (fabs(pivotMnode)<(10^(-precis_in))){
    tmp=0.0;
    for (ss in 1:orderExpand){
      tmp=tmp+((((pivotMnode)*t_in)^(ss-1))/tgamma(ss+1));
    }
    output=1*t_in*tmp;
  }else{
    output=(1/pivotMnode)*(exp(pivotMnode*t_in)-1);
  }
  return output;
}
real geom_series(real pivot2_node2, int precis_in, int orderExpand){
  real tmp;
  real output;
  if (fabs(pivot2_node2-1)<(10^(-precis_in))){
    tmp=0.0;
    for (ss in 1:orderExpand){
      tmp=tmp+((pivot2_node2)^(ss-1));
    }
    output=tmp;
  }else{
    output=1/(1-(pivot2_node2));
  }
  return output;
}
real[] custom_linspace(int n, real start, real end) {
  real result[n];
  real st= (end - start) / (n - 1);
  //st= (end - start) / (n - 1);
  //real step = (end - start) / (n - 1);
  for (i in 1:n) {
    result[i] = start + (i - 1) * st;
  }
  return result;
}
}

// DATA DATA DATA DATA DATA DATA DATA DATA DATA DATA DATA 
data {                                    // Data block
int<lower=0> D_NumGrps[3];  // number of assay timepoints in AAP LAP and IAP subassays
int<lower=1> D_Wf0;        // number of whitefly used (defined in terms of # brought to IAP)
real<lower=0> D_LensA[D_NumGrps[1]];        // timepoints for each subassay 
real<lower=0> D_LensL[D_NumGrps[2]];
real<lower=0> D_LensI[D_NumGrps[3]];
int<lower=0> D_RepsA[D_NumGrps[1]];        // Number of reps for each subassay
int<lower=0> D_RepsL[D_NumGrps[2]]; 
int<lower=0> D_RepsI[D_NumGrps[3]]; 
int<lower=0> D_InfsA[D_NumGrps[1]];        // Number of infected test plants for each subassay
int<lower=0> D_InfsL[D_NumGrps[2]];  
int<lower=0> D_InfsI[D_NumGrps[3]];  
real<lower=-1> D_bgLens[3,3];    // matrix of durations row 1 AAP subassay (col 1 duration acquisition access col 2 duration inoculation access) row 2 IAP subassay NB diagonals are ignored as those durations are varied in each subassay
int<lower=0> D_lsPars[2];        // [2] rough upper-bound on insect lifespan (days); and [1] number of timepoints for each day (for survival discretisation)
int<lower=0> D_ppc;              // ignore, testing mode, defaults to 1 in estimate_virus_parameters_PT
}

// TRANSFORMED DATA TRANSFORMED DATA TRANSFORMED DATA 
transformed data {                                    // Data block

}


// PARAMETERS PARAMETERS PARAMETERS PARAMETERS PARAMETERS 
parameters {                             // Parameters block

real<lower=0> lat[1]; // latent progression rate
real<lower=1/(24*(D_lsPars[1])),upper=(60*6)> bD[1]; // insect mortality - bounded by user supplied rough survival upper bound and default 10s lower bound 
real<lower=0> al[1]; // acquisition rate
real<lower=0> be[1]; // inoculation rate
real<lower=0> mu_lat[1];
}


// TRANSFORMS TRANSFORMS TRANSFORMS TRANSFORMS TRANSFORMS
transformed parameters {
  
  real<lower=0> mu[1]; // insect recovery rate
  real<lower=0,upper=1> DD[1]; // composite parameter (proportion)
  real par1[1];
  real par3[1];
  real par6[1];
  
  real binPam_succA[D_NumGrps[1]];
  real binPam_succL[D_NumGrps[2]];
  real binPam_succI[D_NumGrps[3]];
  
  mu[1]=mu_lat[1]*lat[1];
  DD[1]=be[1]/(mu[1]+be[1]);  
  
  par1[1]=lat[1]-al[1];
  par3[1]=lat[1]-(mu[1]+be[1]);
  par6[1]=mu[1]-al[1];
  
  binPam_succA = rep_array(0.0,D_NumGrps[1]); 
  binPam_succL = rep_array(0.0,D_NumGrps[2]); 
  binPam_succI = rep_array(0.0,D_NumGrps[3]);
  
  {
    int numPts=(D_lsPars[1]*D_lsPars[2])+1; // temporary parameters governing discretised insect survival
    real TT[numPts] = custom_linspace(numPts, 0.0, D_lsPars[1] * 1.0);
    real ls_int_probs[numPts];
    real mn_TT[numPts];
    real likethresh=0.001/D_lsPars[2];
    int count=1;
    
    real termA;
    real termB;
    real temp1;
    real temp2;
    real temp3;
    real temp5;
    real temp6;
    real tempgeo;
    
    real TAAP_assay; // access period durations
    real TLAP_assay;
    real TIAP_assay;
    
    int orderOfExpand=6; // expansion settings
    int precisForSing=30;
    
    real binPam_totsA[D_NumGrps[1],numPts];
    real binPam_totsL[D_NumGrps[2],numPts];
    real binPam_totsI[D_NumGrps[3],numPts];
    
    binPam_totsA = rep_array(0.0,D_NumGrps[1],numPts); 
    binPam_totsL = rep_array(0.0,D_NumGrps[2],numPts); 
    binPam_totsI = rep_array(0.0,D_NumGrps[3],numPts);  
    
    // iterating through discretised survival window until window prob is too small
    while (((exp(-TT[count]*24.0*bD[1])-exp(-TT[count+1]*24.0*bD[1]))>likethresh) && count<(numPts-1)) 
    {      
      ls_int_probs[count]=(exp(-TT[count]*24.0*bD[1])-exp(-TT[count+1]*24.0*bD[1]));
      mn_TT[count]=TT[count+1]*24.0;
      count=count+1;
    }
     // final discretisation window 
    ls_int_probs[count]=(exp(-TT[count]*24.0*bD[1])); // final time window stetches to infinity 
    mn_TT[count]=(TT[count])*24.0;
    
        for (dd in 1:3) //subassays 
        {
          TAAP_assay=0.0;
          TLAP_assay=0.0;
          TIAP_assay=0.0;    
          for (kk in 1:(D_NumGrps[dd])) //timepoint sets in subassays
          {
            for (hh in 1:(count-1)) //discretised windows
            {
              if (dd==1){ // setting the access period durations
                TAAP_assay = D_LensA[kk];
                TLAP_assay = TAAP_assay+D_bgLens[1,2];
                TIAP_assay = TLAP_assay+min([mn_TT[hh],D_bgLens[1,3]]); // IAP insect experiences is bound by the IAP duration by may be shorter (survival less than IAP duration)
              }else if (dd==2){
                TAAP_assay = D_bgLens[2,1];
                TLAP_assay = TAAP_assay+D_LensL[kk];
                TIAP_assay = TLAP_assay+min([mn_TT[hh],D_bgLens[2,3]]);
              }else{
                TAAP_assay = D_bgLens[3,1];
                TLAP_assay = TAAP_assay+D_bgLens[3,2];
                TIAP_assay = TLAP_assay+min([mn_TT[hh],D_LensI[kk]]);
              }
              // expansion terms relating to hypo-exponential distributions (modelled as expansions instead)
              temp1=special_taylor(par1[1], TAAP_assay, precisForSing, orderOfExpand);
              temp2=special_taylor(lat[1], (TIAP_assay-TLAP_assay), precisForSing, orderOfExpand);
              temp3=special_taylor(par3[1], (TIAP_assay-TLAP_assay), precisForSing, orderOfExpand);
              temp5=special_taylor((mu[1]+be[1]), (TLAP_assay-TIAP_assay), precisForSing, orderOfExpand);
              temp6=special_taylor(par6[1], TAAP_assay, precisForSing, orderOfExpand);
              tempgeo=geom_series(mu_lat[1], precisForSing, orderOfExpand);
              termA=tempgeo*temp5*((exp(-(lat[1])*(TLAP_assay))*temp1)-(exp(-(mu[1])*(TLAP_assay))*temp6));
              termB=temp1*DD[1]*exp(-(lat[1])*(TIAP_assay))*(temp2-temp3);
              
              if (dd==1){
                binPam_totsA[kk,hh]=al[1]*((termA*be[1])+(termB*lat[1]))*ls_int_probs[hh];
              }else if (dd==2){
                binPam_totsL[kk,hh]=al[1]*((termA*be[1])+(termB*lat[1]))*ls_int_probs[hh];
              }else{
                binPam_totsI[kk,hh]=al[1]*((termA*be[1])+(termB*lat[1]))*ls_int_probs[hh];
              }
            }
            
            if (dd==1){ 
              binPam_succA[kk]=1-(1-sum(to_vector(binPam_totsA[kk,])))^D_Wf0;
            }else if (dd==2){
              binPam_succL[kk]=1-(1-sum(to_vector(binPam_totsL[kk,])))^D_Wf0;
            }else{
              binPam_succI[kk]=1-(1-sum(to_vector(binPam_totsI[kk,])))^D_Wf0;
            }
          }
        }
  }
}

// MODEL MODEL MODEL MODEL MODEL MODEL MODEL MODEL MODEL 
model {     
 
  // half-normal priors and 1 beta prior for ratio of exponential rates
  lat[1] ~ normal(0,1);   
  al[1] ~ normal(0,1);  
  be[1] ~ normal(0,1); 
  mu_lat[1] ~ beta(1, 5); 
  
  if (!D_ppc)
  {
    for (kk in 1:D_NumGrps[1]) // AAP sub-assay
      target += binomial_lpmf(D_InfsA[kk] | D_RepsA[kk], binPam_succA[kk]);
    
    for (kk in 1:D_NumGrps[2]) // LAP sub-assay
      target += binomial_lpmf(D_InfsL[kk] | D_RepsL[kk], binPam_succL[kk]);
    
    for (kk in 1:D_NumGrps[3])
      target += binomial_lpmf(D_InfsI[kk] | D_RepsI[kk], binPam_succI[kk]);
  }
}
////////////////////////////////////////////////////////////////////  


// GENERATE GENERATE GENERATE GENERATE GENERATE GENERATE 

generated quantities {    

real<lower=0> y_simul_A[D_NumGrps[1]];           
real<lower=0> y_simul_L[D_NumGrps[2]];           
real<lower=0> y_simul_I[D_NumGrps[3]];           

for (kk in 1:D_NumGrps[1]) //D_NumGrps
  y_simul_A[kk] = binomial_rng(D_RepsA[kk], binPam_succA[kk]); // Simulate from likelihood

for (kk in 1:D_NumGrps[2]) //D_NumGrps
  y_simul_L[kk] = binomial_rng(D_RepsL[kk], binPam_succL[kk]); // Simulate from likelihood

for (kk in 1:D_NumGrps[3]) //D_NumGrps
  y_simul_I[kk] = binomial_rng(D_RepsI[kk], binPam_succI[kk]); // Simulate from likelihood

}
