// Description: Probability model of SPT access period assay used in Bayesian analysis, see Donnelly Tankam Gilligan 2024 appendix S2 
// data : access period data NOTE that we assume there are acquisition access and inoculation access components for SPT virus
//
//  stan file is used to estimate viral parameters: alpha beta mu from inputdata 
// the probability model is formulated conditional on discretised insect lab survival 
// i.e. insect mortality is additionally estimated as a nuisance parameter
// user-inputted guess for upper bound on lab insect survival is required (D_lsPars) to structure the discretisation 
// if none is provided estimation will proceed based on upper bound of 50d


functions {
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
int<lower=0> D_NumGrps[2];  // number of assay timepoints in AAP and IAP subassays
int<lower=1> D_Wf0;        // number of whitefly used (defined in terms of no. brought to IAP)
real<lower=0> D_LensA[D_NumGrps[1]];        // timepoints for each subassay 
real<lower=0> D_LensI[D_NumGrps[2]];
int<lower=0> D_RepsA[D_NumGrps[1]];        // Number of reps for each subassay
int<lower=0> D_RepsI[D_NumGrps[2]]; 
int<lower=0> D_InfsA[D_NumGrps[1]];       // Number of infected test plants for each subassay
int<lower=0> D_InfsI[D_NumGrps[2]];  
real<lower=-1> D_bgLens[2,2];            // matrix of durations row 1 AAP subassay (col 1 duration acquisition access col 2 duration inoculation access) row 2 IAP subassay NB diagonals are ignored as those durations are varied in each subassay
int<lower=0> D_lsPars[2];               // [2] rough upper-bound on insect lifespan (days); and [1] number of timepoints for each day (for survival discretisation)
}

// TRANSFORMED DATA TRANSFORMED DATA TRANSFORMED DATA 
transformed data {                                    // Data block

}


// PARAMETERS PARAMETERS PARAMETERS PARAMETERS PARAMETERS 
parameters {                             // Parameters block
// composite parameters that appear in the probability model
real<lower=0> c2[1]; 
real<lower=0> c3[1];
real<lower=0> c1[1]; 
real<lower=1/(24*(D_lsPars[1])),upper=(60*6)> bD[1];   // insect mortality - bounded by user supplied rough survival upper bound and default 10s lower bound 
}

// TRANSFORMS TRANSFORMS TRANSFORMS TRANSFORMS TRANSFORMS
transformed parameters {
  
  real binPam_succA[D_NumGrps[1]];
  real binPam_succI[D_NumGrps[2]];
  binPam_succA = rep_array(0.0,D_NumGrps[1]);  
  binPam_succI = rep_array(0.0,D_NumGrps[2]); 
  
  {  
    int numPts=(D_lsPars[1]*D_lsPars[2])+1; // temporary parameters governing discretised insect survival
    real TT[numPts] = custom_linspace(numPts, 0.0, D_lsPars[1] * 1.0);
    real ls_int_probs[numPts];
    real mn_TT[numPts];
    real likethresh=0.001/D_lsPars[2];
    int count=1;
    
    real pA_T_A[D_NumGrps[1]]; 
    real pA_T_I;
    real pB_T_A[numPts]; 
    real pB_T_I[D_NumGrps[2],numPts]; 
    real binPam_vecsA[D_NumGrps[1],numPts];
    real binPam_vecsI[D_NumGrps[2],numPts];
    real TAAP; 
    real TIAP;
    
    binPam_vecsA = rep_array(0.0,D_NumGrps[1],numPts);
    binPam_vecsI = rep_array(0.0,D_NumGrps[2],numPts); 
    
    // iterating trough discretised survival window until window prob is too small
    while (((exp(-TT[count]*24.0*bD[1])-exp(-TT[count+1]*24.0*bD[1]))>likethresh) && count<(numPts-1)) //D_NumGrps
    {      
      ls_int_probs[count]=(exp(-TT[count]*24.0*bD[1])-exp(-TT[count+1]*24.0*bD[1]));
      mn_TT[count]=TT[count+1]*24.0;
      count=count+1;
    }
    // final discretisation window 
    ls_int_probs[count]=(exp(-TT[count]*24.0*bD[1])); // final time window stetches to infinity
    mn_TT[count]=(TT[count])*24.0;
    
    // separate fitings for AAP and IAP subassays
    // AAP
    for (kk in 1:(D_NumGrps[1]))
    {
      TAAP = D_LensA[kk];
      pA_T_A[kk]=(1-exp(-c2[1]*TAAP));
      for (hh in 1:(count-1)) 
      {
        TIAP = min([mn_TT[hh],D_bgLens[1,2]]);
        pB_T_A[hh]=(1-exp(-c3[1]*TIAP));
        binPam_vecsA[kk,hh]=c1[1]*pA_T_A[kk]*pB_T_A[hh]*ls_int_probs[hh]; 
      }
      binPam_succA[kk]= 1-(1-sum(to_vector(binPam_vecsA[kk,])))^D_Wf0;
    }
    // IAP
    for (ss in 1:(D_NumGrps[2])) 
    {
      TAAP = D_bgLens[2,1];
      pA_T_I=(1-exp(-c2[1]*TAAP));
      for (hh in 1:(count-1)) 
      {
        TIAP = min([mn_TT[hh],D_LensI[ss]]);
        pB_T_I[ss,hh]=(1-exp(-c3[1]*TIAP));
        binPam_vecsI[ss,hh]=c1[1]*pA_T_I*pB_T_I[ss,hh]*ls_int_probs[hh];
      }
      binPam_succI[ss]= 1-(1-sum(to_vector(binPam_vecsI[ss,])))^D_Wf0;
    }
  }
  
}

// MODEL MODEL MODEL MODEL MODEL MODEL MODEL MODEL MODEL 
model {     
  
  // half normal prior distributions (becuase compound parameter declarations restrict to non-zero positives)
  c2[1] ~ normal(0,1); 
  c3[1] ~ normal(0,1);  
  c1[1] ~ normal(0,1);  
  
  // Finally the binommial probability model AAP subassay
  for (kk1 in 1:D_NumGrps[1]) //D_NumGrps  
  target += binomial_lpmf(D_InfsA[kk1] | D_RepsA[kk1],binPam_succA[kk1]);
  
  // Finally the binommial probability model IAP subassay
  for (uu in 1:D_NumGrps[2])
  target += binomial_lpmf(D_InfsI[uu] | D_RepsI[uu], binPam_succI[uu]);
  
}
////////////////////////////////////////////////////////////////////  




// GENERATE GENERATE GENERATE GENERATE GENERATE GENERATE 

generated quantities {      // Generated quantities block. 
// extracting alpha, beta and mu from the composite parameter posteriors 
real<lower=0> albe[1]; // acquisition
real<lower=0> al[1]; // acquisition
real<lower=0> be[1]; // inoculation
real mu[1]; // virus loss
real XX; // composite param

real<lower=0> y_simul_A[D_NumGrps[1]];           
real<lower=0> y_simul_I[D_NumGrps[2]];           

for (kk in 1:D_NumGrps[1]) 
y_simul_A[kk] = binomial_rng(D_RepsA[kk], binPam_succA[kk]); // Simulate from likelihood


for (uu in 1:D_NumGrps[2]) 
y_simul_I[uu] = binomial_rng(D_RepsI[uu], binPam_succI[uu]); // Simulate from likelihood

albe[1]=c1[1]*c2[1]*c3[1];

XX=c2[1]-c3[1]; 

be[1]=(sqrt((XX*XX)+(4*albe[1]))-XX)/2; // nb sqrt will return the positive root ensuring be[1]>0 

al[1]=XX+be[1];

mu[1]=c2[1]-al[1];

}

