//
// Functions to dynmically set priors

// Should probably consolidate to a few generic prior functions
// PriorMatrix, PriorVector, PriorMatrixByCol etc.

real ThetaPriors(matrix theta, matrix priorParam1, matrix priorParam2, vector priorType){
  real out;
  out = 0;
  if(num_elements(priorType) != cols(theta)){
    reject("Number of transformations not equal to number of theta parameters");
  }
  /// may want to investigate other priors. Easily extendable to anything < 2 params.
  for(i in 1:cols(theta)){
    if(priorType[i] == 0){
      // no additional prior information

      out += 0;
    }else if(priorType[i] == 1){
      // uniform prior
      out += uniform_lpdf(theta[,i] | priorParam1[,i], priorParam2[,i]);
    }else if(priorType[i] == 2){
      // a normal prior
      out += normal_lpdf(theta[,i] | priorParam1[,i], priorParam2[,i]);
    }else{
      reject("prior type is not valid");
    }
  }
  return(out);
}

//not implemented yet
real ObsVarPriors(real ObsVar, real priorParam1, real priorParam2, int priorType){
  real out;
  out = 0;
  /// may want to investigate other priors. Easily extendable to anything < 2 params.
  if(priorType == 0){
    // no additional prior information
    reject("Prior not allowed");
    out += 0;
  }else if(priorType == 1){
    // uniform prior
    reject("Prior not allowed");
    out += uniform_lpdf(ObsVar | priorParam1, priorParam2);
  }else if(priorType == 2){
    // a normal prior
    out += normal_lpdf(ObsVar | priorParam1, priorParam2);
  }else if(priorType == 3){
    out += inv_gamma_lpdf(ObsVar | priorParam1, priorParam2);
  }else{
    reject("prior type is not valid");
  }
  return(out);
}

//not implemented yet
real SigmaThetaPriors(vector SigmaTheta, vector priorParam1, vector priorParam2, vector priorType){
  real out;
  out = 0;
  /// may want to investigate other priors. Easily extendable to anything < 2 params.
  if(num_elements(SigmaTheta) != num_elements(priorType)){
    reject("Number of prior types not equal number of Sigma Theta parameters");
  }
  for(i in 1:num_elements(SigmaTheta)){
    if(priorType[i] == 0){
      // no additional prior information
      out += 0;
    }else if(priorType[i] == 1){
      // uniform prior, needs check that lower bound is greater than 0
      out += uniform_lpdf(SigmaTheta[i] | priorParam1[i], priorParam2[i]);
    }else if(priorType[i] == 2){
      // a normal prior truncated at 0
      out += normal_lpdf(SigmaTheta[i]| priorParam1[i], priorParam2[i]) - normal_lccdf(0 | priorParam1[i], priorParam2[i]);
    }else if(priorType[i] == 3){
      out += inv_gamma_lpdf(SigmaTheta[i] | priorParam1[i], priorParam2[i]);
    }else{
      reject("prior type is not valid");
    }
  }

  return(out);
}
// not yet implemented
// just transpose beta if want priors by row
// rows would be covariate level, cols would be by theta parameter level
real BetaPriorCols(matrix betaTheta, matrix priorParam1, matrix priorParam2, vector priorType){
  real out;
  out = 0;
  /// may want to investigate other priors. Easily extendable to anything < 2 params.

  for(i in 1:rows(betaTheta)){
    if(priorType[i] == 0){
      // no additional prior information
      reject("Prior not allowed");
      out += 0;
    }else if(priorType[i] == 1){
      // uniform prior
      reject("Prior not allowed");
      out += uniform_lpdf(betaTheta[,i] | priorParam1[,i], priorParam2[,i]);
    }else if(priorType[i] == 2){
      // a normal prior
      out += normal_lpdf(betaTheta[,i] | priorParam1[,i], priorParam2[,i]);
    }else if(priorType[i] == 3){
      out += inv_gamma_lpdf(betaTheta[,i] | priorParam1[,i], priorParam2[,i]);
    }else{
      reject("prior type is not valid");
    }
  }
  return(out);
}

