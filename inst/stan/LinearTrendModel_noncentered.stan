/// linear trend model

functions{

#include /functions/userDefinedGrowthFunction.stan
#include /functions/userDefinedThetaTransformations.stan
#include /functions/Priors.stan
#include /functions/growthCurveMeans.stan
#include /functions/thetaTransformations.stan

}
data{
  int<lower = 1> N; // Total number of observations. Sum T_u over u.
  int<lower = 1> d; // Num parameters
  // order by year/group then week of the year
  int u[N]; // Index of what year/group observation is from.  presumably ordered.
  real tpt[N]; // tpt for each week presumably ordered.
  int U; //Number of groups/years
  int pTheta; // number of columns in XTheta
  int nrowXTheta; // number of rows in design matrix. Should be = U


  //data
  vector[N] Z;
  matrix[U, pTheta] XTheta;
  int<lower = 0, upper = 1> hasIntercept; // whether intercept term is present
  // 1 is 3 param gompertz; 2 is 4 param logistic
  int meanParameterization;
  // 1 is identity, 2 is exponential, 3 is inverse logit
  vector[d] transformSpecification;
  vector[d] priorThetaSpecification;

  vector[d] priorSigmaThetaSpecification;

  int priorObsVarSpecification;

  // optional
  // likely will be removed
  int<lower = 0, upper = 1> run_estimation; //whether to estimate model or just simulate data

  // Prior Parameters
  // First prior parameter theta. Depends on prior parameterization chosen.
  matrix[U, d] pThetaParam1;
  // Second prior parameter theta. Depends on prior parameterization chosen.
  matrix[U, d] pThetaParam2;
  // vector for means of  betaTheta params
  vector[pTheta*d] pBetaThetaParam1;
  // vector for var of  betaTheta params
  vector[pTheta*d] pBetaThetaParam2;
  // prior on iid variance Z
  vector[2] pObsVar;
  // prior scale for each theta parameter
  vector[d] pSigmaThetaParam1;
  // prior rate for each theta parameter
  vector[d] pSigmaThetaParam2;

}
parameters{
  matrix[U, d] Theta_raw; //Group Mean Params
  vector<lower = 0> [d]SigmaTheta; //diagonal of theta matrix; sd
  matrix[pTheta, d] betaTheta; // beta matrix for covariates in theta mean structure
  real<lower = 0> ObsVar; // sd within group/year is iid
}
transformed parameters{
  matrix[U, d] Theta; //Group Mean Params
  vector[N] mu;
  matrix[U, d] linear_pred_theta;
  matrix[U, d] transTheta;
  //try a matrix beta instead
  linear_pred_theta = XTheta * betaTheta;
  // Theta ~ Normal(X*betaTheta, sqrt(SigmaTheta));
  for(i in 1:U){
    Theta[i] = linear_pred_theta[i] + to_row_vector(sqrt(SigmaTheta)) .* Theta_raw[i];
  }
  transTheta = transformThetaMat(transformSpecification, Theta);

  mu = growthMean(N, u, tpt, transTheta, meanParameterization);

}
model{
  to_vector(betaTheta') ~ normal(pBetaThetaParam1, pBetaThetaParam2);

  target += SigmaThetaPriors(SigmaTheta, pSigmaThetaParam1, pSigmaThetaParam2, priorSigmaThetaSpecification);
  target += ThetaPriors(Theta[(hasIntercept + 1):,],
                        pThetaParam1[(hasIntercept + 1):,],
                        pThetaParam2[(hasIntercept + 1):,],
                        priorThetaSpecification);
  target += ObsVarPriors(ObsVar, pObsVar[1], pObsVar[2], priorObsVarSpecification);


  //yearly Process
  for(i in 1:U){
    Theta_raw[i] ~ normal(0, 1);
  }
  if(run_estimation == 1){
      Z ~ normal(mu,sqrt(ObsVar));
  }
}
generated quantities{
  vector[N] log_lik;

  for(i in 1:N){
    log_lik[i] = normal_lpdf(Z[i] | mu[i], sqrt(ObsVar));
  }
}
