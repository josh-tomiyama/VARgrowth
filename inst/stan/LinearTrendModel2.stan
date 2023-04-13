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
  int tpt[N]; // tpt for each week presumably ordered.
  int U; //Number of groups/years
  int pTheta1; // number of columns in XTheta
  int pTheta2; // number of columns in XTheta
  int pTheta3; // number of columns in XTheta
  int nrowXTheta; // number of rows in design matrix. Should be = U


  //data
  vector[N] Z;
  matrix[U, pTheta1] XTheta1;
  matrix[U, pTheta2] XTheta2;
  matrix[U, pTheta3] XTheta3;
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
  vector[pTheta1] pBetaTheta1Param1;
  // vector for var of  betaTheta params
  vector[pTheta1] pBetaTheta1Param2;
  // vector for means of betaTheta2
  vector[pTheta2] pBetaTheta2Param1;
  // vector for var of  betaTheta params
  vector[pTheta2] pBetaTheta2Param2;
  //vector for means of betaTheta2
  vector[pTheta3l] pBetaTheta3Param1;
  // vector for var of  betaTheta params
  vector[pTheta3] pBetaTheta3Param2;
  // prior on iid variance Z
  vector[2] pObsVar;
  // prior scale for each theta parameter
  vector[d] pSigmaThetaParam1;
  // prior rate for each theta parameter
  vector[d] pSigmaThetaParam2;

}
parameters{
  matrix[U, d] Theta; //Group Mean Params
  vector<lower = 0> [d]SigmaTheta; //diagonal of theta matrix; sd
  vector[pTheta1] betaTheta1; // beta matrix for covariates in theta mean structure
  vector[pTheta2] betaTheta2;
  vector[pTheta3] betaTheta3;
  real<lower = 0> ObsVar; // sd within group/year is iid
}
transformed parameters{
  vector[N] mu;
  vector[U] linear_pred_theta1;
  vector[U] linear_pred_theta2;
  vector[U] linear_pred_theta3;
  matrix[U, d] transTheta = transformThetaMat(transformSpecification, Theta);

  mu = growthMean(N, u, tpt, transTheta, meanParameterization);

  //try a matrix beta instead
  linear_pred_theta1 = XTheta1 * betaTheta1;
  linear_pred_theta2 = XTheta2 * betaTheta2;
  linear_pred_theta3 = XTheta3 * betaTheta3;
}
model{
  // to_vector(betaTheta') ~ normal(pBetaThetaParam1, sqrt(longSigmaTheta) .* pBetaThetaParam2);
  betaTheta1 ~ normal(pBetaTheta1Param1, pBetaTheta1Param2);
  betaTheta2 ~ normal(pBetaTheta2Param1, pBetaTheta2Param2);
  betaTheta3 ~ normal(pBetaTheta3Param1, pBetaTheta3Param2);
  target += SigmaThetaPriors(SigmaTheta, pSigmaThetaParam1, pSigmaThetaParam2,
            priorSigmaThetaSpecification);
  target += ThetaPriors(Theta[(hasIntercept + 1):,],
                        pThetaParam1[(hasIntercept + 1):,],
                        pThetaParam2[(hasIntercept + 1):,],
                        priorThetaSpecification);
  target += ObsVarPriors(ObsVar, pObsVar[1], pObsVar[2], priorObsVarSpecification);


  //yearly Process
  Theta[,1] ~ normal(linear_pred_theta1,
                      sqrt(SigmaTheta[1]));
  Theta[,2] ~ normal(linear_pred_theta2,
                      sqrt(SigmaTheta[2]));
  Theta[,3] ~ normal(linear_pred_theta3,
                      sqrt(SigmaTheta[3]));
  if(run_estimation == 1){
      Z ~ normal(mu, sqrt(ObsVar));
  }
}
generated quantities{
  // These are pointwise log likelihoods
  vector[N] log_lik;
  real theta_samp[U, d, 10000];

  for(i in 1:N){
    log_lik[i] = normal_lpdf(Z[i] | mu[i], sqrt(ObsVar));
//    normal_lpdf(Theta[u[i]] | linear_pred_theta[u[i]], sqrt(SigmaTheta));
  }
}
