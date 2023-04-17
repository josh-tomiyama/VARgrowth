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
  matrix[U, d] Theta; //Group Mean Params
  vector<lower = 0> [d]SigmaTheta; //diagonal of theta matrix; sd
  matrix[pTheta, d] betaTheta; // beta matrix for covariates in theta mean structure
  real<lower = 0> ObsVar; // sd within group/year is iid
}
transformed parameters{
  vector[N] mu;
  matrix[U, d] linear_pred_theta;
  vector[pTheta*d] longSigmaTheta;
  matrix[U, d] transTheta = transformThetaMat(transformSpecification, Theta);

  for(i in 1:pTheta){
    for(j in  1:d){
      longSigmaTheta[j + (i-1)*d] = SigmaTheta[j];
    }
  }

  mu = growthMean(N, u, tpt, transTheta, meanParameterization);

  //try a matrix beta instead
  linear_pred_theta = XTheta * betaTheta;
}
model{
  // to_vector(betaTheta') ~ normal(pBetaThetaParam1, sqrt(longSigmaTheta) .* pBetaThetaParam2);
  to_vector(betaTheta') ~ normal(pBetaThetaParam1, pBetaThetaParam2);

  target += SigmaThetaPriors(SigmaTheta, pSigmaThetaParam1, pSigmaThetaParam2, priorSigmaThetaSpecification);
  target += ThetaPriors(Theta[(hasIntercept + 1):,],
                        pThetaParam1[(hasIntercept + 1):,],
                        pThetaParam2[(hasIntercept + 1):,],
                        priorThetaSpecification);
  target += ObsVarPriors(ObsVar, pObsVar[1], pObsVar[2], priorObsVarSpecification);


  //yearly Process
  for(i in 1:U){
    Theta[i] ~ normal(to_vector(linear_pred_theta[i]),
                      sqrt(SigmaTheta));
  }
  if(run_estimation == 1){
      Z ~ normal(mu,sqrt(ObsVar));
  }
}
generated quantities{
  // These are pointwise log likelihoods
  vector[N] log_lik;
  // real theta_samp[U, d, 10000];
  // matrix[N, 10000] marg_log_lik_samp;
  // vector[N] mu_temp;
  // vector[N] marg_log_lik;

  for(i in 1:N){
    log_lik[i] = normal_lpdf(Z[i] | mu[i], sqrt(ObsVar));
//    normal_lpdf(Theta[u[i]] | linear_pred_theta[u[i]], sqrt(SigmaTheta));
  }

  // for(i in 1:U){
  //   for(j in 1:(10000)){
  //     theta_samp[i, 1:d, j] = normal_rng(to_vector(linear_pred_theta[i]),
  //                             sqrt(SigmaTheta));
  //   }
  // }
  //
  // for(j in 1:(10000)){
  //   mu_temp =
  //   growthMean(N, u, tpt,
  //   transformThetaMat(transformSpecification, to_matrix(theta_samp[,,j])),
  //   meanParameterization);
  //
  //   for(i in 1:N){
  //     marg_log_lik_samp[i, j] = exp(normal_lpdf(Z[i] | mu_temp[i] ,sqrt(ObsVar)));
  //   }
  // }
  //
  // for(i in 1:N){
  //   marg_log_lik[i] = log(mean(marg_log_lik_samp[i, 1:(10000)]));
  // }
}
