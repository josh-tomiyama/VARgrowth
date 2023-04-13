### TODO: Take in full matrix Theta and then subset to first p columns

#' Title
#'
#' @param U
#' @param theta0
#' @param k
#' @param A
#' @param p
#' @param SigmaTheta
#' @param ...
#'
#' @return
#' @export
#' @importFrom mvtnorm rmvnorm
#' @importFrom matrixcalc is.positive.definite
#'
#' @examples
VARSim <- function(U, theta0, k, A, p, SigmaTheta, ...){

  if(is(theta0, "vector")){
    theta0 <- matrix(theta0, ncol = length(theta0), byrow = TRUE)
  }

  if(is(A, "list")){
    Afull <- Reduce(cbind, A)
  }

  if(!is.positive.definite(SigmaTheta)){
    warning("The covariance matrix given is not positive definite")
  }
  d <- ncol(theta0)
  # browser()
  ThetaOut <- matrix(ncol = ncol(theta0), nrow = U + p)
  ThetaOut[1:p,] <- theta0
  eps <- rmvnorm(U, mean = rep(0, d), sigma = SigmaTheta)
  for(i in (p+1):nrow(ThetaOut)){
    ThetaOut[i,] <- k + A %*% matrix(ThetaOut[(i - p):(i-1),])
  }

  list(Theta = ThetaOut,
       theta0 = theta0,
       intercept = k,
       coef = A,
       SigmaTheta = SigmaTheta)
}


#' Title
#'
#' @param U
#' @param XTheta
#' @param betaTheta
#' @param SigmaTheta
#' @param simplify
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
linearTrendSim <- function(U, XTheta, betaTheta, SigmaTheta, simplify = TRUE, seed = NULL, ...){
  # browser()
  if(!is.null(seed)){
    set.seed(seed)
  }
  if(is.vector(betaTheta)){
    d = length(betaTheta)
    betaTheta <- matrix(betaTheta, ncol = d, byrow = TRUE)
  }else{
    d = ncol(betaTheta)
  }
  lin_pred <- XTheta %*% betaTheta



  if(!("matrix" %in% class(SigmaTheta))){
    tryCatch({
      SigmaTheta <- diag(SigmaTheta)
    },
    error = function(cond){
      message(cond)
      stop("\nSigmaTheta not a valid vector/matrix")
    })
  }

  eps <- mvtnorm::rmvnorm(U, rep(0, d),  sigma = SigmaTheta)
  ThetaOut <- lin_pred + eps
  if(simplify){
    list(Theta = ThetaOut)
  }else{
    list(Theta = ThetaOut,
         XTheta = XTheta,
         betaTheta = betaTheta,
         Err = eps,
         SigmaTheta = SigmaTheta)
  }

}

#' Title
#'
#' @param U
#' @param XTheta
#' @param betaTheta
#' @param SigmaTheta
#' @param simplify
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
linearTrendSim2 <- function(U, XTheta, betaTheta, SigmaTheta, simplify = TRUE, seed = NULL, ...){
  # browser()
  if(!is.null(seed)){
    set.seed(seed)
  }
  if(is.vector(betaTheta)){
    d = length(betaTheta)
    betaTheta <- matrix(betaTheta, ncol = d, byrow = TRUE)
  }else{
    d = ncol(betaTheta)
  }
  lin_pred <- mapply(function(XTheta, betaTheta){XTheta %*% betaTheta},
                     XTheta, betaTheta)
  # XTheta %*% betaTheta



  if(!("matrix" %in% class(SigmaTheta))){
    tryCatch({
      SigmaTheta <- diag(SigmaTheta)
    },
    error = function(cond){
      message(cond)
      stop("\nSigmaTheta not a valid vector/matrix")
    })
  }

  eps <- mvtnorm::rmvnorm(U, rep(0, d),  sigma = SigmaTheta)
  ThetaOut <- lin_pred + eps
  if(simplify){
    list(Theta = ThetaOut)
  }else{
    list(Theta = ThetaOut,
         XTheta = XTheta,
         betaTheta = betaTheta,
         Err = eps,
         SigmaTheta = SigmaTheta)
  }

}

#' Title
#'
#' @param U
#' @param XTheta
#' @param betaTheta
#' @param simplify
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
deterministicLinearTrendSim <- function(U, XTheta, betaTheta,
                                        simplify = TRUE, seed = NULL, ...){
  # browser()
  if(!is.null(seed)){
    set.seed(seed)
  }
  if(is.vector(betaTheta)){
    d = length(betaTheta)
    betaTheta <- matrix(betaTheta, ncol = d, byrow = TRUE)
  }else{
    d = ncol(betaTheta)
  }
  lin_pred <- XTheta %*% betaTheta

  ThetaOut <- lin_pred
  if(simplify){
    list(Theta = ThetaOut)
  }else{
    list(Theta = ThetaOut,
         XTheta = XTheta,
         betaTheta = betaTheta)
  }

}

#' Title
#'
#' @param dat
#' @param Theta
#' @param ObsVar
#' @param mean_func
#' @param inv_theta_link
#' @param simplify
#' @param ...
#'
#' @return
#' @export
#' @importFrom dplyr group_by mutate ungroup
#' @importFrom stringr str_detect
#' @importFrom truncnorm rtruncnorm
#'
#' @examples
simObs <- function(dat,
                   Theta,
                   ObsVar,
                   GrowthFunction,
                   ThetaTransform,
                   simplify = TRUE,
                   seed = NULL,
                   ...){

  # browser()
  if(!is.null(seed)){
    set.seed(seed)
  }
  # browser()
  if(!str_detect(GrowthFunction$model, "^trunc_")){
    err <- rnorm(nrow(dat), 0, sqrt(ObsVar))
  }else{
    err <- rtruncnorm(nrow(dat), a = -Inf, b = UpperBound, mean = 0, sd = sqrt(ObsVar))
  }
  ### This may be slower than simply vectorizing the result
  ### would need to repeat the Theta matrix an appropriate amount of times to get the same thing
  meanDefs <- apply(ThetaTransform$inv_theta_link(Theta), 1, GrowthFunction$mean_func)
  names(meanDefs) <- as.character(sort(unique(dat$group)))
  groups <- unique(dat$group)

  dat <- dat %>%
    group_by(group) %>%
    mutate(mu = meanDefs[[unique(as.character(group))]](time)) %>%
    ungroup() %>%
    mutate(err = err,
           obs = mu + err)
  if(simplify){
    list(data = dat)
  }else{
    list(data = dat,
         Theta = Theta,
         transform_theta = ThetaTransform$inv_theta_link(Theta),
         mean_functions = meanDefs,
         ObsVar = ObsVar)
  }
}

#' Title
#'
#' @param ThetaTrend
#' @param GrowthFunction
#'
#' @return
#' @export
#'
#' @examples
PrintParamList <- function(ThetaTrend, GrowthFunction){
  print(setdiff(
    union(args(ThetaTrend$sim_function),
          args(GrowthFunction$sim_function)),
    c("U", "...", "dat", "mean_func", "inv_theta_link")
    )
  )
}


#' Title
#'
#' @param dat
#' @param param_list
#' @param ThetaTrend
#' @param GrowthFunction
#' @param ThetaTransform
#' @param seed
#' @param simplify
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
SimulateData <- function(dat,
                         param_list,
                         ThetaTrend,
                         GrowthFunction,
                         ThetaTransform = NULL,
                         seed = sample.int(.Machine$integer.max, 1),
                         simplify = TRUE, ...){
  if(is.null(ThetaTransform)){
    ThetaTransform <- ThetaTransformations(rep(list(identityTransform),
                                               GrowthFunction$num_params),
                                           inverse = FALSE)
  }
  param_list$seed <- seed
  param_list$simplify <- simplify
  simTheta <- do.call(ThetaTrend$sim_function, c(param_list, ThetaTrend$sim_input))
  param_list$Theta <- simTheta$Theta
  simObs <- do.call(simObs, c(list(dat = dat,
                         # Theta = simTheta$Theta,
                         GrowthFunction = GrowthFunction,
                         ThetaTransform = ThetaTransform),
                    param_list))
  # simObs <- simObs(dat, simTheta$Theta, param_list$ObsVar, GrowthFunction$mean_func, ThetaTransform$inv_theta_link, simplify = simplify, ...)

  if(simplify){
    simObs$dat$obs
  }else{
    list(obs = simObs$dat,
         obsSimParams = simObs,
         thetaSimParams = simTheta)
  }

}
