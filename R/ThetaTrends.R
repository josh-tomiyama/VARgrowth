
#' Title
#'
#' @param num_params
#' @param est_input
#' @param est_function
#' @param sim_function
#' @param stan_input
#' @param calcTheta
#' @param prior_params
#' @param sim_params
#' @param model
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
CreateThetaTrendObj <- function(num_params,
                                est_input,
                                est_function,
                                sim_function,
                                stan_input,
                                calcTheta,
                                prior_params,
                                sim_params,
                                model,
                                ...){

  if(!is.function(est_function)){
    stop("\nest_function is not a function")
  }
  if(!is.function(sim_function)){
    stop("\nsim_function is not a function")
  }
  if(!is.list(stan_input)){
    stop("\nstan_input is not a list")
  }
  if(!is.function(calcTheta)){
    stop("\ncalcTheta is not a function")
  }


  obj <- list(num_params = num_params,
       est_input = est_input,
       est_function = est_function,
       sim_function = sim_function,
       stan_input = stan_input,
       prior_params = prior_params,
       sim_params = sim_params,
       model = model,
       ...)
  class(obj) <- c("ThetaTrend", "list")
  obj

}


is.ThetaTrendObj <- function(obj){
  "ThetaTrend" %in% class(obj)
}

#' Title
#'
#' @param theta
#' @param p
#' @param var_model_control
#'
#' @return
#' @export
#'
#' @examples
VARTrend <- function(theta,  p = 1, var_model_control = list()){
  stop('not done yet')
  var_model_control$p <- p
  CreateThetaTrendObj(est_input = var_model_control,
                      est_function = ,
                      sim_function = VARSim,
                      stan_input = list(),
                      calcTheta = function(x){x},
                      prior_params = c("Theta",
                                       "SigmaTheta",
                                       "coefMatrix",
                                       if(var_model_control$type == 'none'){
                                         NULL
                                       }else{
                                         "intercept"
                                       }),
                      sim_params = c("Theta", "coefMatrix", "SigmaTheta",
                                     if(var_model_control$type == 'none'){
                                       NULL
                                     }else{
                                       "intercept"
                                     }),
                      model = 'VARTrend')
}

#' Specify a Linear model for group trends
#'
#' Use this function to define a linear model mean for the group trends of the form \eqn{g(\theta_d) \sim N(X_\theta \beta_d, \sigma_d)}.
#' \eqn{X_\theta} should be a U x p matrix where U is the number of groups and p is the number of linear predictor parameters.
#' Note that as of now every mean parameter \eqn{\theta} is assumed to have the same design matrix.
#' @param data A VARGrowth data object
#' @param formula a formula specifying how to construct the linear model matrix
#' @param XTheta An optional design matrix with number of rows = number of groups in data
#'
#' @return
#' @export
#' @importFrom stats model.matrix
#'
#' @examples
LinearModelTrend <- function(data, formula = ~ 1, XTheta = NULL){
  # browser()
  if(is.null(XTheta)){
    XTheta <- model.matrix(formula, data)
    XTheta <- XTheta[!duplicated(XTheta),]
  }

  if(is.vector(XTheta) == 1){
    XTheta <- matrix(XTheta,
                   nrow = length(unique(data$group)),
                   ncol = length(XTheta),
                   byrow = TRUE)
  }

  if(nrow(XTheta) != length(unique(data$group))){
    stop("\nDesign matrix different dimension from number of groups")
  }

  CreateThetaTrendObj(num_params = ncol(XTheta),
       est_input = list(XTheta = XTheta),
       est_function = LinearTrendEstFunction,
       sim_input = list(U = length(unique(data$group)),
                        XTheta = XTheta),
       sim_function = linearTrendSim,
       stan_input = list(pTheta = ncol(XTheta),
                         nrowXTheta = nrow(XTheta),
                         XTheta = XTheta),
       calcTheta = function(beta, XTheta){XTheta %*% beta},
       prior_params = c("Theta", "SigmaTheta", "betaTheta"),
       sim_params = c("SigmaTheta", "betaTheta"),
       formula = formula,
       model = 'LinearTrendModel')
}

#' Specify a Linear model for each parameter
#'
#' Use this function to define a linear model mean for the group trends of the form \eqn{g(\theta_d) \sim N(X_\theta_d \beta_d, \sigma_d)}.
#' \eqn{X_\theta} should be a U x p_d matrix where U is the number of groups and p_d is the number of linear predictor parameters for the d-th parameter.
#' Note that as of now every mean parameter \eqn{\theta} is assumed to have the same design matrix.
#' @param data A VARGrowth data object
#' @param formula a list of formulas specifying how to construct the linear model matrix for each parameter
#' @param XTheta An optional list of design matrix with number of rows = number of groups in data
#'
#' @return
#' @export
#'
#' @examples
LinearModelTrend2 <- function(data, formula = ~ 1, XTheta = NULL){

  if(is.null(XTheta)){
    XTheta <- lapply(formula, function(formula, data){
      XTheta_temp <- model.matrix(formula, data)
      XTheta_temp <- XTheta[!duplicated(XTheta_temp),]
    }, data)

  }

  if(is.vector(XTheta) == 1){
    XTheta <- matrix(XTheta,
                     nrow = length(unique(data$group)),
                     ncol = length(XTheta),
                     byrow = TRUE)
  }

  if(nrow(XTheta) != length(unique(data$group))){
    stop("\nDesign matrix different dimension from number of groups")
  }

  CreateThetaTrendObj(num_params = ncol(XTheta),
                      est_input = list(XTheta = XTheta),
                      est_function = LinearTrendEstFunction,
                      sim_input = list(U = length(unique(data$group)),
                                       XTheta = XTheta),
                      sim_function = linearTrendSim,
                      stan_input = list(pTheta = ncol(XTheta),
                                        nrowXTheta = nrow(XTheta),
                                        XTheta = XTheta),
                      calcTheta = function(beta, XTheta){XTheta %*% beta},
                      prior_params = c("Theta", "SigmaTheta", "betaTheta"),
                      sim_params = c("SigmaTheta", "betaTheta"),
                      formula = formula,
                      model = 'LinearTrendModel')
}

#' Specify a Linear model with no variance
#'
#' Use this function to define a linear model mean for the group trends of the form \eqn{g(\theta_d) = X_\theta \beta_d}.
#' \eqn{X_\theta} should be a U x p matrix where U is the number of groups and p is the number of linear predictor parameters.
#' Note that as of now every mean parameter \eqn{\theta} is assumed to have the same design matrix.
#' @param data A VARGrowth data object
#' @param formula a formula specifying how to construct the linear model matrix
#' @param XTheta An optional design matrix with number of rows = number of groups in data
#'
#' @return
#' @export
#'
#' @examples
DeterministicLinearModelTrend <- function(data, formula = ~ 1, XTheta = NULL){

  if(is.null(XTheta)){
    XTheta <- model.matrix(formula, data)
    XTheta <- XTheta[!duplicated(XTheta),]
  }

  if(is.vector(XTheta) == 1){
    XTheta <- matrix(XTheta,
                     nrow = length(unique(data$group)),
                     ncol = length(XTheta),
                     byrow = TRUE)
  }

  if(nrow(XTheta) != length(unique(data$group))){
    stop("\nDesign matrix different dimension from number of groups")
  }

  CreateThetaTrendObj(num_params = ncol(XTheta),
                      est_input = list(XTheta = XTheta),
                      est_function = DeterministicLinearTrendEstFunction,
                      sim_input = list(U = length(unique(data$group)),
                                       XTheta = XTheta),
                      sim_function = deterministicLinearTrendSim,
                      stan_input = list(pTheta = ncol(XTheta),
                                        nrowXTheta = nrow(XTheta),
                                        XTheta = XTheta),
                      calcTheta = function(beta, XTheta){XTheta %*% beta},
                      prior_params = c("betaTheta"),
                      sim_params = c("betaTheta"),
                      formula = formula,
                      model = 'DeterministicLinearTrendModel')
}


#' Estimate Linear trend of thetas
#'
#' @param theta
#' @param XTheta
#' @param ...
#'
#' @return
#' @importFrom  stats lm coef
#' @export
#' @depends stats
#'
#' @examples
LinearTrendEstFunction <- function(theta, XTheta, ...){
  thetaFits <- apply(theta, 2, function(col){ lm(col ~ -1 + XTheta)})
  ### why is it not detecting lm object?
  sThetaFits <- lapply(thetaFits, summary.lm)
  beta <- sapply(thetaFits, coef)
  if(!is.matrix(beta)){
    ## NEED TO GENERALIZE THIS FOR MORE MEAN PARAMS
    beta <- matrix(beta, ncol = 3, byrow = TRUE)
  }
  SigmaTheta <- sapply(sThetaFits, function(sfit){sfit$sigma^2})

  list(fit = thetaFits,
       params = list(beta = beta,
                     SigmaTheta = SigmaTheta)
  )
}

#' Estimate Linear trend of thetas
#'
#' @param theta
#' @param XTheta
#' @param ...
#'
#' @return
#' @importFrom  stats lm coef
#' @export
#' @depends stats
#'
#' @examples
LinearTrendEstFunction2 <- function(theta, XTheta, ...){

  thetaFits <- mapply(function(col, XTheta){
    lm(col ~ -1 + XTheta)
  }, theta, XTheta)
  # thetaFits <- apply(theta, 2, function(col){ lm(col ~ -1 + XTheta)})
  ### why is it not detecting lm object?
  sThetaFits <- lapply(thetaFits, summary.lm)
  beta <- lapply(thetaFits, coef)
  # if(!is.matrix(beta)){
  #   ## NEED TO GENERALIZE THIS FOR MORE MEAN PARAMS
  #   beta <- matrix(beta, ncol = 3, byrow = TRUE)
  # }
  SigmaTheta <- sapply(sThetaFits, function(sfit){sfit$sigma^2})

  list(fit = thetaFits,
       params = list(beta = beta,
                     SigmaTheta = SigmaTheta)
  )
}


#' Estimate Deterministic Linear trend of thetas
#'
#' @param theta
#' @param XTheta
#' @param ...
#'
#' @return
#' @importFrom  stats lm coef
#' @export
#' @depends stats
#'
#' @examples
DeterministicLinearTrendEstFunction <- function(theta, XTheta, ...){
  thetaFits <- apply(theta, 2, function(col){ qr.solve(XTheta, col)})
  list(fit = thetaFits,
       params = list(beta = thetaFits)
  )
}



#' Title
#'
#' @param theta
#' @param XTheta
#' @param var_model_control
#'
#' @return
#' @export
#' @depends vars
#' @examples
VAREstFunction <- function(theta, XTheta, var_model_control){
  browser()
  ### probably want to do intercept model separately?
  fitVAR <- do.call(VAR, c(list(theta), var_model_control))
  sfitVAR <- summary(fitVAR)
  coefMatrix = Acoef(fitVAR0)[[1]]
  list(fit = fitVAR,
       params = list(intercept = if(var_model_control$type == 'none'){
         rep(0, nrow(coefMatrix))
       }else{
         Bcoef(fitVAR)[,ncol(Bcoef(fitVAR))]
       },
       coefMatrix = coefMatrix,
       SigmaTheta = sfitVAR$covres)
  )
}
