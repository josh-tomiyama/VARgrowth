#' Title
#'
#' @param data
#' @param time_nm
#' @param group_nm
#' @param outcome_nm
#'
#' @return
#' @export
#' @importFrom dplyr select arrange everything
#' @examples
VARGrowthData <- function(data, time_nm, group_nm, outcome_nm = NULL){
  g <- data[[group_nm]]
  if(length(unique(g)) < 2){
    stop("\nGrouping variable doesn't have enough levels\n")
  }

  df <- data
  df$group <- data[[group_nm]]
  df$time <- data[[time_nm]]
  if(!group_nm == 'group'){
    df[[group_nm]] <- NULL
  }
  if(!time_nm == 'time'){
    df[[time_nm]] <- NULL
  }
  if(!is.null(outcome_nm)){
    df$outcome = data[[outcome_nm]]
    df <- select(df, time, group, outcome, everything()) %>%
      arrange(group, time)
  }else{
    df <- select(df, time, group, everything()) %>%
      arrange(group, time)
  }
  class(df) <- c("VARGrowthData", "data.frame")
  df
}

is.VARGrowthData <- function(obj){
  "VARGrowthData" %in% class(obj)
}

#' Title
#'
#' @param PriorList
#' @param required_params
#'
#' @return
#' @export
#'
#' @examples
CheckPriorsValid <- function(PriorList, required_params){
  given_params <- Reduce(union,
         lapply(PriorList,
                function(x){
                  x$params
                }))
  if(!all(required_params %in% given_params)){
    stop("\nPriors specified are not for this trend/mean model")
  }
}

#' Title
#'
#' @param data
#' @param ThetaTrend
#' @param GrowthFunction
#' @param ThetaTransform
#' @param PriorList
#'
#' @return
#' @export
#'
#' @examples
VARGrowthInput <- function(data,
                      ThetaTrend,
                      GrowthFunction,
                      ThetaTransform,
                      PriorList){

  CheckPriorsValid(PriorList, union(ThetaTrend$prior_params, GrowthFunction$prior_params))

  stan_input <- c(list(
    # Total number of observations. Sum T_u over u.
    N = nrow(data),
    # Index of what year/group observation is from.  presumably ordered.
    u = as.integer(data$group),
    # tpt for each week presumably ordered.
    tpt = data$time,
    U = length(unique(data$group)),
    # Num parameters
    d = GrowthFunction$num_params,

    # data
    Z = data$outcome,

    run_estimation = as.integer(1),
    hasIntercept = 1
  ),
  ThetaTransform$stan_input,
  ThetaTrend$stan_input,
  GrowthFunction$stan_input,
  Reduce(c,
         lapply(PriorList,
                function(x){
                  x$stan_input
                })
         )

  )
}


#' Title
#'
#' @param data
#' @param ThetaTransformations
#' @param ThetaTrend
#' @param GrowthFunction
#'
#' @return
#' @export
#'
#' @examples
initValEst <- function(data, ThetaTransformations, ThetaTrend, GrowthFunction){
  # browser()
  theta_est <- GrowthFunction$est_function(data)


  # theta_est[["fit"]] <- NULL
  theta <- theta_est$params$theta
  transTheta <- ThetaTransformations$theta_link(theta)

  meanDefs <- apply(transTheta, 1, GrowthFunction$mean_func)
  names(meanDefs) <- as.character(sort(unique(data$group)))
  groups <- unique(data$group)

  dat <- data %>%
    group_by(group) %>%
    mutate(mu = meanDefs[[unique(as.character(group))]](time)) %>%
    ungroup()
  theta_trend_est <- do.call(ThetaTrend$est_function, c(theta = list(transTheta), ThetaTrend$est_input))
  # theta_trend_est[["fit"]] <- NULL
  return(c(theta_est$params, list(transTheta = transTheta), list(mu = dat$mu), theta_trend_est$params))
}

#' Title
#'
#' @param parameter
#' @param perturb_var
#'
#' @return
#' @export
#'
#' @examples
perturb <- function(parameter, perturb_var = 10){
  if(is.data.frame(parameter)){
    parameter <- as.matrix(parameter)
  }

  if(is.vector(parameter) | is.matrix(parameter)){
    eps <- parameter
    eps[1:length(eps)] <- rnorm(length(eps), 0, abs(parameter)/perturb_var)
  }else{
    stop("Not sure what to do with this input type", class(parameter))
  }
  eps + parameter
}



#' Title
#'
#' @param nchains
#' @param init
#' @param seed
#' @param perturb_var
#'
#' @return
#' @export
#'
#' @examples
initVals_multi <- function(nchains, init, seed = NULL, perturb_var = 10, model){

  if(!is.null(seed)){
    set.seed(seed)
  }else{
    message("\nNo seed was set when initializing chains\n")
  }
  ### Will need generalizing past linear trend model
  ### Harder to perturb covariance matrices in this way...
  ### perturb the eigen values maybe?
  ### could use multiple distpatch
  lapply(1:nchains, function(x, model){
    if(model == "LinearTrendModel"){
      list(SigmaTheta = abs(perturb(init$SigmaTheta, perturb_var)),
           betaTheta = perturb(init$beta, perturb_var),
           Theta = perturb(init$transTheta, perturb_var),
           ObsVar = abs(perturb(init$ObsVar, perturb_var))
      )
    }else if(model == "DeterministicLinearTrendModel"){
      list(betaTheta = perturb(init$beta, perturb_var),
           ObsVar = abs(perturb(sqrt(init$ObsVar), perturb_var))
      )
    }else if(model == "BoundedLinearTrendModel"){
      list(SigmaTheta = abs(perturb(sqrt(init$SigmaTheta), perturb_var)),
           betaTheta = perturb(init$beta, perturb_var),
           Theta = perturb(init$transTheta, perturb_var),
           ObsVar = abs(perturb(sqrt(init$ObsVar), perturb_var)),
           UpperBound = init$UpperBound + abs(1, perturb_var)
      )
    }else if(model == "LinearTrendModel_noncentered"){
      list(SigmaTheta = abs(perturb(init$SigmaTheta, perturb_var)),
           betaTheta = perturb(init$beta, perturb_var),
           mu = perturb(init$mu, perturb_var),
           theta_raw = perturb(matrix(0.1, nrow = nrow(init$transTheta),ncol = ncol(init$transTheta)), perturb_var),
           ObsVar = abs(perturb(init$ObsVar, perturb_var))
      )
    }else{
      stop(paste0(model, " does not have an initial value function"))
    }

  }, model = model)
}

#' Title
#'
#' @param data
#' @param ThetaTrend
#' @param GrowthFunction
#' @param PriorList
#' @param ThetaTransform
#' @param nchains
#' @param seed
#' @param perturb_var
#' @param verbose Whether to print information about the model being fit
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
VARGrowth <- function(data,
                      ThetaTrend,
                      GrowthFunction,
                      PriorList,
                      ThetaTransform = NULL,
                      nchains = 4,
                      seed = sample.int(.Machine$integer.max, 1),
                      perturb_var = 10,
                      verbose = TRUE,
                      ...){

  if(is.null(ThetaTransform)){
    ThetaTransform <- ThetaTransformations(rep(list(identityTransform),
                                               GrowthFunction$num_params),
                                           inverse = FALSE)
  }

  if(!is.ThetaTransformation(ThetaTransform)){
    stop("\nInvalid ThetaTransform object")
  }

  if(!is.ObsMeanObj(GrowthFunction)){
    stop("\nInvalid GrowthFunction object")
  }
  if(!is.ThetaTrendObj(ThetaTrend)){
    stop("\nInvalid ThetaTrend object")
  }

  init_est <- initValEst(data, ThetaTransform, ThetaTrend, GrowthFunction)

  out <- list(
    fit = sampling(
      stanmodels[[ThetaTrend$model]],
      data = VARGrowthInput(data,
                            ThetaTrend,
                            GrowthFunction,
                            ThetaTransform,
                            PriorList),
      init =   initVals_multi(nchains, init_est, seed, perturb_var, ThetaTrend$model),
      seed = seed,
      chains = nchains,
      ...
    ),
    data = data,
    ThetaTransform = ThetaTransform,
    ThetaTrend = ThetaTrend,
    GrowthFunction = GrowthFunction,
    PriorList = PriorList
  )
  class(out) <- "VARGrowthFit"
  out

}




DimTheta <- function(data, GrowthFunction){
  c(length(unique(data$group)), GrowthFunction$num_params)
}

DimBeta <- function(GrowthFunction, ThetaTrend){
  c(ThetaTrend$num_params, GrowthFunction$num_params)
}


DimObsVar <- function(GrowthFunction, ThetaTrend = NULL){
  ### Generalize for multivariate normal maybe
  1
}

DimSigmaTheta <- function(GrowthFunction, ThetaTrend = NULL){
  c(GrowthFunction$num_params, GrowthFunction$num_params)
}

