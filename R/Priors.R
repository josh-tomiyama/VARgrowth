
#' Title
#'
#' @return
#' @export
#'
#' @examples
NoPrior <- function(){
  list(params = "Theta",
       prior_idx = 0)
}

#' Title
#'
#' @return
#' @export
#'
#' @examples
UniformPrior <- function(){
  list(params = c("Theta", "UpperBound", "SigmaTheta"),
       prior_idx = 1)
}

#' Title
#'
#' @return
#' @export
#'
#' @examples
NormalPrior <- function(){
  list(params = c("betaTheta", "SigmaTheta", "Theta", "ObsVar"),
       prior_idx = 2)
}

#' Title
#'
#' @return
#' @export
#'
#' @examples
InvGammaPrior <- function(){
  list(params = c("ObsVar", "SigmaTheta"),
       prior_idx = 3)
}


checkValidDistribution <- function(dist, param){
  is_valid <- sapply(dist,
                     function(x){
                       param %in% x$params
                     }
  )

  if(!all(is_valid)){
    stop("\nPrior is not implemented for this parameter")
  }
}


#' Title
#'
#' @param dist
#' @param priorParam1
#' @param priorParam2
#' @param skipFirstTheta
#'
#' @return
#' @export
#'
#' @examples
PriorTheta <- function(dist,
                       priorParam1,
                       priorParam2,
                       skipFirstTheta = 0){

  checkValidDistribution(dist, param = 'Theta')

  if(ncol(priorParam1) != length(dist)){
    stop("\nColumns of Prior Param 1 don't match number of distributions")
  }

  if(ncol(priorParam2) != length(dist)){
    stop("\nColumns of Prior Param 2 don't match number of distributions")
  }

  list(stan_input = list(priorThetaSpecification = sapply(dist, function(x){x$prior_idx}),
                         hasIntercept = skipFirstTheta,
                         pThetaParam1 = priorParam1,
                         pThetaParam2 = priorParam2),
       models = c("linearTrend", "VARTrend"),
       params = "Theta"
       )
}

#' Title
#'
#' @param dist
#' @param priorParam1
#' @param priorParam2
#'
#' @return
#' @export
#'
#' @examples
PriorObsVar <- function(dist, priorParam1, priorParam2){
  checkValidDistribution(dist, param = 'ObsVar')
  list(stan_input = list(priorObsVarSpecification =  sapply(dist, function(x){x$prior_idx}),
                         pObsVar = c(priorParam1, priorParam2)
                         ),
       params = "ObsVar",
       models = c("linearTrend", "VARTrend")
       )
}

#' Title
#'
#' @param dist
#' @param priorParam1
#' @param priorParam2
#'
#' @return
#' @export
#'
#' @examples
PriorBetaTheta <- function(dist, priorParam1, priorParam2){
  checkValidDistribution(dist, param = 'betaTheta')
  list(
    stan_input = list(pBetaThetaParam1 = as.vector(priorParam1),
                      pBetaThetaParam2 = as.vector(priorParam2)),
    models = c("linearTrend"),
    params = "betaTheta"
  )
}

#' Title
#'
#' @param dist
#' @param priorParam1
#' @param priorParam2
#'
#' @return
#' @export
#'
#' @examples
PriorSigmaTheta <- function(dist, priorParam1, priorParam2){
  checkValidDistribution(dist, param = 'SigmaTheta')
  list(
    stan_input = list(priorSigmaThetaSpecification =  sapply(dist, function(x){x$prior_idx}),
                      pSigmaThetaParam1 = priorParam1,
                      pSigmaThetaParam2 = priorParam2),
    models = c("linearTrend", "VARTrend"),
    params = "SigmaTheta"
  )
}


#' Title
#'
#' @param dist
#' @param priorParam1
#' @param priorParam2
#'
#' @return
#' @export
#'
#' @examples
PriorUpperBound <- function(dist, priorParam1, priorParam2){
  checkValidDistribution(dist, param = 'UpperBound')
  list(
    stan_input = list(pUpperBound1 = priorParam1,
                      pUpperBound2 = priorParam2),
    models = c("linearTrend", "VARTrend"),
    params = "SigmaTheta"
  )
}
