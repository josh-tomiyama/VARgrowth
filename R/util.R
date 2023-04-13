#' Convert a VARGrowthFit into a stanfit object
#'
#' @param obj A VARGrowthFit Object
#'
#' @return the stan fit object
#' @export
#'
#' @examples
as.stanfit <- function(obj){
  if(inherits(obj, "VARGrowthFit")){
    obj$fit
  }else{
    stop(paste0("\nDon't know how to convert object of class",
                paste0(class(obj), collapse = ", ")
    )
    )
  }
}


#' Crude check of HMC sampler
#'
#' This function will check the number of divergent transitions is 0, whether bfmi is greater than or equal to 0.3,
#' and whether any transitions met the maximum specified tree depth. If any one of these checks are FALSE, then FALSE is returned.
#' Otherwise, true is returned. Please refer to stan documentation for more information on these diagnostics and how to fix them.
#'
#' @param object an object that can be coerced to a stanfit object
#'
#' @return TRUE if it passes tests, FALSE if not
#' @importFrom rstan get_num_divergent get_bfmi get_num_max_treedepth
#' @export
#'
#' @examples
checkSamplerDiagnostics <- function(object){
  object <- as.stanfit(object)
  get_num_divergent(object) == 0 &
    all(get_bfmi(object) >= 0.3) &
    get_num_max_treedepth(object) == 0
}
