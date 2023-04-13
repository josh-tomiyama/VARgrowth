#' Title
#'
#' @param transform_functions
#' @param inverse
#'
#' @return
#' @export
#'
#' @examples
ThetaTransformations <- function(transform_functions, inverse){
  idx <- sapply(transform_functions, function(f){do.call(f, args = list())$idx})
  out <- list(theta_link = ApplyThetaTransformations(transform_functions, inverse),
       inv_theta_link = ApplyThetaTransformations(transform_functions, !inverse),
       stan_input = list(transformSpecification = idx))
  class(out) <- c("ThetaTransformation", "list")
  out
}

is.ThetaTransformation <- function(obj){
  "ThetaTransformation" %in% class(obj)
}

#' Title
#'
#' @param transform_functions
#' @param inverse
#'
#' @return
#' @export
#'
#' @examples
ApplyThetaTransformations <- function(transform_functions, inverse = FALSE){

  if(length(inverse) == 1){
    inverse <- rep(inverse, length(transform_functions))
  }
  function(theta){
    if(length(transform_functions) != ncol(theta)){
      msg <- paste0("\nNumber transform functions ",
                    length(transform_functions),
                    " do not match dimension of theta ", ncol(theta))
      stop(msg)
    }
    transTheta <- theta
    for(i in 1:length(transform_functions)){
      transTheta[,i] = transform_functions[[i]](theta[,i], inverse[i])$val
    }
    transTheta
  }

}

#' Title
#'
#' @return
#' @export
#'
#' @examples
AvailableThetaTransformations <- function(){
  obnm <- ls("package:VARgrowth")
  obnm[stringr::str_detect(obnm, "Transform$")]
}

#' Title
#'
#' @param x
#' @param inverse
#'
#' @return
#' @export
#'
#' @examples
identityTransform <- function(x = 1, inverse = FALSE){
  list(idx = 1,
       val = x)
}

#' Title
#'
#' @param x
#' @param inverse
#'
#' @return
#' @export
#'
#' @examples
logTransform <- function(x = 1, inverse = FALSE){
  if(!inverse & any(x <= 0)){
    stop("\n values outside of range detected")
  }
  list(idx = 2,
       val =
         if(!inverse){
           log(x)
         }else{
           exp(x)
         }
       )
}

#' Title
#'
#' @param x
#' @param inverse
#'
#' @return
#' @export
#'
#' @examples
logitTransform <- function(x = 0.5, inverse = FALSE){
  if(!inverse & any( x > 1 | x < 0) ){
    stop("\n values outside of range detected")
  }
  list(idx = 3,
       val =
         if(!inverse){
           log(x/(1 - x))
         }else{
            exp(x)/(1 + exp(x))
         }
       )
}
