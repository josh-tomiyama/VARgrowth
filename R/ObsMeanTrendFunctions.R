#' Constructor for an Obs Mean Trend
#'
#' This is not ready for users yet.
#'
#' @param num_params the number of parameters in your mean trend. Do not count any variance parameters.
#' @param est_function A function for estimating the parameters of your model. Must return a named list with "fit" for the fitted object and "params" that extract all parameters estimates (including variance) from the model.
#' @param stan_input
#' @param mean_func A function factory that will take in all the parameters of theta and return the parametric function with argument x for time.
#' @param prior_params Optional vector of the name of parameters that require priors
#' @param sim_params Option vector of names of parameters needed to run a simulation
#' @param model A model name for the parameterization
#' @param param_names Names of the parameters theta (not the variance parameters)
#' @param generator a function to generate simulations from the likelihood
#' @param ... Other meta data to be attached to the object
#'
#' @return An ObsMeanDef object which is a list of all the given arguments.
#' @export
#'
#' @examples
CreateObsMeanObj <- function(num_params, est_function, stan_input,
                             mean_func, prior_params,
                             sim_params, model, param_names = "", ...){
  ### note: can possible extend with an err calculation function to allow for
  ### multivariate normals on the observed process...
  obj <- list(num_params = num_params,
              param_names = param_names,
              est_function = est_function,
              stan_input = stan_input,
              mean_func = mean_func,
              prior_params =prior_params,
              sim_params = c("ObsVar"),
              model = sim_params)
  class(obj) <- c("ObsMeanDef", "list")
  obj
}

is.ObsMeanObj <- function(obj){
  "ObsMeanDef" %in% class(obj)
}

#' Title
#'
#' @return
#' @export
#'
#' @examples
GompertzFunction <- function(dist = 'normal'){
  CreateObsMeanObj(num_params = 3,
       param_names = c("Asymptote",
                       "Offset",
                       "Growth"),
       est_function = GompertzEstFunction_nlme,
       stan_input = list(meanParameterization = 1),
       mean_func = GompertzFunctionBuilder,
       prior_params = c("ObsVar", "Theta"),
       sim_params = c("ObsVar"),
       model = 'Gompertz',
       generator = function(n, param1, param2, ...){rnorm(n, param1, param2)})
}


#' Title
#'
#' @return
#' @export
#'
#' @examples
BoundedGompertzFunction <- function(dist = NormalDist()){
  CreateObsMeanObj(num_params = 3,
                   param_names = c("Asymptote",
                                   "Offset",
                                   "Growth",
                                   "UpperBound"),
                   est_function = BoundedGompertzEstFunction_nlme,
                   stan_input = list(meanParameterization = 1),
                   mean_func = GompertzFunctionBuilder,
                   prior_params = c("ObsVar", "Theta", "UpperBound"),
                   sim_params = c("ObsVar"),
                   model = 'trunc_Gompertz',
                   generator = function(n, param1, param2, ...){rtruncnorm(n, param1, param2)})
}

#' Title
#'
#' @return
#' @export
#'
#' @examples
Logistic4Function <- function(dist = NormalDist()){
  CreateObsMeanObj(num_params = 4,
       param_names = c("Min",
                       "Max",
                       "Inflect",
                       "Scale"),
       est_function = Logistic4EstFunction_nlme,
       stan_input = list(meanParameterization = 2),
       mean_func = Logistic4FunctionBuilder,
       prior_params = c("ObsVar", "Theta"),
       sim_params = c("ObsVar"),
       model = 'Logistic')
}


#' Estimates Initial values for a Gompertz trend using a non-linear mixed effects model
#'
#' @param data
#'
#' @return
#' @export
#' @depends nlme
#' @importFrom nlme nlme groupedData
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#'
#' @examples
GompertzEstFunction_nlme <- function(data){
  gdata <- groupedData(outcome ~ 1 | group, order.groups = FALSE, data = data)
  fitME <- nlme(outcome ~ SSgompertz(time, asymptote, offset, growth),
                random = asymptote + offset + growth ~ 1,
                data = gdata,
                control = nlmeControl(maxIter = 1000, msMaxIter = 1000),
                method = 'ML')
  ### ensure theta is in correct order
  g <- data.frame(group = unique(sort(data$group)))
  cf <- coef(fitME)
  cf <- cbind(rownames(cf), cf)
  colnames(cf)[1] <- "group"
  theta_out <- left_join(g, cf)
  theta_out <- select(theta_out, -group)

  list(fit = fitME,
       params = list(theta = as.matrix(theta_out),
                     ObsVar = fitME$sigma^2)
  )
}

#' Estimates Initial values for a Gompertz trend using a non-linear mixed effects model.
#' Initial values are estimated with the usual non-truncated normal likelihood and bounds are simply 10% above or below the minimum.
#'
#' @param data
#'
#' @return
#' @export
#' @depends nlme
#' @importFrom nlme nlme groupedData
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#'
#' @examples
BoundedGompertzEstFunction_nlme <- function(data){
  gdata <- groupedData(outcome ~ 1 | group, order.groups = FALSE, data = data)
  fitME <- nlme(outcome ~ SSgompertz(time, asymptote, offset, growth),
                random = asymptote + offset + growth ~ 1,
                data = gdata,
                control = nlmeControl(maxIter = 1000, msMaxIter = 1000),
                method = 'ML')
  ### ensure theta is in correct order
  g <- data.frame(group = unique(sort(data$group)))
  cf <- coef(fitME)
  cf <- cbind(rownames(cf), cf)
  colnames(cf)[1] <- "group"
  theta_out <- left_join(g, cf)
  theta_out <- select(theta_out, -group)

  list(fit = fitME,
       params = list(theta = as.matrix(theta_out),
                     ObsVar = fitME$sigma^2,
                     UpperBound = 1.1*max(data$outcome))
  )
}

#' Title
#'
#' @param data
#'
#' @return
#' @importFrom nlme nlmeControl
#' @export
#'
#' @examples
Logistic4EstFunction_nlme <- function(data){
  gdata <- groupedData(outcome ~ 1 | group, order.groups = FALSE, data = data)
  fitME <- nlme(outcome ~ SSfpl(time, min, max, inflect, scale),
                random = min + max  + inflect + scale ~ 1,
                data = data,
                control = nlmeControl(maxIter = 1000, msMaxIter = 1000),
                method = 'ML')
  list(fit = fitME,
       params = list(theta = coef(fitME),
                     ObsVar = fitME$sigma^2)
  )
}

#' Title
#'
#' @param theta_vec
#'
#' @return
#' @export
#'
#' @examples
GompertzFunctionBuilder <- function(theta_vec){

  if(length(theta_vec) > 3){
    warning("\nThis gompertz parameterization has 3 parameters. More were detected.")
  }
  function(x){
    theta_vec[1]*exp(-theta_vec[2]*theta_vec[3]^x)
  }
}

#' Title
#'
#' @param theta_vec
#'
#' @return
#' @export
#'
#' @examples
Logistic4FunctionBuilder <- function(theta_vec){

  if(length(theta_vec) > 4){
    warning("\nThis 4 param logistic parameterization has 4 parameters. More were detected.")
  }
  function(x){
    theta_vec[2] +
      (theta_vec[1] - theta_vec[2]) /
      (exp( (theta_vec[3] - time) / theta_vec[4]) + 1)
  }
}

gompertz_half_point <- function(b2 = 2, b3 = 0.8){
  (log(log(2)) - log(b2))/log(b3)
}

gompertz_temp <- function(x, Asym = 100, b2 = 2, b3 = 0.8){
  Asym*exp(-b2*b3^x)
}

deriv_gompertz <- function(t, Asym = 100, b2 = 2, b3 = 0.8, shift = 0){
  Asym*exp(-b2*b3^t)*(-b2*log(b3)*b3^t) + shift
}
