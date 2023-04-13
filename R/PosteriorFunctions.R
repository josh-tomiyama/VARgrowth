#' Get Posterior Predictive distribution of VARGrowthFit
#'
#' @param object A VARGrowthFit object
#' @param newdata A VARGrowthData object
#' @param ndraws Number of draws to resample from the posterior distribution
#' @param ncores Number of cores to run in parallel. Not used if newdata = NULL
#' @param seed set seed for reproducibility
#' @param ... additional arguments to be passed to posterior::resample_draws
#'
#' @return A matrix of size nrow(newdata) x ndraws
#'
#' @import parallel posterior
#' @export
#' @examples
posterior_predict.VARGrowthFit <- function(object,
                              newdata = NULL,
                              ndraws = NULL,
                              ncores = NULL,
                              seed = sample.int(.Machine$integer.max, 1),
                              ...){
  ### Combine fit into one iteration per row

  set.seed(seed)
  post <- as_draws_matrix(object$fit)
  sim.seeds = sample.int(.Machine$integer.max, nrow(post))

  if(is.null(ndraws)){
    ndraws = nrow(post)
  }
  if(ndraws > nrow(post)){
    warning("\nNumber of draws is larger than posterior samples, using max number")
    ndraws <- nrow(post)
  }
  post <- resample_draws(post, ndraws = ndraws, weights = rep(1, nrow(post)), ...)

  if(is.null(newdata)){
    ### Will need extension if multivariate normal dist used on observed process
    post_mu <- subset_draws(post, "mu")
    post_ObsVar <- subset_draws(post, "ObsVar")
    eps <- apply(post_ObsVar, 1, function(v){rnorm(ncol(post_mu), 0, sqrt(v))})
    post_pred <- post_mu + t(eps)
    post_pred <- as.data.frame(t(post_pred))
    colnames(post_pred) <- paste0("iter_", 1:ncol(post_pred))
    out <- cbind(object$data, as.data.frame(post_pred))
  }else{
    param_list <- CreateParamList(object, post)
    for(i in 1:length(param_list)){
      param_list[[i]]$seed = sim.seeds[i]
    }
    # browser()
    observedGroups <- unique(object$data$group)
    newdataObservedGroups <- dplyr::filter(newdata, group %in% observedGroups)
    newdataUnobservedGroups <- dplyr::filter(newdata, !(group %in% observedGroups))

    if(object$ThetaTrend$model == "LinearTrendModel" |
       object$ThetaTrend$model == "DeterministicLinearTrendModel"
       ){
      if(nrow(newdataObservedGroups) > 0){
        ppObs <- posterior_predict_ObsGroups(object, param_list,
                                             newdataObservedGroups,
                                             ncores)
      } else{
        ppObs <- NULL
      }

      if(nrow(newdataUnobservedGroups) > 0){
        ppUnobs <- posterior_predict_UnobsGroups(object, param_list,
                                                 newdataUnobservedGroups, ncores)
      } else{
        ppUnobs <- NULL
      }

      out <- rbind(ppObs, ppUnobs)
      out <- out %>% arrange(group, time)
    }else{
      stop("\nNot Ready for other types yet")
    }

  }
  class(out) <- c("VARGrowthPosteriorPredictive", class(out))
  out
}

#' Title
#'
#' @param object
#' @param param_list
#' @param newdata
#' @param ncores
#' @param simplify
#'
#' @return
#' @export
#'
#' @examples
posterior_predict_ObsGroups <- function(object,
                                        param_list,
                                        newdata,
                                        ncores = NULL,
                                        simplify = TRUE){
  # browser()
  if(!is.null(ncores)){
    cl <- makeCluster(ncores)
    clusterEvalQ(cl, {library(VARgrowth);library(mvtnorm); library(dplyr)})
    on.exit(stopCluster(cl))
  }
  gf <- object$GrowthFunction
  tt <- object$ThetaTransform
  if(!is.null(ncores)){
    post_pred <- parSapply(cl, param_list,
                           function(pl, gf, tt, newdata){
                             sim <- do.call(simObs,
                                            c(list(dat = newdata,
                                                 GrowthFunction = gf,
                                                 ThetaTransform = tt,
                                                 simplify = simplify),
                                              pl))
                             sim$dat$obs
                           },
                           gf = gf,
                           tt = tt,
                           newdata = newdata)
    colnames(post_pred) <- paste0("iter_", 1:length(param_list))
  }else{
    for(i in 1:length(param_list)){
      post_pred <- matrix(nrow = nrow(newdata), ncol = length(param_list))
      colnames(post_pred) <- paste0("iter_", 1:length(param_list))
      sim <- do.call(simObs,
                     c(list(dat = newdata,
                            GrowthFunction = gf,
                            ThetaTransform = tt,
                            simplify = TRUE),
                       param_list[[i]]))
      post_pred[,i] <- sim
    }
  }
  cbind(newdata, as.data.frame(post_pred))
}


#' Title
#'
#' @param object
#' @param param_list
#' @param newdata
#' @param ncores
#' @param simplify
#'
#' @return
#' @export
#'
#' @examples
posterior_predict_UnobsGroups <- function(object,
                                        param_list,
                                        newdata,
                                        ncores = NULL,
                                        simplify = TRUE){
  # browser()
  if(!is.null(ncores)){
    cl <- parallel::makeCluster(ncores)
    clusterEvalQ(cl, {library(VARgrowth);library(mvtnorm); library(dplyr)})
    on.exit(stopCluster(cl))
  }

  gf <- object$GrowthFunction
  ttrans <- object$ThetaTransform
  if(object$ThetaTrend$model == "LinearModelTrend"){
    ttrend <- LinearModelTrend(newdata, object$ThetaTrend$formula)
  }else if(object$ThetaTrend$model == "DeterministicLinearModelTrend"){
    ttrend <- DeterministicLinearModelTrend(newdata, object$ThetaTrend$formula)
  }else{
    stop("Posterior Predict not updated for this model type")
  }


  if(!is.null(ncores)){

    post_pred <- parSapply(cl, param_list,
                           function(pl, gf, ttrans, ttrend, newdata){
                             pl[["Theta"]] <- NULL
                             sim <- do.call(SimulateData,
                                            list(dat = newdata,
                                                   GrowthFunction = gf,
                                                   ThetaTransform = ttrans,
                                                   ThetaTrend = ttrend,
                                                   simplify = simplify,
                                                   param_list = pl,
                                                   seed = pl$seed))
                             sim
                           },
                           gf = gf,
                           ttrans = ttrans,
                           ttrend = ttrend,
                           newdata = newdata)
    colnames(post_pred) <- paste0("iter_", 1:length(param_list))
  }else{
    for(i in 1:length(param_list)){
      post_pred <- matrix(nrow = nrow(newdata), ncol = length(param_list))
      colnames(post_pred) <- paste0("iter_", 1:length(param_list))
      pl <- param_list[[i]]
      sim <- do.call(SimulateData,
                     list(dat = newdata,
                          GrowthFunction = gf,
                          ThetaTransform = ttrans,
                          ThetaTrend = ttrend,
                          simplify = simplify,
                          param_list = pl,
                          seed = pl$seed)
                     )
      post_pred[,i] <- sim
    }
  }
  cbind(newdata, as.data.frame(post_pred))
}

#' Plot Posterior predictive distribution
#'
#' @param x a VARGrowthPosteriorPredictive object
#' @param y argument is not used
#' @param alpha level to produce (1 - alpha) credible intervals
#' @param estimate a function to produce the point estimate of the posterior predictive distribution
#' @param ... additional arguments that are not used (yet)
#'
#' @return a ggplot object that plots the point estimate and credible interval for each group
#' @importFrom ggplot2 geom_point geom_line facet_wrap
#' @importFrom dplyr rowwise summarize
#' @export
#' @depends nlme stats var
#' @examples
plot.VARGrowthPosteriorPredictive <- function(x, y = NULL, alpha = 0.05, estimate = mean,...){

  spp <- x %>%
    rowwise(group, time) %>%
    summarize(pp_mean = estimate(c_across(starts_with("iter_"))),
              pp_q1 = quantile(c_across(starts_with("iter_")), probs = alpha /2),
              pp_q2 = quantile(c_across(starts_with("iter_")), probs = 1 - alpha/2))

  ggplot(spp, aes(time, pp_mean)) +
    geom_point() +
    geom_line() +
    geom_line(aes(time, pp_q1, color = 'red'), linetype = 'dashed') +
    geom_line(aes(time, pp_q2, color = 'red'), linetype = 'dashed') +
    facet_wrap(~group)
}

#' Title
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' @exportS3Method
print.VARGrowthPosterior <- function(x, ...){

}

#' Title
#'
#' @param object
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
summary.VARGrowthPosterior <- function(object, estimate = mean, alpha = .05, ...){
  spp <- x %>%
    rowwise(group, time) %>%
    summarize(pp_mean = estimate(c_across(starts_with("iter_"))),
              pp_q1 = quantile(c_across(starts_with("iter_")), probs = alpha /2),
              pp_q2 = quantile(c_across(starts_with("iter_")), probs = 1 - alpha/2))
  class(spp) <- c("VARGrowthPosteriorSummary", "data.frame")
  spp
}

#' Title
#'
#' @param x
#' @param y
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot.VARGrowthFit <- function(x, y = NULL, ...){

}



#' Title
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
print.VARGrowthFit <- function(x, ...){

}



#' Title
#'
#' @param object
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
summary.VARGrowthFit <- function(object, ...){
  summary(object$fit)
}


#' Format list of parameters of each posterior sample
#'
#' @param object A VARGrowthFit object
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
CreateParamList <- function(object, posterior, ...){
  params <- union(object$ThetaTrend$prior_params, object$GrowthFunction$prior_params)
  names(params) <- params
  param_list <- lapply(params, FormatParameterOutput, object, posterior)
  out <- list()
  for(i in 1:nrow(posterior)){
    out[[i]] <- lapply(param_list, function(x, idx){x[[idx]]}, idx = i)
  }
  out
}

FormatParameterOutput <- function(param_name, object, posterior){
  posterior_param <- subset_draws(posterior, param_name)
  f <- if(param_name == "betaTheta"){
    function(post_iter, p = object$ThetaTrend$num_params){
      matrix(post_iter, nrow = p)
    }
  }else if(param_name == "Theta"){
    function(post_iter){
      matrix(post_iter, nrow = length(unique(object$data$group)))
    }
  }else if (param_name == "ObsVar"){
    function(post_iter){
      as.vector(post_iter)
    }
  }else if (param_name == "SigmaTheta"){
    function(post_iter){
      as.vector(post_iter)
    }
  }else if (param_name == "UpperBound"){
    function(post_iter){
      as.vector(post_iter)
    }
  }else{
    stop(paste0("\nParam Name", param_name, "not expected"))
  }
  apply(posterior_param, 1, f, simplify=FALSE)
}

dimParams <- function(object, ...){
  UseMethod("dimParams")
}

dimParams.VARGrowthFit <- function(object, ...){
  dimParams(object$GrowthFunction, object$ThetaTrend)
}

dimParams.ThetaTrend <- function(object, GrowthFunction, ...){
  params <- union(object$ThetaTrend$prior_params, object$GrowthFunction$prior_params)
  lapply(params, )
}

dimParams.ObsMeanDef <- function(object, ThetaTrend, ...){

}


# get_post_pred <- function(fit){
#   out2 <- rstan::extract(fit)
#   Z_tilde <- out2$Z_tilde
#   Z_tilde_mean <- apply(Z_tilde, 2, mean)
#   Z_tilde_lwr <- apply(Z_tilde, 2, quantile, prob = .025)
#   Z_tilde_upr <- apply(Z_tilde, 2, quantile, prob = .975)
#   list(mean = Z_tilde_mean,
#        lwr_ci = Z_tilde_lwr,
#        upr_ci = Z_tilde_upr,
#        ci_width = Z_tilde_upr - Z_tilde_lwr)
# }
#
# get_post_theta <- function(fit){
#   out2 <- rstan::extract(fit)
#   Theta <- out2$transTheta
#   Theta_mean <- apply(Theta, 2:3, mean)
#   Theta_sd <- apply(Theta, 2:3, sd)
#   Theta_lwr <- apply(Theta, 2:3, quantile, prob = .025)
#   Theta_upr <- apply(Theta, 2:3, quantile, prob = .975)
#   list(mean = Theta_mean,
#        lwr_ci = Theta_lwr,
#        upr_ci = Theta_upr,
#        ci_width = Theta_upr - Theta_lwr)
# }
#
# get_max_idx <- function(fit){
#   out2 <- rstan::extract(fit)
#   Theta <- out2$Theta
#   max_idx <- apply(Theta, 1,
#                    function(mat){
#                      param <- exp_logit_theta_link(mat)
#                      -log(param[,2])/log(param[,3])
#                    })
#   output <- data.frame(mean = apply(max_idx, 1, mean),
#                        sd = apply(max_idx, 1, sd),
#                        lwr_ci = apply(max_idx, 1, quantile, prob = .025),
#                        upr_ci = apply(max_idx, 1, quantile, prob = .975))
#   output$ci_width = output$upr_ci - output$lwr_ci
#   output
# }
#
# get_pred_max_idx <- function(fit, pp_data){
#   # browser()
#   require(dplyr)
#   out2 <- rstan::extract(fit)
#   Z_tilde <- out2$Z_tilde
#   max_idx <- apply(Z_tilde, 1,
#                    function(row, pp_data){
#                      pp_data$pred = row
#                      max_times <- pp_data %>%
#                        group_by(post_pred_groups) %>%
#                        filter(pred == max(pred)) %>%
#                        arrange(post_pred_groups)
#                      max_times$post_pred_times
#                    }, pp_data = pp_data)
#   output <- data.frame(mean = apply(max_idx, 1, mean),
#                        sd = apply(max_idx, 1, sd),
#                        lwr_ci = apply(max_idx, 1, quantile, prob = .025),
#                        upr_ci = apply(max_idx, 1, quantile, prob = .975))
#   output$ci_width = output$upr_ci - output$lwr_ci
#   output
# }
#
#
# plot_pp_var <- function(pp_dat){
#   pp_dat$pred_mean <- Z_tilde_mean
#   pp_dat$pred_lwr <- Z_tilde_lwr
#   pp_dat$pred_upr <- Z_tilde_upr
#
#   ggplot(pp_dat, aes(x = post_pred_times, y = pred_mean)) +
#     geom_line() +
#     geom_line(aes(x = post_pred_times,
#                   y = pred_lwr,
#                   col = 'red'),
#               linetype = 'dashed') +
#     geom_line(aes(x = post_pred_times,
#                   y = pred_upr,
#                   col = 'red'),
#               linetype = 'dashed') +
#     geom_point(aes(x = post_pred_times, y = obs)) +
#     facet_wrap(~post_pred_groups_nm, scales = 'free_y')
# }
