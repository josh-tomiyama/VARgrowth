CreateLogLikelihood <- function(ll, data, params, ...){
  list(logLik = ll(data, params),
       data = data,
       params = params)
}

#' LogLikelihood calculation for VARGrowthFit Posterior
#'
#' Calculate the log likelihood of the data for each of the posterior samples generated. For Hierarchical models,
#' there is an option whether to return the conditional or marginal (over theta parameters) log likelihood.
#'
#' @param object A VARGrowthFit Object
#' @param type Whether to return the conditional or the marginal likelihood
#' @param ... Additional Args
#'
#' @return S by C by N array where S is the size of the posterior sample, C the number of chains, and N the number of fitted data points
#' @export
#' @importFrom loo extract_log_lik
#'
#' @examples
log_lik.VARGrowthFit <- function(object, type = c("cond", "marg"), ...){
  if(length(type) > 1){
    type <- type[1]
  }

  if(type == 'cond'){
    ll <- loo::extract_log_lik(object$fit, merge_chains = FALSE)
  }else if(type == 'marg'){
    if(object$ThetaTrend$model == "LinearTrendModel"){
      ll <- loo::extract_log_lik(object$fit,
                                 parameter_name = "marg_log_lik",
                                 merge_chains = FALSE)
    }else if(object$ThetaTrend$model == "DeterministicLinearTrendModel"){
      message("The marginalized distribution for the Deterministic Linear Trend model is not meaningful, conditional distribution returned.")
      ll <- loo::extract_log_lik(object$fit, merge_chains = FALSE)
    }else{
      stop("Marginal distribution not implented for this trend model")
    }
    ### should really see if there is a more efficient way to do this.

  }else{
    stop("log likelihood type is not available")
  }

}


#' Conditional Deviance Information Criterion calculation for VARGrowthFit Posterior
#'
#' Calculate the conditional deviance information criterion (cdic) for a VARGrowthFit object. Lower values indicated a better fit.
#'
#' @param object A VARGrowthFit Object
#'
#' @return A VARGrowthCriteria list object with cdic as the criteria, deviance of the posterior mean, mean of the deviances, and effective number of parameters.
#' @export
#' @importFrom posterior as_draws_array
#' @importFrom posterior subset_draws
#' @importFrom posterior summarise_draws
#'
#' @examples
cdic <- function(object){
  ll_pointwise <- log_lik(object, "cond")
  ll <- apply(ll_pointwise, 1:2, sum)
  dev <- -2*ll
  mean_dev <- mean(dev)
  mu <- object$fit %>%
    posterior::as_draws_array() %>%
    subset_draws(variable = "mu") %>%
    summarise_draws("mean")
  var <- object$fit %>%
    posterior::as_draws_array() %>%
    subset_draws(variable = "ObsVar") %>%
    summarise_draws("mean")
  ### calc deviance of posterior mean
  dev_mean <-   -2*sum(dnorm(object$data$outcome,
                          mean = mu$mean,
                          sd = sqrt(var$mean),
                          log = TRUE))

  pd <- mean_dev - dev_mean
  cdic <- dev_mean + 2*pd

  out <- list(Criteria = cdic,
              dev_mean = dev_mean,
              mean_dev = mean_dev,
              pd = pd)
  class(out) <- c("VARGrowthCriteria", "list")
  out
}

#' Marginal Deviance Information Criterion calculation for VARGrowthFit Posterior
#'
#' Calculate the marginal deviance information criterion (mdic) for a VARGrowthFit object. The marginalization takes place over the theta values of the hierarchical model.Lower values indicate a better fit.
#'
#' @param object A VARGrowthFit Object
#'
#' @return A VARGrowthCriteria list object with cdic as the criteria, deviance of the posterior mean, mean of the deviances, and effective number of parameters.
#' @export
#' @importFrom posterior as_draws_array
#' @importFrom posterior subset_draws
#' @importFrom posterior summarise_draws
#'
#' @examples
mdic <- function(object, nsamp = 1000){
# browser()
  if(object$ThetaTrend$model == "LinearTrendModel"){
    gf <- function(theta, t){
      if(is.vector(theta)){
        theta[1]*exp(-theta[2]*theta[3]^t)
      }else if(is.matrix(theta)){
        theta[,1]*exp(-theta[,2] * theta[,3]^t)
      }else{
        stop("\n Theta input not recognized")
      }

    }

    XTheta <- object$ThetaTrend$est_input$XTheta
    betaTheta <- object$fit %>%
      posterior::as_draws_matrix() %>%
      subset_draws(variable = "betaTheta")
    # betaTheta <- matrix(c(beta1), nrow = 1000, ncol = length(c(beta1)), byrow = TRUE)

    linear_pred_theta <- apply(betaTheta, 1, function(beta, XTheta){
      beta_mat <- matrix(beta, nrow = ncol(XTheta))
      linear_pred_theta <- XTheta %*% beta_mat
    },
    XTheta = XTheta,
    simplify = FALSE)

    # betaTheta <- matrix(betaTheta$mean, nrow = ncol(XTheta))

    sigmaTheta <- object$fit %>%
      posterior::as_draws_matrix() %>%
      subset_draws(variable = "SigmaTheta")
    sigmaTheta_list <- lapply(1:nrow(sigmaTheta), function(x, st){st[x,]},
                              st = sigmaTheta)
    # sigmaTheta <- lapply(1:nrow(sigmaTheta), function(idx){sigmaTheta[idx,]})
    # sigmaTheta <- lapply(1:nrow(sigmaTheta), function(idx){c(0.02, 0.02, 0.02)^2})
    var <- object$fit %>%
      posterior::as_draws_matrix() %>%
      subset_draws(variable = "ObsVar")

    theta_samp <- mapply(
           function(lpt, st, nsamp){
             # browser()
             theta_samp <- array(dim = c(dim(lpt), nsamp))
              for(u in 1:nrow(lpt)){
              eps <- sapply(st,
                            function(sigma_sq, nsamp){
                              rnorm(nsamp, 0, sqrt(sigma_sq))
                            },
                            nsamp = nsamp)
              theta_samp[u,,] <- t(lpt[rep(u, nsamp),] + eps)
              }
             theta_samp
          },
          lpt = linear_pred_theta,
          st = sigmaTheta_list,
          nsamp = nsamp,
          SIMPLIFY = FALSE
          )

    marg_dens <- mapply(function(theta, obs_var, object, nsamp){
      # browser()
      gf <- function(theta, t){
        if(is.vector(theta)){
          theta[1]*exp(-theta[2]*theta[3]^t)
        }else if(is.matrix(theta)){
          theta[,1]*exp(-theta[,2] * theta[,3]^t)
        }else{
          stop("\n Theta input not recognized")
        }

      }

      log_pred <- mu <- matrix(data = NA, nrow = nrow(object$data), ncol = nsamp)

      for(j in 1:nsamp){
        mu[,j] <- gf(object$ThetaTransform$inv_theta_link(theta[object$data$group,,j]), object$data$time)
      }

      log_pred <- apply(mu, 2,
                        function(mu, observed, obs_var){
                          dnorm(observed, mu, sqrt(obs_var),FALSE)
                        },
                        observed = object$data$outcome,
                        obs_var = obs_var)
      log_pred_est <- log(apply(log_pred, 1, mean))

    },
    theta = theta_samp,
    obs_var = c(var),
    MoreArgs = list("object" = object,
                    "nsamp" = nsamp)
    )

    ### Calculate marginal by hand
    ll_pointwise <- log_lik(object, "marg")
    ll <- apply(ll_pointwise, 1:2, sum)
    dev <- -2*ll
    mean_dev <- mean(dev)


    # ll_pointwise <- marg_dens





    ### now calculate marginal of mean
    XTheta <- object$ThetaTrend$est_input$XTheta
    betaTheta <- object$fit %>%
      posterior::as_draws_array() %>%
      subset_draws(variable = "betaTheta") %>%
      summarise_draws("mean")

    betaTheta <- matrix(betaTheta$mean, nrow = ncol(XTheta))

    linear_pred_theta <- XTheta %*% betaTheta
    theta_samp <- array(dim = c(dim(linear_pred_theta), nsamp))

    sigmaTheta <- object$fit %>%
      posterior::as_draws_array() %>%
      subset_draws(variable = "SigmaTheta") %>%
      summarise_draws("mean")

    var <- object$fit %>%
      posterior::as_draws_array() %>%
      subset_draws(variable = "ObsVar") %>%
      summarise_draws("mean")

    gf <- function(theta, t){
      if(is.vector(theta)){
        theta[1]*exp(-theta[2]*theta[3]^t)
      }else if(is.matrix(theta)){
        theta[,1]*exp(-theta[,2] * theta[,3]^t)
      }else{
        stop("\n Theta input not recognized")
      }

    }

    for(u in 1:nrow(linear_pred_theta)){
      eps <- sapply(sigmaTheta$mean,
                    function(sigma_sq, nsamp){
                      rnorm(nsamp, 0, sqrt(sigma_sq))
                    },
                   nsamp = nsamp)
      theta_samp[u,,] <- t(linear_pred_theta[rep(u, nsamp),] + eps
)    }

    log_pred <- mu <- matrix(data = NA, nrow = nrow(object$data), ncol = nsamp)

    for(j in 1:nsamp){
      mu[,j] <- gf(object$ThetaTransform$inv_theta_link(theta_samp[object$data$group,,j]), object$data$time)
    }

    log_pred <- apply(mu, 2,
                      function(mu, observed, var){
                        dnorm(observed, mu, sqrt(var),FALSE)
                      },
                      observed = object$data$outcome,
                      var = var$mean)
    log_pred_est <- log(apply(log_pred, 1, mean))
    dev_mean <- -2*sum(log_pred_est)

    pd <- mean_dev - dev_mean
    mdic <- dev_mean + 2*pd

    out <- list(Criteria = mdic,
                dev_mean = dev_mean,
                mean_dev = mean_dev,
                pd = pd)
    class(out) <- c("VARGrowthCriteria", "list")
    out




  }else if(object$ThetaTrend$model == "DeterministicLinearTrendModel"){
    message("The marginalized distribution for the Deterministic Linear Trend model is not meaningful, cdic returned.")
    return(cdic(object))
  }else{
    stop("\nmdic not implented for this trend model type")
  }

}


# can apply loo directly to return object

GompertzConditionalLL <- function(object, stat = PostMean(), ...){
  mu <- object$fit %>%
        posterior::as_draws_array() %>%
        subset_draws(variable = "mu") %>%
        summarise_draws("mean")
  sd <- object$fit %>%
      posterior::as_draws_array() %>%
      subset_draws(variable = "ObsVar") %>%
      summarise_draws("mean")

  sum(dnorm(object$data$outcome,
        mean = mu$mean,
        sd = sqrt(sd$mean),
        log = TRUE))
}

GompertzMarginalLL <- function(object, type = "MonteCarlo"){
  #cubature package for multivariate integration
  #problems: Bounds for the integration?
  #Maybe an approximation would be good...but need to somehow quantify if approximation is good
  mu <- object$fit %>%
    posterior::as_draws_array() %>%
    subset_draws(variable = "linear_pred_theta") %>%
    summarise_draws("mean")
  var <- object$fit %>%
    posterior::as_draws_array() %>%
    subset_draws(variable = "ObsVar") %>%
    summarise_draws("mean")
  thetaVar <- object$fit %>%
    posterior::as_draws_array() %>%
    subset_draws(variable = "SigmaTheta") %>%
    summarise_draws("mean")
}
