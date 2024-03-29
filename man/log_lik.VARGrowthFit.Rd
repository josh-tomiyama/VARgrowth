% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SelectionCriterion.R, R/VARGrowth-package.R
\name{log_lik.VARGrowthFit}
\alias{log_lik.VARGrowthFit}
\alias{log_lik}
\title{LogLikelihood calculation for VARGrowthFit Posterior}
\usage{
\method{log_lik}{VARGrowthFit}(object, type = c("cond", "marg"), ...)
}
\arguments{
\item{object}{A VARGrowthFit Object}

\item{type}{Whether to return the conditional or the marginal likelihood}

\item{...}{Additional Args}
}
\value{
S by C by N array where S is the size of the posterior sample, C the number of chains, and N the number of fitted data points
}
\description{
Calculate the log likelihood of the data for each of the posterior samples generated. For Hierarchical models,
there is an option whether to return the conditional or marginal (over theta parameters) log likelihood.
}
