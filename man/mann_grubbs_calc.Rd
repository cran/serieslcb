% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{mann_grubbs_calc}
\alias{mann_grubbs_calc}
\title{Function to calculate the LCB in the Mann-Grubbs method.}
\usage{
mann_grubbs_calc(s, n, A, alpha)
}
\arguments{
\item{s}{Vector of successes.}

\item{n}{Vector of sample sizes.}

\item{A}{The restricted sum, as caclulated by the mann_grubbs_sum() function.}

\item{alpha}{The significance level; to calculate a 100(1-\eqn{\alpha})\% lower confidence bound.}
}
\value{
The LCB for the Mann-Grubbs method.
}
\description{
Calculate the LCB in the Mann-Grubbs method.
}
