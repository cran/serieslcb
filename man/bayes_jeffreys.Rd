% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{bayes_jeffreys}
\alias{bayes_jeffreys}
\title{Bayesian method (Jeffrey's prior)}
\usage{
bayes_jeffreys(s, n, alpha, MonteCarlo, ...)
}
\arguments{
\item{s}{Vector of successes.}

\item{n}{Vector of sample sizes.}

\item{alpha}{The significance level; to calculate a 100(1-\eqn{\alpha})\% lower confidence bound.}

\item{MonteCarlo}{Number of samples to draw from the posterior distribution for the Monte Carlo estimate.}

\item{...}{Additional arguments to be ignored.}
}
\value{
The 100(1-\eqn{\alpha})\% lower confidence bound.
}
\description{
Calculate a binomial series lower confidence bound using Bayes' method with Jeffrey's prior.
}
\examples{
bayes_jeffreys(s=c(35, 97, 59), n=c(35, 100, 60), alpha=.10, MonteCarlo=1000)
}
