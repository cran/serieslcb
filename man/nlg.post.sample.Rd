% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{nlg.post.sample}
\alias{nlg.post.sample}
\title{Sampling from Posterior of Negative Log Gamma prior and Binomial data.}
\usage{
nlg.post.sample(sample.size, shape, scale, s, n)
}
\arguments{
\item{sample.size}{The number of draws from the posterior distribution.}

\item{shape}{The shape parameter for the NLG prior.}

\item{scale}{The scale parameter for the NLG prior.}

\item{s}{The number of successes for the binomial data (should be a scalar).}

\item{n}{The number of tests for the binomial data (should be a scalar).}
}
\description{
Randomly sample from the posterior distribution resulting from a NLG prior and Binomial data.
}
\examples{
nlg.post.sample(sample.size=50, shape=.2, scale=1, s=29, n=30)
}
