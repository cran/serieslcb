#' Matrix of \emph{p}-vector combinations
#'
#' Calculate a matrix of \emph{p}-vector combinations (component reliabilities) which lie in the specified interval of system reliability. Rows correspond to \emph{p}-vectors and columns correspond to components.
#'
#' @param Rs.int Interval (or single number) of total system reliability.
#' @param m Number of components.
#' @return The \eqn{2^m} by \eqn{m} matrix of \emph{p}-vector combinations.
#' @details Denote Rs.int \eqn{= (R_L, R_U)}. This function calculates all elements of the set \deqn{\Omega ' = \{(p_1, p_2, \dots , p_m): p_i = R_L^{1/m} { or } R_U^{1/m} \forall i \}}
#' @export
#' @examples
#' pm(Rs.int = c(.9, .95), m=3)
pm <- function(Rs.int, m){
  a <- min(Rs.int)^(1/m)
  b <- max(Rs.int)^(1/m)
  if(a==b) out <- rep(a, m)
  if(a!=b){
    list <- vector("list", length=m)
    for(i in 1:m) list[[i]] <- c(a,b)
    p.mat <- expand.grid(list)
    out <- as.matrix(p.mat)
  }
  if(is.vector(out)) out <- t(out) #if only one p-combo, outputs as 1 row matrix
  return(out)
}

#' Bayes' method
#'
#' Calculate a binomial series lower confidence bound using Bayes' method with a Beta prior distribution.
#'
#' @param s Vector of successes.
#' @param n Vector of sample sizes.
#' @param alpha The significance level; to calculate a 100(1-\eqn{\alpha})\% lower confidence bound.
#' @param MonteCarlo Number of samples to draw from the posterior distribution for the Monte Carlo estimate.
#' @param beta.a Shape1 parameter for the Beta prior distribution.
#' @param beta.b Shape2 parameter for the Beta prior distribution.
#' @param ... Additional arguments to be ignored.
#' @return The 100(1-\eqn{\alpha})\% lower confidence bound.
#' @export
#' @examples
#' bayes(s=c(35, 97, 59), n=c(35, 100, 60), alpha=.10, MonteCarlo=1000, beta.a=1, beta.b=1)
bayes <- function(s, n, alpha, MonteCarlo, beta.a, beta.b, ...){
  Rs <- vector(mode='numeric', length=MonteCarlo)
  for(j in 1:MonteCarlo) Rs[j] <- prod(rbeta(n=length(n), shape1=s+beta.a, shape2=n-s+beta.b))
  return(sort(Rs)[alpha*MonteCarlo])
}

#' Bayes' method (Jeffrey's prior)
#'
#' Calculate a binomial series lower confidence bound using Bayes' method with Jeffrey's prior.
#' @param s Vector of successes.
#' @param n Vector of sample sizes.
#' @param alpha The significance level; to calculate a 100(1-\eqn{\alpha})\% lower confidence bound.
#' @param ... Additional arguments to be ignored.
#' @param MonteCarlo Number of samples to draw from the posterior distribution for the Monte Carlo estimate.
#' @return The 100(1-\eqn{\alpha})\% lower confidence bound.
#' @export
#' @examples
#' bayes_jeffreys(s=c(35, 97, 59), n=c(35, 100, 60), alpha=.10, MonteCarlo=1000)
bayes_jeffreys <- function(s, n, alpha, MonteCarlo, ...){
  Rs <- vector(mode='numeric', length=MonteCarlo)
  for(j in 1:MonteCarlo) Rs[j] <- prod(rbeta(n=length(n), shape1=s+.5, shape2=n-s+.5))
  return(sort(Rs)[alpha*MonteCarlo])
}

#' Bayes' method (Uniform prior)
#'
#' Calculate a binomial series lower confidence bound using Bayes' method with a uniform prior distribution.
#'
#' @param s Vector of successes.
#' @param n Vector of sample sizes.
#' @param alpha The significance level; to calculate a 100(1-\eqn{\alpha})\% lower confidence bound.
#' @param MonteCarlo Number of samples to draw from the posterior distribution for the Monte Carlo estimate.
#' @param ... Additional arguments to be ignored.
#' @return The 100(1-\eqn{\alpha})\% lower confidence bound.
#' @export
#' @examples
#' bayes_uniform(s=c(35, 97, 59), n=c(35, 100, 60), alpha=.10, MonteCarlo=1000)
bayes_uniform <- function(s, n, alpha, MonteCarlo, ...){
  Rs <- vector(mode='numeric', length=MonteCarlo)
  for(j in 1:MonteCarlo) Rs[j] <- prod(rbeta(n=length(n), shape1=s+1, shape2=n-s+1))
  return(sort(Rs)[alpha*MonteCarlo])
}


#' Chao-Huwang method
#'
#' Calculate a binomial series lower confidence bound using Chao and Huwang's (1987) method.
#'
#' @param s Vector of successes.
#' @param n Vector of sample sizes.
#' @param alpha The significance level; to calculate a 100(1-\eqn{\alpha})\% lower confidence bound.
#' @param MonteCarlo Number of samples to draw from the posterior distribution for the Monte Carlo estimate.
#' @param ... Additional arguments to be ignored.
#' @return The 100(1-\eqn{\alpha})\% lower confidence bound.
#' @export
#' @examples
#' chao_huwang(s=c(35, 97, 59), n=c(35, 100, 60), alpha=.10, MonteCarlo=1000)
chao_huwang <- function(s, n, alpha, MonteCarlo, ...){
  f <- n-s
  p.hat <- 1-(f+.2)/(n+.2+0)
  Rs <- vector(mode="numeric", length=MonteCarlo)
  for(j in 1:MonteCarlo) Rs[j] <- prod(1-(rbinom(n=length(n), size=n, prob=1-p.hat)+.2)/(n+.2+0))
  return(sort(Rs)[alpha*MonteCarlo])
}

#' Coit's method
#'
#' Calculate a binomial series lower confidence bound using Coit's (1997) method.
#'
#' @param s Vector of successes.
#' @param n Vector of sample sizes.
#' @param alpha The significance level; to calculate a 100(1-\eqn{\alpha})\% lower confidence bound.
#' @param use.backup If TRUE, then a backup.method in the will be used for the methods with calculate LCB = 1 in the case of no failures across all components. If FALSE (default), no backup.method is used.
#' @param backup.method The backup method which is used for the methods which calculate LCB = 1 in the case of zero failures. Use function name.
#' @param ... Additional arguments to be ignored.
#' @return The 100(1-\eqn{\alpha})\% lower confidence bound.
#' @export
#' @examples
#' coit(s=c(35, 97, 59), n=c(35, 100, 60), alpha=.10)
coit <- function(s, n, alpha, use.backup=FALSE, backup.method, ...){
  Ri <- s/n
  var.Ri <- Ri*(1-Ri)/n
  Rs <- prod(Ri)
  var.Rs <- prod(Ri^2 + var.Ri) - prod(Ri^2)
  sig2 <- log(1 + var.Rs/(Rs^2))
  out <- Rs*exp(.5*sig2 - qnorm(1-alpha)*sqrt(sig2))
  if( (sum(s-n)==0) & use.backup  ) out <- backup.method(s=s, n=n, alpha=alpha)
  return(out)
}

#' Easterling's method
#'
#' Calculate a binomial series lower confidence bound using Easterling's (1972) method.
#'
#' @param s Vector of successes.
#' @param n Vector of sample sizes.
#' @param alpha The significance level; to calculate a 100(1-\eqn{\alpha})\% lower confidence bound.
#' @param ... Additional arguments to be ignored.
#' @return The 100(1-\eqn{\alpha})\% lower confidence bound.
#' @export
#' @examples
#' easterling(s=c(35, 97, 59), n=c(35, 100, 60), alpha=.10)
easterling <- function(s, n, alpha, ...){
  p.hat <- s/n
  var.p <- p.hat*(1-p.hat)/n
  for(i in 1:length(n)){
    if(var.p[i]==0){
      fifty.p <- .5^(1/n[i])
      var.p[i] <- fifty.p*(1-fifty.p)/n[i] #lower 50% confidence limit on p_i
    }
  }
  var <- vector(mode='numeric', length=length(n))
  for(j in 1:length(n)) var[j] <- prod(p.hat[-j])^2 * var.p[j]
  Rs <- prod(p.hat)
  n.star <- Rs*(1-Rs)/sum(var)
  #if(n.star==0) n.star <- min(n) #my original implementation
  n.star <- max(n.star, min(n)) #allan's suggested implementation
  return(qbeta(p=alpha, shape1=Rs*n.star, shape2=n.star - Rs*n.star + 1))
}

#' Lindstrom and Madden's method
#'
#' Calculate a binomial series lower confidence bound using Lindstrom and Madden's (1962) method.
#'
#' @param s Vector of successes.
#' @param n Vector of sample sizes.
#' @param alpha The significance level; to calculate a 100(1-\eqn{\alpha})\% lower confidence bound.
#' @param ... Additional arguments to be ignored.
#' @return The 100(1-\eqn{\alpha})\% lower confidence bound.
#' @export
#' @examples
#' lindstrom_madden(s=c(35, 97, 59), n=c(35, 100, 60), alpha=.10)
lindstrom_madden <- function(s, n, alpha, ...){
  s.star <- min(n)*prod(s/n)
  return(qbeta(p=alpha, shape1=s.star, shape2=min(n)-s.star+1))
}

#' Lindstrom and Madden's method with Agresti-Coull
#'
#' Calculate a binomial series lower confidence bound using Agresti-Coull (1998) lower confidence bound calculation in the Lindstrom and Madden's (1962) method.
#'
#' @param s Vector of successes.
#' @param n Vector of sample sizes.
#' @param alpha The significance level; to calculate a 100(1-\eqn{\alpha})\% lower confidence bound.
#' @param ... Additional arguments to be ignored.
#' @return The 100(1-\eqn{\alpha})\% lower confidence bound.
#' @export
#' @examples
#' lindstrom_madden_AC(s=c(35, 97, 59), n=c(35, 100, 60), alpha=.10)
lindstrom_madden_AC <- function(s, n, alpha, ...){
  s.star <- min(n)*prod(s/n)
  nm <- min(n)
  s.star2 <- s.star + .5*qnorm(1-alpha)^2
  nm2 <- nm + qnorm(1-alpha)^2
  Rs <- s.star2/nm2
  return(Rs - qnorm(1-alpha)*sqrt(Rs*(1-Rs)/nm2))
}


#' Normal approximation method
#'
#' Calculate a binomial series lower confidence bound using a normal approximation with MLE estimates.
#'
#' @param s Vector of successes.
#' @param n Vector of sample sizes.
#' @param alpha The significance level; to calculate a 100(1-\eqn{\alpha})\% lower confidence bound.
#' @param use.backup If TRUE, then a backup.method in the will be used for the methods with calculate LCB = 1 in the case of no failures across all components. If FALSE (default), no backup.method is used.
#' @param backup.method The backup method which is used for the methods which calculate LCB = 1 in the case of zero failures. Use function name.
#' @param ... Additional arguments to be ignored.
#' @return The 100(1-\eqn{\alpha})\% lower confidence bound.
#' @export
#' @examples
#' normal_approximation(s=c(35, 97, 59), n=c(35, 100, 60), alpha=.10)
normal_approximation <- function(s, n, alpha, use.backup=FALSE, backup.method, ...){
  phat <- s/n
  var.p <- phat*(1-phat)/n
  var <- vector(mode='numeric', length=length(n))
  for(j in 1:length(n)) var[j] <- prod(phat[-j])^2 * var.p[j]
  out <- prod(phat) + qnorm(alpha)*sqrt(sum(var))
  if( (sum(s-n)==0) & use.backup  ) out <- backup.method(s=s, n=n, alpha=alpha)
  return(out)
}


#' Lagrange multiplier in Madansky's method
#'
#' This function is called in the \code{madansky()} function to solve for the Lagrange multipliers.
#'
#' @param lam The value of the Lagrange multiplier
#' @param s Vector of successes.
#' @param n Vector of sample sizes.
#' @param alpha The significance level; to calculate a 100(1-\eqn{\alpha})\% lower confidence bound.
#' @export
madansky.fun <- function(lam, s, n, alpha){
  return(sum(s*log(1 - lam/s) - n*log(1 - lam/n)) + qchisq(p=1-2*alpha, df=1, lower.tail=TRUE)/2)
}


#' Madansky's method
#'
#' Calculate a binomial series lower confidence bound using Madansky's (1965) method.
#'
#' @param s Vector of successes.
#' @param n Vector of sample sizes.
#' @param alpha The significance level; to calculate a 100(1-\eqn{\alpha})\% lower confidence bound.
#' @param use.backup If TRUE, then a backup.method in the will be used for the methods with calculate LCB = 1 in the case of no failures across all components. If FALSE (default), no backup.method is used.
#' @param backup.method The backup method which is used for the methods which calculate LCB = 1 in the case of zero failures. Use function name.
#' @param ... Additional arguments to be ignored.
#' @return The 100(1-\eqn{\alpha})\% lower confidence bound. Note that if there are zero observed
#'     failures across all components, the output is LCB = 0.
#' @export
#' @examples
#' madansky(s=c(35, 97, 59), n=c(35, 100, 60), alpha=.10)
madansky <- function(s, n, alpha, use.backup=FALSE, backup.method, ...){
  # this returns a one-sided lower LCB only
  # looks like we're doing lindstrom-madden in the case of zero failures across all components...
  out <- 1 #### this is what will be ouput if there are zero failures across all components!
  ss <- NULL
  if(!isTRUE(all.equal(s, n))){
    ind <- !(s==n)
    ss <- s[ind]
    nn <- n[ind]
  }

  if(length(ss)>0){
    if(sum(ss==nn)<length(nn)){
      lam2 <- uniroot(f=madansky.fun, interval=c(0, min(ss)), s=ss, n=nn, alpha=alpha)$root
      out <- prod((ss-lam2)/(nn-lam2))
    }
  }
  if( (sum(s-n)==0) & use.backup  ) out <- backup.method(s=s, n=n, alpha=alpha)
  return(out)
}

#' Myhre and Rennie (modified ML) method
#'
#' Calculate a binomial series lower confidence bound using the Myhre-Rennie (modified ML) method (1986).
#'
#' @param s Vector of successes.
#' @param n Vector of sample sizes.
#' @param alpha The significance level; to calculate a 100(1-\eqn{\alpha})\% lower confidence bound.
#' @param use.backup If TRUE, then a backup.method in the will be used for the methods with calculate LCB = 1 in the case of no failures across all components. If FALSE (default), no backup.method is used.
#' @param backup.method The backup method which is used for the methods which calculate LCB = 1 in the case of zero failures. Use function name.
#' @param ... Additional arguments to be ignored.
#' @return The 100(1-\eqn{\alpha})\% lower confidence bound.
#' @export
#' @examples
#' myhre_rennie1(s=c(35, 97, 59), n=c(35, 100, 60), alpha=.10)
myhre_rennie1 <- function(s, n, alpha, use.backup=FALSE, backup.method, ...){
  p <- s/n
  star <- n*(mean(s) + s)/(mean(s)/mean(p) + n)
  out <- madansky(s=star, n=n, alpha=alpha)
  if( (sum(s-n)==0) & use.backup  ) out <- backup.method(s=s, n=n, alpha=alpha)
  return(out)
}

#' Function of \eqn{\beta} in the Myhre-Rennie 2 method
#'
#' This function is called in \code{myhre_rennie2()} function to solve for the \eqn{\beta} value.
#'
#' @param beta The value of \eqn{\beta}.
#' @param s Vector of successes.
#' @param n Vector of sample sizes.
#' @export
mr.fun <- function(beta, s, n){
  prod((beta + s)/(mean(n) + n)) - prod(s/n)
}

#' Myhre and Rennie (reliability invariant) method
#'
#' Calculate a binomial series lower confidence bound using the Myhre-Rennie (reliability invariant) method (1986).
#'
#' @param s Vector of successes.
#' @param n Vector of sample sizes.
#' @param alpha The significance level; to calculate a 100(1-\eqn{\alpha})\% lower confidence bound.
#' @param use.backup If TRUE, then a backup.method in the will be used for the methods with calculate LCB = 1 in the case of no failures across all components. If FALSE (default), no backup.method is used.
#' @param backup.method The backup method which is used for the methods which calculate LCB = 1 in the case of zero failures. Use function name.
#' @param ... Additional arguments to be ignored.
#' @return The 100(1-\eqn{\alpha})\% lower confidence bound.
#' @export
#' @examples
#' myhre_rennie2(s=c(35, 97, 59), n=c(35, 100, 60), alpha=.10)
myhre_rennie2 <- function(s, n, alpha, use.backup=FALSE, backup.method, ...){
  beta <- ifelse(mr.fun(beta=mean(n), s=s, n=n)==0, mean(n),
                 uniroot(f=mr.fun, interval=c(0, max(n)), s=s, n=n)$root)
  star <- n*( (beta + s)/( mean(n) + n) )
  star <- ifelse(star>s, s, star)
  out <- madansky(s=star, n=n, alpha=alpha)
  if( (sum(s-n)==0) & use.backup  ) out <- backup.method(s=s, n=n, alpha=alpha)
  return(out)
}


#' Function to calculate the restricted sum in the Mann-Grubbs method.
#'
#' Calculate the restricted sum in the Mann-Grubbs method.
#'
#' @param s Vector of successes.
#' @param n Vector of sample sizes.
#' @return The restricted sum.
mann_grubbs_sum <- function(s, n){
  out <- NULL
  A1 <- NULL
  A2 <- NULL #the one to always try

  min.n <- min(n)
  f <- n-s
  ind.z <- which(f==0) #indices of zero failures
  ind.f <- which(f>0) #indices with failures

  if(length(ind.z)==0) out <- min.n*sum(1/n)
  if(length(ind.z)==1){
    ni1 <- n
    ni1[ind.z] <- min.n
    A1 <- min.n*sum(1/ni1)

    ni2 <- n[-ind.z]
    A2 <- min.n*sum(1/ni2)

    out <- c(A1, A2)
  }

  if(length(ind.z)>1){
    ni1 <- n[-ind.z]
    ni1 <- c(ni1, min.n)

    #check if another component with size n.min has a failure
    if( sum(n[ind.f]==min.n)>0 ) ni1 <- n[-ind.z]
    A1 <- min.n*sum(1/ni1)

    ni2 <- n[-ind.z]
    A2 <- min.n*sum(1/ni2)

    out <- c(A1, A2)
  }
  return(out)
}


#' Function to calculate the LCB in the Mann-Grubbs method.
#'
#' Calculate the LCB in the Mann-Grubbs method.
#'
#' @param s Vector of successes.
#' @param n Vector of sample sizes.
#' @param A The restricted sum, as caclulated by the mann_grubbs_sum() function.
#' @param alpha The significance level; to calculate a 100(1-\eqn{\alpha})\% lower confidence bound.
#' @return The LCB for the Mann-Grubbs method.
mann_grubbs_calc <- function(s, n, A, alpha){
  min.n <- min(n)
  phat <- 1 - prod(s/n)
  n0 <- min.n*(1 - .5*phat*phat)*(1 - .5*phat)
  m1 <- .5*(1 + 1/A)/n0  +  phat/(1 - .5*phat)
  v <- .5*(1+1/A)*m1/n0

  m <- m1
  if(m1>.5){
    n.star <- min.n*sum(1/n) / (.5/min.n + .5*sum(1/n))
    r.star <- n.star*phat
    m <- digamma(n.star + 1) - digamma(n.star - r.star)
    v <- -trigamma(n.star + 1) + trigamma(n.star - r.star)
  }

  inside <- (1 - v/((3*m)^2) + qnorm(1-alpha)*sqrt(v)/(3*m))^3
  LCB <- exp(-m * inside)
  return(LCB)
}


#' Mann and Grubb's method
#'
#' Calculate a binomial series lower confidence bound using Mann and Grubb's (1974) method.
#'
#' @param s Vector of successes.
#' @param n Vector of sample sizes.
#' @param alpha The significance level; to calculate a 100(1-\eqn{\alpha})\% lower confidence bound.
#' @param ... Additional arguments to be ignored.
#' @return The 100(1-\eqn{\alpha})\% lower confidence bound.
#' @export
#' @examples
#' mann_grubbs(s=c(35, 97, 59), n=c(35, 100, 60), alpha=.10)
mann_grubbs <- function(s, n, alpha, ...){
  foo <- sapply(X=mann_grubbs_sum(s=s, n=n), FUN=mann_grubbs_calc, s=s, n=n, alpha=alpha)
  return(min(foo))
}

#' Nishime's method
#'
#' Calculate a binomial series lower confidence bound using Nishime's (1959) method.
#'
#' @param s Vector of successes.
#' @param n Vector of sample sizes.
#' @param alpha The significance level; to calculate a 100(1-\eqn{\alpha})\% lower confidence bound.
#' @param ... Additional arguments to be ignored.
#' @return The 100(1-\eqn{\alpha})\% lower confidence bound.
#' @export
#' @examples
#' nishime(s=c(35, 97, 59), n=c(35, 100, 60), alpha=.10)
nishime <- function(s, n, alpha, ...){
  Li <- qbeta(p=alpha, shape1=s, shape2=n-s+1)
  L.prime <- qbeta(p=alpha, shape1=sum(s), shape2=sum(n)-sum(s)+1)
  return( ((L.prime/mean(Li))^length(n)) * prod(Li) )
}

#' Rice and Moore's method
#'
#' Calculate a binomial series lower confidence bound using Rice and Moore's (1983) method.
#'
#' @param s Vector of successes.
#' @param n Vector of sample sizes.
#' @param alpha The significance level; to calculate a 100(1-\eqn{\alpha})\% lower confidence bound.
#' @param MonteCarlo Number of samples to draw from the posterior distribution for the Monte Carlo estimate.
#' @param f.star The number of psuedo-failures to use for a component that exhibits zero observed
#'     failures. The default value is from the log-gamma procedure proposed by Gatliffe (1976), and
#'     is the value used by Rice and Moore.
#' @param ... Additional arguments to be ignored.
#' @return The 100(1-\eqn{\alpha})\% lower confidence bound.
#' @export
#' @examples
#' rice_moore(s=c(35, 97, 59), n=c(35, 100, 60), alpha=.10, MonteCarlo=1000)
rice_moore <- function(s, n, alpha, MonteCarlo, f.star=1.5 - min(n) + .5*sqrt((3-2*min(n))^2 - 4*(min(n)-1)*log(alpha)*qchisq(p=alpha, df=2)), ...){
  f <- n-s
  f <- ifelse(f==0, f.star, f)
  p.hat <- 1-f/n
  sd <- sqrt(p.hat*(1-p.hat)/n)
  Rs <- vector(mode="numeric", length=MonteCarlo)
  for(j in 1:MonteCarlo) Rs[j] <- prod( pmax(0, pmin(1, rnorm(n=length(n), mean=p.hat, sd=sd))) )
  return(sort(Rs)[alpha*MonteCarlo])
}

#' Root Mean Square Error
#'
#' Calculate the root mean squared errors of the LCB's from the true system reliability. A measure of spread.
#'
#' @param LCB Vector of LCB's.
#' @param R The true system reliability .
#' @return The root mean squared error of the LCB's from the true system reliability.
#' @export
rmse.LCB <- function(LCB, R){
  sqrt(sum((LCB - R)^2))
}



#' Launch Shiny App
#'
#' Launches an instance of an R Shiny App, which runs locally on the user's computer.
#'
#' @details If the "Download Histograms" button does not work, it can be fixed by launching the Shiny App on your local browser. This can be done by clicking on "Open in Browser" located at the top of your Shiny App. This seems to be an issue with the Download Handler that Shiny uses.
#' @param MonteCarlo The number of Monte Carlo samples to take. E.g. In a Bayesian method, how many samples to take from a posterior distribution to estimate the lower \eqn{\alpha}-th quantile. The default value is 1000.
#' @param use.backup If TRUE (default), then a backup.method in the will be used for the methods with calculate LCB = 1 in the case of no failures across all components. If FALSE, no backup.method is used.
#' @param backup.method The backup method which is used for the methods which calculate LCB = 1 in the case of zero failures. The default is lindstrom_madden_AC.
#' @export
#' @example
#' launch_app(MonteCarlo=1700)
launch_app <- function(MonteCarlo = 1000, use.backup = TRUE, backup.method = lindstrom_madden_AC){
  shiny::shinyOptions(MonteCarlo = MonteCarlo, use.backup = use.backup, backup.method = backup.method)
  shiny::runApp(appDir = system.file("app.R", package = "serieslcb"),
                display.mode = "normal")
}
