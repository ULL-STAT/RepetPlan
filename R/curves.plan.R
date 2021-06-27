#' Compute the operating characteristic curve of a censored repetitive sampling plan.
#'
#' This function calculates the OC curve that provides the probability of lot acceptance function
#'
#' @param n Sample size
#' @param k Vector (k_r,k_a) of rejection and acceptance constants
#' @param p Single value or vector of nonconforming proportion
#' @param asvar List with the asymptotical variance-covariance matrix of MLE estimators of location and scale lifetime distribution parameters. See `asympt.var()` for more details.
#' @return The probability of lot acceptance of the censored repetitive plan.
#' @export
#' @examples
#' q<- 0.1
#' asvar<-asympt.var(q,"normal")
#' oc.curve(11,c(1.963213, 2.214167),0.067,asvar)
oc.curve<-function(n,k,p,asvar) {

  oc.func<- function(arg1) {
    prob.acc<- probab(n,k[2],arg1,asvar)
    prob.rej<-1-probab(n,k[1],arg1,asvar)
    prob.acc_last<- probab(n,(k[1]+k[2])/2,arg1,asvar)
    prob.acc/(prob.acc+prob.rej)

  }
  oc<-Vectorize(oc.func, vectorize.args='arg1')
  oc(p)
}


#' Compute the average sample number of a censored repetitive sampling plan.
#'
#' This function calculates the ASN of the sampling plan coresponding to a quality level \eqn{p}
#'
#' @param n Sample size
#' @param k Vector \eqn{(k_r,k_a)} of rejection and acceptance constants.
#' @param p Single value or vector of nonconforming proportion
#' @param asvar List with the asymptotical variance-covariance matrix of MLE estimators of location and scale lifetime distribution parameters. See `asympt.var()` for more details.
#' @return The ASN of the sampling plan.
#' @export
#' @examples
#' q<- 0.1
#' asvar<-asympt.var(q,"normal")
#' asn.curve(11,c(1.963213, 2.214167),0.067,asvar)
asn.curve<-function(n,k,p,asvar) {

  func<- function(arg1) {
    prob.acc<- probab(n,k[2],arg1,asvar)
    prob.rej<-1-probab(n,k[1],arg1,asvar)
    n/(prob.acc+prob.rej)
  }
  asn<-Vectorize(func, vectorize.args='arg1')
  asn(p)
}



#' Quantiles function of location and scale distributions
#' @param p Single value or vector of nonconforming proportion
#' @param dist Location and scale distribution name ("normal"= normal, "ev" = extreme-value).
#' Default value is "normal".
#' @noRd
invprob<-function(p,dist) {
  if(dist=="ev") invprob<-log(-log(1-p))
  if(dist=="normal") invprob<-qnorm(p)
  invprob
}

#' Delta function of location and scale distributions
#' @param n Sample size
#' @param k Vector \eqn{(k_r,k_a)} of rejection and acceptance constants.
#' @param p Single value or vector of nonconforming proportion
#' @param asvar List with the asymptotical variance-covariance matrix of MLE estimators of location and scale lifetime distribution parameters. See `asympt.var()` for more details.
#' @noRd
delta <- function(n,k,p,asvar) {
  delta<-sqrt(n)*(invprob(p,asvar$dist)+k)/sqrt(asvar$mat[1,1]-2*k*asvar$mat[1,2]+(k^2)*asvar$mat[2,2])
  delta
}


#' Derivative of delta function of location and scale distributions
#' @param n Sample size
#' @param k Vector \eqn{(k_r,k_a)} of rejection and acceptance constants
#' @param p Single value or vector of nonconforming proportion
#' @param asvar List with the asymptotical variance-covariance matrix of MLE estimators of location and scale lifetime distribution parameters. See `asympt.var()` for more details.
#' @noRd
deriv.delta <- function(n,k,p,asvar) {
  deriv.delta<-sqrt(n)*(1/sqrt(asvar$mat[1,1]-2*k*asvar$mat[1,2]+(k^2)*asvar$mat[2,2]))*
    (1 - 0.5*(invprob(p,asvar$dist)+k)*sqrt(-2*asvar$mat[1,2]+2*k*asvar$mat[2,2])/(asvar$mat[1,1]-2*k*asvar$mat[1,2]+(k^2)*asvar$mat[2,2]))

  deriv.delta
}


#' Probability of lot acceptance for location and scale distributions
#' @param n Sample size
#' @param k Vector \eqn{(k_r,k_a)} of rejection and acceptance constants
#' @param p Single value or vector of nonconforming proportion
#' @param asvar List with the asymptotical variance-covariance matrix of MLE estimators of location and scale lifetime distribution parameters. See `asympt.var()` for more details.
#' @noRd
probab<-function(n,k,p,asvar) {
  return(1-pnorm(delta(n,k,p,asvar)))
}
