#' Obtain the parameters of Generalized Beta (GB) prior according to the information about nonconforming proportion p
#'
#' This function calculates the parameters a and b from GB prior when \eqn{l \le p \le u}. The previous knowledge can
#' consists in the mean and variance of p or the coefficient \eqn{\epsilon} of impartiality.
#'
#'
#' @param p Vector of acceptance and rejection quality levels, \eqn{p_alpha} and \eqn{p_beta}
#' @param l Lower limit of \eqn{p}
#' @param l Upper limit of \eqn{p}
#' @param know_p List containing the mean and variance of p or the coefficient \eqn{\epsilon} of impartiality.
#' In the case of impartiality, then Pr(\eqn{p \le p_alpha})=Pr(\eqn{p \ge p_beta}).
#' @return A data.frame with the following variables of the parameters of GB prior
#' \itemize{
#'  \item{"p_alpha": }{Acceptable quality level (AQL)}
#'  \item{"p_beta": }{Rejectable quality level (RQL)}
#'  \item{"a": }{Parameter a}
#'  \item{"b": }{Parameter b}
#'  \item{"l": }{Lower limit of \eqn{p}}
#'  \item{"u": }{Upper limit of \eqn{p}}
#'  \item{"mean_p": }{Mean of \eqn{p}}
#'  \item{"var_p": }{Variance of \eqn{p}}
#'  \item{"cdf_left": }{Pr(\eqn{p \le p_alpha}) }
#'  \item{"cdf_right": }{Pr(\eqn{p \ge p_beta}) }
#' }
#' @export
#' @examples
#' p<-c(0.00654, 0.0426)
#' l<- p[1]/5
#' u<- p[2]+(p[1]-l)
#'
#' # knowledge of mean and variance of p distribution
#' know_p<-list(mean_p=p[1],var_p=((p[2]-p[1])/4)^2)
#' beta.params(p,l,u, know_p)
#'
#' # knowledge of epsilon-impartiality of p distribution
#' know_p<-list(epsilon=0.2)
#' beta.params(p,l,u, know_p)
beta.params<-function(p,l,u,know_p) {
  if( ("mean_p" %in% names(know_p)) &
      ("var_p" %in% names(know_p)) ) {

    # Objective function (variance of p)
    var_prior <- function(a,l,u,mean_p,var_p) {
      b =  a*(u - l)/(mean_p - l)-a
      cond<-a*b*(u - l)^2/((a+b)^2*(a + b + 1))-var_p
      cond
    }

    sol<-nleqslv::nleqslv(1, var_prior, control=list(ftol=1e-7, allowSingular=TRUE),
                          jacobian=TRUE, method="Newton",
                          l=l,u=u,mean_p=know_p$mean_p,var_p=know_p$var_p)

    dat<-data.frame(p_alfa=p[1],p_beta=p[2],
                    a=sol$x,b=(sol$x*(u -l)/(know_p$mean_p -l)) - sol$x,
                    l=l,u=u,
                    mean_p=know_p$mean_p,var_p=know_p$var_p,
                    cdf_left=gentrunc.beta.cdf(p[1],sol$x,(sol$x*(u -l)/(know_p$mean_p -l)) - sol$x,l,u),
                    cdf_right=1-gentrunc.beta.cdf(p[2],sol$x,(sol$x*(u -l)/(know_p$mean_p -l)) - sol$x,l,u))

  }

  if( ("epsilon" %in% names(know_p))  ) {

    beta_mean<-function(a,b,l,u) {
      l+a*(u-l)/(a+b)
    }

    beta_var<-function(a,b,l,u) {
      a*b*(u-l)^2/((a+b)^2*(a+b+1))
    }

    # Objective function
    prob_tails <- function(x,p,l,u,epsilon) {
      x1 = x[1]; x2 = x[2]
      y = numeric(2)
      y[1]<-ifelse(x1<0|x2<0,1000,gentrunc.beta.cdf(p[1],x1,x2,l,u)-epsilon)
      y[2]<-ifelse(x1<0|x2<0,1000,gentrunc.beta.cdf(p[2],x1,x2,l,u)-(1-epsilon))
      y
    }
    sol<-nleqslv::nleqslv(c(1,1), prob_tails, control=list(ftol=.0001, allowSingular=TRUE),
                          jacobian=TRUE, method="Newton",
                          p=p,l=l,u=u,epsilon=know_p$epsilon)
    dat<-data.frame(p_alfa=p[1],p_beta=p[2],
                    a=sol$x[1],b=sol$x[2],
                    l=l,u=u,
                    mean_p=beta_mean(sol$x[1],sol$x[2],l,u),
                    var_p=beta_var(sol$x[1],sol$x[2],l,u),
                    cdf_left=gentrunc.beta.cdf(p[1],sol$x[1],sol$x[2],l,u),
                    cdf_right=1-gentrunc.beta.cdf(p[2],sol$x[1],sol$x[2],l,u))
  }
  dat
}


#' Probability density function of Generalized Beta (GB) prior
#' @param p Nonconforming proportion
#' @param a Parameter \eqn{a} of GB Prior
#' @param b Parameter \eqn{b} of GB Prior
#' @param l Lower limit of \eqn{p}
#' @param u Upper limit of \eqn{p}
#' @noRd
gentrunc.beta.pdf<-function(p,a,b,l=0,u=1) {
  h<-dbeta((p-l)/(u-l),a,b)/(u-l)
  return(h)
}

#' Probability distribution function of Generalized Beta (GB) prior
#' @param p Nonconforming proportion
#' @param a Parameter \eqn{a} of GB Prior
#' @param b Parameter \eqn{b} of GB Prior
#' @param l Lower limit of \eqn{p}
#' @param u Upper limit of \eqn{p}
#' @noRd
gentrunc.beta.cdf<-function(p,a,b,l=0,u=1) {
  h<-pbeta((p-l)/(u-l),a,b)
  return(h)
}
