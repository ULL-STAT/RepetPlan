#' Determine the design of a censored repetitive sampling plan using a GB prior to model the knowledge about \eqn{p}.
#'
#' This function computes the design of the censored repetitive sampling plan using a model prior of \eqn{p} and expected sampling risks
#' given the requirements of the producer and consumer about maximum allowable risks and quality levels
#'
#' @import dplyr
#' @param risks Vector of producer and consumer maximum sampling risks
#' @param p Vector of acceptance and rejection quality levels
#' @param asvar List with the asymptotical variance-covariance matrix of MLE estimators of location and scale lifetime distribution parameters. See `asympt.var()` for more details.
#' @param beta.pars Data frame with the GB parameters according to the available knowledge about \eqn{p}. See `beta.params()` for more details.
#' @return A data.frame with the following variables of the censored repetitive design of the sampling plans.
#' \itemize{
#'  \item{"q": }{Censoring degree}
#'  \item{"n": }{Sample size}
#'  \item{"n_low": }{Lower bound of sample size}
#'  \item{"n_up": }{Upper bound of sample size}
#'  \item{"ka": }{Acceptance constant}
#'  \item{"kr": }{Rejection constant}
#'  \item{"termcd": }{Termination code of Newton-Raphson optimization}
#'  \item{"message": }{Final message of Newton-Raphson optimization}
#'  \item{"p_alpha": }{Acceptable quality level (AQL)}
#'  \item{"p_beta": }{Rejectable quality level (RQL)}
#'  \item{"a": }{Parameter a}
#'  \item{"b": }{Parameter b}
#'  \item{"l": }{Lower limit of \eqn{p}}
#'  \item{"u": }{Upper limit of \eqn{p}}
#'  \item{"mean_p": }{Mean of \eqn{p}}
#'  \item{"var_p": }{Variance of \eqn{p}}
#'  \item{"dist": }{Distribution name}
#'  \item{"alpha": }{Maximum producer's risk}
#'  \item{"beta": }{Maximum consumer's risk}
#'  \item{"asn_alpha": }{ASN at AQL}
#'  \item{"asn_beta": }{ASN at RQL }
#'  \item{"asn_avg": }{Average of ASN at AQL and ASN at RQL}
#'  \item{"easn": }{Expected ASN}
#'  \item{"p_asn_max": }{Quality level in which ASN maximizes}
#'  \item{"asn_max": }{Maximum ASN over the quality level}
#' }
#' @export
#' @examples
#' risks<-c(0.05,0.10)
#' p<-c(0.00654, 0.0426)
#' q<- 0.1
#' asvar<-asympt.var(q,"normal")
#' l<- p[1]/5
#' u<- p[2]+(p[1]-l)
#'
#' # GB parameters for a knowledge of mean and variance of p distribution
#' know_p<-list(mean_p=p[1],var_p=((p[2]-p[1])/4)^2)
#' beta.parms<-beta.params(p,l,u, know_p)
#'
#' designs<-repGBprior.plan(risks,p,asvar, beta.parms)
#'
#' optimal.design<-designs %>% group_by(q,dist,p_alpha,p_beta) %>%
#'                 filter( (abs(alpha-risks[1])<1e-05) & (abs(risks[2]-beta)<1e-05) & (termcd==1)) %>%
#'                 group_by(q,p_alpha,p_beta,a,b,l,u,dist) %>%
#'                 mutate(easn_min=min(easn)) %>%
#'                 slice(which.min(easn)) %>% as.data.frame()
#'
repGBprior.plan<-function(risks,p,asvar,beta.pars,pg_bar=TRUE) {

  # function to compute lower and upper bounds of n
  bounds.n <-function(risks,p,asvar,pg_bar=FALSE) {
    tol<-1e-05
    start<-rep.plan(risks,p,asvar)

    start<-start %>% filter( (abs(alpha-risks[1])<tol) & (abs(risks[2]-beta)<tol) & (termcd==1))
    start_low<-start %>% slice(which.min(n)) %>% as.data.frame()
    start_up<-start %>% slice(which.max(n)) %>% as.data.frame()
    return(list(low=start_low[,c("n","ka","kr")],up=start_up[,c("n","ka","kr")]))
  }
  n.up<-bounds.n(risks,c(p[1],p[2]),asvar)$up
  n.low<-if(beta.pars$l>0 & beta.pars$u<1) bounds.n(risks,c(beta.pars$l,beta.pars$u),asvar)$low else data.frame(n=2,ka=NaN,kr=NaN)

  # function to compute expected sampling risks
  expected.risk<-function(n,k,p,beta.pars,asvar,type) {
    tol <- 1e-10
    func_prod<- function(arg1) {
      integrate(function(t) (1-oc.curve(n,k,t,asvar))*gentrunc.beta.pdf(t,beta.pars$a,beta.pars$b,beta.pars$l,beta.pars$u),
                lower = beta.pars$l+tol, upper = arg1)$value/gentrunc.beta.cdf(arg1,beta.pars$a,beta.pars$b,beta.pars$l,beta.pars$u)
    }
    risk_prod<-Vectorize(func_prod, vectorize.args='arg1')

    func_cons<- function(arg1) {
      integrate(function(t) oc.curve(n,k,t,asvar)*gentrunc.beta.pdf(t,beta.pars$a,beta.pars$b,beta.pars$l,beta.pars$u),
                lower = arg1, upper = beta.pars$u)$value/(1-gentrunc.beta.cdf(arg1,beta.pars$a,beta.pars$b,beta.pars$l,beta.pars$u))
    }
    risk_cons<-Vectorize(func_cons, vectorize.args='arg1')

    switch(type,
           prod=risk_prod(p),
           cons=risk_cons(p))
  }


  # function to compute ka, kr verifying the maximum expected sampling risks
  risks.cond <- function(x,n,risks,p,beta.pars,asvar) {
    cond<-c(expected.risk(n,c(x[1], x[2]),p[1],beta.pars,asvar,type="prod"),
            expected.risk(n,c(x[1], x[2]),p[2],beta.pars,asvar,type="cons")) - c(risks[1],risks[2])
    cond
  }

  #function to compute the expected ASN
  easn.curve<-function(n,k,beta.pars,asvar) {
    easn<-tryCatch({
      integrate(function(t) asn.curve(n,k,t,asvar)*gentrunc.beta.pdf(t,beta.pars$a,beta.pars$b,beta.pars$l,beta.pars$u),
                lower = beta.pars$l, upper = beta.pars$u)$value
    },
    error=function(e) return(NA)
    )
    return(easn)
  }



  #iterative search of censored repetitive plans
  found<-FALSE
  des<-list()
  n<-n.up$n

  pb <- if(pg_bar) {
    cat("repGBprior.plan:\n")
    txtProgressBar(min = 0, max = n, style = 3)
  } else NA


  while(!found) {
    kaprox<- kval.aprox(n,risks,p,asvar)
    k <- c(kaprox$kr,kaprox$ka)

    sol<-nleqslv::nleqslv(k, risks.cond, control=list(ftol=.0001, allowSingular=TRUE),jacobian=TRUE,
                          method="Newton",n=n,risks=risks,p=p,beta.pars=beta.pars,asvar=asvar)

    if( (sol$x[1]<sol$x[2]) & (n>1) ) {
      des<-append(des,list(data.frame(q=asvar$q,n=n,n_low=n.low$n,n_up=n.up$n,
                                      kr=sol$x[1],ka=sol$x[2],
                                      termcd=sol$termcd,message=sol$message,
                                      p_alpha=p[1],p_beta=p[2],

                                      a=beta.pars$a,b=beta.pars$b,l=beta.pars$l,u=beta.pars$u,
                                      mean_p=beta.pars$mean,
                                      var_p=beta.pars$var_p,
                                      dist=asvar$dist)))
    }
    n<-n-1
    if(n<2) found<-TRUE
    if(pg_bar) setTxtProgressBar(pb, n.up$n-n)

  }
  if(pg_bar) close(pb)

  des<-do.call(rbind,des) %>% rowwise() %>%
    mutate(
      alpha=expected.risk(n,c(kr, ka),p_alpha,beta.pars,asvar,type="prod"),
      beta=expected.risk(n,c(kr, ka),p_beta,beta.pars,asvar,type="cons"),
      asn_alpha=asn.curve(n,c(kr,ka),p_alpha,asvar),
      asn_beta=asn.curve(n,c(kr,ka),p_beta,asvar),
      asn_avg=sum(asn_alpha,asn_beta)*0.5,
      easn=easn.curve(n,c(kr,ka),beta.pars,asvar)) %>%
    mutate(p_asn_max=asn.max(n,c(kr,ka),c(p_alpha,p_beta),asvar),
           asn_max=asn.curve(n,c(kr,ka),p_asn_max,asvar)) %>%
    as.data.frame()

  return(des)
}
