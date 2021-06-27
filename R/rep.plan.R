#' Determine the designs of conventional censored repetitive sampling plans.
#'
#' This function computes the designs of censored repetitive sampling plans using conventional sampling risks
#' given the requirements of the producer and consumer about maximum allowable risks and quality levels
#'
#' @import dplyr
#' @param risks Vector of producer and consumer maximum sampling risks
#' @param p Vector of acceptance and rejection quality levels
#' @param asvar List with the asymptotical variance-covariance matrix of MLE estimators of location and scale lifetime distribution parameters. See `asympt.var()` for more details.
#' @return A data.frame with the following variables of the censored repetitive design of the sampling plans.
#' \itemize{
#'  \item{"q": }{Censoring degree}
#'  \item{"n": }{Sample size}
#'  \item{"ka": }{Acceptance constant}
#'  \item{"kr": }{Rejection constant}
#'  \item{"termcd": }{Termination code of Newton-Raphson optimization}
#'  \item{"message": }{Final message of Newton-Raphson optimization}
#'  \item{"p_alpha": }{Acceptable quality level (AQL)}
#'  \item{"p_beta": }{Rejectable quality level (RQL)}
#'  \item{"dist": }{Distribution name}
#'  \item{"alpha": }{Maximum producer's risk}
#'  \item{"beta": }{Maximum consumer's risk}
#'  \item{"asn_alpha": }{ASN at AQL}
#'  \item{"asn_beta": }{ASN at RQL }
#'  \item{"asn_avg": }{Average of ASN at AQL and ASN at RQL}
#'  \item{"p_asn_max": }{Quality level in which ASN maximizes}
#'  \item{"asn_max": }{Maximum ASN over the quality level}
#' }
#' @export rep.plan
#' @examples
#' risks<-c(0.05,0.10)
#' p<-c(0.00654, 0.0426)
#' q<- 0.1
#' asvar<-asympt.var(q,"normal")
#' designs<-rep.plan(risks,p,asvar)
#'
#' optimal.design<-designs %>% group_by(q,dist,p_alpha,p_beta) %>%
#'                 filter( (abs(alpha-risks[1])<1e-05) & (abs(risks[2]-beta)<1e-05) & (termcd==1)) %>%
#'                 slice(which.min(asn_avg)) %>% arrange(q,p_alpha,p_beta) %>% as.data.frame()
#'
rep.plan<-function(risks,p,asvar,pg_bar=TRUE) {


  # function to compute ka, kr verifying the maximum sampling risks
  risks.cond <- function(x,n,risks,p,asvar) { # Objective function
    kr = x[1]; ka = x[2]
    p_alpha = p[1]; p_beta = p[2]
    alpha=risks[1]; beta=risks[2]
    cond_risk = numeric(2)
    cond<-oc.curve(n,c(kr,ka),c(p_alpha,p_beta),asvar)-c(1-alpha,beta)
    cond_risk[1]<-cond[1]*ifelse(abs(cond[1])>-1e-4 ,1000,1)
    cond_risk[2]<-cond[2]*ifelse(abs(cond[2])>-1e-4 ,1000,1)
    cond_risk
  }


  # iterative search of censored repetitive plans
  start<-single.plan(c(risks[1],risks[2]),c(p[1],p[2]),asvar)
  found<-FALSE
  des<-list()
  nstart <- round(start$n)+10 # small increment of the initial value of n of the single sampling plan
  n<-nstart
  pb <- if(pg_bar) {
    cat("rep.plan:\n")
    txtProgressBar(min = 0, max = n, style = 3)
    } else NA

  while(!found) {
    kaprox<- kval.aprox(n,risks,p,asvar)
    k <- c(kaprox$kr,kaprox$ka)

    sol<-nleqslv::nleqslv(k, risks.cond, control=list(ftol=.0001, allowSingular=TRUE),jacobian=TRUE,
                 method="Newton",n=n,risks=risks,p=p,asvar=asvar)

    if( (sol$x[1]<sol$x[2]) & (n>1) ) {
      des<-append(des,list(data.frame(q=asvar$q,n=n,kr=sol$x[1],ka=sol$x[2],
                                      termcd=sol$termcd,message=sol$message,
                                      p_alpha=p[1],p_beta=p[2],dist=asvar$dist)))
    }
    n<-n-1
    if(n<2) found<-TRUE
    if(pg_bar) setTxtProgressBar(pb, nstart-n)

  }
  if(pg_bar) close(pb)


  des<-do.call(rbind,des) %>% rowwise() %>%
    mutate(alpha=1-oc.curve(n,c(kr,ka),p_alpha,asvar),
           beta=oc.curve(n,c(kr,ka),p_beta,asvar),
           asn_alpha=asn.curve(n,c(kr,ka),p_alpha,asvar),
           asn_beta=asn.curve(n,c(kr,ka),p_beta,asvar),
           asn_avg=sum(asn_alpha,asn_beta)*0.5) %>%
    mutate(p_asn_max=asn.max(n,c(kr,ka),c(p_alpha,p_beta),asvar),
           asn_max=asn.curve(n,c(kr,ka),p_asn_max,asvar)) %>%
    as.data.frame()
  return(des)
}



#' Design of single sampling plan
#' @param risks Vector of producer and consumer maximum sampling risks
#' @param p Vector of acceptance and rejection quality levels
#' @param asvar List with the asymptotical variance-covariance matrix of MLE estimators of location and scale lifetime distribution parameters. See `asympt.var()` for more details.
#' @noRd
single.plan<-function(risks,p,asvar) {
  alpha=risks[1]; beta=risks[2]
  p_alpha = p[1]; p_beta = p[2]
  k<- (qnorm(alpha)*invprob(p_beta,asvar$dist)+qnorm(beta)*invprob(p_alpha,asvar$dist))/(-(qnorm(alpha)+qnorm(beta)))
  n <- (asvar$mat[1,1]-2*k*asvar$mat[1,2]+(k^2)*asvar$mat[2,2])*(qnorm(alpha)+qnorm(beta))^2/(invprob(p_beta,asvar$dist)-invprob(p_alpha,asvar$dist))^2
  return(list(n=n,k=k))
}


#' Compute maximum of ASN over p
#' @param n Sample size
#' @param k Vector \eqn{(k_r,k_a)} of rejection and acceptance constants.
#' @param p Vector of acceptance and rejection quality levels
#' @param asvar List with the asymptotical variance-covariance matrix of MLE estimators of location and scale lifetime distribution parameters. See `asympt.var()` for more details.
#' @noRd
asn.max <-function(n,k,p,asvar) {
  op<-optimize(function(p) asn.curve(n,k,p,asvar), c(p[1],p[2]), maximum = TRUE)
  return(op$maximum)
}


#' Compute start values of ka and kr in Newton method
#' @param n Sample size
#' @param risks Vector of producer and consumer maximum sampling risks
#' @param p Vector of acceptance and rejection quality levels
#' @param asvar List with the asymptotical variance-covariance matrix of MLE estimators of location and scale lifetime distribution parameters. See `asympt.var()` for more details.
#' @noRd
kval.aprox<-function(n,risks,p,asvar) {

  k0<-(qnorm(1-risks[1])*invprob(p[2],asvar$dist)+qnorm(1-risks[2])*invprob(p[1],asvar$dist))/
    (-qnorm(1-risks[1])-qnorm(1-risks[2]))
  k0 <- max(k0,asvar$mat[1,2]/asvar$mat[2,2]+1e-6)

  coefs1<-c(pnorm(delta(n,k0,p[1],asvar)) , dnorm(delta(n,k0,p[1],asvar))*deriv.delta(n,k0,p[1],asvar) )
  coefs2<-c(pnorm(delta(n,k0,p[2],asvar)) , dnorm(delta(n,k0,p[2],asvar))*deriv.delta(n,k0,p[2],asvar) )

  ka<- k0 + ((1-risks[1])*risks[2])/(risks[2]-(1-risks[1]))*( (coefs2[1] -(1-risks[2]))/(coefs2[2]*risks[2]) - (coefs1[1] -risks[1])/(coefs1[2]*(1-risks[1])) )

  kr<- k0 - (risks[1]*risks[2])/(risks[2]-(1-risks[1]))*( (coefs2[1] -(1-risks[2]))/(coefs2[2]*risks[2]) - (coefs1[1] -risks[1])/(coefs1[2]*(1-risks[1])) )+
    (risks[1]-coefs1[1])/(coefs1[2]*(1-risks[1]))

  return( data.frame(k0,kr,ka) )
}


#' Compute the asymptotical variance-covariance (AVAR) matrix of MLEs of location and scale parameters
#'
#' This function calculates the AVAR matrix of the maximum likelihood estimators of the location and scale
#' parameters of normal and extreme value distributions.
#'
#' @param q Censoring degree (\eqn{0 \le q < 1})
#' @param dist Location and scale distribution name ("normal"= normal, "ev" = extreme-value).
#' Default value is "normal".
#' @return A list composed by q (censoring degree), dist (distribution) and mat (AVAR matrix of parameters estimators computed as the inverse of the Fisher information matrix).
#' @export
#' @examples
#' asympt.var(0.5)
asympt.var <- function(q,dist="normal") {
  k<-matrix(rep(NA,4),nrow=2,ncol=2)

  if(dist=="ev") {
    if(q==0) {
      k[1,1:2]<-c(1,1-(-digamma(1)) )
      k[2,1:2]<-c(k[1,2],pi^2/6+(1-(-digamma(1)))^2 )
    }
    else {
      v<-function(i,j,q) {
        integrate( function(x) x^i*(log(x)^j*exp(-x)), 0, -log(q))$value
      }

      k[1,1:2]<-c(1-q, v(1,1,q)-log(-log(q))*q*log(q) )
      k[2,1:2]<-c(k[1,2], q-1-2*v(0,1,q)+2*v(1,2,q)+2*k[1,2]+(k[1,2]-v(1,1,q))*log(-log(q))  )
    }
  }

  if(dist=="normal") {
    if(q==0) {
      k[1,1:2]<-c(1,0 )
      k[2,1:2]<-c(k[1,2],2 )
    }
    else {
      lambda<- -qnorm(q)-dnorm(qnorm(q))/q
      k[1,1:2]<-c(1-q-dnorm(qnorm(q))*lambda, dnorm(qnorm(q))*(qnorm(q)*lambda-1) )
      k[2,1:2]<-c(k[1,2], 2*(1-q)- qnorm(q)*k[1,2] )
    }
  }

  var<-matrix(c(k[2,2],-k[1,2],-k[1,2],k[1,1]),nrow=2,ncol=2)/(k[1,1]*k[2,2]-k[1,2]^2)
  return(list(q=q,dist=dist,mat=var))

}
