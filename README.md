README
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

# RepetPlan

<!-- badges: start -->

[![License:
MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://cran.r-project.org/web/licenses/MIT)
[![GitHub
release](https://img.shields.io/github/release/ULL-STAT/RepetPlan.svg)](https://gitHub.com/ULL-STAT/RepetPlan/releases/)
[![Github all
releases](https://img.shields.io/github/downloads/ULL-STAT/RepetPlan/total.svg)](https://gitHub.com/ULL-STAT/RepetPlan/releases/)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5035753.svg)](https://doi.org/10.5281/zenodo.5035753)
<!-- badges: end -->

RepetPlan is an R package developed to obtain failured-censored
repetitive group sampling plans.

# Important note

THE GITHUB REPO IS ANONYMIZED UNTIL FINISHING REVIEWING PROCESS OF PAPER. THE R CODE IS ALSO AVAILABLE IN ZENODO REPOSITORY FOR INSTALLATION (https://zenodo.org/badge/latestdoi/380837603). 

## Installation

To install the current version of the code from GitHub:

``` r
if(!require(devtools)){install.packages("devtools")}  #install if needed
devtools::install_github("ULL-STAT/RepetPlan")
```

## Load and Help

To load the RepetPlan package:

``` r
# load library and dependant libraries 
library(RepetPlan)
```

To see all available functions in the package use the command below

``` r
# To get index of help on all functions
help(package="RepetPlan")
```

## Examples

### Design of repetitive group sampling plans using conventional sampling risks

Suppose that *T* represents a lifetime variable and *X* = log (*T*)
follows a log-location and scale distribution. This is an example which
shows how to determine the designs of **conventional** censored
repetitive sampling plans for the given requirements of maximum risks
and quality levels

``` r
risks<-c(0.05,0.10)     #vector of producer and consumer maximum sampling risks
p<-c(0.00654, 0.0426)   #vector of acceptance and rejection quality levels
q<- 0.1                 #censoring degree
asvar<-asympt.var(q,"normal")    #asymptotical variance-covariance matrix of MLE estimators of location and scale paramters
designs<-rep.plan(risks,p,asvar) #designs satisfying the previous requirements
```

The first designs returned by the function *rep.plan()* are
<table class="table table-striped table-hover table-condensed" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:right;">
q
</th>
<th style="text-align:right;">
n
</th>
<th style="text-align:right;">
kr
</th>
<th style="text-align:right;">
ka
</th>
<th style="text-align:right;">
termcd
</th>
<th style="text-align:left;">
message
</th>
<th style="text-align:right;">
p_alpha
</th>
<th style="text-align:right;">
p_beta
</th>
<th style="text-align:left;">
dist
</th>
<th style="text-align:right;">
alpha
</th>
<th style="text-align:right;">
beta
</th>
<th style="text-align:right;">
asn_alpha
</th>
<th style="text-align:right;">
asn_beta
</th>
<th style="text-align:right;">
asn_avg
</th>
<th style="text-align:right;">
p_asn_max
</th>
<th style="text-align:right;">
asn_max
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
49
</td>
<td style="text-align:right;">
2.053225
</td>
<td style="text-align:right;">
2.055370
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
Function criterion near zero
</td>
<td style="text-align:right;">
0.00654
</td>
<td style="text-align:right;">
0.0426
</td>
<td style="text-align:left;">
normal
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
49.04823
</td>
<td style="text-align:right;">
49.06290
</td>
<td style="text-align:right;">
49.05557
</td>
<td style="text-align:right;">
0.0188643
</td>
<td style="text-align:right;">
49.16223
</td>
</tr>
<tr>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
48
</td>
<td style="text-align:right;">
2.048891
</td>
<td style="text-align:right;">
2.060242
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
Function criterion near zero
</td>
<td style="text-align:right;">
0.00654
</td>
<td style="text-align:right;">
0.0426
</td>
<td style="text-align:left;">
normal
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
48.25590
</td>
<td style="text-align:right;">
48.32952
</td>
<td style="text-align:right;">
48.29271
</td>
<td style="text-align:right;">
0.0188293
</td>
<td style="text-align:right;">
48.84398
</td>
</tr>
<tr>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
47
</td>
<td style="text-align:right;">
2.044397
</td>
<td style="text-align:right;">
2.065329
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
Function criterion near zero
</td>
<td style="text-align:right;">
0.00654
</td>
<td style="text-align:right;">
0.0426
</td>
<td style="text-align:left;">
normal
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
47.47311
</td>
<td style="text-align:right;">
47.60146
</td>
<td style="text-align:right;">
47.53728
</td>
<td style="text-align:right;">
0.0187927
</td>
<td style="text-align:right;">
48.52982
</td>
</tr>
<tr>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
46
</td>
<td style="text-align:right;">
2.039733
</td>
<td style="text-align:right;">
2.070647
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
Function criterion near zero
</td>
<td style="text-align:right;">
0.00654
</td>
<td style="text-align:right;">
0.0426
</td>
<td style="text-align:left;">
normal
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
46.70045
</td>
<td style="text-align:right;">
46.87892
</td>
<td style="text-align:right;">
46.78969
</td>
<td style="text-align:right;">
0.0187547
</td>
<td style="text-align:right;">
48.22016
</td>
</tr>
<tr>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
45
</td>
<td style="text-align:right;">
2.034888
</td>
<td style="text-align:right;">
2.076213
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
Function criterion near zero
</td>
<td style="text-align:right;">
0.00654
</td>
<td style="text-align:right;">
0.0426
</td>
<td style="text-align:left;">
normal
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
45.93857
</td>
<td style="text-align:right;">
46.16215
</td>
<td style="text-align:right;">
46.05036
</td>
<td style="text-align:right;">
0.0187151
</td>
<td style="text-align:right;">
47.91550
</td>
</tr>
<tr>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
44
</td>
<td style="text-align:right;">
2.029850
</td>
<td style="text-align:right;">
2.082044
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
Function criterion near zero
</td>
<td style="text-align:right;">
0.00654
</td>
<td style="text-align:right;">
0.0426
</td>
<td style="text-align:left;">
normal
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
45.18813
</td>
<td style="text-align:right;">
45.45141
</td>
<td style="text-align:right;">
45.31977
</td>
<td style="text-align:right;">
0.0186738
</td>
<td style="text-align:right;">
47.61636
</td>
</tr>
</tbody>
</table>

The *ASN*<sub>*avg*</sub>-optimal design can be obtained as

``` r
optimal.design<-designs %>% group_by(q,dist,p_alpha,p_beta) %>%
  filter( (abs(alpha-risks[1])<1e-05) & (abs(risks[2]-beta)<1e-05) & (termcd==1)) %>%
  slice(which.min(asn_avg)) %>% arrange(q,p_alpha,p_beta) %>% as.data.frame()
```

<table class="table table-striped table-hover table-condensed" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:right;">
q
</th>
<th style="text-align:right;">
n
</th>
<th style="text-align:right;">
kr
</th>
<th style="text-align:right;">
ka
</th>
<th style="text-align:right;">
termcd
</th>
<th style="text-align:left;">
message
</th>
<th style="text-align:right;">
p_alpha
</th>
<th style="text-align:right;">
p_beta
</th>
<th style="text-align:left;">
dist
</th>
<th style="text-align:right;">
alpha
</th>
<th style="text-align:right;">
beta
</th>
<th style="text-align:right;">
asn_alpha
</th>
<th style="text-align:right;">
asn_beta
</th>
<th style="text-align:right;">
asn_avg
</th>
<th style="text-align:right;">
p_asn_max
</th>
<th style="text-align:right;">
asn_max
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
21
</td>
<td style="text-align:right;">
1.799825
</td>
<td style="text-align:right;">
2.3974
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
Function criterion near zero
</td>
<td style="text-align:right;">
0.00654
</td>
<td style="text-align:right;">
0.0426
</td>
<td style="text-align:left;">
normal
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.1000001
</td>
<td style="text-align:right;">
34.71139
</td>
<td style="text-align:right;">
32.27409
</td>
<td style="text-align:right;">
33.49274
</td>
<td style="text-align:right;">
0.0168694
</td>
<td style="text-align:right;">
46.02262
</td>
</tr>
</tbody>
</table>

### Design of repetitive group sampling plans using a generalized beta (GB) prior model of *p* and expected sampling risks

In this case, the censored repetitive sampling plans can be determined
when a *GB* prior is assumed and there is a knowledge about the mean
and variance of *p*. For given requirements of maximum expected risks
and quality levels, the sampling plans are

``` r
risks<-c(0.05,0.10)     #vector of producer and consumer maximum sampling risks
p<-c(0.00654, 0.0426)   #vector of acceptance and rejection quality levels
q<- 0.1                 #censoring degree
asvar<-asympt.var(q,"normal")    #asymptotical variance-covariance matrix of MLE estimators of location and scale paramters
l<- p[1]/5              #lower limit of p
u<- p[2]+(p[1]-l)       #upper limit of p

# GB parameters for a knowledge of mean and variance of p distribution
know_p<-list(mean_p=p[1],var_p=((p[2]-p[1])/4)^2)
beta.parms<-beta.params(p,l,u, know_p)

designs<-repGBprior.plan(risks,p,asvar, beta.parms)
```

Then, the function *repGBprior.plan()* returns these designs. The first
plans are
<table class="table table-striped table-hover table-condensed" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:right;">
q
</th>
<th style="text-align:right;">
n
</th>
<th style="text-align:right;">
n_low
</th>
<th style="text-align:right;">
n_up
</th>
<th style="text-align:right;">
kr
</th>
<th style="text-align:right;">
ka
</th>
<th style="text-align:right;">
termcd
</th>
<th style="text-align:left;">
message
</th>
<th style="text-align:right;">
p_alpha
</th>
<th style="text-align:right;">
p_beta
</th>
<th style="text-align:right;">
a
</th>
<th style="text-align:right;">
b
</th>
<th style="text-align:right;">
l
</th>
<th style="text-align:right;">
u
</th>
<th style="text-align:right;">
mean_p
</th>
<th style="text-align:right;">
var_p
</th>
<th style="text-align:left;">
dist
</th>
<th style="text-align:right;">
alpha
</th>
<th style="text-align:right;">
beta
</th>
<th style="text-align:right;">
asn_alpha
</th>
<th style="text-align:right;">
asn_beta
</th>
<th style="text-align:right;">
asn_avg
</th>
<th style="text-align:right;">
easn
</th>
<th style="text-align:right;">
p_asn_max
</th>
<th style="text-align:right;">
asn_max
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
24
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
49
</td>
<td style="text-align:right;">
2.200833
</td>
<td style="text-align:right;">
2.200877
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
Function criterion near zero
</td>
<td style="text-align:right;">
0.00654
</td>
<td style="text-align:right;">
0.0426
</td>
<td style="text-align:right;">
0.1862234
</td>
<td style="text-align:right;">
1.469713
</td>
<td style="text-align:right;">
0.001308
</td>
<td style="text-align:right;">
0.047832
</td>
<td style="text-align:right;">
0.00654
</td>
<td style="text-align:right;">
8.13e-05
</td>
<td style="text-align:left;">
normal
</td>
<td style="text-align:right;">
0.0500619
</td>
<td style="text-align:right;">
0.1000161
</td>
<td style="text-align:right;">
24.00090
</td>
<td style="text-align:right;">
24.00042
</td>
<td style="text-align:right;">
24.00066
</td>
<td style="text-align:right;">
24.00039
</td>
<td style="text-align:right;">
0.0121556
</td>
<td style="text-align:right;">
24.00107
</td>
</tr>
<tr>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
23
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
49
</td>
<td style="text-align:right;">
2.187946
</td>
<td style="text-align:right;">
2.216942
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
Function criterion near zero
</td>
<td style="text-align:right;">
0.00654
</td>
<td style="text-align:right;">
0.0426
</td>
<td style="text-align:right;">
0.1862234
</td>
<td style="text-align:right;">
1.469713
</td>
<td style="text-align:right;">
0.001308
</td>
<td style="text-align:right;">
0.047832
</td>
<td style="text-align:right;">
0.00654
</td>
<td style="text-align:right;">
8.13e-05
</td>
<td style="text-align:left;">
normal
</td>
<td style="text-align:right;">
0.0500908
</td>
<td style="text-align:right;">
0.1000238
</td>
<td style="text-align:right;">
23.58511
</td>
<td style="text-align:right;">
23.27324
</td>
<td style="text-align:right;">
23.42917
</td>
<td style="text-align:right;">
23.30382
</td>
<td style="text-align:right;">
0.0120410
</td>
<td style="text-align:right;">
23.69109
</td>
</tr>
<tr>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
22
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
49
</td>
<td style="text-align:right;">
2.173528
</td>
<td style="text-align:right;">
2.234751
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
Function criterion near zero
</td>
<td style="text-align:right;">
0.00654
</td>
<td style="text-align:right;">
0.0426
</td>
<td style="text-align:right;">
0.1862234
</td>
<td style="text-align:right;">
1.469713
</td>
<td style="text-align:right;">
0.001308
</td>
<td style="text-align:right;">
0.047832
</td>
<td style="text-align:right;">
0.00654
</td>
<td style="text-align:right;">
8.13e-05
</td>
<td style="text-align:left;">
normal
</td>
<td style="text-align:right;">
0.0500003
</td>
<td style="text-align:right;">
0.1000001
</td>
<td style="text-align:right;">
23.20269
</td>
<td style="text-align:right;">
22.56201
</td>
<td style="text-align:right;">
22.88235
</td>
<td style="text-align:right;">
22.63308
</td>
<td style="text-align:right;">
0.0119251
</td>
<td style="text-align:right;">
23.40886
</td>
</tr>
<tr>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
21
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
49
</td>
<td style="text-align:right;">
2.158403
</td>
<td style="text-align:right;">
2.254219
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
Function criterion near zero
</td>
<td style="text-align:right;">
0.00654
</td>
<td style="text-align:right;">
0.0426
</td>
<td style="text-align:right;">
0.1862234
</td>
<td style="text-align:right;">
1.469713
</td>
<td style="text-align:right;">
0.001308
</td>
<td style="text-align:right;">
0.047832
</td>
<td style="text-align:right;">
0.00654
</td>
<td style="text-align:right;">
8.13e-05
</td>
<td style="text-align:left;">
normal
</td>
<td style="text-align:right;">
0.0500319
</td>
<td style="text-align:right;">
0.1000146
</td>
<td style="text-align:right;">
22.82986
</td>
<td style="text-align:right;">
21.85480
</td>
<td style="text-align:right;">
22.34233
</td>
<td style="text-align:right;">
21.97611
</td>
<td style="text-align:right;">
0.0118011
</td>
<td style="text-align:right;">
23.12469
</td>
</tr>
<tr>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
20
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
49
</td>
<td style="text-align:right;">
2.142068
</td>
<td style="text-align:right;">
2.275593
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
Function criterion near zero
</td>
<td style="text-align:right;">
0.00654
</td>
<td style="text-align:right;">
0.0426
</td>
<td style="text-align:right;">
0.1862234
</td>
<td style="text-align:right;">
1.469713
</td>
<td style="text-align:right;">
0.001308
</td>
<td style="text-align:right;">
0.047832
</td>
<td style="text-align:right;">
0.00654
</td>
<td style="text-align:right;">
8.13e-05
</td>
<td style="text-align:left;">
normal
</td>
<td style="text-align:right;">
0.0500924
</td>
<td style="text-align:right;">
0.1000943
</td>
<td style="text-align:right;">
22.47555
</td>
<td style="text-align:right;">
21.15569
</td>
<td style="text-align:right;">
21.81562
</td>
<td style="text-align:right;">
21.34039
</td>
<td style="text-align:right;">
0.0116246
</td>
<td style="text-align:right;">
22.84847
</td>
</tr>
<tr>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
49
</td>
<td style="text-align:right;">
2.123503
</td>
<td style="text-align:right;">
2.300221
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
Function criterion near zero
</td>
<td style="text-align:right;">
0.00654
</td>
<td style="text-align:right;">
0.0426
</td>
<td style="text-align:right;">
0.1862234
</td>
<td style="text-align:right;">
1.469713
</td>
<td style="text-align:right;">
0.001308
</td>
<td style="text-align:right;">
0.047832
</td>
<td style="text-align:right;">
0.00654
</td>
<td style="text-align:right;">
8.13e-05
</td>
<td style="text-align:left;">
normal
</td>
<td style="text-align:right;">
0.0500073
</td>
<td style="text-align:right;">
0.1000022
</td>
<td style="text-align:right;">
22.18107
</td>
<td style="text-align:right;">
20.48217
</td>
<td style="text-align:right;">
21.33162
</td>
<td style="text-align:right;">
20.74946
</td>
<td style="text-align:right;">
0.0114847
</td>
<td style="text-align:right;">
22.62697
</td>
</tr>
</tbody>
</table>

and the *EASN*-optimal design is obtained as

``` r
optimal.design<-designs %>% group_by(q,dist,p_alpha,p_beta) %>%
                 filter( (abs(alpha-risks[1])<1e-05) & 
                           (abs(risks[2]-beta)<1e-05) & (termcd==1)) %>%
                 group_by(q,p_alpha,p_beta,a,b,l,u,dist) %>%
                 mutate(easn_min=min(easn)) %>%
                 slice(which.min(easn)) %>% as.data.frame()
```

<table class="table table-striped table-hover table-condensed" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:right;">
q
</th>
<th style="text-align:right;">
n
</th>
<th style="text-align:right;">
n_low
</th>
<th style="text-align:right;">
n_up
</th>
<th style="text-align:right;">
kr
</th>
<th style="text-align:right;">
ka
</th>
<th style="text-align:right;">
termcd
</th>
<th style="text-align:left;">
message
</th>
<th style="text-align:right;">
p_alpha
</th>
<th style="text-align:right;">
p_beta
</th>
<th style="text-align:right;">
a
</th>
<th style="text-align:right;">
b
</th>
<th style="text-align:right;">
l
</th>
<th style="text-align:right;">
u
</th>
<th style="text-align:right;">
mean_p
</th>
<th style="text-align:right;">
var_p
</th>
<th style="text-align:left;">
dist
</th>
<th style="text-align:right;">
alpha
</th>
<th style="text-align:right;">
beta
</th>
<th style="text-align:right;">
asn_alpha
</th>
<th style="text-align:right;">
asn_beta
</th>
<th style="text-align:right;">
asn_avg
</th>
<th style="text-align:right;">
easn
</th>
<th style="text-align:right;">
p_asn_max
</th>
<th style="text-align:right;">
asn_max
</th>
<th style="text-align:right;">
easn_min
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
49
</td>
<td style="text-align:right;">
1.921685
</td>
<td style="text-align:right;">
2.611877
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
Function criterion near zero
</td>
<td style="text-align:right;">
0.00654
</td>
<td style="text-align:right;">
0.0426
</td>
<td style="text-align:right;">
0.1862234
</td>
<td style="text-align:right;">
1.469713
</td>
<td style="text-align:right;">
0.001308
</td>
<td style="text-align:right;">
0.047832
</td>
<td style="text-align:right;">
0.00654
</td>
<td style="text-align:right;">
8.13e-05
</td>
<td style="text-align:left;">
normal
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.1000001
</td>
<td style="text-align:right;">
21.7853
</td>
<td style="text-align:right;">
16.33357
</td>
<td style="text-align:right;">
19.05943
</td>
<td style="text-align:right;">
18.15345
</td>
<td style="text-align:right;">
0.0099061
</td>
<td style="text-align:right;">
22.38735
</td>
<td style="text-align:right;">
18.15345
</td>
</tr>
</tbody>
</table>
