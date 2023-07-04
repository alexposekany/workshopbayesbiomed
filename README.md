# workshopbayesbiomed

You will find the R Markdown File for the Workshop under Presentation_Bayes_Biomed.Rmd 

To install relevant packages use 

```
list.of.packages <- c('ggplot2','dplyr',"tidyr",'lattice',"Pareto","HDInterval","rstanarm","rstan","rjags","bayess","MCMCpack","MASS","data.table","bayesreg","boot","MCMCvis","Rgraphviz","graph","invgamma")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(Pareto)
library(HDInterval)
library(rstanarm)
library(rjags)
library(bayess)
library(MCMCpack)
library(INLA)
library(data.table)
library(MASS)
library(bayesreg)
library(boot)
library(MCMCvis)
library(dplyr)
library(invgamma)
```

## Covid Example

We recalculate the estimation of the prevalence of Covid19 in spring 2020. Samples from 1279 persons were analysed with PCR testing procedures. Out of all those not a single randomly selected person was tested positively. This obviously breaks standard testing mechanisms for extimating the proportion of infected person in Austria. 

However, additional information is available from similar tests in Germany which had a comparable behaviour of the spread of the disease at that time. In the same time span 4 positive cases out of 4068 had been found.

###  Prior Choices
```
alpha=1.4
beta=407.5
xrate=seq(0,0.01,length.out=300)
prior=dbeta(xrate, shape1=alpha, shape2=beta)
plot(xrate, prior, 
     xlab="proportion of COVID positive people",
     ylab="density",
     main="Beta prior based on German study",
     type="l",
     col="red",
     lwd=2)
```

### Posterior Choices

```
alphaPost=alpha
betaPost=beta+1279
posterior=dbeta(xrate, alphaPost, betaPost)
plot(xrate, posterior, 
     xlab="proportion of COVID positive people",
     ylab="density",
     main="Beta posterior and prior based on German study",
     type="l",
     col="blue",
     lwd=2)

lines(xrate, prior, col="red")
legend("topright", c("posterior", "prior"), col=c("blue","red"),lwd=4)
```
### Plotting Posteriors

```

posteriorMean=alphaPost/(alphaPost+betaPost)
posteriorMode=(alphaPost-1)/(alphaPost+betaPost-2)
posteriorMedian=(alphaPost-1/3)/(alphaPost+betaPost-2/3)
posteriorhdi=hdi(qbeta, 0.95, shape1=alphaPost, shape2=betaPost) 

HDIlower=as.double(posteriorhdi["lower"])
HDIupper=as.double(posteriorhdi["upper"])

knitr::kable(cbind(c("Mean","Mode","Median","lower 95%-HPD boundary",
                     "upper 95%-HPD boundary"),
                   c(posteriorMean, 
                     posteriorMode, 
                     posteriorMedian,
                     HDIlower,
                     HDIupper)),
             col.names=c("Estimator","Posterior Point Estimate"),
             caption="Point esimators for the posterior distribution based on
             the prior that uses German COVID-19 cases")

posterior=dbeta(xrate, alphaPost, betaPost)
plot(xrate, posterior, 
     xlab="proportion of COVID positive people",
     ylab="density",
     main="Beta posterior and HDI based on German study",
     type="l",
     col="blue",
     lwd=2)

plotSeq=seq(HDIlower,HDIupper,length.out=100)
polygon(x=c(HDIlower,
            plotSeq,
            HDIupper),
        y=c(0,
            dbeta(plotSeq,
                  shape1=alphaPost,
                  shape2=betaPost),
            0),col="#FF990022")

legend("topright", c("posterior", "HPD interval 95%"), col=c("blue","#FF990022"),lwd=4)
```

## Bayesian Regression - Practical How to

```
library(bayesreg)
bayesregDNase<-bayesreg(density~I(conc^(1/2)),data=DNase)
summary(bayesregDNase)

lmDNase<-lm(density~I(conc^(1/2)),data=DNase)
summary(lmDNase)

newconc<-seq(min(DNase$conc),max(DNase$conc),length.out=200)
plot((DNase$conc),(DNase$density))
lines(newconc,predict(lmDNase,newdata = data.frame(conc=newconc)),col="red",lwd=6)
lines(newconc,predict(bayesregDNase,newdata = data.frame(conc=newconc)),col="darkblue",lwd=3)
legend("bottomright",legend = c("Bayes","Classical"),col=c("darkblue","red"),lwd = 5)
```
