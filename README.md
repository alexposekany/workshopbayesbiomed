# Workshop on Bayesian Methods in Biology and Medicine

You will find the R Markdown File for the Workshop under Presentation_Bayes_Biomed.Rmd 
The Video Recording of this course is available on [YouTube](https://youtu.be/h_pEORzO8VE) as part of a lecture series on [Statistical Computing and computer intensive methods](https://www.youtube.com/playlist?list=PLM2jwuiJI2BNQvbNBCTu5xdsJc978vTwX).

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

 Compare the Bayesian against the frequentist results for dataset DNase and model lm(density~I(conc)^(1/2),data=DNase)

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

# Algortihms 

We first adapt the example from [Montana University](https://math.montana.edu/ahoegh/teaching/stat532/labs/STANDemo.html) for our body size example.

```
## Simulate Data
set.seed(45725)
N = 5000
mu = 175
sigma = 14
sigmasq = sigma^2
y = rnorm(N, mean = mu, sd = sigma)

## set up priors
mu.0 <- 0
tausq.0 <- 100
nu.0 <- .5
sigmasq.0 <- 1
rate.param <- nu.0 * sigmasq.0 / 2
shape.param <- nu.0 / 2
```

## Example - Body Sizes - Gibbs Sampler

```
### initialize vectors and set starting values 
num.sims <- 10000
mu.samples <- rep(0, num.sims)
sigmasq.samples <- rep(1, num.sims)

mean.y <- mean(y)
nu.n <- nu.0 + N

for (iter in 2:num.sims){
  # sample theta from full conditional
  mu.n <- (mu.0 / tausq.0 + N * mean.y / sigmasq.samples[iter - 1]) / (1 / tausq.0 + N / sigmasq.samples[iter - 1] )
  tausq.n <- 1 / (1/tausq.0 + N / sigmasq.samples[iter - 1])
  mu.samples[iter - 1] <- rnorm(1,mu.n,sqrt(tausq.n))
  
  # sample (1/sigma.sq) from full conditional
  sigmasq.n.theta <- 1/nu.n*(nu.0*sigmasq.0 + sum((y - mu.samples[iter])^2))
  sigmasq.samples[iter] <- rinvgamma(1,shape = nu.n/2, rate = nu.n*sigmasq.n.theta/2)
}
# compute posterior mean and quantiles
mean(mu.samples)
quantile(mu.samples, probs = c(.025, .975))

mean(sigmasq.samples)
quantile((sigmasq.samples)^(1/2), probs = c(.025, .975))

par(mfrow=c(1,2))
# plot marginal posterior of mu
hist(mu.samples,xlab=expression(mu),main=expression('Marginal Posterior of ' ~ mu),probability=T,xlim = c(160,190),breaks = seq(0,200,by=5),ylim=c(0,0.16))
lines(density(mu.samples))
abline(v=mu,col='red',lwd=2)
# plot marginal posterior of sigmasq
hist((sigmasq.samples)^(1/2),xlab=expression(sigma[2]),main=expression('Marginal Posterior of ' ~ sigma[2]),probability=T,xlim = c(160,190),breaks = seq(0,200,by=5),ylim=c(0,0.16))
lines(density((sigmasq.samples)^(1/2)))
abline(v=13.35^2,col='red',lwd=2)

par(mfrow=c(1,2))
# plot trace plots
plot(mu.samples,type='l',ylab=expression(mu), main=expression('Trace plot for ' ~ mu))
abline(h=mu,lwd=2,col='red')
plot(sigmasq.samples^(1/2),type='l',ylab=expression(sigma[2]), main=expression('Trace plot for ' ~ sigma[2]))
abline(h=13.35^2,lwd=2,col='red')
```
Rerun the example with a more informative prior to remove the bias of the mean

```
## set up priors
mu.0 <- 170
tausq.0 <- 100
nu.0 <- .5
sigmasq.0 <- 1
rate.param <- nu.0 * sigmasq.0 / 2
shape.param <- nu.0 / 2

library(invgamma)
### initialize vectors and set starting values 
num.sims <- 10000
mu.samples <- rep(0, num.sims)
sigmasq.samples <- rep(1, num.sims)

mean.y <- mean(y)
nu.n <- nu.0 + N

for (iter in 2:num.sims){
  # sample theta from full conditional
  mu.n <- (mu.0 / tausq.0 + N * mean.y / sigmasq.samples[iter - 1]) / (1 / tausq.0 + N / sigmasq.samples[iter - 1] )
  tausq.n <- 1 / (1/tausq.0 + N / sigmasq.samples[iter - 1])
  mu.samples[iter - 1] <- rnorm(1,mu.n,sqrt(tausq.n))
  
  # sample (1/sigma.sq) from full conditional
  sigmasq.n.theta <- 1/nu.n*(nu.0*sigmasq.0 + sum((y - mu.samples[iter])^2))
  sigmasq.samples[iter] <- rinvgamma(1,shape = nu.n/2, rate = nu.n*sigmasq.n.theta/2)
}
# compute posterior mean and quantiles
mean(mu.samples)
quantile(mu.samples, probs = c(.025, .975))

mean(sigmasq.samples)
quantile((sigmasq.samples)^(1/2), probs = c(.025, .975))

par(mfrow=c(1,2))
# plot marginal posterior of mu
hist(mu.samples,xlab=expression(mu),main=expression('Marginal Posterior of ' ~ mu),probability=T,xlim = c(160,190),breaks = seq(0,200,by=5),ylim=c(0,0.16))
lines(density(mu.samples))
abline(v=mu,col='red',lwd=2)
# plot marginal posterior of sigmasq
hist((sigmasq.samples)^(1/2),xlab=expression(sigma[2]),main=expression('Marginal Posterior of ' ~ sigma[2]),probability=T,xlim = c(160,190),breaks = seq(0,200,by=5),ylim=c(0,0.16))
lines(density((sigmasq.samples)^(1/2)))
abline(v=13.35^2,col='red',lwd=2)

par(mfrow=c(1,2))
# plot trace plots
plot(mu.samples,type='l',ylab=expression(mu), main=expression('Trace plot for ' ~ mu))
abline(h=mu,lwd=2,col='red')
plot(sigmasq.samples^(1/2),type='l',ylab=expression(sigma[2]), main=expression('Trace plot for ' ~ sigma[2]))
abline(h=13.35^2,lwd=2,col='red')

```

## Metropolis-Hastings Sampler

```
# Run Metropolis Algorithm
num.mcmc <- 15000
stepsize.mu <- 2
stepsize.sigmasq <- .75
acceptratio.mu <- acceptratio.sigmasq <- rep(0,num.mcmc)
mu.samplesmh <- rep(170,num.mcmc)
sigmasq.samplesmh <- sigmasq.star <- rep(3, num.mcmc)

for (iter in 2:num.mcmc){
  # mu
  mu.star <- mu.samplesmh[iter] + rnorm(1,0,stepsize.mu)
  log.p.current <- sum(dnorm(y, mean = mu.samplesmh[iter - 1],
                             sd = sqrt(sigmasq.samplesmh[iter - 1]), log=T)) +
  dnorm(mu.samplesmh[iter - 1],mean = mu.0, sd = sqrt(tausq.0), log=T)
  log.p.star <- sum(dnorm(y, mean = mu.star, sd = sqrt(sigmasq.samplesmh[iter - 1]), log=T)) +
    dnorm(mu.star, mean = mu.0, sd = sqrt(tausq.0), log=T)

  log.r <- log.p.star - log.p.current

  if (log(runif(1)) < log.r){
    mu.samplesmh[iter] <- mu.star
    acceptratio.mu[iter] <- 1
  } else{
    mu.samplesmh[iter] <- mu.samplesmh[iter - 1]
  }
  # sigma
  sigmasq.star[iter] <- sigmasq.samplesmh[iter-1] + rnorm(1,0,stepsize.sigmasq)
  log.p.current <- sum(dnorm(y, mean = mu.samplesmh[iter], 
                             sd = sqrt(sigmasq.samplesmh[iter - 1]), log=T)) + 
     dinvgamma(sigmasq.samplesmh[iter - 1], rate = rate.param, shape = shape.param, log=T)
  log.p.star <- sum(dnorm(y, mean = mu.samplesmh[iter], sd = sqrt(sigmasq.star[iter]), log=T)) +
    dinvgamma(sigmasq.star[iter], rate = rate.param, shape = shape.param, log=T)

  log.r <- log.p.star - log.p.current
  if (log(runif(1)) < log.r){
    sigmasq.samplesmh[iter] <- sigmasq.star[iter]
    acceptratio.sigmasq[iter] <- 1
  } else{
    sigmasq.samplesmh[iter] <- sigmasq.samplesmh[iter - 1]
  }
}

paste("Mean Acceptance rate of mu:",mean(acceptratio.mu))
paste("Mean Acceptance rate of sigmasq:",mean(acceptratio.sigmasq))

paste("Posterior Mean of mu:",round(mean(mu.samplesmh),digits=2))
paste("Posterior 95% HPDI of mu:",paste(round(quantile(mu.samplesmh, probs = c(.025, .975)),digits=2),collapse = ", "))

paste("Posterior Mean of sigmasq:",round(mean((sigmasq.samplesmh)^(1/2)),digits=2))
paste("Posterior 95% HPDI of sigmasq:",paste(round(quantile(sigmasq.samplesmh^(1/2), probs = c(.025, .975)),digits=2),collapse = ", "))

par(mfrow=c(1,2))
plot(mu.samplesmh,type='l')
abline(h=mu,lwd=2,col='red')

plot((sigmasq.samplesmh)^(1/2),type='l')
abline(h=sigmasq,lwd=2,col='red')
```

## JAGS

```

library(rjags)
# Define the model:
modelString = "
model {
    for (i in 1:N) {
        y[i] ~ dnorm(mu, tau.sq)
    }
    mu ~ dnorm(0, 1/ 100)
    tau.sq  ~ dgamma(.005, .005)
  sigmasq <- 1 / tau.sq
}
" 
writeLines( modelString , con="TEMPmodel.txt" )

jags <- jags.model('TEMPmodel.txt',
                   data = list('y' = y,
                               'N' = N),
                   n.chains = 4,
                   n.adapt = 100)
 


codaSamples = coda.samples( jags , variable.names=c("mu","tau.sq", 'sigmasq') ,
                            n.iter=1000)

summary(codaSamples)
plot(codaSamples)
```

## RStan

```
library(rstan)
set.seed(11122017)
N = 5000
mu = 170
sigmasq = 10
y = rnorm(N, mu, sqrt(sigmasq))
num.mcmc <- 1000

stancode <- 'data {
  int<lower=0> N; 
  vector[N] y;
}
parameters {
  real<lower=0> sigmasq;
  real mu; 
} 
model {
  mu ~ normal(0, sqrt(100));
  sigmasq ~ inv_gamma(.005, .005);
  y ~ normal(mu, sqrt(sigmasq));
}' 

stan_data = list(N=N, y=y) # data passed to stan 
# set up the model
stan_model = stan(model_code = stancode, data = stan_data, chains = 1)
stanfit = stan(fit = stan_model, data = stan_data,
               iter=num.mcmc) # run the model print(stanfit,digits=2)
print(stanfit)
```

# Practical Example for Logistic Regression with Pima Indian Data 

```
library(MASS)
summary(Pima.tr2)
## put histograms on the diagonal
panel.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}
## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y,use="pairwise.complete.obs",method="spearman"))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r,col=ifelse(r>=0.6,"red",ifelse(r>=0.4,"orange","black")))
}
pairs(MASS::Pima.tr2, lower.panel = panel.smooth, upper.panel = panel.cor,
diag.panel = panel.hist, las=1)
```
## Reference Analysis (Classical Logistic Regression)

```
summary(logitPima<-glm(formula = type ~ . - npreg - skin, family = "binomial",
data = Pima.tr2))
PredictTrain <- round(predict(logitPima,newdata = Pima.tr2, type = "response"),digits = 0)
table(PredictTrain,Pima.tr2$type)
```

You can obtain a more detailed [Base Line Reference Logistic Regression of Pima Indian Diabetes Data (Kaggle)](https://www.kaggle.com/code/ksp585/pima-indian-diabetes-logistic-regression-with-r) on Kaggle.

## RStan 

```
library(rstanarm)
fit <- stan_glm(type ~  glu + bp  + bmi + ped + age,
                data = MASS::Pima.te,
                family = binomial(),
                prior_intercept = normal(0, 10),
                prior = normal(0, 2.5),
                prior_aux = cauchy(0, 2.5),
                chains = 4,
                iter = 2000,
                seed = 12345)
print(fit)
fit$stanfit
plot(fit)
```
You can obtain a more detailed [Base Line Reference Logistic Regression of Pima Indian Diabetes Data with rstanarm (Kaggle)](https://www.kaggle.com/code/avehtari/bayesian-logistic-regression-with-rstanarm) on Kaggle. 

## INLA

```
library(INLA)
library(data.table)

library(MASS)
# Load the Pima Indian diabetes dataset from the mlbench package
data(Pima.tr2, package = "MASS")

# Split the data into training and testing sets using a 70/30 split
set.seed(1234)
train <- sample(nrow(Pima.tr2), 0.7 * nrow(Pima.tr2))
test <- setdiff(1:nrow(Pima.tr2), train)

y <- as.matrix(as.numeric(Pima.tr2[train, 8])-1)
x <- as.matrix(Pima.tr2[train, -8])

formula <- y ~ x

model <- inla(formula, family = "binomial", data = list(y = y, x = x), control.compute = list(dic = TRUE))

summary(model)
```




