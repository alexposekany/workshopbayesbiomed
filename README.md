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
