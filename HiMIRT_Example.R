# Regularized Multidimensional Item Response Theory Modeling 
# Examples

# load all necessary packages
library(MASS) # mvrnorm() used
library(progress) # used for progress bar
library(clue) # solve_LSAP() from it is used for finding A_best
library(GPArotation) # GPFoblq() used for geomin rotation

# source the r script for functions
source("./HiMIRT_Functions.R") # remember to modify the path if needed

## Notation
### A: loading matrix
### d: difficulty vector
### Theta: latent score matrix



#######################################
#####-----CJMLE with rotation-----#####
#######################################


J=300;tau=10;N=J*tau;K=10
# J: Number of Items;
# N: sample size, N=tau*J
# K: Number of factors

A_true <- generate_A(J,K,zero_prob=.7,diag=F)
d <- runif(J,-2,2)
Data <- simuMIRTData(J,K,N,rho=0,A_true,d)
# rho: correlation between factors (0 for orthogonal)
# zero_prob: controls the proportion of zeros in the loading matrix
Y <- Data$Y

## run CJMLE with rotation
result <- HiMMIRT(Y,K,max_iter=200,tol=1e-4,rotation=T,orthogonal=TRUE)
#View(result$A)

## hard thresholding
result$A[abs(result$A)<.4] <- 0

A_best <- ref.perp(result$A,A_true) ## account for sign and column indeterminacy

# check sensitivity
sens <- sum(ifelse(A_best!=0&A_true!=0, 1, 0))/sum(ifelse(A_true!=0,1,0))
# check specificity
spec <- sum(ifelse(A_best==0&A_true==0, 1, 0))/sum(ifelse(A_true==0,1,0))

View(A_best)
View(A_true)

### Note: rotation might sometimes not converge, especially when the dataset is
### large. Hence, you can change the setting to the most complicated and repeat 
### the algorithm until the warning message that indicates non-convergence 
### occurs, then calculate sens and spec, which is expected to be not optimal.



################################
#####-----lasso_HiMIRT-----#####
################################


# Non-orthogonal case for example
J=300;tau=10;N=J*tau;K=10

A_true <- generate_A(J,K,zero_prob=.7,diag=F)
d <- runif(J,-2,2)
Data <- simuMIRTData(J,K,N,rho=0.3,zero_prob=0.7)
Y <- Data$Y

## Initialize by less accurate CJMLE with rotation
Params <- HiMMIRT(Y,K,max_iter=200,tol=1e-3,rotation=T,orthogonal=F)

system.time(
  result_list <- autolasso_HiMIRT(Y,K,Params,lambda_range=c(0,300),
                                  orthogonal=F,progress_shown = T)
  # remember to change the range of lambda if the dataset size is changed
)

showResults(result_list)

### The model chosen by BIC is presented in red points. It prefers more
### complex models a bit and chooses the sub-optimal.

## pick the result by the previous step of inspection
picked_result <- result_list[[18]] # the best number may change

## refitting by the selected design matrix
Lambda <- ifelse(picked_result$A==0,0,1)
best_result <- CHiMMIRT(Y,K,Lambda,Params,max_iter=200,orthogonal=F)

A_best <- ref.perp(best_result$A,A_true)

# check sensitivity
sens <- sum(ifelse(A_best!=0&A_true!=0, 1, 0))/sum(ifelse(A_true!=0,1,0))
# check specificity
spec <- sum(ifelse(A_best==0&A_true==0, 1, 0))/sum(ifelse(A_true==0,1,0))
# check loading error
loading_error <- sum((A_best-A_true)^2)/(J*K)



################################
#####-----card_HiMIRT-----######
################################


J=300;tau=10;N=J*tau;K=10

A_true <- generate_A(J,K,zero_prob=.7,diag=F)
d <- runif(J,-2,2)
Data <- simuMIRTData(J,K,N,rho=0,zero_prob=0.7)
Y <- Data$Y

system.time(
  result <- card_HiMIRT(Y,K,random_start=F,max_iter=200,
                        orthogonal=T,kappaC=sum(Data$A!=0))
)

A_best <- ref.perp(result$A,Data$A)
A_true <- Data$A

sens <- sum(ifelse(A_best!=0&A_true!=0, 1, 0))/sum(ifelse(A_true!=0,1,0))
spec <- sum(ifelse(A_best==0&A_true==0, 1, 0))/sum(ifelse(A_true==0,1,0))
loading_error <- sum((A_best-A_true)^2)/(J*K)

