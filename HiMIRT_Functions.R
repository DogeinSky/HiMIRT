# Regularized Multidimensional Item Response Theory Modeling in a High-dimensional Setting
# All Functions

# Author: Travis Yang
# Last modified: 03-12-2025

#################################################
#####-----Functions for generating data-----#####
#################################################

# generate_q: generate vector q that indicates the sparse pattern of a_j
#          K: the length of vector q=the number of factors
#  zero_prob: the probability of a loading compnent to be zero

generate_q <- function(K,zero_prob){
  first <- rbinom(K-1, 1, 1-zero_prob)
  last <- ifelse(sum(first)==0, 1, rbinom(1,1,1-zero_prob))
  c(first, last)
}


# generate_A: generate loading matrix A
# diag: whether it should be in the diagonal form

generate_A <- function(J,K,zero_prob,diag=F){
  A <- t(vapply(1:J,
                function(j) generate_q(K,zero_prob)*runif(K, 0.5, 2.5),
                numeric(K)))
  
  if (diag) {
    A[1:K,] <- diag(K)*runif(K,0.5,2.5)
  }
  
  return(A)
}

# simuMIRTData: generate data for MIRT
#          rho: the correlation between factors; rho=0 means orthogonal
# ***The output contains***:
#   - The response matrix: Y
#   - The true difficulty vector: d
#   - The true loading matrix: A
#   - The true latent matrix: Theta

simuMIRTData <- function(J,K,N,rho,A,d){
  cov_mat <- matrix(rho, nrow = K, ncol = K)
  diag(cov_mat) <- 1
  
  Theta <- mvrnorm(n = N, mu = rep(0,K), Sigma = cov_mat)
  
  if (rho==0) {
    Theta=scale(svd(Theta)$u)
  } else {
    Theta=scale(Theta)
  }
  
  prob.mat <- plogis(d+A%*%t(Theta))
  Y <- matrix(rbinom(nrow(prob.mat)*ncol(prob.mat), 1, prob.mat),
              nrow=nrow(prob.mat), byrow=FALSE)
  list(Y=t(Y), d=d, A=A, Theta=Theta)
}

# ref.perp: function for reflection and permutation (align A_hat with A_true)
#    A_hat: the estimated loading matrix
#    A_true: the true loading matrix to align

ref.perp <- function(A_hat,A_true){
  K=ncol(A_hat)
  matrix_data=t(vapply(1:K, function(i){
    colSums(ifelse(A_hat[,i]!=0 & A_true!=0, 1,0))},
    numeric(K)))
  assignment <- solve_LSAP(t(matrix_data), maximum = TRUE)
  row_indices <- as.integer(assignment)
  A_best <- A_hat[,row_indices]
  
  corsign <- sign(diag(cor(A_true, A_best)))
  A_best <- A_best %*% diag(corsign)
  return(A_best)
}

###################################################
#####-----Functions for basic computation-----#####
###################################################

# logl: calculate the negative log-likelihood value

logl <- function(X){
  return(-sum(log(X)))
}

# compute_H: calculate the matrix H

compute_H <- function(Q, mat_prod) {
  return(1 / (1 + exp(-Q*mat_prod)) - 1) 
}

# compute_Z: calculate the matrix Z

compute_Z <- function(mat_prod, Q, H){
  return(mat_prod-4*Q*H)
}

# soft_threshold: soft-threshold operator
#         lambda: the penalty parameter for l1 penalty
#            rho: the pre-specified penality parameter for augmented lagrangian term
#              x: the element to impose soft-threshold on

soft_threshold <- function(lambda, rho, x){
  x_abs <- abs(x)
  shrinkage_factor <- pmax(0, 1 - (lambda / (rho * x_abs)))
  return(shrinkage_factor*x)
}

# get_Residuals: calculate residuals for ADMM iterations
#         A_old: matrix A from the previous iteration

getResiduals <- function(A,A_tilde,A_old,U){
  r=sum((A-t(A_tilde))^2)/sum(A^2)
  s=sum((A-A_old)^2)/sum(U^2)
  
  return(c(r,s))
}

# grad_Ad: calculate the gradient with respect to A and d

grad_Ad <- function(Q,X,Theta){
  coef_mat=Q*(1-X)
  return(-t(coef_mat)%*%Theta)
}

# proj_C: projection operator
# Params: The list of parameters containing Theta, A, d
#      C: Normally set to 5\sqrt{K}

proj_C <- function(Params, C, K){
  Theta=Params$Theta
  Ad=cbind(Params$A,Params$d)
  
  vec_1= which((apply(Theta^2,1,sum)+1)>C^2)
  vec_2= which(apply(Ad^2,1,sum)>C^2)
  
  if (length(vec_1)!=0){
    Theta[vec_1,] <- t(vapply(vec_1, function(x) sqrt((C^2-1)/sum(Theta[x,]^2)) * Theta[x,], numeric(K)))
  }
  
  if (length(vec_2)!=0){
    Ad[vec_2,] <- t(vapply(vec_2, function(x) C/sqrt(sum(Ad[x,]^2))*Ad[x,], numeric(K+1)))
  }
  
  Params$Theta <- Theta
  Params$A <- Ad[,1:K]
  Params$d <- Ad[,K+1]
  return(Params)
}

########################################
#####-----Esimating MIRT by MM-----#####
########################################

#     HiMMIRT: Estimate High-dimensional MIRT models by MM algorithm
#    max_iter: A pre-specified maximum number of iteration
#  orthogonal: whether orthogonality is assumed or not
# ***The output contains***:
#   - The estimated difficulty vector: d
#   - The estimated loading matrix: A
#   - The estimated latent matrix: Theta

HiMMIRT <- function(Y,K,max_iter,tol=1e-4,rotation=FALSE,orthogonal=c(TRUE,FALSE),progress_shown=TRUE){
  N=nrow(Y)
  J=ncol(Y)
  
  if (any(rowSums(Y) == 0 | rowSums(Y) == ncol(Y)) ||
      any(colSums(Y) == 0 | colSums(Y) == nrow(Y))){
    projection=TRUE
  } else {
    projection=FALSE
  }
  
  Q=2*Y-1
  
  Params <- paramsINIT(Y, K)
  mat_prod=cbind(Params$Theta,1)%*%t(cbind(Params$A,Params$d))
  X=stats::plogis(Q*(mat_prod))
  initial_logl_value=logl(X)
  print(sprintf("The initial logl value is %.2f", initial_logl_value))
  
  if (progress_shown){
    pb <- progress::progress_bar$new(
      format = "Convergence: [:bar] :progress%  | log-likelihood: :loss | Iter: :iter",
      total = 100,
      clear = FALSE,
      width = 80
    )
  }
  
  logl_values=vector(mode="numeric", length=max_iter)
  
  for (i in 1:max_iter){
    A=Params$A
    d=Params$d
    
    H=compute_H(Q,mat_prod)
    Z=compute_Z(mat_prod,Q,H)
    Z_D=t(t(Z)-Params$d)
    
    Params$Theta=Z_D%*%A%*%solve(t(A)%*%A)
    Theta=cbind(Params$Theta,1)
    
    Ad=t(Z)%*%Theta%*%solve(crossprod(Theta))
    Params$A=Ad[,1:K]
    Params$d=Ad[,K+1]
    
    if (projection) {
      Params=proj_C(Params,C=sqrt(5)*K,K)
    }
    
    mat_prod=cbind(Params$Theta,1)%*%t(cbind(Params$A,Params$d))
    X=stats::plogis(Q*mat_prod)
    logl_values[i]=logl(X)
    
    relative_change=ifelse(i==1,
                           (initial_logl_value-logl_values[i])/initial_logl_value,
                           (logl_values[i-1]-logl_values[i])/logl_values[i-1]
    )
    
    if (progress_shown){
      if (i==1){
        max_diff <- relative_change
      }
      progress_ratio <- min(1, abs(max_diff - relative_change) / (max_diff - tol))
      pb$update(progress_ratio, tokens = list(
        loss = sprintf("%.3f", logl_values[i]),
        iter = i,
        progress=round(progress_ratio*100)
      ))
      Sys.sleep(0.02) # leave some time to show the progress bar
    }
    
    if (relative_change < tol){
      break
    }
    
    if (relative_change < 0){
      warning("increase in the objective")
      break
    }
  }
  cat("\u2714 Algorithm finished!!!\n")
  
  result <- standardiseMIRT(Params,orthogonal=T,admm=F)
  
  if (K!=1 & rotation){
    rotate_list <- GPArotation::GPFoblq(result$A,method="geomin")
    result$A <- rotate_list$loadings # use geomin rotation
    rotation_M <- rotate_list$Th
    result$Theta <- result$Theta %*% rotation_M
  }
  return(result) # standardisation would be different for admm iterations
}

#   CHiMMIRT: Esitmate confirmatory MIRT models by MM algorithms (remove shrinkage)
#     Lambda: the design matrix \Lambda
#     Params: the starting value; use result from JMLE by MM usually helps

CHiMMIRT <- function(Y,K,Lambda,Params,max_iter,tol=1e-4,orthogonal=c(TRUE,FALSE)){
  N=nrow(Y)
  J=ncol(Y)
  Q=2*Y-1
  
  Params$A=ifelse(Lambda==0,0,Params$A)
  mat_prod=cbind(Params$Theta,1)%*%t(cbind(Params$A,Params$d))
  X=stats::plogis(Q*(mat_prod))
  initial_logl_value=logl(X)
  #print(sprintf("The initial logl value is %.2f", initial_logl_value))
  
  logl_values=vector(mode="numeric", length=max_iter)
  
  for (i in 1:max_iter){
    A=Params$A
    d=Params$d
    
    H=compute_H(Q,mat_prod)
    Z=compute_Z(mat_prod,Q,H)
    Z_D=t(t(Z)-Params$d)
    
    if (orthogonal) {
      svd_result <- svd(1/sqrt(N)*Z_D%*%A)
      Params$Theta <- sqrt(N)*svd_result$u%*%t(svd_result$v)
      Theta=cbind(Params$Theta,1)
    } else {
      Params$Theta=Z_D%*%A%*%solve(t(A)%*%A)
      Theta=cbind(Params$Theta,1)
    }
    
    Ad=t(Z)%*%Theta%*%solve(crossprod(Theta))
    Params$A=Ad[,1:K]
    Params$A=Params$A*Lambda # Update as normal, but only keep the non-zeros
    Params$d=Ad[,K+1]
    
    mat_prod=cbind(Params$Theta,1)%*%t(cbind(Params$A,Params$d))
    X=stats::plogis(Q*mat_prod)
    logl_values[i]=logl(X)
    
    relative_change=ifelse(i==1,
                           (initial_logl_value-logl_values[i])/initial_logl_value,
                           (logl_values[i-1]-logl_values[i])/logl_values[i-1]
    )
    
    if (relative_change < tol){
      break
    }
    
    if (relative_change < 0){
      warning("increase in the objective")
      break
    }
  }
  
  if (orthogonal==F){
    Params=standardiseMIRT(Params,orthogonal=F,admm=F)
  }
  
  Params$logl_value=logl_values[i]
  
  return(Params)
}

#   paramsINIT: Initialize a starting set of parameters
#      epsilon: A pre-specified value, set to .01 by default
# Reference:
# Chen, Y., Li, X., & Zhang, S. (2019). Joint maximum likelihood estimation 
# for high-dimensional exploratory item factor analysis. 
# Psychometrika, 84(1), 124-146.
# ***Note***: X (appeared in the referenced text) is replaced as Q in the function

paramsINIT <- function(Y,K,epsilon= .01){
  N=nrow(Y)
  p=1
  Params_start <- list()
  
  Q <- 2*Y-1
  Q.SVD <- svd(Q)
  U.Q <- Q.SVD$u
  V.Q <- Q.SVD$v
  Sigma.Q <- diag(Q.SVD$d)
  J <- sum(Q.SVD$d>=2*sqrt(N*p))
  
  Q_tilde <- U.Q[,1:J]%*%Sigma.Q[1:J,1:J]%*%t(V.Q[,1:J])
  
  M <- suppressWarnings(
    ifelse(Q_tilde < -1+epsilon, qlogis(epsilon),
           ifelse(Q_tilde > 1-epsilon, qlogis(1-epsilon),
                  qlogis(0.5*(Q_tilde+1)))))
  
  Params_start$d <- apply(M, 2, sum)/N
  
  M_tilde <- t(t(M)-Params_start$d)
  
  M.SVD <- svd(M_tilde)
  U.M <- M.SVD$u
  V.M <- M.SVD$v
  Sigma.M <- diag(M.SVD$d)
  
  Params_start$Theta <- U.M[,1:K,drop=F]*sqrt(N)
  Params_start$A <- V.M[,1:K,drop=F]%*%Sigma.M[1:K,1:K,drop=F]/sqrt(N)
  
  return(Params_start)
}

#   standardiseMIRT: standardisation of the parameter set

standardiseMIRT <- function(Params,orthogonal=c(TRUE,FALSE),admm=c(TRUE,FALSE)){
  if (orthogonal) {
    if (admm) {
      A_tilde=Params$A_tilde
      d=Params$d
      Theta=Params$Theta
      
      N=nrow(Theta)
      Theta_center <- t(t(Theta)-colMeans(Theta))
      Params$d <- d+rowMeans(t(A_tilde)%*%t(Theta))
      svd_result <- svd(Theta_center)
      Params$Theta <- sqrt(N)*svd_result$u
      Params$A_tilde <- t(t(A_tilde)%*%svd_result$v%*%diag(svd_result$d)/sqrt(N))
    } else {
      A=Params$A
      d=Params$d
      Theta=Params$Theta
      
      N=nrow(Theta)
      Theta_center <- t(t(Theta)-colMeans(Theta))
      Params$d <- d+rowMeans(A%*%t(Theta))
      svd_result <- svd(Theta_center)
      Params$Theta <- sqrt(N)*svd_result$u
      Params$A <- A%*%svd_result$v%*%diag(svd_result$d)/sqrt(N)
    }
  } else {
    if (admm){
      A_tilde=Params$A_tilde
      d=Params$d
      Theta=Params$Theta
      
      Theta_mean <- colMeans(Theta)
      Theta_sd <- apply(Theta, 2, sd)
      Params$Theta <- scale(Params$Theta)
      Params$d <- drop(d+Theta_mean %*% A_tilde)
      Params$A_tilde <- A_tilde*Theta_sd
    } else{
      A=Params$A
      d=Params$d
      Theta=Params$Theta
      
      Theta_mean <- colMeans(Theta)
      Theta_sd <- apply(Theta, 2, sd)
      Params$Theta <- scale(Params$Theta)
      Params$d <- drop(d+Theta_mean %*% t(A))
      Params$A <- t(t(A)*Theta_sd)
    }}
  
  return(Params)
}

######################################################
#####-----Regularized HiMIRT (lasso version)-----#####
######################################################

#   lasso_HiMIRT: estimate regularized HDMIRT models via l1 penalty

lasso_HiMIRT <- function(Y,K,Params=NULL,max_iter,lambda,tol=1e-4,orthogonal,progress_shown=TRUE){
  N=nrow(Y)
  J=ncol(Y)
  Q=2*Y-1
  
  if (is.null(Params)){
    cat("Starting HiMMIRT!!!\n")
    Params <- HiMMIRT(Y,K,max_iter=200,tol=1e-3,rotation=TRUE,orthogonal=T,progress_shown=progress_shown) # tol=1e-3: less accurate
  }
  Params$A_tilde=t(Params$A)
  
  mat_prod=cbind(Params$Theta,1)%*%t(cbind(t(Params$A_tilde),Params$d))
  initial_loss_value=logl(stats::plogis(Q*(mat_prod)))+lambda*sum(abs(Params$A))
  cat("Starting lasso_HiMIRT\n")
  #print(sprintf("The initial logl value is %.2f", loss_values[2]))
  
  if (progress_shown){
    pb <- progress::progress_bar$new(
      format = "Convergence: [:bar] :progress%  | loss value: :loss | Iter: :iter",
      total = 100,  # We’ll manually control this
      clear = FALSE,
      width = 80
    )
  }
  
  loss_values=vector(mode="numeric", max_iter)
  
  for (i in 1:max_iter){
    H=compute_H(Q,mat_prod)
    Params$d=colMeans(t(t(-4*Q*H)+Params$d))
    mat_prod=cbind(Params$Theta,1)%*%t(cbind(t(Params$A_tilde),Params$d))
    #print(sprintf("After %d iterations (d) the value is %.2f",
    #m-1, logl(plogis(Q*(mat_prod)))+lambda*sum(abs(Params$A))))
    
    H=compute_H(Q,mat_prod)
    Z=compute_Z(mat_prod,Q,H)
    Z_D=t(t(Z)-Params$d)
    
    #print(sprintf("kappa: %.2f", kappa(Params$A_tilde%*%t(Params$A_tilde))))
    if (orthogonal) {
      svd_result <- svd(1/sqrt(N)*Z_D%*%t(Params$A_tilde))
      Params$Theta <- sqrt(N)*svd_result$u%*%t(svd_result$v)
    } else {
      Params$Theta <- Z_D%*%t(Params$A_tilde)%*%solve(Params$A_tilde%*%t(Params$A_tilde))
      
      Params <- standardiseMIRT(Params,orthogonal=F,admm=T)
      Z_D=t(t(Z)-Params$d)
    }
    
    mat_prod=cbind(Params$Theta,1)%*%t(cbind(t(Params$A_tilde),Params$d))
    #print(sprintf("After %d iterations (Theta) the value is %.2f",
    #m-1, logl(plogis(Q*(mat_prod)))+lambda*sum(abs(Params$A))))
    
    Params <- subADMM(Q,Z_D,Params,lambda,20)
    
    mat_prod=cbind(Params$Theta,1)%*%t(cbind(t(Params$A_tilde),Params$d))
    loss_values[i]=logl(stats::plogis(Q*(mat_prod)))+lambda*sum(abs(Params$A))
    #print(sprintf("After %d iterations the value is %.2f",
    #m-2, loss_values[m]))
    
    relative_change <- ifelse(i==1,
                              (initial_loss_value-loss_values[i])/initial_loss_value,
                              (loss_values[i-1]-loss_values[i])/loss_values[i-1]
    )
    
    if (progress_shown){
      if (i==1){
        max_diff <- relative_change
      }
      progress_ratio <- min(1,  abs(max_diff - relative_change) / (max_diff - tol))
      pb$update(progress_ratio, tokens = list(
        loss = sprintf("%.3f", loss_values[i]),
        iter = i,
        progress = round(progress_ratio * 100)
      ))
    }
    
    if (relative_change < tol){
      break
    }
  }
  cat("\u2714 Algorithm finished!!!\n")
  
  Params$logl_value <- logl(stats::plogis(Q*(mat_prod)))
  return(Params)
}

#   subADMM: function for solving the subproblem by ADMM iterations

subADMM <- function(Q,Z_D,Params,lambda,N_iter){
  J=nrow(Params$A)
  K=ncol(Params$A)
  A=Params$A
  A_tilde=Params$A_tilde
  Theta=Params$Theta
  U=matrix(0,nrow=J,ncol=K)
  
  G=crossprod(Theta)/4
  rho=sum(diag(G))/K
  L=t(chol(G+rho*diag(1,nrow(G)),pivot=F))
  M=t(Theta)%*%Z_D/4
  
  #loss_values=vector(mode="numeric", N_iter)
  #mat_prod=cbind(Theta,1)%*%t(cbind(t(A_tilde),Params$d))
  #loss_values[1]=logl(plogis(Q*(mat_prod)))+lambda*sum(abs(Params$A))+rho*sum((A-t(A_tilde)+U)^2)/2
  #print(sprintf("A-update: The initial logl value is %.2f", loss_values[1]))
  
  relative_res=c(1,1)
  
  for (i in 1:N_iter){
    A_tilde <- forwardsolve(L,forwardsolve(L, M+rho*t(A+U)),transpose = T)
    #mat_prod=cbind(Theta,1)%*%t(cbind(t(A_tilde),Params$d))
    #print(sprintf("A-update: After %d iterations (A_tilde) the value is %.2f",
    #m, logl(plogis(Q*(mat_prod)))+lambda*sum(abs(A))+rho*sum((A-t(A_tilde)+U)^2)/2))
    
    A_old=A
    A=soft_threshold(lambda,rho,t(A_tilde)-U)
    
    U=U+A-t(A_tilde)
    
    relative_res=getResiduals(A,A_tilde,A_old,U)
    
    #mat_prod=cbind(Theta,1)%*%t(cbind(t(A_tilde),Params$d))
    #loss_values[m]=logl(plogis(Q*(mat_prod)))+lambda*sum(abs(A))+rho*sum((A-t(A_tilde)+U)^2)/2
    #print(sprintf("A-update: After %d iterations the value is %.2f",
    #m-1, loss_values[m]))
    
    if (all(relative_res<1e-4)){
      break
    }
  }
  
  Params$A=A
  Params$A_tilde=A_tilde
  Params$U=U
  
  return(Params)
}

# autolasso_HiMIRT:

autolasso_HiMIRT <- function(Y,K,Params=NULL,lambda_range,try_out=25,orthogonal,progress_shown=TRUE){
  Q=2*Y-1
  
  lambda_values=seq(lambda_range[2],lambda_range[1],length.out=try_out)
  cat("fitting the most sparse model!\n")
  
  if (is.null(Params)){
    result <- lasso_HiMIRT(Y,K,max_iter=200,lambda=lambda_values[1],
                           orthogonal=orthogonal,progress_shown=progress_shown)
  } else {
    result <- lasso_HiMIRT(Y,K,Params,max_iter=200,lambda=lambda_values[1],
                           orthogonal=orthogonal,progress_shown=progress_shown)
  }
  
  Lambda=ifelse(result$A==0,0,1)
  cat("removing shrinkage!\n")
  result <- CHiMMIRT(Y,K,Lambda,result,max_iter=200,orthogonal=orthogonal)
  
  mat_prod=cbind(result$Theta,1)%*%t(cbind(result$A,result$d))
  X=stats::plogis(Q*mat_prod)
  deriv.A <- grad_Ad(Q, X, result$Theta)
  
  result_set <- vector("list",length(lambda_values))
  names(result_set) <- lambda_values
  result_set[[1]] <- result
  
  if (progress_shown){
    pb <- progress::progress_bar$new(
      format = "Progress: [:bar] :progress% | lambda: :lambda",
      total = 100,
      clear = FALSE,
      width = 80
    )
  }
  
  cat("trying out all lambda values given\n")
  
  for (i in 2:length(lambda_values)){
    lambda=lambda_values[i]
    threshold=2*lambda-lambda_values[i-1]
    discarded_mat=abs(deriv.A)>=threshold
    
    if (sum(discarded_mat)!=0){
      Lambda=ifelse(discarded_mat|result$A!=0,1,0)
      result=CHiMMIRT(Y,K,Lambda,result,max_iter=200,orthogonal=orthogonal)
      
      mat_prod=cbind(result$Theta,1)%*%t(cbind(result$A,result$d))
      X=stats::plogis(Q*mat_prod)
      deriv.A <- grad_Ad(Q, X, result$Theta)
      
      check_mat=ifelse(abs(deriv.A)>lambda & result$A==0, 1, 0)
      
      while (sum(check_mat)>0) {
        Lambda[check_mat==1]=1
        result=CHiMMIRT(Y,K,Lambda,result,max_iter=200,orthogonal=orthogonal)
        
        mat_prod=cbind(result$Theta,1)%*%t(cbind(result$A,result$d))
        X=stats::plogis(Q*mat_prod)
        deriv.A <- grad_Ad(Q, X, result$Theta)
        check_mat=ifelse(abs(deriv.A)>lambda & result$A==0, 1, 0)
      }
    }
    
    result_set[[i]]=result
    if (progress_shown){
      progress_ratio <- (i - 1) / (length(lambda_values) - 1)
      pb$update(progress_ratio, tokens = list(
        progress=round(progress_ratio * 100),
        lambda=round(lambda)
      ))
      
      Sys.sleep(0.02)
    }
  }
  
  result_set$lambda_values=lambda_values
  cat("\u2714 Whole procedure finished!!!\n")
  
  return(result_set)
}

############################################################
#####-----Regularized HiMIRT (Cardinality version)-----#####
############################################################

card_HiMIRT <- function(Y,K,Params=NULL,random_start=c(TRUE,FALSE),tol=1e-4,max_iter,kappaC,orthogonal=c(TRUE,FALSE),progress_shown=TRUE){
  N=nrow(Y)
  J=ncol(Y)
  Q=2*Y-1
  kappaC <- J*K-kappaC
  
  if (is.null(Params)){
    cat("Starting HiMMIRT!!!\n")
    Params <- HiMMIRT(Y,K,max_iter=200,tol=1e-3,rotation=TRUE,orthogonal=T,progress_shown=progress_shown)
  }
  
  if (random_start){
    a_max <- max(abs(Params$A))
    Params$A <- a_max*matrix(stats::runif(J*K,-1,1),nrow=J)
  }
  
  mat_prod=cbind(Params$Theta,1)%*%t(cbind(Params$A,Params$d))
  X=stats::plogis(Q*(mat_prod))
  H=compute_H(Q,mat_prod)
  Z=compute_Z(mat_prod,Q,H)
  initial_logl_value=logl(X)
  #print(sprintf("The initial logl value is %.2f", logl_values[1]))
  cat("Starting card_HiMIRT\n")
  
  if (progress_shown){
    pb <- progress::progress_bar$new(
      format = "Convergence: [:bar] :progress%  | log-likelihood: :loss | Iter: :iter",
      total = 100,
      clear = FALSE,
      width = 80
    )
  }
  
  logl_values=vector(mode="numeric",length=max_iter)
  
  for (i in 1:max_iter){
    Params$d=colMeans(t(t(-4*Q*H)+Params$d))
    
    Z_D=t(t(Z)-Params$d)
    
    if (orthogonal){
      S=t(Z_D)%*%Z_D/N
      EVD_result=eigen(t(Params$A)%*%S%*%Params$A)
      L=EVD_result$vectors
      Lambda_inv=diag(1/sqrt(EVD_result$values))
      B_right=Params$A%*%L%*%Lambda_inv%*%t(L)
      B=S%*%B_right
      Params$Theta=Z_D%*%B_right
      
      Params$A <- cardize(Params,B,kappaC)
    } else {
      Params$Theta <- Z_D%*%Params$A%*%solve(t(Params$A)%*%Params$A)
      
      Params <- standardiseMIRT(Params,orthogonal=F,admm=F)
      Z_D=t(t(Z)-Params$d)
      
      Theta_squared <- crossprod(Params$Theta)
      alpha=max(eigen(Theta_squared)$values)
      
      for (sub_i in 1:20){
        MM_loss_before <- sum((Z_D-Params$Theta%*%t(Params$A))^2)
        B=t(t(Params$A)-(Theta_squared%*%t(Params$A)-t(Params$Theta)%*%Z_D)/alpha)
        
        Params$A <- cardize(Params,B,kappaC)
        MM_loss_after <- sum((Z_D-Params$Theta%*%t(Params$A))^2)
        
        relative_change_sub <- ifelse(sub_i==1,
                                      10,
                                      (MM_loss_before-MM_loss_after)/MM_loss_before)
        
        if (relative_change_sub<1e-5){
          break
        }
      }
    }
    
    mat_prod=cbind(Params$Theta,1)%*%t(cbind(Params$A,Params$d))
    X=stats::plogis(Q*(mat_prod))
    H=compute_H(Q,mat_prod)
    Z=compute_Z(mat_prod,Q,H)
    logl_values[i]=logl(X)
    #print(sprintf("After %d iterations the value is %.5f",
    #m-2, logl_values[m]))
    
    relative_change=ifelse(i==1,
                           10,
                           (logl_values[i-1]-logl_values[i])/logl_values[i-1])
    
    if (progress_shown){
      if (i==1){
        max_diff <- relative_change
      }
      progress_ratio <- min(1, abs(max_diff - relative_change) / (max_diff - tol))
      pb$update(progress_ratio, tokens = list(
        loss = sprintf("%.3f", logl_values[i]),
        iter = i,
        progress=round(progress_ratio*100)
      ))
    }
    
    if (relative_change < tol){
      break
    }
    
    if (relative_change < 0){
      warning("increase in the objective")
      break
    }
  }
  
  cat("\u2714 Algorithm finished!!!\n")
  Params$value=logl_values[i]
  return(Params)
}

cardize <- function(Params,B,kappaC){
  if (length(kappaC)==1){
    B_rank=order(abs(B))
    Params$A <- B
    Params$A[B_rank[1:kappaC]] <- 0
  }else{
    Params$A <- B
    for (i in 1:K){
      B_rank <- order(abs(B[,i]))
      Params$A[B_rank[1:kappaC[i]],i] <- 0
    }
  }
  
  return(Params$A)
}

# Data output
showResults <- function(result_list){
  df <- data.frame(
    lambda=result_list$lambda_values,
    prob_diff=NA,
    sens=NA,
    spec=NA,
    load_diff=NA,
    BIC=NA
  )
  
  df[2:6] <- t(vapply(result_list[1:length(result_list$lambda_values)], function(result) {
    A_hat <- result$A
    Theta_hat <- result$Theta
    d_hat <- result$d
    
    A_true <- Data$A
    d_true <- Data$d
    Theta_true <- Data$Theta
    
    prob_diff <- sum((plogis(d_hat+A_hat%*%t(Theta_hat))-
                        plogis(d_true+A_true%*%t(Theta_true)))^2)/(N*J)
    
    A_best=ref.perp(A_hat,A_true)
    
    sens <- sum(ifelse(A_best!=0&A_true!=0, 1, 0))/sum(ifelse(A_true!=0,1,0))
    spec <- sum(ifelse(A_best==0&A_true==0, 1, 0))/sum(ifelse(A_true==0,1,0))
    
    load_diff <- sum((A_best-A_true)^2)/(J*K)
    
    BIC = 2*result$logl_value+(sum(result$A!=0)+J+N*K)*log(N)
    return(c(prob_diff,sens,spec,load_diff,BIC))
  }, numeric(5)))
  
  pt_col <- rep("black", length(result_list$lambda_values))
  pt_col[which.min(df$BIC)] <- "red"
  
  plot(df$lambda,df$prob_diff,col=pt_col,pch=16)
  lines(df$lambda,df$prob_diff)
  
  plot(df$lambda,df$sens,col=pt_col,pch=16)
  lines(df$lambda,df$sens)
  
  plot(df$lambda,df$spec,col=pt_col,pch=16)
  lines(df$lambda,df$spec)
  
  plot(df$lambda,df$load_diff,col=pt_col,pch=16)
  lines(df$lambda,df$load_diff)
  
  plot(df$lambda,df$BIC,col=pt_col,pch=16)
  lines(df$lambda,df$BIC)
}

calResults <- function(result_list,Data,N,J,K){
  df <- data.frame(
    lambda=result_list$lambda_values,
    prob_diff=NA,
    sens=NA,
    spec=NA,
    load_diff=NA,
    card_diff=NA,
    BIC=NA
  )
  
  df[2:7] <- t(vapply(result_list[1:length(result_list$lambda_values)], function(result) {
    A_hat <- result$A
    Theta_hat <- result$Theta
    d_hat <- result$d
    
    A_true <- Data$A
    d_true <- Data$d
    Theta_true <- Data$Theta
    
    prob_diff <- sum((plogis(d_hat+A_hat%*%t(Theta_hat))-
                        plogis(d_true+A_true%*%t(Theta_true)))^2)/(N*J)
    
    A_best=ref.perp(A_hat,A_true)
    
    sens <- sum(ifelse(A_best!=0&A_true!=0, 1, 0))/sum(ifelse(A_true!=0,1,0))
    spec <- sum(ifelse(A_best==0&A_true==0, 1, 0))/sum(ifelse(A_true==0,1,0))
    
    load_diff <- sum((A_best-A_true)^2)/(J*K)
    card_diff <- abs(sum(A_best!=0)-sum(A_true!=0))
    
    BIC = 2*result$logl_value+sum(result$A!=0)*J*log(N)
    return(c(prob_diff,sens,spec,load_diff,card_diff,BIC))
  }, numeric(6)))
  
  return(df)
}
