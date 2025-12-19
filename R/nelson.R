## ===============================
## Nelson–Siegel related functions
## ===============================

library(INLA) 
library(rSPDE)
library(fmesher)
library(inlabru)
library(Matrix)

precision.ar1 = function(dates, rho, prec){
  N <- max(dates)
  Q = Matrix(0, N, N)
  diag(Q) = rep(1 + rho^2, nrow(Q))
  for (i in 1:(N-1)) {
    Q[i, i+1] = -rho
    Q[i+1, i] = -rho
  }
  Q[1,1] = 1
  Q[N,N] = 1
  return(prec * Q)
}

build_nelsonsiegel_model <- function(dates, maturity , lambda){
  Q1_temp = precision.ar1(dates, 1, 1)
  Q2_temp = precision.ar1(dates, 1, 1)
  Q3_temp = precision.ar1(dates, 1, 1)
  
  Q_graph = bdiag(Q1_temp,Q2_temp,Q3_temp)
  
  tmat <- length(maturity)
  N <- max(dates) - min(dates) + 1
  
  A_mat <- function(dates, lambda, tmat){
    N <- max(dates) - min(dates) + 1
    p1 = rep(1, tmat)
    p2 = (1 - exp(-lambda * maturity))/(lambda * maturity)
    p3 = (1 - exp(-lambda * maturity))/(lambda * maturity) - exp(-lambda * maturity)
    A_1 = Matrix(diag(1, N, N))
    A_2 = Matrix(diag(p2[1], N, N))
    A_3 = Matrix(diag(p3[1], N, N))
    for (i in 2:tmat){
      A_1 = rbind(A_1, diag(p1[i], N,N))
      A_2 = rbind(A_2, diag(p2[i], N,N))
      A_3 = rbind(A_3, diag(p3[i], N,N))
    }
    return (cbind(A_1,A_2,A_3))
  }
  
  
  'inla.rgeneric.nelsonsiegel.model' <- function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const","log.prior", "quit"), theta = NULL) {
    interpret.theta <- function() {
      return(
        list(prec_1 = exp(theta[1L]),
             rho_1 = 1 / (1 + exp(-theta[2L])),
             prec_2 = exp(theta[3L]),
             rho_2 = 1 / (1 + exp(-theta[4L])),
             prec_3 = exp(theta[5L]),
             rho_3 = 1 / (1 + exp(-theta[6L])))
      )
    }
    
    Q <- function() {
      require(Matrix)
      N = max(dates) - min(dates) + 1
      param <- interpret.theta()
      precision.ar1 = function(N, rho, prec){
        Q = Matrix(0, N, N)
        diag(Q) = rep(1 + rho^2, nrow(Q))
        for (i in 1:(N-1)) {
          Q[i, i+1] = -rho
          Q[i+1, i] = -rho
        }
        Q[1,1] = 1
        Q[N,N] = 1
        return(prec * Q)
      }
      
      N = max(dates) - min(dates) + 1
      rho_1 <- param$rho_1
      prec_1 <- param$prec_1
      rho_2 <- param$rho_2
      prec_2 <- param$prec_2
      rho_3 <- param$rho_3
      prec_3 <- param$prec_3
      
      Q1 = precision.ar1(N, rho_1, prec_1)
      Q2 = precision.ar1(N, rho_2, prec_2)
      Q3 = precision.ar1(N, rho_3, prec_3)
      
      Q_all = bdiag(Q1,Q2,Q3)
      return(Q_all)
    }
    
    
    graph <- function(){
      require(Matrix)
      return(Q_graph)
    }
    
    
    mu = function(){
      return(numeric(0))
    }
    
    
    log.norm.const <- function() {
      return(numeric(0))
    }
    
    
    log.prior <- function() {
      param = interpret.theta()
      res <- dgamma(param$prec_1, 1, 5e-05, log = TRUE) + log(param$prec_1) +
        log(1) + log(param$rho_1) + log(1 - param$rho_1)
      res <- res + dgamma(param$prec_2, 1, 5e-05, log = TRUE) + log(param$prec_2) +
        log(1) + log(param$rho_2) + log(1 - param$rho_2)
      res <- res + dgamma(param$prec_3, 1, 5e-05, log = TRUE) + log(param$prec_3) +
        log(1) + log(param$rho_3) + log(1 - param$rho_3)
      return(res)
    }
    
    
    initial <- function() {
      return(rep(0,6))
    }
    
    
    quit <- function() {
      return(invisible())
    }
    
    res <- do.call(match.arg(cmd), args = list())
    return(res)
    
  }
  
  nelsonsiegel_model <- inla.rgeneric.define(inla.rgeneric.nelsonsiegel.model, dates = dates, Q_graph = Q_graph, N = N)
  nelsonsiegel_model$n_Q <- nrow(Q_graph)
  nelsonsiegel_model$dates <- dates
  nelsonsiegel_model$tmat <- tmat
  nelsonsiegel_model$A_mat <- A_mat
  nelsonsiegel_model$lambda <- lambda
  nelsonsiegel_model$Q <- Q_graph
  nelsonsiegel_model$N <- N
  class(nelsonsiegel_model) <- c("nelsonsiegel_model", class(nelsonsiegel_model))
  return(nelsonsiegel_model)
}

register_nelsonsiegel_inlabru <- function(){
  library(inlabru)
  
  bru_get_mapper.nelsonsiegel_model <- function(model, ...) {
    mapper <- list(model = model)
    inlabru::bru_mapper_define(mapper, new_class = "bru_mapper_nelsonsiegel_model")
  }
  
  ibm_n.bru_mapper_nelsonsiegel_model <- function(mapper, ...) {
    model <- mapper[["model"]]
    return(model$n_Q) # return nrow(Q)
  }
  
  ibm_values.bru_mapper_nelsonsiegel_model <- function(mapper, ...) {
    seq_len(inlabru::ibm_n(mapper)) # can use the same
  }
  
  ibm_jacobian.bru_mapper_nelsonsiegel_model <- function(mapper, input, ...) {
    model <- mapper[["model"]]
    dates <- model$dates
    tmat <- model$tmat
    lambda <- model$lambda
    rank_index <- as.numeric(factor(input[, 2], levels = sort(unique(input[, 2]))))
    input <- cbind(input, rank_index)
    slice <- c(input[, 1] + length(dates) * (input[, 3]-1))
    A_tmp <- model$A_mat(dates,lambda,tmat)[slice,]
    return(A_tmp)
  }
  .S3method("bru_get_mapper", "nelsonsiegel_model", bru_get_mapper.nelsonsiegel_model)
  .S3method("ibm_n", "bru_mapper_nelsonsiegel_model", ibm_n.bru_mapper_nelsonsiegel_model)
  .S3method("ibm_values", "bru_mapper_nelsonsiegel_model", ibm_values.bru_mapper_nelsonsiegel_model)
  .S3method("ibm_jacobian", "bru_mapper_nelsonsiegel_model", ibm_jacobian.bru_mapper_nelsonsiegel_model)
}

