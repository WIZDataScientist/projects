#   -----------------------------------------------------------------------
# logit link function
#   -----------------------------------------------------------------------

inverse_logit <- function(value){
  exp(value) / (1 + exp(value))
}

#   -----------------------------------------------------------------------
# Probability to be in category j given input variables/ covariates for person i, xi
# P(Y_i=j|xi)=P(Y_i<=j|xi)−P(Y_i<=j−1|x_i)
# where P(Y_i<=j|x_i)=exp(a_j+b^Txi)(1+exp(a_j+b^Tx_i))P(Y_i<=j|x_i)=exp(a_j+b^Tx_i)(1+exp(a_j+b^Tx_i))
#   -----------------------------------------------------------------------
prob_j         <- function(a, j, b, x_i, min_obs_resp, max_obs_resp){
  
  # probability of being in category <= j
  prob         <- function(z) 
    return(inverse_logit(a[z] + x_i %*% b))
  
  
  # return probability of being in category j is conditioned by the function "log_prob"
  if(j == min_obs_resp)
    return(prob(j))
  else if(j == max_obs_resp)
    return(1 - prob(j-1))
  else
    return(prob(j) - prob(j-1))
  
}

#   -----------------------------------------------------------------------
# The loglikelihood function
#   -----------------------------------------------------------------------
log_likelihood <- function(theta, M, C, N, y, x, min_obs_resp, max_obs_resp){
  # browser()
  # defining parameters
  b <- theta[seq_len(M)]     # parameters 
  a <- theta[M + seq_len(C)] # intercept parameters
  
  lik_contribution <- function(i)
    log(prob_j(a, j = y[i], b, x[i, ], min_obs_resp, max_obs_resp))
  
  # sum for each individual to give loglikelihood value 
  return(-sum(vapply(seq_len(N), lik_contribution, numeric(1))))
}


#   -----------------------------------------------------------------------
# Evaluate and maximize likelihood function
# Maximizing log likelihood conditioned on intercept variables -infty = a0 < a1< a2< a3 < a4 < a5 = infty (the cumulative logit condition)
#   -----------------------------------------------------------------------

eval_likelihood <- function(response, data, theta_init = NULL, outer.iterations = NULL){
  # all sorted ordinal categories
  obs_resp <- unique(data[, response]) %>% sort
  
  # max and min of ordinal categories for response variable
  min_obs_resp <- min(obs_resp)
  max_obs_resp <- max(obs_resp)
  
  # index of col for response variable in data set
  which_response <- which(colnames(data) == response)
  
  # No for data set
  N <- nrow(data)           # no of individuals
  M <- ncol(data) - 1       # no of features
  C <- length(obs_resp) - 1 # no of ordinal categories in responsevariable
  
  # convert data to matrix and define design matrix and output variable
  matrix_data <- as.matrix(data)
  
  x <- matrix_data[, -which_response]
  y <- matrix_data[, which_response]
  
  log_likelihood_optim <- function(theta)
    log_likelihood(theta, M, C, N, y, x, min_obs_resp, max_obs_resp)
  
  # conditions for the mathematical optimization problem
  Amat           <- diag(M + C)
  index_param    <- seq_len(M + C)
  forall_a_param <- index_param[ - seq_len(M + 1)]
  
  for(i in forall_a_param){
    Amat[i,i-1] <- -1
  }
  
  # the condition is Amat >= bvec
  bvec <- c(rep(-10000, M), -10, rep(1e-6, C - 1))
  
  # auto set initial param values
  if(is.null(theta_init)){
    tempini <-  rnorm(M)
    tempini[tempini == 0] <- rnorm(M)
    theta_init <- c(tempini, 1, 2, 3, 4)
  }
  
  # auto set iterations
  if(is.null(outer.iterations))
    outer.iterations <- 1000
  
  
  output <- constrOptim(theta = theta_init, log_likelihood_optim, NULL, Amat, bvec, outer.iterations)
  
  return(list(alpha = output$par[-seq_len(M)], beta = output$par[seq_len(M)]))
}