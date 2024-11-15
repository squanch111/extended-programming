# Contribution:
# Team members:
#   Member1 (S2702040)
#   Member2 (S2663848)
#   Member3 (S2678971)
#
# Team Member 1(35%):
#   Responsible for implementing and refining the "LMMprof()" function.
#   Implements the Cholesky decomposition for calculating beta_hat and ensuring
#   the numerical stability of the code.
#
# Team Member 2:(35%)
#   Responsible for debugging, testing, and ensuring that the overall code meets
#   the task requirements.
#   Ensures code readability and maintainability, and makes code logic easier to understand.
#
# Team Member 3(30%):
#   Responsible for the main framework design of the code, including the writing
#   of the core functions "lmm()" and "LMMsetup()".
#   Optimized parameter initialization after many tests to make the model more stable.



# Overview:
# This code implements a linear mixed model (LMM) framework to estimate the effects
# of both fixed and random factors on response data.This implementation provides
# an alternative to the "lmer" function from the "lme4" package, allowing users to
# estimate model parameters, through maximum likelihood estimation.
#
# The code consists of three main functions:
# 1. "lmm()": The main function for fitting the linear mixed model to data by setting
#             up model matrices and optimizing the negative log-likelihood function
#             using the “optim” function
#             
# 2. "LMMsetup()": Prepares the design matrices for fixed effects (X) and random
#             effects (Z) and other parameters based on the model formula 
#
# 3. "LMMprof()": The likelihood function that calculates the negative log-likelihood
#             value of the model based on the current parameters
#
# Comparison with "lmer":
#   Requires the "lme4" package for linear mixed model comparisons.
#   To validate the implementation, this code also compares results from the custom
# linear mixed model function with the "lmer" function from the "lme4" package.
# It verifies that both fixed effects and random effects obtained align closely
# with those from "lmer"
### Note: The theta values output by the lmm function are in log-standard-deviation form,
### so a transformation may be needed when comparing with the parameters from "lmer".
### (Using "exp(theta)" in lmm function can get the same value as "lmer")



# lmm() function
# Purpose:
#   This function estimates the fixed and random effects of a linear mixed model
#   using maximum likelihood estimation. The function performs optimization on the
#   negative log-likelihood function, using the "optim" function to find the optimal 
#   values for variance parameters (theta). The resulting model provides estimates
#   for both the fixed effects (beta) and the random effect (theta)
# Inputs:
#   form: Formula for the fixed effects
#   dat: Data frame containing the data for the model
#   ref: List of random effect groupings
# Outputs:
#   A list containing the estimated fixed effects (beta_hat) and optimized
#   random effects (theta)

lmm = function(form, dat, ref = list()) {
  # Set up matrices and initial parameters needed for optimization
  setup = LMMsetup(form, dat, ref)
  
  # Optimize the negative log-likelihood to find the best parameters (theta)
  # LMMprof is the Function to be minimized; In this code, the function is the
  # negative log-likelihood function defined in LMMprof() function.
  opt_result <- optim(par = setup$init_theta, fn = LMMprof, setup = setup, method = "BFGS")
  
  # Recalculate LMMprof to extract beta_hat using the optimized theta
  lmm_result <- LMMprof(opt_result$par, setup)
  beta_hat <- attr(lmm_result, "beta_hat") # Extract beta_hat
  
  # Return estimated fixed effects (beta) and optimized random effects (theta)
  list(beta_hat = beta_hat, theta = opt_result$par)
}



# LMMsetup() function
# Purpose:
#   Creates the design matrix X for fixed effects and Z for random effects.
#   Extract response value y and initialize theta with log-transformed variance
# Inputs:
#   form: Formula specifying the fixed effects
#   dat: Data frame containing the data
#   ref: List of variable names to be used as random effects
# Outputs:
#   A list with X, y, Z, and initial theta values

LMMsetup = function(form, dat, ref) {
  # Fixed effects design matrix (X) and response variable (y)
  X <- model.matrix(form, data = dat)
  y <- model.response(model.frame(form, data = dat))
  
  # Initialize Z (random effects design matrix) and theta structure
  Z_blocks <- list()
  init_theta <- c()  #  store the initial variance parameters for each random effect block
  block_sizes <- c() #  storing the size of each random effect block
  
  # Generate random effects design matrix (Z)
  if (length(ref) > 0) {
    # Create a list of model matrices for each grouping in 'ref'
    # Each element of 'ref' specifies a group of variables for random effects.
    # Use interaction terms to construct random effect design matrices.
    Z_blocks <- lapply(ref, function(vars) {
      # Define a formula for interactions between the specified variables in 'vars'
      # For each grouping, create a model matrix without an intercept using "-1"
      interaction_formula <- as.formula(paste("~", paste(vars, collapse = ":"), "-1"))
      # Construct the model matrix for the current grouping using the interaction formula
      Z_block <- model.matrix(interaction_formula, data = dat)
      
      # Gets the size of the current random effect block (ncol)
      block_size <- ncol(Z_block)
      # Add the current block size to the block_sizes vector
      block_sizes <<- append(block_sizes, block_size)
      # Assign an initial variance parameter (log(sigma)) to the block
      init_theta <<- c(log(sd(y)) , rep(log(10),length(block_sizes)))
      return(Z_block)
    })
    # Combine each individual random effect matrix to create matrix 'Z'
    # This matrix 'Z' will contain each random effect specified in 'ref'.
    Z <- do.call(cbind, Z_blocks)
  } 
  else {
    # If there are no random effects, Z will be empty.
    Z <- matrix(0, nrow = nrow(dat), ncol = 0) 
    init_theta <- log(sd(y)) # Initialize only sigma component
  }
  
  # Return a list with X, y, Z, initial theta values and block_sizes.
  list(X = X, y = y, Z = Z, init_theta = init_theta, block_sizes = block_sizes)
}



# LMMprof() function
# Purpose:
#   Calculates the log-likelihood based on theta values by using matrix decomposition
#   to avoid direct inversion and improve computational stability.
# Inputs:
#   theta: Vector containing random effect parameters
#   setup: List of model matrices and parameters set up in LMMsetup() function
# Outputs:
#   Negative log-likelihood with an attribute for the estimated fixed effects (beta_hat)

LMMprof <- function(theta, setup) {
  # Get X, Z, y and block_sizes
  X <- setup$X
  Z <- setup$Z
  y <- setup$y
  block_sizes <- setup$block_sizes
  # Define the number of samples and number of random effects
  n <- length(y)  # sample size
  # number of random effects (columns in the random effects design matrix Z)
  p <- ncol(Z)    
  
  # Extract sigma^2
  sigma2 <- exp(2 * theta[1]) 
  
  #when colomns of Z is bigger than 0, which means there is random effect
  if(p > 0){ 
    # Construct the vector for random effect variances (psi_vector) based on block 
    # sizes in Z. Each element of theta corresponds to a block of random effects in Z
    # For each block in Z, replicate the corresponding variance value to match the 
    # number of columns in that block.
    psi_vector <- unlist(lapply(seq_along(block_sizes), function(i) {
      # Since theta represents the log standard deviations of the random effect variances, 
      # we apply exp(2 * theta) here to transform theta into the form of variances.
      rep(exp(2 * theta[i + 1]), block_sizes[i])
    }))
    
    # Convert psi_vector into a diagonal matrix (psi_diag), which represents the
    # variance-covariance matrix for the random effects,
    psi_diag <- diag(psi_vector)
    
    # Perform QR decomposition on Z
    QR_Z <- qr(Z)
    R <- qr.R(QR_Z)
    
    # Construct the inverse covariance matrix M without including the Q matrix
    # M is a block diagonal matrix that will be used to represent the inverse of
    # the covariance matrix for the random effects and residuals.
    # The matrix M is constructed using two main blocks: M_upper and M_lower.
    
    # Compute the Cholesky decomposition of the matrix (R*psi_diag*R^t + sigma2*I_p),
    # which represents the covariance of the random effects in the transformed space.
    # U is an upper triangular matrix from the Cholesky decomposition.
    U <- chol(R %*% psi_diag %*% t(R) + diag(sigma2, p))
    
    # Construct M_upper, the inverse of (R*psi_diag*R^t + sigma2*I_p)
    # M_upper is calculated by inverting U via backsolve
    M_upper <- backsolve(U, backsolve(U, diag(p), transpose = TRUE))
    
    # Construct M_lower, which represents the inverse of the residual variance for
    # the remaining n - p observations.
    # M_lower is a (n - p) x (n - p) diagonal matrix with entries 1 / sigma2.
    # This part corresponds to the inverse of the covariance for the residuals in 
    # the lower right block.
    M_lower <- diag(1 / sigma2, n - p)
    
    # Construct the block diagonal matrix M using cbind and rbind
    #   Here, cbind is used to horizontally bind the M_upper block, creating the top
    #   row of blocks for M.
    #   Similarly, cbind combines a zero matrix of size (n - p) x p with M_lower,
    #   creating the bottom row of blocks for M.
    #   Finally, rbind vertically binds these two rows, resulting in a block diagonal
    #   matrix M of size n x n.
    M <- rbind(cbind(M_upper, matrix(0, p, n - p)), cbind(matrix(0, n - p, p), M_lower))
    
    # Calculate Q^t*X and Q^t*y
    Qt_X <- qr.qty(QR_Z, X)
    Qt_Y <- qr.qty(QR_Z, y)
    
    # Calculate W*X and W*y, where W = (Z * psi_theta * Z^t + sigma^2 * I)^(-1)
    W_X <- qr.qy(QR_Z, M %*% Qt_X)
    W_y <- qr.qy(QR_Z, M %*% Qt_Y)
    
    # Calculate X^t*W*X and X^t*W*y
    Xt_W_X <- t(X) %*% W_X
    Xt_W_y <- t(X) %*% W_y
    
    # Beta_hat = (X^t*W*X)^(-1) * X^t*W*y
    # Calculate beta_hat using the Cholesky decomposition
    L <- chol(Xt_W_X)
    beta_hat <- backsolve(L, forwardsolve(t(L), Xt_W_y))
    
    # using the Cholesky decomposition to compute the log-determinant of the 
    # covariance matrix W.
    log_det_M <-  2*sum(log(diag(chol(R %*% psi_diag %*% t(R) + diag(sigma2, p))))) + (n - p) * log(sigma2)
    
    # Calculate the residual_term for negative log_likelihood
    residual_term <- t(y - X %*% beta_hat) %*% (W_y - W_X %*% beta_hat)
  }
  
  # when there are no random effects
  else{
    # W is simply the inverse of the residual variance matrix
    W <- diag(1/sigma2, n)
    
    # Since there are no random effects, beta_hat is computed directly using W.
    Xt_W_X <- t(X) %*% W %*% X  # (X^T * W * X)
    Xt_W_y <- t(X) %*% W %*% y  # (X^T * W * y)
    
    # Beta_hat = (X^t*W*X)^(-1) * X^t*W*y
    # Calculate beta_hat using the Cholesky decomposition
    L <- chol(Xt_W_X)
    beta_hat <- backsolve(L, forwardsolve(t(L), Xt_W_y))
    
    # Calculate the log-determinant of W
    log_det_M <- n * log(sigma2)
    
    # Calculate residual term for negative log-likelihood
    residual_term <- t(y - X %*% beta_hat) %*% W %*% (y - X %*% beta_hat)
  }
  
  # Calculate the negative log-likelihood 
  neg_log_likelihood <- 0.5 * (residual_term + log_det_M)
  
  # Attach the estimated fixed effects beta_hat as an attribute to the negative
  # log-likelihood
  attr(neg_log_likelihood, "beta_hat") <- beta_hat
  
  return(neg_log_likelihood)
}

