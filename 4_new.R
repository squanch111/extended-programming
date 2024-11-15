
# main function used for linear mixed structure
lmm = function(form,dat,ref = list()){
  # Set the matrix and parameters required for the model
  setup = LMMsetup(form, dat, ref)
  # Optimize negative log likelihood with optim to get the best theta and beta
  opt_result <- optim(par = setup$init_theta, fn = LMMprof, setup = setup, method = "BFGS")
  # Returns maximum likelihood estimates, including beta and theta
  # 重新计算 LMMprof 以提取 beta_hat
  lmm_result <- LMMprof(opt_result$par, setup)
  beta_hat <- attr(lmm_result, "beta_hat")  # 提取 beta_hat 属性
  neg_log_likelihood <-attr(lmm_result, "neg_log_likelihood")
  theta = opt_result$par[2] #随机效应的theta
  print( opt_result) #optim的结果
  std <- exp(2*theta) #随即效应的方差
  list(beta_hat = beta_hat,theta = std,neg_log_likelihood=neg_log_likelihood)
}



# Create the design matrix(X,Z) and other necessary elements(with initialize theta)
LMMsetup = function(form,dat,ref) {
  # set up variable X and response Y 
  X <- model.matrix(form, data = dat)
  y <- model.response(model.frame(form, data = dat))
  # set up random effects Z with c("z","x","w3"), and combine model.matrix(˜ z:x:w3-1,dat)
  if (length(ref)>0){
    Z_blocks <- lapply(ref, function(vars) {
      interaction_formula <- as.formula(paste("~", paste(vars, collapse = ":"), "-1"))
      model.matrix(interaction_formula, data = dat)
    })
    Z <- do.call(cbind, Z_blocks)
    
    # initialize theta (Includes log sigma and logarithmic standard deviation of random effect variance)
    init_theta <- c(log(sd(y)), 0)
  }
  else{
    #set Z with object matrix but in 0
    Z <- matrix(0, nrow = nrow(dat), ncol = 0)
    # initial only contain one element
    init_theta <- log(sd(y))
    
  }
  print(dim(Z))
  #print(init_theta)

  list(X = X, y = y, Z = Z, init_theta = init_theta)
}



# LMMprof 函数：计算负对数似然和 beta 值
LMMprof <- function(theta, setup) {
  X <- setup$X
  Z <- setup$Z
  y <- setup$y
  # 提取 sigma^2 和随机效应方差参数 psi_diag
  sigma2 <- exp(2 * theta[1])
  psi <- exp(2 * theta[2])  # 计算随机效应方差(同方差)
  #print(psi)
  #psi_diag <- exp(2 * theta)
  # QR 分解 Z 矩阵
  QR_Z <- qr(Z)
  R <- qr.R(QR_Z)
  p <- ncol(R) 
  #print(p)
  psi_diag <- rep(psi,p)
  #print(diag(psi_diag))
  # 构造 (R psi R^T + I_p * sigma^2) 矩阵
  R_psi_Rt <- R %*% diag(psi_diag) %*% t(R) + diag(sigma2, nrow = p)
  
  # 使用 Cholesky 分解 R_psi_Rt 以避免显式求逆
  chol_R_psi_Rt <- chol(R_psi_Rt)
  
  # 使用 qr.qty 将 X 和 y 映射到 Q^T 空间中
  QtX <- qr.qty(QR_Z, X)  # 相当于 Q^T %*% X
  Qty <- qr.qty(QR_Z, y)  # 相当于 Q^T %*% y
  # 取前 p 行计算 WX_upper 和 Wy_upper
  WX_upper <- forwardsolve(t(chol_R_psi_Rt), QtX[1:p, ], upper.tri = FALSE)
  Wy_upper <- forwardsolve(t(chol_R_psi_Rt), Qty[1:p], upper.tri = FALSE)
  
  # 后 n-p 行直接除以 sigma^2
  WX_lower <- QtX[(p + 1):nrow(QtX), ] / sigma2
  Wy_lower <- Qty[(p + 1):length(Qty)] / sigma2
  
  # 合并上下部分
  WX_combined <- rbind(WX_upper, WX_lower)
  Wy_combined <- c(Wy_upper, Wy_lower)
  
  # 将结果映射回原空间
  WX <- qr.qy(QR_Z, WX_combined)
  Wy <- qr.qy(QR_Z, Wy_combined)
  
  # 计算 beta 的估计值
  XtWX <- crossprod(WX)
  XtWy <- crossprod(WX, Wy)
  chol_XtWX <- chol(XtWX)  # Cholesky分解 XtWX
  
  beta_hat <- backsolve(chol_XtWX, forwardsolve(t(chol_XtWX), XtWy, upper.tri = FALSE), upper.tri = TRUE)
  #print(beta_hat)
  # 计算残差和负对数似然
  resid <- y - X %*% beta_hat #原始残差
  #print(resid)
  # 使用 qr.qty 投影残差到 Q^T 空间
  Q_resid <- qr.qty(QR_Z, resid)
  
  # 处理前 p 行和后 n-p 行
  W_resid_upper <- backsolve(chol_R_psi_Rt, forwardsolve(t(chol_R_psi_Rt), Q_resid[1:ncol(R)], upper.tri = FALSE), upper.tri = TRUE)
  W_resid_lower <- Q_resid[(ncol(R) + 1):length(Q_resid)] / sigma2  # 直接处理下部分
  
  # 合并上下部分残差
  W_resid <- c(W_resid_upper, W_resid_lower) #得到加权残差
  
  # 还原残差到原空间
  #W_resid <- qr.qy(QR_Z, W_resid)
  
  #?
  neg_log_likelihood <- sum(W_resid^2) / (2*sigma2)
 
  # 返回带 beta 估计值的负对数似然
  attr(neg_log_likelihood, "beta_hat") <- beta_hat
  attr(neg_log_likelihood, "neg_log_likelihood") <- neg_log_likelihood
  return(neg_log_likelihood)
}


# 加载必要的包和数据集
library(nlme)
library(lme4)

# 测试数据
data("Machines", package = "nlme")

# 使用自定义的 lmm 函数进行拟合
lmm_result <- lmm(score ~ Machine, Machines, list("Worker", c("Worker", "Machine")))

# 使用 lmer 函数拟合相同的模型
lmer_result <- lmer(score ~ Machine + (1 | Worker) + (1 | Worker:Machine), data = Machines, REML = FALSE)

# 输出两个模型的参数估计值，比较结果
cat("lmm 函数的参数估计：\n")
print(lmm_result) #lmm结果
cat("\nlmer 函数的参数估计：\n")
print((lmer_result))    # lmer结果










