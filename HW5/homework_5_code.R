

set.seed(8172013)
library(mvtnorm)
expit <- function(x) {
  expit <- exp(x)/(1+exp(x))
  return(expit)}

# Number of Datasets
S <- 100

# Sample Size
n <- 200

# Parameter and Standard Error Estimates 
parm.est <- matrix(0, nrow = S, ncol = 5)
se.est <- matrix(0, nrow = S, ncol = 5)
colnames(parm.est) <- colnames(se.est) <- paste("Method", letters[1:5], sep = ".")

for (q in 1:S) {
  X_123 <- rmvnorm(n, mean = rep(0, 3), 
                   sigma = matrix(c(1, 0.6, 0.36, 0.6, 1, 0.6, 0.36, 0.6, 1), 
                                  nrow = 3, ncol = 3))
  X1 <- X_123[, 1]
  X2 <- X_123[, 2]
  X3 <- X_123[, 3]
  X4 <- rbinom(n, 1, expit(0.5*X1 + 0.5*X2 - 0.5*X3))
  X5 <- rbinom(n, 1, expit(0.5*X1 - 0.5*X2 + 0.5*X3 - 0.5*X4))
  Y <- 2 + X1 + X2 + X3 + 2*X4 + 2*X5 + rnorm(n, mean = 0, sd = 3)
  
  # Add in some missing data
  X1_obs <- rbinom(n, 1, expit(1 + 0.25*X2 + 0.25*X3 + 0.125*X4 + 0.125*X5))
  X4_obs <- rbinom(n, 1, expit(1 + 0.25*X1 + 0.25*X2 + 0.25*X3 +  0.125*X5))
  X1 <- ifelse(X1_obs == 1, X1, NA)
  X4 <- ifelse(X4_obs == 1, X4, NA)
  final_data <- data.frame(Y, X1, X2, X3, X4, X5)
  
  # Method A
  method.A <- summary(lm( Y ~ X1 + X2 + X3 + X4 + X5, 
                          data = final_data))$coefficients
  method.A.coef <- method.A[2, 1]
  method.A.se <- method.A[2, 2]
  
  # Method B
  X1_impB <- ifelse(is.na(X1) == TRUE, mean(X1, na.rm = TRUE), X1)
  X4_impB <- ifelse(is.na(X4) == TRUE, mean(X4, na.rm = TRUE), X4)
  final_data <- data.frame(final_data, X1_impB, X4_impB)
  method.B <- summary(lm( Y ~ X1_impB + X2 + X3 + X4_impB + X5, 
                          data = final_data))$coefficients
  method.B.coef <- method.B[2, 1]
  method.B.se <- method.B[2, 2]
  
  # Method C
  X1_mod <- lm(X1 ~ X2 + X3 + X5, data = final_data)
  X4_mod <- glm(X4 ~ X2 + X3 + X5, data = final_data, family = "binomial")
  X1_impC <- ifelse(is.na(X1) == TRUE, predict(X1_mod), X1)
  X4_impC <- ifelse(is.na(X4) == TRUE, predict(X4_mod, type = "response"), X4)
  final_data <- data.frame(final_data, X1_impC, X4_impC)
  method.C <- summary(lm( Y ~ X1_impC + X2 + X3 + X4_impC + X5, 
                          data = final_data))$coefficients
  method.C.coef <- method.C[2, 1]
  method.C.se <- method.C[2, 2]
  
  # Method D
  impD <- mice(final_data[, 1:6], maxit = 20, m = 1, print = FALSE)
  method.D <- summary(lm( Y ~ X1 + X2 + X3 + X4 + X5, 
                          data = complete(impD)))$coefficients
  method.D.coef <- method.D[2, 1]
  method.D.se <- method.D[2, 2]
  
  # Method E
  impE <- mice(final_data[, 1:6], maxit = 20, m = 20, print = FALSE)
  method.E <- summary(pool(with(impE, lm( Y ~ X1 + X2 + X3 + X4 + X5))))
  method.E.coef <- method.E[2, 2]
  method.E.se <- method.E[2, 3]
  
  # All methods
  parm.est[q, ] <- c(method.A.coef, method.B.coef, 
                     method.C.coef, method.D.coef, method.E.coef)
  se.est[q, ] <- c(method.A.se, method.B.se, 
                   method.C.se, method.D.se, method.E.se)
  
  print(q)
}

