

set.seed(8172013)
library(mvtnorm)
expit <- function(x) {
  expit <- exp(x)/(1+exp(x))
  return(expit)}

# Number of Datasets
S <- 100

# Sample Size
N <- seq(from = 100, to = 500, by = 100)

# Intercepts: control amount of missingness 
INT <- seq(from = 0.5, to = 2, by = .5)



# Parameter and Standard Error Estimates 

# in the original code we have 5 columns for 5 methods: 
  # we need to add N, INT, so add two more columns 
  # also I want to track how much X1 and X4 is missing, so add two more columns 

  # nrow also changes, according to the M x N counting rule! 

empty_mat <- matrix(0, nrow = S * length(N) * length(INT), ncol = 9)

parm.est <- empty_mat
se.est <- empty_mat

colnames(parm.est) <- colnames(se.est) <- 
  c(paste("Method", letters[1:5], sep = "."), "N", "INT", "X1_p_missing", "X4_p_missing")

## S <- controls the number of iterations we use to create a sampling distribution 
## N <- sample size for the data: currently, we have n = 200. We should probably try 100, 200, ... , 1000 
## INT <- intercept
    ## currently: X1_obs <- rbinom(n, 1, expit(1 + 0.25*X2 + 0.25*X3 + 0.125*X4 + 0.125*X5))
    ## thus, intercept is currently set to 1, try to use 0.5 to 2 with a .25 increment 

iter = 1

for (q in 1:S) {
  for(i in 1:length(N)){
    for(j in 1:length(INT)){
      
      n <- N[i]
      
      X_123 <- rmvnorm(n, mean = rep(0, 3), 
                       sigma = matrix(c(1, 0.6, 0.36, 0.6, 1, 0.6, 0.36, 0.6, 1),  # this one is supposed to be a 3x3 matrix 
                                                                                  # to include all variances and covariances 
                                      nrow = 3, ncol = 3))
      X1 <- X_123[, 1] # store individual variables 
      X2 <- X_123[, 2]
      X3 <- X_123[, 3]
      X4 <- rbinom(n, 1, expit(0.5*X1 + 0.5*X2 - 0.5*X3)) 
      X5 <- rbinom(n, 1, expit(0.5*X1 - 0.5*X2 + 0.5*X3 - 0.5*X4))
      Y <- 2 + X1 + X2 + X3 + 2*X4 + 2*X5 + rnorm(n, mean = 0, sd = 3) # true coefficients are here. They are supposed to be 1 and 2
      
      # Add in some missing data
      X1_obs <- rbinom(n, 1, expit(INT[j] + 0.25*X2 + 0.25*X3 + 0.125*X4 + 0.125*X5))
      X1_p_missing <- (n - sum(X1_obs))/sum(X1_obs)
      
      X4_obs <- rbinom(n, 1, expit(INT[j] + 0.25*X1 + 0.25*X2 + 0.25*X3 +  0.125*X5))
      X4_p_missing <- (n - sum(X4_obs))/sum(X4_obs)
      
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
      

      
      parm.est[iter, ] <- c(method.A.coef, method.B.coef, 
                         method.C.coef, method.D.coef, method.E.coef, n, INT[j], X1_p_missing, X4_p_missing)
      
      se.est[iter, ] <- c(method.A.coef, method.B.coef, 
                         method.C.coef, method.D.coef, method.E.coef, n, INT[j], X1_p_missing, X4_p_missing)
      
      print(paste0(iter, " of ", S * length(N) * length(INT)))
      iter <- iter + 1
    }
  }
}

#   summary(parm.est)

#   summary(se.est)

# save estimates to help wtih rendering HW. Simulation takes a few minutes to run 
write.csv(parm.est, "./HW5/problem 2 parameter estimates extra credit.csv")
write.csv(se.est, "./HW5/problem 2 errors estimates extra credit.csv")
