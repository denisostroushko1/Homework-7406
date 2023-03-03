
rm(list = ls())

load("/Users/denisostroushko/Desktop/UofM MS/MS 2023 - 1 Spring/PUBH 7406/Homework-7406/In class bootstrapp/oregon_hie_data.rda")
head(ohie_data)

library(tidyverse)
library(boot)
library(bootstrap)

# 1. 
# Use a non-parametric bootstrap to obtain standard error estimates for the difference
# in the treatment effect when controlling for and not controlling for survey wave and
# household size.

fit1_true <- glm(doc_num_mod_12m ~ treatment + wave_survey12m +
            hhsize * wave_survey12m,
            data = ohie_data,
            family = poisson()
            )

fit2_true <- glm(doc_num_mod_12m ~ treatment,
            data = ohie_data,
            family = poisson()
            )


set.seed(1637894634)
N_iter <- 2000
nrow_data <- nrow(ohie_data)

boot_f <- 
  function(iter, data){
    
    res <- 
      data.frame(
        iter_number = seq(from = 1, to = iter, by = 1), 
        T1 = NA, 
        T2 = NA
      )
    
    for(i in 1:N_iter){
      
      if(i %% 100 == 0){print(paste("On iter:", i))}
      
      resampled <- data[sample(rownames(data), size = nrow(data), replace =  T),  ]
      
      fit1 <- glm(doc_num_mod_12m ~ treatment + wave_survey12m +
            hhsize * wave_survey12m,
            data = resampled,
            family = poisson()
            )

      fit2 <- glm(doc_num_mod_12m ~ treatment,
                  data = resampled,
                  family = poisson()
                  )
      
      res$T1[i] <- summary(fit1)$coefficients[rownames(summary(fit1)$coefficients) == "treatment"][1]
      res$T2[i] <- summary(fit2)$coefficients[rownames(summary(fit2)$coefficients) == "treatment"][1]
      
    }
    
    res$diff <- with(res, T1 - T2)
    
    return(res)
  }

boot_results <- boot_f(iter = N_iter, data = ohie_data)

# 2. 
# Make a visualization of the bootstrap sampling distribution of this estimated difference in effect.


summary(boot_results$diff)
sd(boot_results$diff)

ggplot(data = boot_results, 
       aes(x = diff)) + 
  geom_histogram()

# 3. 
# Construct a confidence interval for this difference using the 
# i) normal approximation method, 
# ii) percentile bootstrap method, 
# iii) and basic bootstrap method. Is the difference significant?

# (I)
lower_bound <- mean(boot_results$diff) - qt(0.975, df = nrow(ohie_data)) * sd(boot_results$diff)
upper_bound <- mean(boot_results$diff) + qt(0.975, df = nrow(ohie_data)) * sd(boot_results$diff)

print(paste0("Mean: ", round(mean(boot_results$diff), 4), " | Normal approximation C.I.: (" ,
             round(lower_bound, 4), ", ", round(upper_bound, 4), ")"))

# (II)
lower_bound_pct <- quantile(boot_results$diff, .025) 
upper_bound_pct <- quantile(boot_results$diff, .975) 

print(paste0("Mean: ", round(mean(boot_results$diff), 4), " | Quantile C.I.: (" ,
             round(lower_bound_pct, 4), ", ", round(upper_bound_pct, 4), ")"))

# (III)

# summary of this approach: take the difference of booted 
# take quartile of the difference 
# then take away true difference from each quartile 

# true difference 
true_fit_1 <- summary(fit1_true)$coefficients[rownames(summary(fit1_true)$coefficients) == "treatment"][1]
true_fit_2 <- summary(fit2_true)$coefficients[rownames(summary(fit2_true)$coefficients) == "treatment"][1]

true_diff <- true_fit_1 - true_fit_2

quantiles <- quantile(boot_results$diff - true_diff, probs = c(0.975, 0.025))
true_diff - quantiles

## Bias estimation: 
mean(boot_results$diff) - true_diff

#############
# 

mean(boot_results$T1)
sd(boot_results$T1)

model_res <- round(c(summary(fit1_true)$coefficients[rownames(summary(fit1_true)$coefficients) == "treatment"]),4)[c(1,2)]

names(model_res) <- c("Est.", "Error")

print(paste0("True coefficient: ", model_res[1], ", with standard error: ", model_res[2]))
print(paste0("Bootleg coefficient: ", round(mean(boot_results$T1), 4), 
             ", with standard error: ", round(sd(boot_results$T1), 4)))

