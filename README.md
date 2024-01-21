# Homework-7406

Homework for Spring 2023 PUBH 7406 at the University of Minnesota. I found a lot of cool ways to present tables with summary statistics, as well as test and models output. I also 
found a lot of useful packages that interact with `ggplot` mainly. 

Also, textbook and lecture materials cover way more intro-level regression methods and theory, and I wish to use this file as a table of  contents. 

Last update: `r Sys.Date()`

By default, I use `kable`, `kableExtra`, `tidyverse`, and `ggplot2` a lot, so I will only be making notes of new packages 
I find and use for the analysis. 

# HW1

* Statistical Concepts: 
  + Two-way ANOVA models. Preliminary data analysis using Box and Line plots, as well as tables rendered with kable 
  + Statistical definition of two-way ANOVA (also refereed to as cell-mean) models and all relevant assumptions 
  + Interpretation of sum-of-squares estimates in the context of such models 
  + Residual diagnostics, similar to topics in previous course [PUBH 7405](https://github.com/denisostroushko1/Homework-7405)
  + Bonferroni adjustment for multiple comparisons
  
* Implementation via R code
  + base R code for `lm` models 
  + use of `emmeans` package combined with fitted `lm` object to statistically compare model fitted values 
  + use of `pwpp` function with an `emmeans` object to produce a plot of model estimates 
  + Use of `emmeans` package to perform Bonferroni adjustment without coding extra manipulations of p-values. 
    - Other adjustments are also available in the `emmeans` package 
  
* Packages: 
  + `emmeans`: useful guides and small examples for this package are [here](https://yuzar-blog.netlify.app/posts/2022-11-29-emmeans/) and [here](https://aosmith.rbind.io/2019/03/25/getting-started-with-emmeans/)

# HW2

* Statistical Concepts: 

  + Binomial family GLM with linear, log, and logit link (logistic regression) examples. Interpretation of coefficients and 
    quantities. 
  + Probability rate of change at given values of independent predictor. 
  + Use of the Delta Method to obtain standard errors for odds ration using coefficients from the logistic regression 
  + Hypothesis testing for logistic regression coefficients using both Wald test and Deviance test 
    - Wald test uses standard error estimate from a fitted glm 
    - Deviance test uses the drop in the deviance statistic between models with and without a predictor of interest 
  + Variance calculation for a linear predictor using basic theory of statistics and probability 
  
* Implementation via R code
  + All calculations are done with the use of functions and base R code 
  + Delta methods functions use variance-covariance matrices of fitted glm objects 

# HW3

* Statistical Concepts: 
  + Logistic regression models with interaction terms. Estimation of conditional odds ratios. Interpretation of conditional 
    odds ratios arising from the interaction terms. 
  + More hypothesis testing for model terms in logistic regression 
  + Holm-Bonferroni Adjustment of p-values for control of the Family Wise Error Rate (FWER)
  + Benjamini-Hochberg Adjustment of p-value to control for the False Discovery Rate 
    - continuation of the Bonferroni and Tukey Adjustment applications 
  + Deviance test for the logistic regression goodness of fit. 
  
* Implementation via R code
  + more base R code for calculation of standard errors using the Delta method
  + More applications of `emmeans` for contrasts and pairwise mean comparisons 
  + Examples of p-value adjustment using `p.adjust` for Holm and BH adjustments 
  + Examples of deviance test in R using `anova` function and two `glm` objects 
  + Implementation of the Delta Methods and Bootstrap methods to calculate the ratio of odds ratios. 
    - Comparison of results to validate the proper application of the delta method 


# HW4

* Statistical Concepts: 

  + Poisson family GLM (Poisson Regression) application. Compassion of variance of coefficients using Poisson, Quasipoisson, and Bootstrap methods. Illustration that bootstrapping coefficients tends to produce more optimistic (smaller) standard errors 
  + Estimation of confidence intervals using normal approximation and quantile method from bootstrap sampling distributions. 
    - comparison with the quasipoisson confidence intervals 
  + The Beverton-Holt Model for stable population levels. Application of bootstrap method to find variance of the stable population estimator. 
  + 
  
* Implementation via R code
  + Use of `for` loops and base R code to implement a bootstrap simulation. 

# HW5

* Statistical Concepts: 

  + Multiple Imputation for missing data and application of the Rubin rules 
  + use of pooled standard error estimates from the multiple imputation data sets to conduct statistical inference 
  + First look at the simulation study approach to compare imputation methods 
  + Simulation study aims to compare imputation methods in terms of bias and MSE
  
* Implementation via R code
  + use of `mice` function and package to leverage predictive mean matching, Bayesian linear models, nearest neighbors, etc... 
    to imputed data 
  + use `plot` for imputation results to assess converge of the imputation algorithms 
  + 

* Packages: 
  + `mice`: main package to easily implement multiple imputation using various modeling methods. Package is well documented and 
    has examples in [this book](https://stefvanbuuren.name/fimd/)

# HW6

* Statistical Concepts: 

  + First look at the use of Lasso shrinkage for logistic regression model. 
  + Variable selection for the predictive classifier using LASSO model. 
  + Basic examples of K (10)-fold cross validation to understand possible variability of classifier performance. AUC is used as an example metric
  + Use of calibration plots to understand predicted probabilities 
  
* Implementation via R code
  + `val.prob` function applied to a vector of predicted values to construct a calibration plot 
  + extracting information from `roc` to construct pretty ROC curve plots using ggplot
