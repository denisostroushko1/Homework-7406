
# List of all packages that we need for the 7405 assignments and exams 

library(MASS)
require(tidyverse) # require instead of library to make sure that other packages do not overwrite tidyverse packages 
library(kableExtra)
library(readxl)
library(gridExtra)
library(ggeffects)
library(mltools) # one hot encoding outside of caret package 
library(data.table) # need this for mltools to work 
library(olsrr) # a better package for stepwise regression 
library(DescTools)
library(car)
library(broom) # For converting models into data frame
library(eulerr)        # For creating Euler and Venn diagrams
library(MatchIt)
library(caret)
library(survminer)
require(survival)
library(flexsurv)
library(emmeans)