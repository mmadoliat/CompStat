rm(list = ls())
library(PermAlgo);library(mvtnorm);library(survival);library(ggfortify);library(gridExtra)
library(grid);library(ggplot2) ;library(fdapace);library(mgcv);library(tidyverse)  
library(fda);library(refund);library(ggpubr);library(tidyfun)
library(dplyr);library(tibble);library(tidyr);library(refundr)
library(purrr);library(forcats); library(funData); library(MFPCA)


######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## Data generation ######## ######## ######## ######## ######## ######## 
######## ######## ######## ######## ######## ######## ######## ######## ######## 
NNN=300
Results=CR=matrix(0,NNN,3)
Ctime=c()
#for(kk in 1:NNN){
kk=1
  set.seed(kk)
  nsample <- 100 # number of indivduals
  # real values of parameters
  gamma <- rep(-0.5,6)
  sigma <- 0.5
  alpha <- c(1,-1)
  rho <- 0.5
  Sigma <- (1-rho)*diag(14)+rho
  Real=c(alpha,gamma)
  # gap between longitudinal measurements
  gapLongi <- 0.05
  gap <- 0.01
  # observed times for longitudinal marker
  followup <- 2
  t <- seq(from = 0, to = followup, by = gap)
  timesLongi <- t[which(round(t - round(t / gapLongi, 0) * gapLongi, 6) == 0)] # visit times
  time <- rep(t, nsample)
  # max. number of individual measurements
  nmesindiv <- followup / gap + 1
  # max. total number of longitudinal measurements
  nmesy <- nmesindiv * nsample
  #  individual id for longitudinal
  idY <- rep(1:nsample, each = nmesindiv)
  
  ###############
  
  b <- rmvnorm(nsample, rep(0, 14), Sigma)
  b1 <- rep(b[, 1], each = nmesindiv) # random intercept Y1
  b2 <- rep(b[, 2], each = nmesindiv) # random slope Y1
  b3 <- rep(b[, 3], each = nmesindiv) # random slope Y1
  
  b4 <- rep(b[, 4], each = nmesindiv) # random intercept Y1
  b5 <- rep(b[, 5], each = nmesindiv) # random slope Y1
  b6 <- rep(b[, 6], each = nmesindiv) # random slope Y1
  
  
  b7 <- rep(b[, 7], each = nmesindiv) # random intercept Y1
  b8 <- rep(b[, 8], each = nmesindiv) # random slope Y1
  b9 <- rep(b[, 9], each = nmesindiv) # random slope Y1
  
  
  b10 <- rep(b[, 10], each = nmesindiv) # random intercept Y1
  b11 <- rep(b[, 11], each = nmesindiv) # random slope Y1
  b12 <- rep(b[, 12], each = nmesindiv) # random slope Y1
  
  
  b13 <- rep(b[, 13], each = nmesindiv) # random intercept Y1
  b14 <- rep(b[, 14], each = nmesindiv) # random slope Y1

  x1 <- rnorm(nsample, 0, 1)
  x2 <- rbinom(nsample, 1, 0.6)
  X1 <- rep(x1, each = nmesindiv)
  X2 <- rep(x2, each = nmesindiv)
  w1 <- rnorm(nsample, 0, 1)
  W1<- rep(w1, each = nmesindiv)
  w2 <- rbinom(nsample, 1, 0.5)
  W2<- rep(w2, each = nmesindiv)
  
  # linear predictor
  # ضرایب تابعی واقعی
  true_beta0 <- function(t) {
    2-2*t 
    # (t-followup/2)
  }
  
  true_beta1 <- function(t) {
    (t-followup/2)^2
  }
  
  true_beta2 <- function(t) {
    (t-followup/2)^3
    
    
  }
  
  true_beta0(t)
  true_beta1(t)
  true_beta2(t)
  
  Beta0 <- rep(true_beta0(t), nsample)
  Beta1 <- rep(true_beta1(t), nsample)
  Beta2 <- rep(true_beta2(t), nsample)
  eta1 <- Beta0+Beta1 * X1 + Beta2* X2 +(b1+b2*time+b3*time^2)
  Y1 <- rnorm(nmesy, eta1, sqrt(sigma))##### Polynomial 
  
  # linear predictor
  # ضرایب تابعی واقعی
  true_beta0 <- function(t) {
    ;1-2*t
  }
  
  true_beta1 <- function(t) {
    ;1
  }
  
  true_beta2 <- function(t) {
    ;1
  }
  
  
  true_beta0(t)
  true_beta1(t)
  true_beta2(t)
  
  Beta0 <- rep(true_beta0(t), nsample)
  Beta1 <- rep(true_beta1(t), length(t)*nsample)
  Beta2 <- rep(true_beta2(t), length(t)*nsample)
  
  eta2 <- Beta0+Beta1 * X1 + Beta2* X2 +(b4+b5*time)
  Y2 <- rnorm(nmesy, eta2, sqrt(sigma))##### Liear
  ###################################################
  true_beta0 <- function(t) {
    1-t^2
  }
  
  true_beta1 <- function(t) {
    0.5+6*t-4*t^2-6*t^3  }
  
  true_beta2 <- function(t) {
    1+3*t-2*t^2-6*t^3
  }
  

  true_beta0(t)
  true_beta1(t)
  true_beta2(t)
  
  Beta0 <- rep(true_beta0(t), nsample)
  Beta1 <- rep(true_beta1(t), nsample)
  Beta2 <- rep(true_beta2(t), nsample)
  
  
  
  eta3 <- Beta0+Beta1 * X1 + Beta2* X2 +(b6+b7*time)#*Gammat

  Y3 <- rnorm(nmesy, eta3, sqrt(sigma))
  #########################################
  
  # linear predictor
  # ضرایب تابعی واقعی
  true_beta0 <- function(t) {
    1-1*t 
    # (t-followup/2)
  }
  
  true_beta1 <- function(t) {
    (t-followup/2)^2
  }
  
  true_beta2 <- function(t) {
    (t-followup/2)^3
    
    
  }
  
  true_beta0(t)
  true_beta1(t)
  true_beta2(t)
  
  Beta0 <- rep(true_beta0(t), nsample)
  Beta1 <- rep(true_beta1(t), nsample)
  Beta2 <- rep(true_beta2(t), nsample)
  eta4 <- Beta0+Beta1 * X1 + Beta2* X2 +(b8+b9*time+b10*time^2)
  Y4 <- rnorm(nmesy, eta4, sqrt(sigma))##### Polynomial 
  
  # linear predictor
  # ضرایب تابعی واقعی
  true_beta0 <- function(t) {
    ;1-3*t
  }
  
  true_beta1 <- function(t) {
    ;1
  }
  
  true_beta2 <- function(t) {
    ;1
  }
  
  
  true_beta0(t)
  true_beta1(t)
  true_beta2(t)
  
  Beta0 <- rep(true_beta0(t), nsample)
  Beta1 <- rep(true_beta1(t), length(t)*nsample)
  Beta2 <- rep(true_beta2(t), length(t)*nsample)
  
  eta5 <- Beta0+Beta1 * X1 + Beta2* X2 +(b11+b12*time)
  Y5 <- rnorm(nmesy, eta5, sqrt(sigma))##### Liear
  ###################################################
  true_beta0 <- function(t) {
    1-t^2
  }
  
  true_beta1 <- function(t) {
    0.5+6*t-4*t^2-6*t^3  }
  
  true_beta2 <- function(t) {
    1+3*t-2*t^2-6*t^3
  }
  
  
  true_beta0(t)
  true_beta1(t)
  true_beta2(t)
  
  Beta0 <- rep(true_beta0(t), nsample)
  Beta1 <- rep(true_beta1(t), nsample)
  Beta2 <- rep(true_beta2(t), nsample)
  
  eta6 <- Beta0+Beta1 * X1 + Beta2* X2 +(b13+b14*time)#*Gammat
  Y6 <- rnorm(nmesy, eta6, sqrt(sigma))
  #########################################
  
  # Permutation algorithm to generate survival times
  Data_permu <- permalgorithm(nsample, nmesindiv,
                              Xmat = cbind(eta1,eta2,eta3,eta4,eta5,eta6,W1,W2),
                              eventRandom = round(rexp(nsample, .0001) + 1, 0),
                              censorRandom = runif(nsample, followup, nmesindiv),
                              XmatNames = c("eta1","eta2","eta3","eta4","eta5","eta6","W1","W2"),
                              betas = c(gamma,alpha)
  )
  
  
  # extract last line for each id (= death/censoring time)
  Data_permu2 <- Data_permu[c(which(diff(Data_permu[, "Id"]) == 1), dim(Data_permu)[1]), c("Id", "Event", "Stop")]
  Data_permu2$survtime <- t[Data_permu2$Stop + 1] # survtime
  surv.data <- Data_permu2[, c("Id", "survtime", "Event")]
  surv.data$w1 <- w1
  surv.data$w2 <- w2
  Data_permu$time <- t[Data_permu$Start + 1] # measurements times of the biomarker
  Data_permu$Uid <- paste(Data_permu$Id, Data_permu$time) # unique identifier to match covariates and observed biomarker values
  long.data3 <- merge(Data_permu[, c("Uid", "Id", "time")], cbind("Uid" = paste(idY, time), X1, X2, Y1,Y2,Y3,Y4,Y5,Y6), by = c("Uid"))
  long.data <- sapply(long.data3[long.data3$time %in% timesLongi, -1], as.numeric)
  long.data <- as.data.frame(long.data[order(long.data[, "Id"], long.data[, "time"]), ])
  colnames(long.data) <- c("id", "obstime", "x1", "x2", "Y1", "Y2", "Y3", "Y4", "Y5", "Y6")
  colnames(surv.data) <- c("id", "survtime", "death", "w1", "w2")
  
  
  head(long.data)
  head(surv.data)
  mean(table(long.data$id))
  median(table(long.data$id))
  min(table(long.data$id))
  max(table(long.data$id))
  long.data$id=as.character(long.data$id)
  
  long.data1 <- cbind("Y1",long.data[, c("id", "obstime", "x1", "x2", "Y1")])
  long.data2 <- cbind("Y2",long.data[, c("id", "obstime", "x1", "x2", "Y2")])
  long.data3 <- cbind("Y3",long.data[, c("id", "obstime", "x1", "x2", "Y3")])
  long.data4 <- cbind("Y4",long.data[, c("id", "obstime", "x1", "x2", "Y4")])
  long.data5 <- cbind("Y5",long.data[, c("id", "obstime", "x1", "x2", "Y5")])
  long.data6 <- cbind("Y6",long.data[, c("id", "obstime", "x1", "x2", "Y6")])
  
  colnames(long.data1)=colnames(long.data2)= 
    colnames(long.data3)=colnames(long.data4)= 
    colnames(long.data5)=colnames(long.data6)= c("dim","id", "obstime", "x1", "x2", "Y")

  long_data=rbind(long.data1,long.data2,long.data3,long.data4,
  long.data5,long.data6)
   
  
  colnames(long_data)=c("Variable","id","obstime","x1","x2","Value")
  head(long_data)
  
  names(long_data)
  
  ######## ######## ######## ######## ######## ######## ######## ######## ######## 
  ######## ######## ######## ######## ######## ######## ######## ######## ######## 
  ######## ######## ######## ######## ######## ######## ######## ######## ######## 
  ######## ######## ######## ######## ######## ######## ######## ######## ######## 
  ######## ######## ######## ######## ######## ######## ######## ######## ######## 
  ######## ######## ######## ######## ######## ######## ######## ######## ######## 
  
  
  
  
  library(data.table)
  
  
  longitudinal_long <- as.data.table(long_data)
  head(longitudinal_long)
  
 
  # Create Spaghetti Plot
  ggplot(longitudinal_long, aes(x = obstime, y = Value, group = id, color = Variable)) +
    geom_line(alpha = 0.4, size = 0.7) +  # Light individual trajectories
    facet_wrap(~Variable, scales = "free_y") +  # Separate plots for Y1, Y2, Y3
    theme_minimal() +
    labs(title = "Spaghetti Plot of Longitudinal Data",
         subtitle = "Each line represents an individual's trajectory",
         x = "Time Point",
         y = "Observed Value",
         color = "Variable") +
    theme(text = element_text(size = 14),
          legend.position = "none")  # Optional: Remove legend for cleaner visualization
  
  # Step 1: Convert Each Response Variable into Functional Data
  time_points <- unique(longitudinal_long$obstime)
  
  
  
  library(data.table)
  
  # Ensure your data is a data.table
  setDT(longitudinal_long)
  
  # Get unique combinations of id and Variable
  id_var_combos <- unique(longitudinal_long[, .(id, Variable, x1, x2)])
  
  # Get all unique time points
  time_points <- unique(longitudinal_long$obstime)
  
  # Create full grid of all combinations
  full_grid <- CJ(id = id_var_combos$id,
                  Variable = id_var_combos$Variable,
                  obstime = time_points,
                  unique = TRUE)
  
  # Merge back x1 and x2 (assuming these don't change with obstime for each id-variable combo)
  full_data <- merge(full_grid, id_var_combos, by = c("id", "Variable"), all.x = TRUE)
  
  # Merge with original data to get Value
  extended_longitudinal <- merge(full_data, longitudinal_long,
                                 by = c("id", "Variable", "obstime", "x1", "x2"),
                                 all.x = TRUE)
  
  # Now extended_longitudinal contains NAs for missing Value observations
  head(extended_longitudinal)
  
  
  
  # Create Spaghetti Plot
  ggplot(extended_longitudinal, aes(x = obstime, y = Value, group = id, color = Variable)) +
    geom_line(alpha = 0.4, size = 0.7) +  # Light individual trajectories
    facet_wrap(~Variable, scales = "free_y") +  # Separate plots for Y1, Y2, Y3
    theme_minimal() +
    labs(title = "Spaghetti Plot of Longitudinal Data",
         subtitle = "Each line represents an individual's trajectory, with missing values handled",
         x = "Time Point",
         y = "Observed Value",
         color = "Variable") +
    theme(text = element_text(size = 14),
          legend.position = "none")  # Optional: Remove legend for cleaner visualization
  
  
  
  
  ids <- unique(extended_longitudinal[obstime==0.05 & is.na(Value)]$id)
  ids
  long_data=extended_longitudinal[!(id%in%ids)]

  # Function to convert data into funData object (handles missing values)
  convert_to_funData <- function(df, variable) {
    df_var <- df %>% filter(Variable == variable)
    data_matrix <- matrix(df_var$Value, nrow = length(unique(df_var$id)), byrow = TRUE)
    funData(argvals = time_points, X = data_matrix)
  }
  
  fd_Y1 <- convert_to_funData(long_data, "Y1")
  fd_Y2 <- convert_to_funData(long_data, "Y2")
  fd_Y3 <- convert_to_funData(long_data, "Y3")
  fd_Y4 <- convert_to_funData(long_data, "Y4")
  fd_Y5 <- convert_to_funData(long_data, "Y5")
  fd_Y6 <- convert_to_funData(long_data, "Y6")
  
  # Step 2: Create Multi-Functional Data Object
  mfd <- multiFunData(list(fd_Y1, fd_Y2, fd_Y3,fd_Y4, fd_Y5, fd_Y6))
  plot(mfd)
  
  # Step 3: Perform Multivariate Functional PCA (MFPCA)
  mfpca_results <- MFPCA(mfd, M = 5,
                         uniExpansions = list(list(type = "uFPCA"),list(type = "uFPCA"),
                                              list(type = "uFPCA"),list(type = "uFPCA"),
                                              list(type = "uFPCA"),
                                              list(type = "uFPCA")))  # Extract 3 principal components
  
  
  pred <- predict(mfpca_results) # default reconstructs data used for the MFPCA fit
  
  # Step 5: Convert Predictions Back to Data Frame
  predicted_data <- data.frame(
    Time = rep(time_points, times = nsample-length(ids)),
    id = rep((1:nsample)[-as.numeric(ids)], each = length(time_points)),
    Y1_pred = as.vector(t(pred[[1]]@X)),                                          #Adding "t(.)
    Y2_pred = as.vector(t(pred[[2]]@X)),                                          #Adding "t(.)
    Y3_pred = as.vector(t(pred[[3]]@X)),                                          #Adding "t(.)
    Y4_pred = as.vector(t(pred[[4]]@X)),                                          #Adding "t(.)
    Y5_pred = as.vector(t(pred[[5]]@X)),                                          #Adding "t(.)
    Y6_pred = as.vector(t(pred[[6]]@X))                                           #Adding "t(.)
  )
  
  
  
  # Step 6: Visualize Predictions vs Observed Data
  predicted_long <- predicted_data %>%
    pivot_longer(cols = starts_with("Y"), names_to = "Variable", values_to = "Predicted_Value")
  
  
  
  
  head(predicted_long)
  # Create Spaghetti Plot
  ggplot(predicted_long, aes(x = Time, y = Predicted_Value, group = id, color = Variable)) +
    geom_line(alpha = 0.4, size = 0.7) +  # Light individual trajectories
    geom_point(alpha = 0.5, size = 0.6) +  # Show observed points
    facet_wrap(~Variable, scales = "free_y") +  # Separate plots for Y1_pred, Y2_pred, Y3_pred
    theme_minimal() +
    labs(title = "Spaghetti Plot of Predicted Responses Using MFPCA",
         subtitle = "Dashed lines represent individual predicted trajectories",
         x = "Time Point",
         y = "Predicted Value",
         color = "Variable") +
    theme(text = element_text(size = 14),
          legend.position = "none")  # Optional: remove legend for cleaner visualization
  
  
  