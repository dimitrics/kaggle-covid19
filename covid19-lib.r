
glm.transform <- function(gfit1, alpha = 0.05, logodds_str = "(exp(beta_i)-1)*100", ci = F, stdout = T){

  pvals = coef(summary(gfit1))[,4] #get p-vaues
  
  results = data.frame( variable = as.character(names(pvals)) )
  results[, logodds_str] = NA
  sig_coefs = coefficients(gfit1)[pvals <= alpha]
  results[, logodds_str][pvals <= alpha] = round( (exp(sig_coefs)-1)*100, 3)
  results$p_value = round(pvals, 3)
  results$beta_i = coefficients(gfit1)
  
  if(ci == T){
    ci_df = confint(gfit1, level = 1-alpha, )
    results = cbind(results, ci_df)
  }
  
  if(stdout){
    print(head(results))
  } 
  
  return(results)
} #glm.transform(gfit)


glm.test <- function(gfit1, alpha = 0.05, stdout = T){
  pvals = coef(summary(gfit1))[,4] #get p-vaues
  
  #GLM method: 1 - pchisq(scaled_deviance, df=2)
  #df: no. of observations - no. of model parameters (betas)
  
  #Likelihood ratio test (Omnibus - must reject H0: Mi = M0):
  pval_m0 = 1 - pchisq(gfit1$null.deviance - gfit1$deviance, gfit1$df.null - gfit1$df.residual) #get test pvalue (must reject)
  
  #Deviance test: (must not reject H0: Mi = Ms)    #https://thestatsgeek.com/2014/04/26/deviance-goodness-of-fit-test-for-poisson-regression/
  pval_ms = 1 - pchisq(gfit1$deviance, df = gfit1$df.residual)
  
  if(stdout){
    cat("\n---- Coefficient p-values ----\n")
    print(pvals)
    cat("\n Significance:\n")
    print(pvals <= alpha)
    cat("------ Goodness of fit ---------\n")
    cat("Likelihood Ratio Test pvalue (must reject):", pval_m0, "\n")
    cat("Deviance Test pvalue: (must not reject)", pval_ms, "\n")
    cat("--------------------------------\n")
  }  
} #glm.test(gfit)


glm.plot <- function(gfit1, alpha = 0.05, logodds_str = "(exp(beta_i)-1)*100", plot = T, stdout = T){
  pvals = coef(summary(gfit1))[,4] #get p-vaues
  
  #Get significant coefficients and corresponding log-odds:
  results = glm.transform(gfit1, stdout = stdout)
  
  #using the formula: (odds_ratio - 1)*100
  if(plot){
    barplot(results[,2], names.arg = as.character(results$variable), 
            ylab = logodds_str, ylim=range(min(results[,2], na.rm = T), max(results[,2], na.rm = T)), cex.names=0.9, 
            col = "grey", main = "Percentage change in response probability")
  }
  #sum(results[,logodds_str][!is.na(results[,logodds_str])]>0)
} #glm.plot(gfit)



myprecision <- function(y_pred, y_actual){
  precision = sum((y_pred==1) & (y_actual==1), na.rm = T) / sum(y_pred==1, na.rm = T) #precision
  if(length(y_pred) != length(y_actual)){
    warning("myprecision(): vectors to compare are not of the same length.")
  }
  return(precision)
}

myrecall <- function(y_pred, y_actual){
  recall = sum((y_pred==1) & (y_actual==1), na.rm = T) / sum(y_actual==1, na.rm = T) #recall
  if(length(y_pred) != length(y_actual)){
    warning("myrecall(): vectors to compare are not of the same length.")
  }
  return(recall)
} 


bagged_glm <- function(dat, gfit, type = c("glm", "lm"), family = c("binomial", "poisson", "gaussian"), class_balancing = F, alpha = 0.05, bag_pct = 0.70, bag_size = NULL, n_bags = 1000, ocv = F){
  #-------------------------
  # Bagged inference:
  
  #Initialise bootstrap parameters:
  if(is.null(bag_size)){
    bag_size = round(bag_pct*nrow(dat))   #override bag_pct
  }
  gfit_list = list()
  res_list = list()
  coef_bags_df = as.data.frame(matrix(0, n_bags, length(coefficients(gfit))))
  names(coef_bags_df) = names(coefficients(gfit))
  if(ocv){
    precision_bag_vec = numeric(n_bags)
    recall_bag_vec = numeric(n_bags)
  }
  
  #Class separation (for the case of class rebalancing)
  response_name = strsplit(as.character(gfit$formula), '~')[[2]]
  inds_1_vec = dat$y==1    #paste("dat$", response_name, sep = "")
  dat_1s = dat[inds_1_vec, ]
  dat_0s = dat[!inds_1_vec, ]
  
  #Bagging loop:
  for(i in 1:n_bags){
    
    #Create bag of training data:
    if(class_balancing == F){ #Unbalanced classes (random subsample, assumes representative training sample)
      ind_bag = sample(x = 1:nrow(dat_0s), size = bag_size, replace = F)
      dat_bag = dat[ind_bag, ]
      if(ocv){dat_out_of_bag = dat[setdiff(1:nrow(dat), ind_bag), ]}
    }else{ #Balance classes (subsample 50-50, assumes class 1 is undersampled)
      # Draw as many 0s as there are 1s:
      ind_bag = sample(x = 1:nrow(dat_0s), size = sum(inds_1_vec), replace = F)
      dat_bag = rbind(dat_0s[ind_bag, ], dat_1s)
      if(ocv){dat_out_of_bag = dat[setdiff(1:nrow(dat), c(ind_bag, which(dat$y==1))), ]}
    }
    
    #Train model on current bag data:
    gfit_list[[i]] = glm(gfit$formula, data = dat_bag, binomial)
    
    #Predict on out-of-bag data: (out-of-bag cross validation)
    if(ocv){
      y_pred = round(predict(gfit_list[[i]], dat_out_of_bag, type = "response"))
      precision_bag_vec[i] = myprecision(y_pred = y_pred, y_actual = dat_out_of_bag$y) #sum((y_pred==1) & (dat_out_of_bag$y==1), na.rm = T) / sum(y_pred==1, na.rm = T)
      recall_bag_vec[i] = myrecall(y_pred = y_pred, y_actual = dat_out_of_bag$y) #sum((y_pred==1) & (dat_out_of_bag$y==1), na.rm = T) / sum(dat_out_of_bag$y==1, na.rm = T)
    }
    
    #Store model diagnostics:
    res_bag = karma.glm_diagnose(gfit_list[[i]], alpha = alpha, stdout = F, plot = F)
    res_list[[i]] = res_bag
    
    #Store model coefficients:
    coef_bags_df[i,] = res_bag$beta_i
  }
  
  bagged_results = list()
  bagged_results$coef_bags_df = coef_bags_df
  if(ocv){
    bagged_results$precision_bag_vec = precision_bag_vec
    bagged_results$recall_bag_vec = recall_bag_vec
  }else{
    bagged_results$precision_bag_vec = NA
    bagged_results$recall_bag_vec = NA    
  }
  class(bagged_results) = "boot.glm"
  return(bagged_results)
}


bagged_cv <- function(dat, formula = "", method = c("glm", "rf"), class_balancing = F, bag_pct = 0.70, bag_size = NULL, n_bags = 1000, stdout = T){
  #-------------------------
  # Bagged cross validation for binary classification (response name must be 'y'):
  
  #Initialise bootstrap parameters:
  if(is.null(bag_size)){
    bag_size = round(bag_pct*nrow(dat))   #override bag_pct
  }
  gfit_list = list()
  precision_bag_vec = numeric(n_bags)
  recall_bag_vec = numeric(n_bags)
  
  #Class separation (for the case of class rebalancing)
  #response_name = strsplit(as.character(formula), '~')[[2]]
  inds_1_vec = dat$y==1    #eval(parse(text=paste("dat$", response_name, sep = "")))
  dat_1s = dat[inds_1_vec, ]
  dat_0s = dat[!inds_1_vec, ]
  
  #Bagging loop:
  for(i in 1:n_bags){
    
    #Create bag of training data:
    if(class_balancing == F){ #Unbalanced classes (random subsample, assumes representative training sample)
      ind_bag = sample(x = 1:nrow(dat_0s), size = bag_size, replace = F)
      dat_bag = dat[ind_bag, ]
      dat_out_of_bag = dat[setdiff(1:nrow(dat), ind_bag), ]
    }else{ #Balance classes (subsample 50-50, assumes class 1 is undersampled)
      # Draw as many 0s as there are 1s:
      ind_bag = sample(x = 1:nrow(dat_0s), size = sum(inds_1_vec), replace = F)
      dat_bag = rbind(dat_0s[ind_bag, ], dat_1s)
      dat_out_of_bag = dat[setdiff(1:nrow(dat), c(ind_bag, which(dat$y==1))), ]
    }
    
    #Train model on current bag data:
    if(method == "glm"){
      gfit_list[[i]] = glm(formula, data = dat_bag, binomial)
    }else if(method == "rf"){
      gfit_list[[i]] = randomForest(eval(parse(text=formula)), data = dat_bag)
    }
    
    #Predict on out-of-bag data: (out-of-bag cross validation)
    if(method == "glm"){
      y_pred = round(predict(gfit_list[[i]], dat_out_of_bag, type = "response"))
    }else if(method == "rf"){
      y_pred = predict(gfit_list[[i]], dat_out_of_bag)
    }
    precision_bag_vec[i] = myprecision(y_pred = y_pred, y_actual = dat_out_of_bag$y) #sum((y_pred==1) & (dat_out_of_bag$y==1), na.rm = T) / sum(y_pred==1, na.rm = T)
    recall_bag_vec[i] = myrecall(y_pred = y_pred, y_actual = dat_out_of_bag$y) #sum((y_pred==1) & (dat_out_of_bag$y==1), na.rm = T) / sum(dat_out_of_bag$y==1, na.rm = T)
  }
  
  bagged_cv = list()
  bagged_cv$precision_bag_vec = precision_bag_vec
  bagged_cv$recall_bag_vec = recall_bag_vec
  
  if(stdout == T){
    cat("Bagged precision: ", mean(bagged_cv$precision_bag_vec, na.rm = T), "\n")
    cat("Bagged recall:    ", mean(bagged_cv$recall_bag_vec, na.rm = T), "\n")
  }

  class(bagged_cv) = "boot.cv"
  return(bagged_cv)
}


boot.dist <- function(Y, alpha){
  # Create bootstramp samples, bootstrap statistics, and a bootstrap distribution:
  Y_boot = numeric(length(Y)) #bootstrap sample
  theta_boot_dist = numeric(length(Y_boot)) #bootstrap distribution of bootstrap statisic (mean)
  
  #Generate Bootstrap distribution:
  for(i in 1:length(theta_boot_dist)){
    # Create bootstrap sample from Y: (sample random observations with replacement)
    for(j in 1:length(Y_boot)){    #NB: could also be done in one go with: Y_boot = sample(Y, 1000, replace = T)
      Y_boot[j] = sample(Y, 1)  #sample a random observation
    }
    theta_boot_dist[i] = quantile(Y_boot, probs = 1-alpha) #calculate bootstrap statistic
    #theta_boot_dist[i] = quantile(sort(Y_boot), probs = 1-alpha/2) - quantile(sort(Y_boot), probs = alpha/2)  #Bayesian credible interval method
  }
  return(theta_boot_dist) #returns bootstrap disribution of theta
} #usage: boot_dist = boot.dist(Y, 0.05)


bootstrap_ci <- function(data, bag_n = 10000, boot_n = 10000, alpha = 0.05){
  n = length(data)
  x_bar = mean(data)
  boot_sample = replicate( boot_n, mean(sample(data, n, rep = T)) - x_bar )
  #boot_sample = replicate( boot_n, mean(sample(data, n, rep = T)))
  boot_qt = quantile(boot_sample, c(1-alpha/2, alpha/2))
  boot_ci = x_bar - boot_qt
  
  boot_dist = replicate( boot_n, sample(data, n, rep = T))
  
  boot_obj = list()
  class(boot_obj) = "bootstrap_ci"
  boot_obj$boot_ci = boot_ci
  boot_obj$boot_dist = boot_dist
  return(boot_obj)
}


split_data <- function(dat0, train_pct = 0.70){
  train_n <- floor(train_pct * nrow(dat0))
  train_ind <- sample(seq_len(nrow(dat0)), size = train_n)
  dat_train <- dat0[train_ind, ]
  dat_test <- dat0[-train_ind, ]
  
  split_obj = list()
  split_obj$train_ind = train_ind
  split_obj$dat_train = dat_train
  split_obj$dat_test = dat_test
  class(split_obj) = "split_obj"
  return(split_obj)
} #usage: split_data(dat, train_pct = 0.70)


validate_model <- function(model, dat_train, dat_test, response = c("class", "score"), threshold = 0.5, stdout = T){
  #(classification scenario)
  #Training set:
  y_train = dat_train$y
  y_hat_train = predict(model)   #training prediction
  precision_train = myprecision(y_pred = y_hat_train, y_actual = y_train)    #training error
  recall_train = myrecall(y_pred = y_hat_train, y_actual = y_train)  
  
  #Test set:
  y_test = dat_test$y
  y_hat_test = predict(model, dat_test)    #test prediction
  if(response == "score"){
    y_hat_test[y_hat_test > threshold] = 1
    y_hat_test[y_hat_test <= threshold] = 0
  }
  precision_test = myprecision(y_pred = y_hat_test, y_actual = y_test)
  recall_test = myrecall(y_pred = y_hat_test, y_actual = y_test)
  
  if(stdout == T){
    cat("test precision: ", precision_test, "\n")
    cat("test recall:    ", recall_test, "\n")
  }
  
  valid_obj = list()
  valid_obj$precision_train = precision_train
  valid_obj$precision_test = precision_test
  valid_obj$recall_train = recall_train
  valid_obj$recall_test = recall_test
  class(valid_obj) = "valid_obj"
  return(valid_obj)
}

#cv.kfold()





