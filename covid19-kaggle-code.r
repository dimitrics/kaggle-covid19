#Installing missing packages:
#install.packages('mice')
library(mice)
library(VIM)
library(jtools)
library(boot)
library(randomForest)
library(ggplot2)
library(gridExtra)
options(repr.plot.width = 14, repr.plot.height = 8)

#Load data:
dat0 = read.csv('C:/Users/dtziotis002/Documents/PROJECTS/Coronavirus/data/novel-corona-virus-2019-dataset/COVID19_line_list_data.csv', header = T, sep = ',')
names(dat0)


dat1 = dat0[c("location", "country", "gender", "age", "reporting.date", "symptom_onset", "hosp_visit_date", "visiting.Wuhan", "from.Wuhan", "death", "recovered", "symptom")]
print(head(dat1, 5))

#Create response:
dat1$y = 0
dat1$y[dat1$death==1] = 1   #binarize mortality
dat1$y = as.factor(dat1$y)


table(dat1$y)

# Plot missing values:
mice_plot <- aggr(dat1, col=c('green','red'),
                  numbers=TRUE, sortVars=TRUE,
                  labels=names(dat1), cex.axis=.7,
                  gap=3, ylab=c("Missing data","Pattern"))


cat("Percentage of NAs in hosp_visit_date: ", sum(is.na(dat1$hosp_visit_date))/nrow(dat1), "\n")
cat("Percentage of NAs in symptom_onset: ", sum(is.na(dat1$symptom_onset))/nrow(dat1), "\n")
barplot(table(dat1$country))

dat1$from.Wuhan = as.factor(dat1$from.Wuhan)

#Create time variables:
dat1$hosp_visit_date = as.character(dat1$hosp_visit_date)
dat1$symptom_onset = as.character(dat1$symptom_onset)
dat1$hosp_visit_date = as.Date(dat1$hosp_visit_date, format="%m/%d/%y")
dat1$symptom_onset = as.Date(dat1$symptom_onset, format="%m/%d/%y")

#Create onset time in days:
dat1$days_before_hosp = as.numeric(as.Date(dat1$hosp_visit_date, format="%m/%d/%y") - as.Date(dat1$symptom_onset, format="%m/%d/%y"))
dat1$days_before_hosp[dat1$days_before_hosp < 0] = NA

# Loading custom functions used in the analysis:
source('~/PROJECTS/Coronavirus/covid19-lib.r')


#----------------------------------------------------------------------------------
# Model 0: gender + age
#----------------------------------------------------------------------------------

#Train model:
gfit0 = glm(y~gender+age, data = dat1, binomial)

#Summarize and inspect:
summ(gfit0, confint = TRUE, pvals = TRUE, digits = 3)




# Comparing current model to the null model (LTR):
pval_m0 = 1 - pchisq(gfit0$null.deviance - gfit0$deviance, gfit0$df.null - gfit0$df.residual)
round(pval_m0, 3)


# Comparing current model to the saturated model (Deviance test):
pval_ms = 1 - pchisq(gfit0$deviance, df = gfit0$df.residual)
round(pval_ms, 3)


#Asymptotic distribution of model coefficients (model 0)
plot_summs(gfit0, scale = TRUE, plot.distributions = TRUE)


#Make prediction: #https://stats.stackexchange.com/questions/2335/output-of-logistic-model-in-r
x1 = data.frame(gender="male", age=1:90)
plot(predict(gfit0, x1, type = "response"))   #same as: exp(predict(gfit2,x1))/(1+exp(predict(gfit2,x1))) because: y_hat = exp( log(Odds)) ) = exp( log(p/(1+p)) ) <- transformation done by type="response"
x2 = data.frame(gender="female", age=1:90)
points(predict(gfit0, x2, type = "response"), col="red")


#Create imputed dataset using 'mice':
dat1_obj <- mice(dat1[c("gender", "age", "days_before_hosp", "from.Wuhan", "country")], m=5, maxit = 50, method = 'pmm', seed = 400)
dat1_imp = complete(dat1_obj)
dat1_imp$y = dat1$y


#Train models on imputed data:
gfit1a = glm(y~gender+age, data = dat1_imp, binomial)
gfit1b = glm(y~gender+age+days_before_hosp, data = dat1_imp, binomial)
gfit1c = glm(y~gender+age+days_before_hosp+from.Wuhan, data = dat1_imp, binomial)

#Compare fits:
export_summs(gfit1a, gfit1b, gfit1c, scale = F, error_format = "[{conf.low}, {conf.high}]", digits = 3, model.names = c("model 1a", "model 1b", "model 1c"))
plot_summs(gfit1a, gfit1b, gfit1c, scale = TRUE, plot.distributions = F, inner_ci_level = .95, model.names = c("model 1a", "model 1b", "model 1c"))


# Create reduced dataset for model 2 (filter out rows with NA values):
not_na_mask = !is.na(dat1$y) & !is.na(dat1$days_before_hosp) & !is.na(dat1$age) & !is.na(dat1$gender) & dat1$days_before_hosp >= 0
dat2 = dat1[not_na_mask, ]
table(dat2$y)


#----------------------------------------------------------------------------------
# Model 2: gender + age + time to hospitalization (on reduced dataset)
#----------------------------------------------------------------------------------

#Train model:
gfit2a = glm(y~gender+age, data = dat2, binomial)
gfit2b = glm(y~gender+age+days_before_hosp, data = dat2, binomial)
gfit2c = glm(y~gender+age+days_before_hosp+from.Wuhan, data = dat2, binomial)

#Merge model summaries:
export_summs(gfit2a, gfit2b, gfit2c, scale = F, error_format = "[{conf.low}, {conf.high}]", digits = 3,  model.names = c("model 2a", "model 2b", "model 2c"))

#Compare asymptotic distributions of coefficients:
plot_summs(gfit2a, gfit2b, gfit2c, scale = TRUE, plot.distributions = F, inner_ci_level = .95, model.names = c("model 2a", "model 2b", "model 2c"))

#Summarize and inspect:
#summ(gfit2_3, confint = TRUE, pvals = TRUE, digits = 3)


#Merge summaries of two models (imputed vs. reduced)
export_summs(gfit1c, gfit2c, scale = TRUE, error_format = "[{conf.low}, {conf.high}]", digits = 3, model.names = c("model 1c", "model 2c"))

#Compare asymptotic distributions of coefficients:
p1 = plot_summs(gfit1c, gfit2c, scale = TRUE, plot.distributions = T, model.names = c("model 1c", "model 2c"))
p2 = plot_summs(gfit1c, gfit2c, scale = TRUE, plot.distributions = F, inner_ci_level = .95, model.names = c("model 1c", "model 2c"))
grid.arrange(p1, p2, ncol=2)


#--------------------------------------------------------------------
#Inspecting proportion of mortality rate in Wuhan natives:
#--------------------------------------------------------------------
#Create dataframe:
tmp = data.frame(rbind(
  table(dat1_imp$y[dat1_imp$from.Wuhan==1]),
  table(dat1_imp$y[dat1_imp$from.Wuhan==0])
))
names(tmp) = c("total", "deaths")
tmp$death_rate = round(tmp$deaths/tmp$total, 3)
tmp$from_wuhan = as.factor(c(1,0))

#Compare proportions (deaths per cases between groups):
se1 = sqrt(tmp$death_rate[1]*(1-tmp$death_rate[1])/tmp$total[1])  #standard errors of proportions
se2 = sqrt(tmp$death_rate[2]*(1-tmp$death_rate[2])/tmp$total[2])
tmp$prop_se = round(c(se1, se2), 4)
print(tmp)
print(prop.test(x = tmp$deaths, n = tmp$total, alternative = "greater"))
# Barplot of proportions:
ggplot(tmp, aes( y=death_rate, x=from_wuhan)) + 
  geom_bar(position="dodge", stat="identity", width=0.4, color="black", fill="cyan", alpha=.2) + 
  geom_errorbar(aes(ymin=death_rate - prop_se, ymax=death_rate + prop_se), width=.1, position=position_dodge(.9))


#Compare average number of days gone by before hospitalization between both groups:
d1 = dat1_imp$days_before_hosp[dat1_imp$from.Wuhan==1]
d2 = dat1_imp$days_before_hosp[dat1_imp$from.Wuhan==0]
sem1 = t.test(d1)$stderr
sem2 = t.test(d1)$stderr
tmp$avg_days = c(mean(d1), mean(d2))
tmp$mean_se = c(sem1, sem2)
print(tmp)
t.test(d1, d2, alternative = "greater")
#Barplot:
b1 = ggplot(tmp, aes( y=avg_days, x=from_wuhan, fill=from_wuhan)) + 
  geom_bar(position="dodge", stat="identity", width = .4, alpha=1) + 
  geom_errorbar(aes(ymin=avg_days - mean_se, ymax=avg_days + mean_se), width=.1, position=position_dodge(.9)) 
#Boxplot:
df = data.frame(days = c(d1, d2), from_wuhan = as.factor(c(numeric(length(d1))+1, numeric(length(d2)) )) )
b2 = ggplot(df, aes( y=days, x=from_wuhan, fill=from_wuhan)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=T) +
  stat_summary(fun=mean, geom="point", shape=23, size=4)
grid.arrange(b1, b2, ncol=2)



dat = dat1_imp

dat$gdp = 0
#levels(dat$country)
gdp_per_country = c(2182, 16091, 54799, 55171, 51991, 50904, 5004, 52144, 
                    20984, 29207, 14800, 49548, 48640, 55306, 66527, 9027,
                    17832, 40337, 41582, 46827, 67891, 15599, 34567, 3550,
                    10094, 30820, 105689, 46452, 43007, 14509, 55989, 67558,
                    57214, 21361, 70441, 48169, 67426, 8677)
for(i in 1:length(gdp_per_country)){
  country = levels(dat$country)[i]
  country_inds = dat$country == country
  dat$gdp[country_inds] = gdp_per_country[i]
}



#----------------------------------------------------------------
# Model 4: Adding GDP (PPP)
#----------------------------------------------------------------

#Fit model with GDP:
gfit1d = glm(y~gender+age+days_before_hosp+from.Wuhan+gdp, data = dat, binomial)

#Compare model with GDP with one without GDP:

#Merge model summaries:
export_summs(gfit1c, gfit1d, scale = F, error_format = "[{conf.low}, {conf.high}]", model.names = c("without GDP", "with GDP"))

#Compare symtotic distributions:
f1 = plot_summs(gfit1c, gfit1d, scale = TRUE, plot.distributions = TRUE, model.names = c("without GDP", "with GDP"))
f2 = plot_summs(gfit1c, gfit2c, scale = TRUE, plot.distributions = F, inner_ci_level = .95, model.names = c("without GDP", "with GDP"))
grid.arrange(f1, f2, ncol=2)

#Final model summary:
summ(gfit1d, scale = F, plot.distributions = TRUE, inner_ci_level = .9, digits = 3)
gfit = gfit1d   #rename final model


#Mean GDP per group (from Wuhan):
d3 = dat$gdp[dat$from.Wuhan==1]
d4 = dat$gdp[dat$from.Wuhan==0]
t.test(d3, d4)
sem3 = t.test(d3)$stderr
sem4 = t.test(d4)$stderr
tmp$avg_gdp = c(mean(d3), mean(d4))
tmp$mean_se_gdp = c(sem3, sem4)

#Barplot:
f3 = ggplot(tmp, aes( y=avg_gdp, x=from_wuhan, fill=from_wuhan)) + geom_bar(position="dodge", stat="identity", width = .5) + geom_errorbar(aes(ymin=avg_gdp-mean_se_gdp, ymax=avg_gdp+mean_se_gdp), width=.1, position=position_dodge(.9)) 
#Boxplot:
df = data.frame(days = c(d1, d2), gdp = c(d3, d4), from_wuhan = as.factor(c(numeric(length(d1))+1, numeric(length(d2)) )) )
f4 = ggplot(df, aes( y=gdp, x=from_wuhan, fill=from_wuhan)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) + 
  stat_summary(fun=mean, geom="point", shape=23, size=4)
grid.arrange(f3, f4, ncol=2)


#----------------------------------------------------------------------------------
# Coefficient stability (final model):
#----------------------------------------------------------------------------------
bagged_results = bagged_glm(dat = dat, gfit = gfit, alpha = 0.05, class_balancing = F, bag_pct = 0.80, ocv = T)
coef_bags_df = bagged_results$coef_bags_df


#Bootstrap means of transformed coefficients:
cat(
  "\n Bootstrap means of back-transformed coefficients:\n",
  "Intercept: ", mean(coef_bags_df$`(Intercept)`), 
  "Gender: ", mean(coef_bags_df$gendermale),    #ta odds ths pithanothtas thanatou auksanonta kata X% stous antres ap'oti stis gynaikes
  "Age: ", mean(coef_bags_df$age),     #ta odds ths pithanothtas thanatou auksanonta kata X% otan h hlikia auksanetai kata mia monada
  "Days before hospitalization: ", mean(coef_bags_df$days_before_hosp),
  "From Wuhan: ", mean(coef_bags_df$from.Wuhan1),
  "\n"
)
#Compare with asymptotic values: (merge with table above)
coefficients(gfit)
confint(gfit)

# #Plot histograms of transformed coefficient empirical distributions:
# h0 = qplot(coef_bags_df$"(Intercept)", geom="histogram", bins = 10, xlab = "Intercept", fill=I("blue"), col=I("black"), alpha=I(.2))
# h1 = qplot(coef_bags_df$gendermale, geom="histogram", bins = 10, xlab = "gender", fill=I("blue"), col=I("black"), alpha=I(.2))
# h2 = qplot(coef_bags_df$age, geom="histogram", bins = 10, xlab = "age", fill=I("blue"), col=I("black"), alpha=I(.2))
# h3 = qplot(coef_bags_df$days_before_hosp, geom="histogram", bins = 10, xlab = "Days before hosp", fill=I("blue"), col=I("black"), alpha=I(.2))
# h4 = qplot(coef_bags_df$from.Wuhan, geom="histogram", bins = 10, xlab = "from Wuhan", fill=I("blue"), col=I("black"), alpha=I(.2))
# h5 = qplot(coef_bags_df$gdp, geom="histogram", bins = 10, xlab = "GDP", fill=I("blue"), col=I("black"), alpha=I(.2))
# h6 = qplot(bagged_results$precision_bag_vec, geom="histogram", bins = 10, xlab = "Precision", fill=I("blue"), col=I("black"), alpha=I(.2))
# h7 = qplot(bagged_results$recall_bag_vec, geom="histogram", bins = 10, xlab = "Recall", fill=I("blue"), col=I("black"), alpha=I(.2))
# grid.arrange(h0, h1, h2, h3, h4, h5, h6, h7, ncol=4, nrow=2)


# Test null distribution directly: H0: m2-m1 = 0
gg_vec = list()
for(i in 1:length(bagged_results$coef_bags_df)){
  d = bagged_results$coef_bags_df[[i]] - coefficients(gfit)[i]
  tt = quantile(d, c(.975, .025))  #tt = t.test(d)
  #hist(d)
  df2 = data.frame(d=d)
  gg_vec[[i]] =  ggplot(df2, aes(x=d)) +
    geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
    geom_density(alpha=.2) + labs(x = names(gfit$coefficients[i])) +
    geom_vline(aes(xintercept=0, color="red"), linetype="dashed") +
    geom_vline(aes(xintercept=tt[1], color="green"), linetype="dashed") +
    geom_vline(aes(xintercept=tt[2], color="green"), linetype="dashed")
}
grid.arrange(gg_vec[[1]], gg_vec[[2]], gg_vec[[3]], gg_vec[[4]], gg_vec[[5]], gg_vec[[6]], ncol=2, nrow=3)


# Test null distribution directly: H0: m2-m1 = 0 [ggplot doesn't plot correctly on Jupyter]
par(mfcol=c(3,2))
for(i in 1:length(bagged_results$coef_bags_df)){
  d = bagged_results$coef_bags_df[[i]] - coefficients(gfit)[i]
  qtt = quantile(d, c(.975, .025))  
  att = t.test(d)
  hist(d, col="lightblue", main = names(coefficients(gfit)[i]), breaks = 30)
  abline(v = 0, col="blue", lwd=3, lty=2)
  abline(v = qtt[1], col="red", lwd=3, lty=2)
  abline(v = qtt[2], col="red", lwd=3, lty=2)
  #abline(v = att[1], col="red", lwd=3, lty=2)
  #abline(v = att[2], col="red", lwd=3, lty=2)
}
dev.off()

bootstrap_ci <- function(data, bag_n = 10000, boot_n = 10000){
  n = length(data)
  x_bar = mean(data)
  boot_sample = replicate( boot_n, mean(sample(data, n, rep = T)) - x_bar )
  boot_qt = quantile(boot_sample, c(.975, .025))
  boot_ci = x_bar - boot_qt
  
  boot_dist = replicate( boot_n, sample(data, n, rep = T))
  
  boot_obj = list()
  class(boot_obj) = "bootstrap_ci"
  boot_obj$boot_ci = boot_ci
  boot_obj$boot_dist = boot_dist
  return(boot_obj)
}


# #Compare overlapping confidence intervals between asymptotic and computational inference (using ggplot)
gg_vec2 = list()
for(i in 1:length(bagged_results$coef_bags_df)){
  #Simulate bootstrap confidence intervals:
  x = bagged_results$coef_bags_df[[i]]
  boot_obj2 = bootstrap_ci(x)
  boot_ci2 = boot_obj2$boot_ci
  
  #Create df:
  d1 = x  #boot_obj2$boot_dist
  d2 = rnorm(n = length(d1), 
             mean = coefficients(gfit)[2], 
             sd = sqrt(diag(vcov(gfit)))[2] )
  groups = factor(rep( c("computational", "asymptotic"), c(length(d1), length(d2)) ))
  df = data.frame(observations = c(d1, d2), inference_method = groups)
  mu = data.frame(inference_method = c("computational", "asymptotic"), avg = c(mean(df$observations[df$inference_method=="computational"]), mean(df$observations[df$inference_method=="asymptotic"])))
  mu$ci_lo = c(boot_ci2[1], confint(gfit)[2,1])
  mu$ci_up = c(boot_ci2[2], confint(gfit)[2,2])
  
  # Plot with color by groups:
  gg_vec2[[i]] = ggplot(df, aes(x=observations, color=inference_method, fill=inference_method)) + 
    geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
    geom_density(alpha=.2) + labs(x = names(gfit$coefficients[i])) +
    #geom_vline(data=mu, aes(xintercept=avg, color=inference_method), linetype="dashed") +
    geom_vline(data=mu, aes(xintercept=ci_lo, color=inference_method), linetype="dashed") +
    geom_vline(data=mu, aes(xintercept=ci_up, color=inference_method), linetype="dashed")
}
grid.arrange(gg_vec2[[1]], gg_vec2[[2]], gg_vec2[[3]], gg_vec2[[4]], gg_vec2[[5]], gg_vec2[[6]], ncol=2, nrow=3)



#Compare overlapping confidence intervals between asymptotic and computational inference (using hist)
par(mfcol=c(3,2))
for(i in 1:length(bagged_results$coef_bags_df)){
  #Simulate bootstrap confidence intervals:
  x = bagged_results$coef_bags_df[[i]]
  boot_obj2 = bootstrap_ci(x)
  boot_ci2 = boot_obj2$boot_ci    
  
  #Create df:
  d1 = x
  d2 = rnorm(n = length(d1), 
             mean = coefficients(gfit)[2], 
             sd = sqrt(diag(vcov(gfit)))[2] )
  groups = factor(rep( c("computational", "asymptotic"), c(length(d1), length(d2)) ))
  df = data.frame(observations = c(d1, d2), inference_method = groups)
  mu = data.frame(inference_method = c("computational", "asymptotic"), avg = c(mean(df$observations[df$inference_method=="computational"]), mean(df$observations[df$inference_method=="asymptotic"])))
  mu$ci_lo = c(boot_ci2[1], confint(gfit)[2,1])
  mu$ci_up = c(boot_ci2[2], confint(gfit)[2,2])
  
  # Plot with color by groups:
  hist(df$observations, col="lightblue", main = names(coefficients(gfit)[i]), breaks = 30)
  #hist(df$observations[df$inference_method=="asymptotic"], col="lightgreen", main = names(coefficients(gfit)[i]), breaks = 30, add = T)
  abline(v = mu$ci_lo[1], col="red", lwd=3, lty=2)
  abline(v = mu$ci_up[1], col="red", lwd=3, lty=2)
  abline(v = mu$ci_lo[2], col="blue", lwd=3, lty=2)
  abline(v = mu$ci_up[2], col="blue", lwd=3, lty=2)
}
dev.off()


# ### Bootstrap confidence intervals:
# samplemean <- function(x, d) {
#   return(mean(x[d]))
# }
# bootobject = boot(data = coef_bags_df$gendermale, statistic = samplemean, R = 1000)  #do the same for variance
# #hist(bootobject)
# boot_ci = boot.ci(bootobject, conf=0.95, type="all")   #<- use this to prove that the coefficient distributions are or aren't different 
# boot_ci$percent


#----------------------------------------------------------------
# Prediction (compare rf and logreg predictions)
#----------------------------------------------------------------

# 70-30 split:
train_obj = split_data(dat, train_pct = 0.70)
dat_train = train_obj$dat_train
dat_test = train_obj$dat_test

# Train random forest for regression:
rf = randomForest(y~gender+age+days_before_hosp+from.Wuhan+gdp, data = dat_train, ntree = 500)
#print(rf)
#varImpPlot(rf)
# Train logreg:
gfitp = glm(y~gender+age+days_before_hosp+from.Wuhan+gdp, data = dat_train, binomial)

#----------------------
# Validate both models:
cat("\nRandom Forest test set accuracy:\n")
rf_valid_obj = validate_model(model = rf, dat_train, dat_test, response = c("class"))
cat("\nBinary logistic regression test set accuracy:\n")
gfit_valid_obj = validate_model(model = gfitp, dat_train, dat_test, response = c("score"))



#-----------------------------
# Bagged CV on both models:  (probably redudant for RF)
rf_bagged_cv = bagged_cv(dat = dat, formula = "y~gender+age+days_before_hosp+from.Wuhan+gdp", method = "rf", class_balancing = F, bag_pct = 0.70, n_bags = 1000, stdout = T)
gfit_bagged_cv = bagged_cv(dat = dat, formula = "y~gender+age+days_before_hosp+from.Wuhan+gdp", method = "glm", class_balancing = F, bag_pct = 0.70, n_bags = 1000, stdout = T)

#Out-of-bag cross validation (for GLM):
cat(
  "\n Bagged classification accuracy (GLM) \n",
  "Precision: ", mean(gfit_bagged_cv$precision_bag_vec, na.rm = T),
  "\n Recall:    ", mean(gfit_bagged_cv$recall_bag_vec, na.rm = T),
  "\n"
)

# #Plot OCV:
# par(mfrow = c(1,2))
# hist(bagged_results$precision_bag_vec)
# hist(bagged_results$recall_bag_vec)
# dev.off()


# #Out-of-bag error (for RF)
# round(1-rf$err.rate[nrow(rf$err.rate),1], 4)   #the i-th element is the (OOB) error rate for all trees up to the i-th.

# #Plot OBE:
# par(mfrow = c(1,2))
# plot(rf$err.rate[,1], type = "l")
# plot(rf)
# dev.off()


# Test null (ggplot): H0: m2-m1 = 0, where m2 is 70-30 RF precision/recall an m1 is bagged GLM precition/recall
gg_vec = list()
dp = gfit_bagged_cv$precision_bag_vec - rf_bagged_cv$precision_bag_vec
dr = gfit_bagged_cv$recall_bag_vec - rf_bagged_cv$recall_bag_vec
tt_p = quantile(dp, c(.975, .025), na.rm = T)  #tt = t.test(d)
tt_r = quantile(dr, c(.975, .025), na.rm = T)  #tt = t.test(d)
#hist(d)
#abline(v = tt[1], col="red", lwd=3, lty=2)
#abline(v = tt[2], col="red", lwd=3, lty=2)
# or with ggplot:
df3 = data.frame(dp=dp, dr=dr)
gg_vec[[1]] =  ggplot(df3, aes(x=dp)) +
  geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
  #gg_vec[[i]] = qplot(d, geom="histogram", bins = 10) +
  geom_density(alpha=.2) +
  geom_vline(aes(xintercept=0, color="red"), linetype="dashed") +
  geom_vline(aes(xintercept=tt[1], color="green"), linetype="dashed") +
  geom_vline(aes(xintercept=tt[2], color="green"), linetype="dashed")
gg_vec[[2]] =  ggplot(df3, aes(x=dr)) +
  geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
  #gg_vec[[i]] = qplot(d, geom="histogram", bins = 10) +
  geom_density(alpha=.2) +
  geom_vline(aes(xintercept=0, color="red"), linetype="dashed") +
  geom_vline(aes(xintercept=tt[1], color="green"), linetype="dashed") +
  geom_vline(aes(xintercept=tt[2], color="green"), linetype="dashed")
grid.arrange(gg_vec[[1]], gg_vec[[2]], ncol=2, nrow=1)



# Test null (hist): H0: m2-m1 = 0, where m2 is 70-30 RF precision/recall an m1 is bagged GLM precition/recall
gg_vec = list()
dp = gfit_bagged_cv$precision_bag_vec - rf_bagged_cv$precision_bag_vec
dr = gfit_bagged_cv$recall_bag_vec - rf_bagged_cv$recall_bag_vec
tt_p = quantile(dp, c(.975, .025), na.rm = T)  #tt = t.test(d)
tt_r = quantile(dr, c(.975, .025), na.rm = T)  #tt = t.test(d)

par(mfcol=c(1,2))
# Plot with color by groups:
hist(df3$dp, col="lightblue", main = "Bagged precision", breaks = 30)
abline(v = 0, col="blue", lwd=3, lty=2)
abline(v = tt_p[1], col="red", lwd=3, lty=2)
abline(v = tt_p[2], col="red", lwd=3, lty=2)
# Plot with color by groups:
hist(df3$dr, col="lightblue", main = "Bagged recall", breaks = 30)
abline(v = 0, col="blue", lwd=3, lty=2)
abline(v = tt_r[1], col="red", lwd=3, lty=2)
abline(v = tt_r[2], col="red", lwd=3, lty=2)
dev.off()


glm.transform(gfit, alpha = .10, logodds_str = "exp(beta)%", ci = T, stdout = F)


glm.plot(gfit, alpha = 0.10, logodds_str = "exp(beta)%", stdout = F)






