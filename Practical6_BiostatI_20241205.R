# Practical 6 : 2025/12/05


#### envir ####
getwd() # know where you are
# setwd(file.path("...", "...")) # setup a working directory
# list.files()


#### pkg I will need ####
library(readxl)
library(writexl)
library(dplyr)
library(gtsummary)
library(car)
library(pROC)

#### 6.1- PG6 ####

# 1) Load the PG6 dataset into R. 
mydata <- read_excel("PG6.xlsx")

# 2) Format the data: Check or correct the formats of the variables. 

str(mydata)
# tibble [31 × 4] (S3: tbl_df/tbl/data.frame)
# $ id_p  : num [1:31] 1 2 3 4 5 6 7 8 9 10 ...
# $ weight: num [1:31] 4.17 5.58 5.18 6.11 4.5 4.61 5.17 4.53 5.33 5.14 ...
# $ group : chr [1:31] "ctrl" "ctrl" "ctrl" "ctrl" ...
# $ site  : chr [1:31] "Site1" "Site2" "Site1" "Site2" ...

# mydata$id_p <- as.character(mydata$id_p) # can be character (just information?)
mydata$group <- as.factor(mydata$group)
levels(mydata$group) # "ctrl" "trt1" "trt2"
# so ctrl is well my reference, and trt1 and 2 will be compared to this ref. 
mydata$site <- as.factor(mydata$site)

# we can ignor id column, not a variable to analyse, just a meta information...
mydata$id_p <- NULL

str(mydata)
# tibble [31 × 3] (S3: tbl_df/tbl/data.frame)
# $ weight: num [1:31] 4.17 5.58 5.18 6.11 4.5 4.61 5.17 4.53 5.33 5.14 ...
# $ group : Factor w/ 3 levels "ctrl","trt1",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ site  : Factor w/ 2 levels "Site1","Site2": 1 2 1 2 1 2 1 2 1 2 ...

# 3) Describe the dataset as a whole and by treatment condition. 
# several possibilities ... 

summary(mydata)
#     weight       group       site   
# Min.   :3.590   ctrl:10   Site1:15  
# 1st Qu.:4.550   trt1:10   Site2:16  
# Median :5.155   trt2:11             
# Mean   :5.073                       
# 3rd Qu.:5.530                       
# Max.   :6.310                       
# NA's   :1                           
table(mydata$group, useNA = "ifany")
# ctrl trt1 trt2 
# 10   10   11 
prop.table(table(mydata$group))
# ctrl      trt1      trt2 
# 0.3225806 0.3225806 0.3548387 

# library(gtsummary)
gtsummary::tbl_summary(
  mydata,
  statistic = list(all_continuous() ~ "{mean} +/-{sd} {median} ({IQR})")
)

# by group
# library(dplyr)
mydata %>% 
  group_by(group) %>% 
  summarise(
    mean = mean(weight, na.rm = TRUE),
    # ...
    n_NA = sum(is.na(weight)), 
    n = length(weight)
  )
# # A tibble: 3 × 4
#   group  mean  n_NA     n
#   <fct> <dbl> <int> <int>
# 1 ctrl   5.03     0    10
# 2 trt1   4.66     0    10
# 3 trt2   5.53     1    11

gtsummary::tbl_summary(
  mydata, 
  by = group, 
  statistic = list(all_continuous() ~ "{mean} +/-{sd} {median} ({IQR})")
)

aggregate(weight ~ group, data = mydata, FUN = length, na.action = identity)
aggregate(weight ~ group, data = mydata, FUN = mean) # by default na.action = na.omit
aggregate(weight ~ group, data = mydata, FUN = sd)
aggregate(
  weight ~ group, data = mydata, FUN = function(x) sum(is.na(x)),
  na.action = identity
)


### save the desc : 

# save it
globaldesc <- as_tibble(
  gtsummary::tbl_summary(
    mydata,
    statistic = list(all_continuous() ~ "{mean} +/-{sd} {median} ({IQR})")
  )
)
groupdesc <-  as_tibble(
  gtsummary::tbl_summary(
    mydata, 
    by = group, 
    statistic = list(all_continuous() ~ "{mean} +/-{sd} {median} ({IQR})")
  )
)

writexl::write_xlsx(list(globaldesc, groupdesc), "gtsummary_GP.xlsx")


# 4) Clean the data: If missing values exist,
#use na.omit() to remove unusable rows. Example of use:

mydata_clean <- na.omit(mydata)  # Removes rows with NA
  
# 5) Visualize the relationships between variables with appropriate graphs.

# png(...)
boxplot(
  weight ~ group, data = mydata_clean,
  main = "Plant weight by treatment group",
  xlab = "Treatment group", ylab = "Weight"
)
means <- aggregate(weight ~ group, data = mydata_clean, FUN = mean)
points(1:3, means$weight, col = "red")
# dev.off()

# Optional: with ggplot2
# library(ggplot2)
# ggplot(mydata_clean, aes(x = group, y = weight)) + geom_boxplot()

# 6) Perform an ANOVA test to check whether there is a significant difference between the means of weight according to the different groups (group).
anova_result <- aov(weight ~ group, data = mydata_clean)
summary(anova_result)
#               Df Sum Sq Mean Sq F value Pr(>F)  
# group        2  3.766  1.8832   4.846 0.0159 *
# Residuals   27 10.492  0.3886                 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

saveRDS(anova_result, file = "anova_result.rds")

# 7) Validate the conditions for applying the ANOVA test. 
plot(anova_result)

# no pb of equal var, look good, but :
plot(anova_result, 2)
## qqplot looks to have pb... let's test :
shapiro.test(anova_result$residuals)
# Shapiro-Wilk normality test
# data:  anova_result$residuals
# W = 0.96607, p-value = 0.4379

# Conclusion : finally it is ok (p>0.05), normal residual, good. 

# 8) Post-hoc test: If the ANOVA is significant, perform a post-hoc test (Tukey's test) to identify which pairs of groups are significantly different.
TukeyHSD(anova_result)

# $group
#             diff        lwr       upr     p adj
# trt1-ctrl -0.371 -1.0622161 0.3202161 0.3908711
# trt2-ctrl  0.494 -0.1972161 1.1852161 0.1979960
# trt2-trt1  0.865  0.1737839 1.5562161 0.0120064

## trt 2 is signif different from other 2

# 9) The second scientific question in this exercise is whether plant growth
# differs depending on the treatment condition and the site.

# a) Perform a two-factor ANOVA test to analyze the interactions of these two variables,
# as follow :
anova_two_factor <- aov(weight ~ group * site, data = mydata_clean)
summary (anova_two_factor)
#               Df Sum Sq Mean Sq F value Pr(>F)  
# group        2  3.766  1.8832   5.401 0.0116 *
# site         1  0.356  0.3564   1.022 0.3221  
# group:site   2  1.768  0.8839   2.535 0.1003  

# b) Interpret the results: Is there an interaction between group and site? 
# If the two-factor ANOVA shows a significant interaction, 
# it means that the effect of group depends on the level of site.

# here as you can see, site has no signif effect, and site across group neither. 

# c)
boxplot(weight ~ site, data = mydata_clean)
boxplot(weight ~ group * site, data = mydata_clean)


#### 6.2- GLM log mtcars ####

# 1) From the mtcars dataset, create a binary variable named low_mpg:
data(mtcars)
median_mpg <- median(mtcars$mpg)
mtcars$low_mpg <- ifelse(mtcars$mpg <= median_mpg, 1, 0)
mtcars$low_mpg <- factor(mtcars$low_mpg)
  
# 2) The mtcars dataset contains many variables, but today we will only focus on 
# low_mpg explained by the following parameters:
# disp(Displacement), qsec (acceleration) and am (automatic transmission = binary variable). Create a dataset named “mycars” that will only contain these four columns of interest (and all rows). 
# For the rest of this exercise, use only the “mycars” dataset.
mycars <- mtcars[, c("low_mpg", "disp", "qsec", "am")]

# 3) Check your data format and describe. 
summary(mycars)
mycars$am <- as.factor(mycars$am)

gtsummary::tbl_summary(
  mycars, 
  by = low_mpg, 
  statistic = list(all_continuous() ~ "{mean} +/-{sd} {median} ({IQR})")
)

# 4) Check graphically the link between the “low_mpg” binary status and the 
# variables of interest. 

with(mycars, table(low_mpg, am))
boxplot(qsec ~ low_mpg, data = mycars)
boxplot(disp ~ low_mpg, data = mycars)


# 5) Create the regression model that will explain the “low_mpg” binary
#status based on the 3 parameters of interest: disp, qsec and am.

glm_mod <- glm(
  low_mpg ~ disp + qsec + am, 
  data = mycars, 
  family = binomial(link = "logit")
)

summary(glm_mod)
# Coefficients:
#                 Estimate Std. Error z value Pr(>|z|)
# (Intercept)  40.12968   29.53570   1.359    0.174
# disp          0.01632    0.01696   0.962    0.336
# qsec         -2.21039    1.45759  -1.516    0.129
# am1         -10.46683    6.62252  -1.580    0.114

# 6) Are the variables significant? Does this seem consistent with the 
# initial observations from the box plots?

# no var signif.. strange according the boxplot,
# we could expect a strong link between disp and low_mpg! 

# 7) Checking for multicollinearity: 

# library(car)
vif(glm_mod)

# Identify a variable with a high VIF (> 5 or > 10).
#     disp     qsec       am 
# 1.159923 8.436442 8.017619
## qsec has a high VIF

cor.test(mycars$disp, mycars$qsec)
# Pearson's product-moment correlation
# data:  mycars$disp and mycars$qsec
# t = -2.6363, df = 30, p-value = 0.01314
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.6796151 -0.1001493
# sample estimates:
#        cor 
# -0.4336979 

## Correlation between disp and qsec is a problem here !

# Adjust a model without this variable.
glm_mod_correct <- glm(
  low_mpg ~ disp + am, 
  data = mycars, 
  family = binomial(link = "logit")
)

summary(glm_mod_correct)
#             Estimate Std. Error z value Pr(>|z|)   
# (Intercept) -5.73431    2.40390  -2.385  0.01706 * 
# disp         0.03199    0.01227   2.608  0.00911 **
# am1         -2.27128    1.91286  -1.187  0.23508   


# 8) Prediction and confusion matrix
# a) Obtain the predicted probabilities:
mycars$prob <- predict(glm_mod_correct, type = "response")
# b) Choose a threshold (default: 0.5)
mycars$pred_class <- ifelse(mycars$prob > 0.5, 1, 0)
# c) Display the confusion matrix: 

table(Pred = mycars$pred_class, Obs = mycars$low_mpg)
with(mycars, table(pred_class, low_mpg))

# d) Calculate accuracy, sensitivity (Se), specificity (Sp). 
# Accuracy = (TP + TN) / (TP + TN + FP + FN) # correct classification rate
# Se = TP / (TP + FN)
# Sp = TN / (TN + FP)
# Where T is for TRUE, F is for FALSE, P for Positive, and N for Negative. 

TP <- sum(mycars$low_mpg %in% 1 & mycars$pred_class %in% 1)
TP
TN <- sum(mycars$low_mpg %in% 0 & mycars$pred_class %in% 0)
TN
FP <- sum(mycars$low_mpg %in% 0 & mycars$pred_class %in% 1)
FP
FN <- sum(mycars$low_mpg %in% 1 & mycars$pred_class %in% 0)
FN

AccuracyRate <- (TP + TN) / (TP + TN + FP + FN) 
AccuracyRate # 0.90625

Se <- TP / (TP + FN)
Sp <- TN / (TN + FP)

Se ; Sp
# [1] 0.8823529
# [1] 0.9333333

# e) save csv 
write.csv(x = mycars, file = "mycars_prediction.csv", row.names = FALSE)

# 9) ROC Curve
# The choice of threshold is, by definition, a parameter that can be varied. Se and Sp depend on the threshold s.
# At each threshold, we obtain a point (1-Sp(s), Se(s)) on the ROC curve.
# Let's try to plot this curve with the pROC package:
# library(pROC)
roc_obj <- roc(mycars$low_mpg, mycars$prob)
plot(roc_obj, col = "blue")
auc(roc_obj)

# a) What is the AUC of the model? 
# Area under the curve: 0.9843
# b) Does the model have acceptable discriminatory power?
# 0.5: poor
# ~0.7: acceptable
# ~0.8-0.9: good
# 0.9: excellent ==> Excellent ! 


# 10) also look at 
library(performance)

performance(glm_mod_correct)

performance::check_model(glm_mod_correct)
# automatic selection of appropriate check according the nature of the model 



performance::performance_roc(glm_mod_correct)
# AUC: 98.43%

performance::check_collinearity(glm_mod) # previous model with pb...
performance::check_collinearity(glm_mod_correct) # all good now

performance::check_predictions(glm_mod_correct)
?performance::performance_accuracy()
# to go further , estimate accuracy with cross validation (if enought data), or bootstrap otherwise
performance::performance_accuracy(glm_mod_correct, method = "boot")

# and many other functions... 



## Poisson reg ##

reg_poisson <- glm(NB_CHUTE ~ BRAS + SCORE_COMBINE + SEXE + AGE, data = base_count, family = "poisson")
summary(reg_poisson)
performance::check_overdispersion(reg_poisson)
# we test hyp about dispersion
# https://easystats.github.io/performance/reference/check_overdispersion.html 

# also possible :
# AER::dispersiontest(reg1,trafo=NULL) # library(AER) # pour dispersiontest


