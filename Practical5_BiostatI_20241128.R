# Practical 5 : 2025/11/28

# LM, compare modeles, GLM

#### envir ####
getwd() # know where you are
# setwd(file.path("...", "...")) # setup a working directory
# list.files()

#### 5.1- trees ####

# data(datasets::trees)
head(trees)
pairs(trees)

str(trees)

# Volume = β_0 + β_1 × Girth + β_2 × Height + ε
reg <- lm(
  formula = Volume ~ Girth + Height, 
  data = trees
)
summary(reg)

## 1) 
# n=31, see str, also can get it via df (d.d.l.)
# or 
nrow(trees) # get the number of lines. 
# and
# 2 variables used to explain the variable Volume
# but 3 parameters are estimated :
# B0 (for intercept), B1 the effect of Girth and B2 the effect of Height

## 2) beta_2 is for Height = 0.3393
coef(reg) # returns coefficents
coef(reg)[3]  # returns directly Beta2 
# we can also do :
res <- summary(reg) # and  
res$coefficients

## 3) 
# See summary(reg) 
# Highly signif yes because p-valeur <0.001 => symbol ***
  
## 4) 
confint(reg, level = 0.95)
# or more precisly : 
confint(reg, level = 0.95)[3, ]

## 5) 
# see `summary(reg)`
# "Multiple R-squared:  0.948,	Adjusted R-squared:  0.9442"
# We can also get it via
res <- summary(reg)
res$r.squared # and 
res$adj.r.squared

## 6)
# See summary(reg)
# "F-statistic:   255" and "p-value: < 2.2e-16
# H0 : "all betas = 0" vs H1 : “At least one coefficient not null”.
# here we reject H0, it means the model is relevant 

# FYI: The Fisher statistic can be found here too :
res$fstatistic
# With the function `pf`, We could re-compute the exact p-value of Fisher 
pf(
  q = res$fstatistic[["value"]],
  df1 = res$fstatistic[["numdf"]],
  df2 = res$fstatistic[["dendf"]],
  lower.tail=FALSE
)

## 7) prediction 
## 7-a) 
## with R function: 
predict(reg, data.frame(Girth = 8.3, Height = 70))
# predict(reg)[1]
# Note that this is precisely the first observation in the dataset,
# and predict is based by default on the dataset used for regression. 

## Manually : prediction is actually y_x = b0 * 1 + b1 * Girth + b2 * Height
sum(reg$coefficients * c(const = 1, Girth = 8.3, Height = 70))

## 7-b) 
predict(
  reg, 
  data.frame(Girth = 8.3, Height = 70),
  interval = "confidence"
)

## 8) Graph diag

e <- residuals(reg)
plot(e)
abline(h = 0, col = "red")
## ou directly plot on the reg object to get the 4 diag graphs : 
par(mfrow = c(2, 2))
plot(reg)
## To focus on the 1st graph, it seems to have a problem : 
par(mfrow = c(1,1)); plot(reg, 1)
## a structure is well visible,
# We do NOT have symmetry around the axis y = 0. 
# One might assume that residues are dependent

# to check normality : QQplot 
par(mfrow = c(1,1)); plot(reg, 2)
par(mfrow = c(1,1)); plot(reg, 4)

## To go further about independence of errors : 
par(mfrow = c(1, 2))
acf(e)
pacf(e)

par(mfrow = c(1, 1)) ## get back params graphic 

# The goal is to test for temporal independence or the absence of autocorrelation 
# among the residuals. If the model has accurately captured the structure present 
# in the data, then the residuals should be “like noise” — no correlation between 
# a value at time t and the value at time t + lag.
# ACF plots the total correlation between the residual series and itself at 
# different lags. PACF plots the partial correlation—that is, the direct
# correlation at one lag, after “removing” the influence of intermediate lags.
# If the residuals are “good,” we expect these correlations to be (statistically) 
# close to zero for all lags > 0 - this is the behavior of white noise (i.i.d.).

# # how to interpret graphs
# If the majority of spikes for all lags are within the confidence intervals 
# =>  no significant correlation =>
#   the residuals are compatible with white noise.
# 
# In other words, there is no visible autocorrelation structure, 
# which means that your model has not left any unmodeled “memory” in the residuals.
# 
# In this case, this would support the hypothesis of error independence—one 
# of the classic assumptions of linear models.


## 9) multi colinearity 
cor(x = trees$Girth, y = trees$Height)^2
# to compare with
summary(reg)$r.squared

# cor returns 0.26, far from R2 = 0.948.
# So according the Klein's rule,
# no linear link between Girth and Height. good !

## 10) Extrem values

## 10-a) cook d
cook <- cooks.distance(reg)
cook[cook>1]

# no distance > 1, 
# however the graphs `plot(e)` and `plot(reg, 2)` (above) 
# shows a point a little far from others... (num 31).
# so we can check (10-B)

## 10-b) 
summary(influence.measures(reg))

# The answer is mixed because, using only Cook's distance 
# as a criterion, no values were outliers. However, 
# we can still see that observation 31 is “influential”
# (meaning that it can “drive” the signal, i.e., direct
# (or pull up or down) our coefficients—the slope of the 
# regression line).
# It would be a good idea to look more closely at whether 
# observation “31” stands out “too much” from the others 
# in the data graphs (e.g., “pairs”) and perhaps remove it 
# from our model.


#### 5.2 - new dataset : compare models ####

## 1) read
library(readxl)
mydata <- readxl::read_excel(
  file.path("../Desktop/cp biostat I", "PS5_trees_data.xlsx")
)
## 2) check 
str(mydata)

##3) anova between 2 lm models 
reg1 = lm(Volume ~ Girth + Height + X3 + X4, data = mydata)
reg2 = lm(Volume ~ Girth + Height, data = mydata)
anova(reg1, reg2)

## 4) 
# Reg1 is the alternative model tested, Reg2 is the model null or the base. 

# We read the p-value = 0.8403.
# If p-value > 0.05, then the variables studied do not contribute significantly
# to the model.
# Here, the best model is model reg2: 
# meaning that X3 and X4 do not contribute significantly to the model.

## 5) AIC BIC

message("reg1")
message("AIC = ", AIC(reg1))
message("BIC = ", BIC(reg1))

message("reg2")
message("AIC = ", AIC(reg2))
message("BIC = ", BIC(reg2))

# No, the Reg2 model still has the smallest values.
# Specifically, the AIC criterion does not seem to be particularly helpful 
#(almost identical values for reg1 and reg2).
# However, the BIC criterion is much smaller for reg2. Based on this information, 
# reg2 is a preferable model to reg1.



#### 5.3 - GLM ####

## 1) 
T2Ddata <- data.frame(
  weight = c(
    35.9, 38.3, 55.7, 41.7, 43.2, 49.1, 45, 45.3, 46.1, 46.9, 48.1, 
    48.9, 49.2, 51.2, 56.4, 51.7, 51.8, 52.6, 52.9, 51.3, 53.7, 55, 
    55.4, 55.8, 58, 58.7, 60.3, 61.1, 61.5, 63.1
  ), 
  CC = factor(
    x = c(rep("0", 15), rep("1", 15)),
    levels = c("0", "1"), 
    labels = c("CTRL", "CAS")
  ) 
)
str(T2Ddata)

## 2) viz 
plot(T2Ddata$weight, T2Ddata$CC)
boxplot(T2Ddata$weight, T2Ddata$CC)
boxplot(T2Ddata$weight ~ T2Ddata$CC) # better to get the labels

## 3) glm
reg <- glm(
  formula = CC ~ weight,
  data = T2Ddata,
  family = binomial(link = "logit")
)
summary(reg)

## 4) OR
OR <- exp(coef(reg)) # store it 
OR # show it
OR > 1 # check it > 1

# The OR for weight is > 1, so an increase of one unit in weight leads
# to an increase in the probability that {CC = 1} will occur, i.e.,
# according to our data, an increase in weight increases the risk of 
# developing diabetes.

## check
plot(reg)
plot(reg, 1)
# horizontal line, good
plot(reg, 2)
# std deviance residual looks ~ Normal
# may be pb with obs 3 and 15
plot(reg, 3)
# may be pb with obs 3 and 15
plot(reg, 4)
# ~ extrem obs but real... 
