# Practical 4 : 2025/11/21

# Study of juice yield from three apple varieties

#### envir ####
# getwd()
# setwd(file.path("...", "..."))
# list.files()

#### 4.2- Dataset introduction ####

## data in group
pommes_by_group <- data.frame(
  Golden = c(48,46,52,50),
  Delicious = c(52,50,49,49), 
  Jonagold = c(53,51,55,57)
)

## data long format
pommes <- data.frame(
  rendement = c(48,46,52,50,52,50,49,49,53,51,55,57),
  variete = factor(rep(c("Golden","Delicious","Jonagold"), rep(4,3)))
)

# test 
variete_facteur <- gl(n = 3, k = 4, label = c("Golden","Delicious","Jonagold"))
help("gl")


#### 4.3 - Understand the studied population before testing. ####

str(pommes)
## 'data.frame':	12 obs. of  2 variables:
##  $ rendement: num  48 46 52 50 52 50 49 49 53 51 ...
##  $ variete  : Factor w/ 3 levels "Delicious","Golden",..: 2 2 2 2 1 1 1 1 3 3 ...
is.data.frame(pommes)
## [1] TRUE
is.numeric(pommes$rendement)
## [1] TRUE
is.factor(pommes$variete) 
## [1] TRUE
class(pommes)
## [1] "data.frame"
class(pommes$rendement)
## [1] "numeric"
class(pommes$variete)
## [1] "factor"


?class 

## viz 

?boxplot

## base / stats

means <- aggregate(rendement ~ variete, data = pommes, FUN = mean)
boxplot(formula = rendement ~ variete, data = pommes, main = "Boite à moustaches")
points(1:3, means$rendement, col = "red")

## Since the columns of the apples object have the correct class , 
## the plot function adapts and also returns the boxplot graph.  
## try 
plot(formula = rendement ~ variete, data = pommes)

# Boxplots are displayed for the observations corresponding to each level (modality) 
# of the “variete” factor. If the boxplots are shifted, we can suspect an effect of the factor.

summary(pommes$rendement)
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   46.00   49.00   50.50   51.00   52.25   57.00
summary(pommes$variete) # summary.factor(pommes$variete)
## Delicious    Golden  Jonagold 
##         4         4         4
table(pommes$variete)
## 
## Delicious    Golden  Jonagold 
##         4         4         4


## what about missing data ?
# summary will specify the number of missing observations in the NA's column.
# table only displays the number of missing observations if the useNA option is specified
 # with the value “ifany” or “always.”
# To obtain this answer, we can either read the documentation for 
?table
# or directly experiment on a small test vector table(c(1, 2, 3, NA)).


## sum up data by group

## base / stat
n_group <- aggregate(rendement ~ variete, data = pommes, FUN = length)
mean_group <- aggregate(rendement ~ variete, data = pommes, FUN = mean)
sd_group <- aggregate(rendement ~ variete, data = pommes, FUN = sd)
na_group <- aggregate(rendement ~ variete, data = pommes, FUN = function(x) sum(is.na(x)))

## example with tapply
mean_tapply <- tapply(X = pommes$rendement, INDEX = pommes$variete, FUN = mean)

## NB : 
# It is important to keep in mind the default options for these basic functions when applying them.
# A good practice would be to make them explicit each time they are used. 
# For example, pay attention to the option na.rm = FALSE, 
# which indicates that missing values will not be removed by default. 
# Here, the calculation of the mean would fail if there were missing values.


## check normality 
## With a density plot
plot(
  density(pommes$rendement), 
  col = "black", 
  lty = 1, 
  main = "Distribution"
)
mr = mean(pommes$rendement, na.rm = TRUE)
sdr = sd(pommes$rendement, na.rm = TRUE)
x_norm = seq(-4,4,length=100) * sdr + mr
lines(
  x = x_norm, 
  y = dnorm(x = x_norm, mean = mr, sd = sdr), 
  col = "red", lty = 2
)

## with histogram
hist(x = pommes$rendement, main = "Histogram")

##  with a normal qqplot
qqnorm(pommes$rendement)
qqline(pommes$rendement)


shapiro.test(x = pommes$rendement)
##  Shapiro-Wilk normality test
## 
## data:  pommes$rendement
## W = 0.97643, p-value = 0.9653

# We recall the null hypothesis H0: the variable is normally distributed.
# If the p-value is less than a chosen alpha level (e.g., 0.05), 
# then the null hypothesis is rejected.
# To assume normality of the residuals, 
# it is therefore necessary to obtain a p-value > 0.05. 
# (Here p-value = 0.96, we do not reject H0)


#### 4.4 - Finally, we test it. ####

?stats::aov
?stats::lm
?stats::anova

# ?aov said :
  # "Fit an analysis of variance model by a call to lm for each stratum."

# ?lm said :
  # lm is used to fit linear models. 
  # It can be used to carry out regression, 
  # single stratum analysis of variance and analysis of covariance 
  # (although aov may provide a more convenient interface for these).

# ?anova said :
  # Compute analysis of variance (or deviance) tables for one or more
  # fitted model objects.

##  so ? 
# The aov function allows us to perform variance analysis,
# and we can see that it calls the lm function for each level 
# (the modalities of our factor).
# Therefore, lm is also suitable for answering our question, 
# since we want to perform variance analysis.
# We can see that lm(rendement ~ variete, data = pommes) is therefore
# equivalent to aov(formula = rendement~variete, data = pommes) in our case.

# However, the output of the aov function is more convenient
# for answering the ANOVA question.
# The anova function allows us to test the significance of the predictors,
# so anova(lm(rendement ~ variete, data = pommes)) would also answer our question.

# N.B.: FYI, the anova function also allows us to test the addition of a predictor 
# between a “null” model and an alternative model, for example 
# (comparison between two or more models), 
# but this point will not be discussed here.

## Perform the ANOVA 

my_anova <- aov(formula = rendement~variete, data = pommes)
my_anova
## Call:
##    aov(formula = rendement ~ variete, data = pommes)
## 
## Terms:
##                 variete Residuals
## Sum of Squares       56        46
## Deg. of Freedom       2         9
## 
## Residual standard error: 2.260777
## Estimated effects may be unbalanced

# Actually, “my_anova” is a list with 13 components, 
# the names of which can be found by typing names(my_anova).
# Tip: to access the elements of a list, 
# you can start typing my_anova$ in R and use auto-complete with the “TAB” key 
# on your keyboard. R will then suggest the names known in the list.

summary(my_anova)
##             Df Sum Sq Mean Sq F value Pr(>F)  
## variete      2     56  28.000   5.478 0.0278 *
## Residuals    9     46   5.111                 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


## conclusion

# Therefore, taking alpha=0.05, since p.value = 0.02778 < 0.05,
# we reject the null hypothesis H0:(m1=m2=m3).
# The varieties have different juice yields, and in saying this, 
# we only have a 5 in 100 chance of being wrong. 
# 
# We conclude that there is a significant variety effect.
# 
# N.B.: However, if we take alpha=0.01, we do not reject H0. 
# The effect is not very significant.


## Perform the diagnosis

par(mfrow=c(2, 2));
plot(my_anova)
par(mfrow=c(1,1))

## diag : 
# + No correlation between residuals: 
#   The residuals do not appear to be correlated with each other 
# (Residuals vs. Fitted) since the red line is horizontal and at Y = 0. 
# The value of the residuals does not appear to depend on the treatment 
# since they are all centered around Y = 0 (Residuals vs. Factor Levels).
# 
# + To check the normality of the residuals, we look at the Normal QQplot.
# Here, the points are well distributed along the line, 
# which means that the residuals are distributed according to a normal
# distribution. The fact that the points are centered on 0 (on the y-axis)
# shows that their mean is equal to 0.
# 
# + The assumption of homogeneity of variances, i.e.,
# the assumption that the residuals have a constant variance,
# can be evaluated using the “Scale-Location” graph. The graphical
# method consists of plotting the standardized residuals against the 
# predicted values (the means of the different factors). Here we can 
# see that the dispersions of the residuals (their vertical deviations) 
# relative to each treatment modality are broadly identical,
# so the assumption of homogeneity of residuals is accepted.

# Visually, these graphs appear to be within the norm. 
# The stochastic assumptions are therefore acceptable,
# and we can therefore “trust” the results of the ANOVA test.

# Another solution: We can also access the model residuals, i.e.,
# my_anova$residuals, and then test their normality, autocorrelation, etc.

## Test post-hoc

?TukeyHSD

TukeyHSD(x = my_anova)
##   Tukey multiple comparisons of means
##     95% family-wise confidence level
## 
## Fit: aov(formula = rendement ~ variete, data = pommes)
## 
## $variete
##                    diff        lwr     upr     p adj
## Golden-Delicious     -1 -5.4633295 3.46333 0.8101561
## Jonagold-Delicious    4 -0.4633295 8.46333 0.0784642
## Jonagold-Golden       5  0.5366705 9.46333 0.0296317
par(cex.axis=0.4)
plot(TukeyHSD(my_anova))

## NB : 
# This test provides us with adjusted p-values in the “p adj” column. 
# Multiple test correction has already been taken into account.

## It is the Jonagold variety that has, on average,
# a different yield from the others.

## to try 
pairwise.t.test(x = pommes$rendement, g = pommes$variete, p.adj = "bonf")
pairwise.t.test(x = pommes$rendement, g = pommes$variete, p.adj = "holm")

