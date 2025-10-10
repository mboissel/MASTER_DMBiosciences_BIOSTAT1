## Cours M1 Biostat  ##
## La Catho Univ ##

## Mathilde Boissel 


#### basic introduction ####

1 + 1
1 - 1
10/2
10*3
10%%2 # modulo

# store the result with arrow "<-"
x <- 1 + 1
y <- 2 * 2
z <- x - y
vec <- c(x, y, z)
mylist <- list(vec, x)
mylist

# you can store with "=" but not recommand to avoid confusion with logic test "=="

# logic test
x == 2
x %in% vec
vec %in% x
TRUE & TRUE
TRUE & FALSE
FALSE & FALSE
TRUE | TRUE
TRUE | FALSE
FALSE | FALSE

is.numeric(x)
is.character(x)
is.character(x)
is.factor(x)

?iris # from : # library(datasets)
head(iris)

is.data.frame(iris)
is.list(iris) # a data frame is a list of vector

# access to an element of a list : "$"
iris$Sepal.Length
is.vector(iris$Sepal.Length)
is.numric(iris$Sepal.Length)

# also possible with `[[`
iris[["Species"]]
is.factor(iris[["Species"]])
is.character(iris[["Species"]])

my_tab <- data.frame(
  let = LETTERS[1:10], 
  num = 11:20, 
  other = as.factor(c(rep("blue", 5), rep("red", 5)))
)
my_tab
str(my_tab)
# 'data.frame':	10 obs. of  3 variables:
#   $ let  : chr  "A" "B" "C" "D" ...
# $ num  : int  11 12 13 14 15 16 17 18 19 20
# $ other: Factor w/ 2 levels "blue","red": 1 1 1 1 1 2 2 2 2 2




#### Set up envir ####
# load pkg
library(data.table)
library(ggplot2)

#### 1.1 – ELEMENTS OF CALCULUS ####

##### 1.1 - D	In practice: R #####
iris # from datasets
dim(iris) ; nrow(iris) ; ncol(iris) ;
head(iris)
summary(iris)
str(iris)
summary(iris$Sepal.Length)
mean(iris$Sepal.Length) ; var(iris$Sepal.Length) ; sd(iris$Sepal.Length) ;
quantile(x = iris$Sepal.Length, probs = seq(0, 1, 0.25))
length(iris$Sepal.Length) 
table(iris$Species) ; prop.table(table(iris$Species)) ;

## more code to test, if wanted: 

# open help panel with question mark before a function name:
?str

# Simulate a simple dataset
set.seed(42) # set the seed to reproduce the random generation.
x <- rnorm(100, mean = 50, sd = 10) # 100 values from a normal distribution

# Basic descriptive stats
mean_x <- mean(x)
sd_x <- sd(x) # Standard Deviation
se_x <- sd_x / sqrt(length(x)) # Standard Error

# Display the results
cat("Mean:", round(mean_x, 2), "\n")
cat("Standard Deviation (SD):", round(sd_x, 2), "\n")
cat("Standard Error (SE):", round(se_x, 2), "\n")

# Test basic graph
plot(density(x))

# Test ggplot style
library(ggplot2)

df <- data.frame(x = x)
ggplot(data = df, mapping = aes(x = x)) +
  geom_density() 


# Visualize Mean, SD and SE difference !
ggplot(data = df, mapping = aes(x = x)) +
  geom_histogram(binwidth = 2, fill = "skyblue", color = "black") +
  geom_vline(xintercept = mean_x, color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = mean_x + sd_x, color = "orange", linetype = "dotted", size = 1) +
  geom_vline(xintercept = mean_x - sd_x, color = "orange", linetype = "dotted", size = 1) +
  geom_vline(xintercept = mean_x + se_x, color = "green", linetype = "dotdash", size = 1) +
  geom_vline(xintercept = mean_x - se_x, color = "green", linetype = "dotdash", size = 1) +
  labs(
    title = "Histogram with SD and SE around the Mean",
    subtitle = "Red = Mean | Orange = ± SD | Green = ± SE",
    x = "x", y = "Frequency"
  ) + 
  theme_light()

#### 2.2 - PROBABILITY THEORY ####

##### 2.2 - A PMF – Discrete Variables ##### 

# Probability Mass Function (PMF) – Discrete Variables
x <- 0:1
probs <- c(0.5, 0.5)
barplot(probs, names.arg = x, main = "PMF of fair coin", ylab = "P(X=x)")

##### 2.2 - B PDF - Continuous Variables ##### 

# Probability Density Function (PDF) – Continuous Variables



#### 2.3 - USUAL PROBABILITY DISTRIBUTIONS ####

##### 2.3-A-2 Binomial #####

# Simulate and plot a binomial distribution
x <- 0:10
probs <- dbinom(x, size = 10, prob = 0.5)
barplot(probs, names.arg = x, main = "Binomial(n=10, p=0.5)", ylab = "P(X = x)")

## real world example binomial
n = 150
p = 0.1
mu = n*p # 15
sigma = sqrt(n*p*(1-p)) # 3.67
# visualiser la majorité de la distribution
cat("[", mu - 3 * sigma, " ; ", mu + 3 * sigma, "]")
# [ 3.977296  ;  26.0227 ]
## La plupart des probabilités se concentrent entre 3 et 26 patients ayant un effet secondaire.
# Define the range of possible outcomes (e.g., 0 to 40 patients with side effects)
x <- 0:40
# Binomial probabilities for each possible number of patients with the side effect
probs <- dbinom(x, size = 150, prob = 0.1)
# Plot the distribution
barplot(
  probs,
  names.arg = x,
  main = "Binomial(n = 150, p = 0.1)",
  xlab = "Number of patients with side effect",
  ylab = "P(X = x)",
  las = 2,  # rotate x-axis labels for readability
  cex.names = 0.7
)  # adjust label size
## ou Tracer la fonction de masse (densité) de la loi binomiale
# avec type = "h"
plot(
  x, probs,
  type = "h", # trace des segments verticaux (comme des tiges) : adapté aux lois discrètes
  lwd = 2, col = "blue",
  main = "Binomial PMF: n = 150, p = 0.1",
  xlab = "Number of patients with side effects",
  ylab = "Probability"
)
points(x, probs, pch = 16, col = "blue")  # Ajouter les points en haut des lignes


##### 2.3-A-3 Poisson #####

# Simulate and plot a poisson distribution
x <- 0:15
probs <- dpois(x, lambda = 4)
barplot(probs, names.arg = x, main = "Poisson(λ=4)", ylab = "P(X = x)")

## real world example poisson

lambda <- 4 # Define the rate (average number of events)
x <- 0:15 # Range of values to plot (from 0 to 15 events/day) 
probs <- dpois(x, lambda)
# Plot the probability mass function # type = "h"
plot(
  x, probs, type = "h", lwd = 2, col = "darkgreen",
  main = "Poisson Distribution (λ = 4)",
  xlab = "Number of asthma-related visits per day",
  ylab = "P(X = x)"
)
points(x, probs, pch = 16, col = "darkgreen")  # Add points at the top

##### 2.3-A-4 NB #####

# Simulate and plot a NB distribution
x <- 0:20
probs <- dnbinom(x, size = 5, prob = 0.4)
barplot(probs, names.arg = x, main = "Negative Binomial(r=5, p=0.4)", ylab = "P(X = x)")

## real world example NB

# Parameters for the Negative Binomial
mu <- 2       # Mean number of infections
size <- 1.5   # Dispersion parameter (lower = more overdispersion)
x <- 0:15 # Range of possible values for the count of infections
probs <- dnbinom(x, size = size, mu = mu)
# Plot the probability mass function
plot(x, probs, type = "h", lwd = 2, col = "steelblue",
     main = "Negative Binomial Distribution (μ = 2, size = 1.5)",
     xlab = "Number of infections per patient",
     ylab = "P(X = x)")
points(x, probs, pch = 16, col = "steelblue")


##### 2.3-B-1 Norm #####

# Simulate and plot a normal distribution
x <- seq(-4, 4, length = 100)
y <- dnorm(x, mean = 0, sd = 1)
plot(x, y, type = "l", main = "Standard Normal Distribution", ylab = "Density")

library(ggplot2)
ggplot(data = data.frame(x = c(-6, 6)), aes(x)) +
  stat_function(aes(color = "mu = 0"), fun = dnorm, n = 101, args = list(mean = 0, sd = 1)) + 
  stat_function(aes(color = "mu = 1"), fun = dnorm, n = 101, args = list(mean = 1, sd = 1)) + 
  stat_function(aes(color = "mu = 2"), fun = dnorm, n = 101, args = list(mean = 2, sd = 1)) + 
  ylab("Density") +
  labs(title = "Normal Distributions N(mu, s = 1)") + 
  # scale_y_continuous(breaks = NULL) + 
  scale_color_manual(
    name = "Means", 
    values = c("mu = 0" = "firebrick", "mu = 1" = "dodgerblue", "mu = 2" = "green")
  ) + 
  theme_light()

ggplot(data = data.frame(x = c(-6, 6)), aes(x)) +
  stat_function(aes(color = "s = 1"), fun = dnorm, n = 101, args = list(mean = 0, sd = 1)) + 
  stat_function(aes(color = "s = 2"), fun = dnorm, n = 101, args = list(mean = 0, sd = 2)) + 
  stat_function(aes(color = "s = 2.5"), fun = dnorm, n = 101, args = list(mean = 0, sd = 2.5)) + 
  ylab("Density") +
  labs(title = "Normal Distributions N(mu = 0, s)") + 
  # scale_y_continuous(breaks = NULL) + 
  scale_color_manual(
    name = "SD (s)", 
    values = c("s = 1" = "firebrick", "s = 2" = "dodgerblue", "s = 2.5" = "green")
  ) + 
  theme_light()


##### 2.3-B-2 Exp #####

# Simulate and plot a exp distribution
x <- seq(0, 10, length = 100)
y <- dexp(x, rate = 0.5)
plot(x, y, type = "l", main = "Exponential(λ=0.5)", ylab = "Density")

##### 2.3-B-3 Gamma #####

# Simulate and plot a gamma distribution
x <- seq(0, 20, length = 100)
y <- dgamma(x, shape = 2, rate = 0.5)
plot(x, y, type = "l", main = "Gamma(α=2, β=0.5)", ylab = "Density")

##### 2.3-B-4 Unif #####

# Simulate and plot a uniform distribution
x <- seq(0, 1, length = 100)
y <- dunif(x, min = 0, max = 1)
plot(x, y, type = "l", main = "Uniform(0, 1)", ylab = "Density")


#### 3. DATA OBS ####

#### 3.1 – DESCRIPTIVE STATISTICS ####

##### ExA #####

collection1 <- c(170, 190, 201, 204, 209, 250, 254, 257, 280)
length(collection1)
summary(collection1)
mean(collection1)
median(collection1)
var(collection1)
sqrt(var(collection1))
sd(collection1)
quantile(collection1)

##### ExB #####

collection2 <- c(170, 190, 201, 204, 209, 250, 254, 257, 280, 148)
sort(collection2)
dput(sort(collection2))

collection2 <- c(148, 170, 190, 201, 204, 209, 250, 254, 257, 280)
length(collection2)
summary(collection2)
mean(collection2)
median(collection2)
var(collection2)
sqrt(var(collection2))
sd(collection2)
quantile(collection2)

##### ExC #####
# copy the value from word :
collection3 <- c("Yes	Yes	No	No	No	Yes	No	No	No	No")
collection3
# but we want a vect with sep elements : 
strsplit(x = collection3, split = "\t")
collection3 <- strsplit(x = collection3, split = "\t")[[1]]

table(collection3)
table(collection3, useNA = "ifany")

##### ExD #####
Travel_out_FR <- c("Yes
No
No
No
Yes
No")
Travel_out_FR
strsplit(Travel_out_FR, split = "\n")
Travel_out_FR <- strsplit(Travel_out_FR, split = "\n")[[1]]
table(Travel_out_FR)

Where <- c("Italia



Spain

")
Where
Where <- strsplit(Where, split = "\n")[[1]]
table(Where, useNA = "always")

is.na(Where) # TRUE / FALSE
Where %in% "" # FALSE  TRUE  TRUE  TRUE FALSE  TRUE
Where[Where %in% ""] <- "France"
# if it was NA : 
# Where[is.an(Where)] <- "France"
table(Where, useNA = "always")

# if we want it as NA : 
Where[Where %in% "France"] <- NA
table(Where, useNA = "always")


##### ExE ####
boxplot(...)

IQR(...)


#### 3.2 – EXPLORATORY DATA ANALYSIS ####

set.seed(42)
data <- rnorm(1000, mean = 10, sd = 2)

##### 3.2-A-1- hist and bins #####

# Base R
hist(
  data, breaks = 30, col = "lightblue",
  main = "Histogram", xlab = "Values"
)
# or
plot(density(data), col = "firebrick", lwd = 2, main = "Density")
# ggplot2
library(ggplot2)
ggplot(data.frame(x = data), aes(x)) +
  geom_histogram(
    aes(y = ..density..), bins = 30,
    fill = "skyblue", color = "white"
  ) +
  geom_density(color = "firebrick2", size = 1) +
  ggtitle("Histogram + Density") + 
  theme_minimal()


##### 3.2-A-2- boxplot #####

# Base R
boxplot(data, main = "Boxplot")
points(mean(data), col = "red", pch = 19) # Add mean

# ggplot2
ggplot(data.frame(x = data), aes(x = "", y = x)) +
  geom_boxplot(fill = "lightgreen") +
  stat_summary(
    fun = mean, geom = "point", shape = 20, color = "red", size = 3
  ) +
  ggtitle("Boxplot with Mean") + 
  xlab("") +
  theme_light()

##### 3.2-A-3- Pie chart // bar plot #####

set.seed(92)
comm_method <- sample(
  c("Phone", "Email", "In-person", "Text"), 
  size = 2000, replace = TRUE, 
  prob = c(0.28, 0.27, 0.23, 0.22)  # Very close proportions
)

# Create a frequency table
table_comm <- table(comm_method)
df_comm <- as.data.frame(table_comm)
colnames(df_comm) <- c("method", "count")

# Bar plot (ordered)
library(ggplot2)
ggplot(df_comm, aes(x = reorder(method, -count), y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Preferred Communication Method",
       x = "Method", y = "Number of Participants") +
  theme_light()


# Trap : Compare with a Pie Chart 
# Pie chart (harder to interpret small differences)
pie(
  df_comm$count,
  labels = df_comm$method, 
  main = "Pie Chart of Communication Methods"
)


##### 3.2-B-1- Scatter plot #####

set.seed(42)
x <- rnorm(100)
y <- 2 * x + rnorm(100)
# Base R
plot(x, y, main = "Scatter Plot")
abline(lm(y ~ x), col = "red")
# ggplot2
df <- data.frame(x, y)
ggplot(df, aes(x, y)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  ggtitle("Scatter Plot with Linear Trend") +
   theme_minimal()



##### 3.2-B-2- Confounders #####

set.seed(123)
library(ggplot2)

n <- 1000
q <- rbinom(n, size = 1, prob = .35)
x <- 2 * q + rnorm(n)
y <- -3 * q + rnorm(n)
confounder_data <- data.frame(x, y, q = as.factor(q))

p1 <- ggplot(data = confounder_data, aes(x, y)) +
  geom_point(alpha = .2) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  ggtitle("Not adjusting for `q` : biased") + 
  theme_light()

p2 <- ggplot(data = confounder_data, aes(x, y, color = q)) +
  geom_point(alpha = .2) +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("Adjusting for `q` : unbiased") + 
  theme_light()

p1
p2

##### 3.2-B-3- Heatmap Correlation #####

library(ggplot2)
library(reshape2)
library(data.table)

df <- mtcars[, 1:6]
corr <- round(cor(df), 2)
## with reshape
melted <- reshape2::melt(corr) 
head(melted,10)
## with data.table
melted_df <- data.table::melt(
  data = as.data.table(corr, keep.rownames = TRUE), 
  id.vars = "rn", 
  measure.vars = names(mtcars)[1:6],
  variable.name = "cn", 
  value.name = "valueofcorrelation"
)
melted_df
# viz heatmap = geom_tile
ggplot(melted_df, aes(rn, cn, fill = valueofcorrelation)) +
  geom_tile() +
  geom_text(aes(label = valueofcorrelation)) +
  scale_fill_gradient2(
    "Correlation",
    low = "blue", high = "red", mid = "white"
  ) +
  ylab(NULL) + xlab(NULL) + 
  ggtitle("Correlation Heatmap") + 
  theme_minimal()

# idea : Corr test, show heap map of pvalues 
# --here todo

##### 3.2-B-4- quali vs quanti #####

# Example: blood pressure by treatment group
df <- data.frame(
  treatment = rep(c("A", "B", "C"), each = 50),
  bp = c(rnorm(50, 120), rnorm(50, 130), rnorm(50, 125))
)

ggplot(df, mapping = aes(treatment, bp)) +
  geom_boxplot(fill = "lightblue") +
  stat_summary(
    fun = mean, geom = "point", shape = 20, color = "red", size = 3
  ) +
  ggtitle("BP by Treatment Group") +
  theme_minimal()

##### 3.2-B-5- spagetti plot #####

# Example: blood pressure over time
set.seed(42)
n_patient = 10
df <- data.frame(
  subjectid = as.factor(rep(1:n_patient, 4)),
  visits = as.factor(rep(c("V1", "V2", "V3", "V4"), each = n_patient)),
  bp = c(
    rnorm(n = n_patient, mean = 120, sd = 3), rnorm(n_patient, 130),
    rnorm(n_patient, 135), rnorm(n_patient, 145)
  )
)
df_bp <- data.frame(
  visits = as.factor(c("V1", "V2", "V3", "V4")), 
  bp_mean = tapply(X = df$bp, INDEX = df$visits, FUN = mean),
    # aggregate(formula = bp ~ visits, data = df, FUN = mean) # same
  bp_sd = tapply(X = df$bp, INDEX = df$visits, FUN = sd)
)

ggplot(data = df, mapping = aes(x = visits, y = bp, color = subjectid)) +
  geom_line(data = df, aes(group = subjectid), size = 0.5) + 
  geom_point(size = 1) + 
  geom_point(
    data = df_bp, mapping = aes(x = visits, y = bp_mean),
    shape = 4, inherit.aes = FALSE
  ) + 
  geom_errorbar(
    data = df_bp, 
    mapping = aes(
      x = visits,
      ymin = bp_mean - 2*bp_sd, ymax = bp_mean + 2*bp_sd
    ), 
    width = 0.2,
    inherit.aes = FALSE
  ) + 
  ggtitle("Evolution of Blood Pressure Over Time") +
  ylab("Blood Pressure (mmHg)") + 
  xlab("Visits") +
  theme_light() +
  theme(
    axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 12)
  )


##### 3.2-B-5- spagetti plot, labels #####

library(data.table)
library(ggplot2)
library(ggrepel) # for geom_text_repel
set.seed(42)
# data
propIV_by_service <- data.table::rbindlist(l = list(
  data.table::setDT(expand.grid(
    date_M = c("M1"),
    service = c("Emergency SV", "Rumathology SV", "Emergency SP", "Rumathology SP", "Geriatry SP")
  ))[
      , prop_IV := runif(5, min = 50, max = 100)
  ], 
  data.table::setDT(expand.grid(
    date_M = c("M3", "M6", "M12"),
    service = c("Emergency SV", "Rumathology SV", "Emergency SP", "Rumathology SP", "Geriatry SP")
  ))[
    , date_M := as.factor(date_M)
  ][
    , prop_IV := c(
      runif(2*3, min = 50, max = 100),
      runif(3*3, min = 0, max = 52)
    )
  ] 
))[, date_M := as.factor(date_M)]
prop_ends <- propIV_by_service[
  , .SD[which.max(date_M)],
  by = service
]
## gg 
evolution_service <- ggplot(
  data = propIV_by_service,
  mapping = aes(
    x = date_M, y = prop_IV,
    color = service, group = service
  )
) +
  geom_line() +
  geom_point() +
  geom_vline(
    xintercept = 1.5,
    linetype="dashed", 
    color = "aquamarine", size = 1
  ) + 
  scale_y_continuous(limits = c(0,100)) +
  scale_x_discrete(expand = expansion(add = c(0, 2))) + # ajouter une marge à droite
  coord_cartesian(clip = "off") + # ne pas tronquer les labels
  geom_text_repel(
    data = prop_ends,
    aes(label = service),
    hjust = -1,
    direction = "y", 
    nudge_x = 0.5,
    segment.alpha = 0.3,
    size = 3.2, 
    force = 1.5
  ) +
  xlab("Date Follow Up") + 
  ylab("IV proportion (%)") + 
  ggtitle("Evolution of IV proportion (%), by service" ) +
  theme_light() +
  guides(color = "none")

evolution_service


#### 4.3 – TESTING – PARAMETRIC OR NON-PARAMETRIC STATISTICS ####

set.seed(42)
data <- data.frame(age = rnorm(n = 42, mean = 25, sd = 3.3))

##### 4.3-A Graph #####
# A) How to Check for Normality: Graphical Method

qqnorm(data$age, main = "Q-Q Plot") ; qqline(data$age, col = "red", lwd = 2) ;
ggpubr::ggqqplot(data$age, title = "ggpubr Q-Q Plot of Age")

hist(data$age, col = "skyblue", main = "Histogram of Age", xlab = "Values", border = "white")
plot(density(data$age), main = "Density Plot", lwd = 2) 

boxplot(data$age, main = "Boxplot of Age")

##### 4.3-B Test #####
# B) How to Check for Normality: Computational Method

shapiro.test(data$age)
ks.test(data$age, "pnorm", mean = mean(data$age), sd = sd(data$age))
print(nortest::lillie.test(data$age)) 
print(nortest::ad.test(data$age))
print(nortest::cvm.test(data$age))
jarque.bera.test(data$age)

##### 4.3-C Trouble Shooting #####
# C) Your Data Fails the Normality Test: Trouble Shooting

Y_transformed <- log(Y)   # valid if Y > 0
Y_transformed <- log(Y + 1)   # if some values of Y are zero
Y_transformed <- sqrt(Y + c)   # where c is a constant ensuring positivity
Y_transformed <- 1 / Y   # don’t use if Y has 0s.

#### 4.4 – BIVARIATE TESTING ####

##### 2-sample T-test #####

t.test(BP ~ Group, data = df, var.equal = TRUE)    # Student t-test
t.test(BP ~ Group, data = df, var.equal = FALSE)   # Welch’s correction

##### One-way ANOVA #####

aov(Sepal.Length ~ Species, data = iris)

fit <- aov(Cholesterol ~ Diet, data = df)
summary(fit)

#####	Wilcoxon-Mann-Whitney test  #####

wilcox.test(VPSAH ~ Group, data = df)

##### Kruskal wallis #####

kruskal.test(Cholesterol ~ Diet, data = df)

##### Chi2 #####

tbl <- table(df$Smoking, df$Gender)
chisq.test(tbl)

##### fisher exact #####

tbl <- table(df$advers_event, df$ttt)
fisher.test(tbl)   # valid for small or unbalanced samples

##### Correlation test #####

cor.test(df$BMI, df$BP, method = "pearson")
cor.test(df$BMI, df$BP, method = "spearman")

##### ICC #####
# warning : aggrement is not correlation : 

# Packages for ICC
library(irr)

# Simulate Teacher A's grades
set.seed(42)
teacherA <- round(runif(30, min = 0, max = 20), 1)
# Teacher B always adds +2 points
teacherB <- teacherA + 2

# Combine into a data frame
grades <- data.frame(teacherA, teacherB)

# 1. Pearson Correlation
correlation <- cor(
  grades$teacherA, grades$teacherB, method = "pearson"
)
print(paste("Pearson correlation:", correlation))
# [1] "Pearson correlation: 1"

# 2. Intraclass Correlation Coefficient
icc_result <- irr::icc(
  grades, model = "twoway", unit = "single", type = "agreement"
)
print(icc_result)
# Single Score Intraclass Correlation
# 
# Model: twoway 
# Type : agreement 
# 
# Subjects = 30 
# Raters = 2 
# ICC(A,1) = 0.945
# 
# F-Test, H0: r0 = 0 ; H1: r0 > 0 
# F(29,1) = Inf , p = 0 
# 
# 95%-Confidence Interval for ICC Population Values:
#   0.017 < ICC < 0.99

# 3. Visualization
## basic graph
plot(
  grades$teacherA, grades$teacherB, 
  main = "Teacher A vs Teacher B",
  xlab = "Teacher A", ylab = "Teacher B", pch = 19, col = "blue"
)
abline(0, 1, col = "red", lwd = 2, lty = 2) # Perfect agreement line
## ggplot2 graph
library(ggplot2)
ggplot(grades, aes(x = teacherA, y = teacherB)) +
  geom_point(color = "blue", size = 2, shape = 21) +
  geom_abline(
    intercept = 0, slope = 1, color = "firebrick",
    linetype = "dashed", size = 1
  ) +
  labs(
    title = "Perfect Correlation but Imperfect ICC", 
    subtitle = "Teacher A vs Teacher B",
    x = "Teacher A",
    y = "Teacher B"
  ) +
  theme_minimal(base_size = 14)

## alt pkg for ICC : psych
library(psych)
icc_result <- psych::ICC(
  x = grades, missing = TRUE, alpha = 0.05
) ## but in this example, return an error : 
# Erreur dans eval_f(x, ...) : Downdated VtV is not positive definite
## That error with ICC() in psych happens when you only provide two 
## raters and the data is perfectly collinear 
#" (like my +2 constant shift). 
## The covariance matrix becomes singular (not positive definite).

# fix this example, for demonstration : 
library(psych)
teacherB <- teacherA + 2 + rnorm(30, 0, 0.01) # almost perfect +2, tiny noise
grades <- data.frame(teacherA, teacherB)
psych::ICC(grades, missing = TRUE, alpha = 0.05)


#### 4.5 – PAIRED TESTING ####

##### Paired t-test #####

t.test(tumor_exp, normal_exp, paired = TRUE)
t.test(Before, After, paired = TRUE)

## fig 
# ggpaired(ToothGrowth, x = "supp", y = "len",
#          color = "supp", line.color = "gray", line.size = 0.4,
#          palette = "jco")+
#   stat_compare_means(paired = TRUE)
## https://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/

##### Wilcoxon signed rank test #####

wilcox.test(Before, After, paired = TRUE)

##### Anova repeated measures #####

aov_result <- aov(cholesterol ~ time + Error(patient/time), data = data)
summary(aov_result)

##### Friedman test #####

friedman.test(pain ~ treatment | patient, data = data)


##### Quade #####

# Quade test
quade.test(y ~ group | block, data = df)
# Skillings-Mack test
library(Skillings.Mack)
skillings.mack(y ~ group | subject, data = df)

#### 4.6 – MULTIPLE TESTING CORRECTIONS ####

p.adjust(p_values, method = "bonferroni")
p.adjust(p_values, method = "holm")
p.adjust(p_values, method = "BH")
 
#### step 5 ####

#### Regression ####

reg <- lm(Y ~ X1 + X2 + X3)
reg ## returns the estimated coefficients beta
## To extract a specific coefficient estimate, e.g beta_2:
reg$coef[3]
# To obtain the fitted (predicted mean) values of Y, based on the observed X1, X2 et X3 
predict(reg) # (ou fitted(reg))
# To predict the mean value of Y for a new observation (X1, X2, X3) = (1.2, 2.2, 6):
predict(reg, data.frame(X1 = 1.2, X2 = 2.2, X3 = 6))
# If beta_0 has no meaningful interpretation in your model, remove it with  "-1": 
reg = lm(Y ~ X1 + X2 + X3 - 1)

# R2
summary(reg)
# conf int : 
confint(reg, level = 0.95)
predict(reg, data.frame(X1 = 1.2, X2 = 2.2, X3 = 6), interval = "confidence")

# validation 

par(mfrow = c(2, 2)); plot(reg, 1:4)

plot(w) # or pairs(w) 
#  par exemple : pairs(~ Y + X1 + X4).
car::scatterplotMatrix(w)

# transformation ex 
 reg = lm(log(Y) ~ sqrt(X1) + exp(X2) + I(X3^4))

# check residuals 
residuals(reg)
rstandard(reg)

acf(residuals(reg))
pacf(residuals(reg))

library(lawstat)
Box.test(residuals(reg), type = "Ljung")

# interaction 

# To include an interaction term between X1 and X2, you can write:
reg = lm(Y ~ X1 * X2, data = df)
# This automatically expands to:
reg = lm(Y ~ X1 + X2 + X1:X2, data = df)
# If you only want the interaction term without main effects (rare, but sometimes used), use:
lm(Y ~ X1:X2, data = df)

interaction.plot(df$X1, df$X2, df$Y)

# multicollinerarity 

cor(df[, c("X1", "X2", "X3")])

library(car)
vif(reg)

#### 5.3 GLM ####

res <- glm(y ~ x, data = df)
summary(res)
plot(res)

## glm gaussian
# res <- glm(x ~ factor1 + factor2 + ..., family = gaussian)
res <- glm(y ~ factor1 + factor2, family = gaussian(link = "identity"))
summary(res)
plot(res)
anova(res, test = "F")

# glm poisson 
# res <- glm(y_count~factor1+factor2, family = poisson)
res <- glm(y_count ~ factor1+factor2, family = poisson(link = "log"))
summary(res)
plot(res)
anova(res, test = "Chisq")

res <- glm(y_count ~ factor1 + factor2, family = quasipoisson(link = "log"))

# glm binom
curve(log(x/(1-x)), 0, 1)

res <- glm(cbind(y1, y2) ~ factor1 + factor2, family = binomial(link = "logit"))
summary(res)
plot(res)
anova(res, test = "Chisq")

res <- glm(y_binary ~ factor1 + factor2, family = binomial(link = "logit"))
summary(res)
plot(res)
anova(res, test = "Chisq")

#### 5.4  Model Selection ####

## Metrics
mean((residuals(reg))^2) # MSE
AIC(model1, model2, model3)
BIC(model1, model2, model3)

## CV 80/20
set.seed(123)
train_index <- sample(1:nrow(df), 0.8 * nrow(df))
train <- df[train_index, ]
test <- df[-train_index, ]

model <- lm(Y ~ X1 + X2, data = train)
pred <- predict(model, newdata = test)
mean((test$Y - pred)^2) # MSE on test data

## K-fold CV 
library(caret)
train_control <- trainControl(method = "cv", number = 10)
train(Y ~ X1 + X2, data = df, method = "lm", trControl = train_control)

## LOOCV
train_control <- trainControl(method = "LOOCV")
train(Y ~ X1 + X2, data = df, method = "lm", trControl = train_control)

