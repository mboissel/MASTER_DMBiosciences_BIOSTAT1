## Cours M1 Biostat  ##
## La Catho Univ ##

## Mathilde Boissel 

#### Set up envir ####
library(data.table)
library(ggplot2)

#### 1.1 – ELEMENTS OF CALCULUS ####

##### 1.1 - D	In practice: R #####

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

fit <- aov(Cholesterol ~ Diet, data = df)
summary(fit)

#####	Wilcoxon-Mann-Whitney test  #####

wilcox.test(VPSAH ~ Group, data = df)

##### Chi2 #####

tbl <- table(df$Smoking, df$Gender)
chisq.test(tbl)

##### fisher exact #####

tbl <- table(df$advers_event, df$ttt)
fisher.test(tbl)   # valid for small or unbalanced samples

##### Correlation test #####

cor.test(df$BMI, df$BP, method = "pearson")
cor.test(df$BMI, df$BP, method = "spearman")

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
 




#### Practical session ####

# https://bookdown.org/ael/rexplor/chap1.html

