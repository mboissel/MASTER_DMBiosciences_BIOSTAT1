#### Global cmd ####


##### Help #####

# get help in R via :
# help(<nom.fonct>)
# or the question mark
# ?<nom.fonct>

# also ok to get help via google, many documentation on line.


##### Setup your environment #####

R.version # check R version 
.libPaths() # check where your pkg are installed 
getwd() # check where "is living" = working directory
# setwd("path/your/folder") # set your working directory = project folder path. 
# setwd("path\your\folder") # warning about / or \ in windows, mac or linux...
# to avoid the error on / or \ use :
# file.path("path", "to", "your", "folder", "myfile") # will select the good separator.

list.files() # list file available in your folder 

dir.create(
  "myresults", # create the folder myresults
  showWarnings = FALSE, recursive = TRUE
)

##### Calculus & basic cmd #####

2 + 3        # addition
8 - 12       # subtraction
5 ^ 2        # power (5 squared)
10 / 3       # division
14 * 25      # multiply 
sqrt(16)     # square root
# usual operator priority : 
exp((7 * 3 + 12 / 2 - 7^2 + sqrt(4) - abs(log(0.05))) / (2^3 + 1))


# store i.e. save your date with "<-" 
a <- 5                       # assign numeric value
b <- "hello world"           # character string
c <- TRUE                    # logical (TRUE/FALSE)

# show it, type the name of the object in the consol. 
a
print(a)                     # same = show content of object

# check class = type
class(a); class(b); class(c) # check type of each object
# is function ...
is.data.frame()
is.numeric()
is.factor() 

# use the object in computation...
a^2
sqrt(a)

# modify the object by re-saving the modification
a <- 1
a <- a + 5
a # now a is 6, no more 1.

# copy the object
y <- a # y is a copy of a.
y
a <- 101 # a is modify, not y.
cat("y =", y, "and a =", a)

## Maths function available in R
sqrt(4)      # racine carrée de 4
log(4)       # logarithme népérien (base e) de 4
exp(4)       # exponentielle de 4
abs(-4)      # valeur absolue de -4
sin(4)       # sinus de 4
factorial(4) # factorielle de 4, càd. 4!
choose(4, 2) # combinaisons de 2 éléments parmi 4, càd. 4!/2!(4-2)!
?choose # read the documentation in help panel

# some objets/values ever defined  
NULL  # an empty object: « nothing » !
Inf   # + infinity
-Inf  # - infinity
pi    # pi = 3.14...
letters     # alphabet
LETTERS
month.name  

##### Play with iris #####

data(iris)           # load dataset (optional), 
# this one is nativly available in R.
?iris                # help # documation in {datasets} r pkg
?datasets::iris      # same, you can specify pkg:: before. 
head(iris)           # show first rows
head(iris, 10)       # show first 10 rows
str(iris)            # structure of dataset
summary(iris)        # summary statistics
class(iris)          # type of object (data.frame)

names(iris) # get the name
dput(names(iris)) # get the name in vector format, easier to copy and paste.   


# a dataset is a collection of i line and j column, 
# select the line i before the "," 
# and the column j after the ","

# Select the first 10 rows
iris[1:10, ]
# Select the first column, via its number
iris[, 1]

# select line and column in same time ! 
dfiris2 <- iris[
  iris$Species %in% c("setosa", "virginica"), # lines
  c("Petal.Length", "Species") # columns 
]

# Select sepal length column via its name
iris$Sepal.Length
# Filter rows where Species = “setosa”
iris_setosa <- subset(iris, Species == "setosa")
# Or using indexing
iris[iris$Species == "setosa", ]
iris[iris$Species %in% "setosa", ]

# store your selection
mydf <- iris[iris$Species %in% "setosa", ]

# select data without NA = complet
df_withoutna <- dt_select[complete.cases(dt_select), ]
# same :
# df_withoutna <- dt_select[
#   !is.na(dt_select$bill_length_mm) & !is.na(dt_select$sex),
# ]
# also see na.omit()

##### Play with mtcars #####

?datasets::mtcars

##### Play with penguins #####

# install.packages("palmerpenguins")
library(palmerpenguins)
data("penguins")
dput(names(penguins)) # easier to copy and paste names

##### Play with trees #####
data(datasets::trees)
head(trees)

##### Save object #####

## as R obj ##

# Save an object to a file
saveRDS(mydf, file = "01-mydf.rds")
rm(mydf) # remove or hard clean : Session > Restart R
# Restore the object
mydf <- readRDS(file = "01-mydf.rds")

## Save multiple objects
# save(mydf, a, file = "data.RData")
# # Load all objects again
# load("data.RData")
# ls() # list anything in you envir.

## Save the entire workspace
# save.image(file = "my_work_space.RData")
# # Restore your workspace
# load("my_work_space.RData")
# ls()

## as excel file ##

# install.packages("writexl")
library(writexl)

write_xlsx(
  list(
    sheet1 = one_dataframe,
    sheet2 = anotherone_dataframe
  ),
  "myexcelsaved.xlsx"
)



#### Read data ####

## read flat text file (txt) ##
tab1a <- read.table("yourfile.txt", header = TRUE)

## read flat text file csv  ##
## comma "," separator or ";"
tab2a <- read.csv("english_format.csv")
tab2b <- read.csv2("french_format.csv")

## read Excel ##
# install.packages("readxl")
library(readxl)
tab3a <- read_excel("dataset.xlsx")
tab3b <- read_excel("dataset.xlsx", sheet = 2)
tab3b <- read_excel("dataset.xlsx", sheet = "put the name") # safer
readxl::excel_sheets("dataset.xlsx") # list sheet names

## read a zone :
tab4a <- read_excel("experiment.xlsx", range = "F20:I25")


#### Format ####

# format factor and numeric

tab1a$Species <- as.factor(tab1a$Species)
tab1a$Species <- factor(tab1a$Species, levels, labels...)
# so we can put specific order
nlevels() # check number of levels  
# either with 
tab1a <- droplevels(tab1a) # remove unused levels


tab1a$Sepal.Length <- as.numeric(tab1a$Sepal.Length)

# format a complete dataset : 
pommes <- data.frame(
  rendement = c(48,46,52,50,52,50,49,49,53,51,55,57), # numeric
  variete = factor(rep(c("Golden","Delicious","Jonagold"), rep(4,3)))
)
# test 
variete_facteur <- gl(n = 3, k = 4, label = c("Golden","Delicious","Jonagold"))
help("gl")

#### Describe ####

## a numeric YYY variable ##
summary(tab$YYY)
mean(tab$YYY)
median(tab$YYY)
sd(tab$YYY)
quantile(tab$YYY)
quantile(tab$YYY, 0.25) # only Q1 wanted = 25%
sum(is.na(tab$YYY))
length(tab$YYY)  # n

nrow(tab) # nlines

# format a summary 
# broom::tidy(summary(tab$YYY))

# round to not keep many digits 
round(0.25848551, 3) # 0.258

## a factor XXX factor ##

table(tab$XXX)
table(tab$XXX, useNA = "ifany")
table(tab$XXX, useNA = "always")
prop.table(table(tab$XXX)) * 100
sum(is.na(tab$XXX))
length(tab$XXX)   # n

## summary by group ofvariete
means <- aggregate(rendement ~ variete, data = pommes, FUN = mean)
mean_tapply <- tapply(
  X = pommes$rendement, INDEX = pommes$variete, FUN = mean
)

#### Tidyverse ####

# install.packages("tidyverse")
library(tidyverse)

# %>% is call pipe, it helps to pass the result of a line into the one under.
  
# describe all numerics
tab3b %>%
  select(where(is.numeric)) %>%
  pivot_longer(everything()) %>%
  group_by(name) %>%
  summarise(
    n = sum(!is.na(value)),
    na = sum(is.na(value)),
    mean = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    q1 = quantile(value, 0.25, na.rm = TRUE),
    q3 = quantile(value, 0.75, na.rm = TRUE)
  )


# describe all factors

tab3b %>%
  select(where(is.factor)) %>%
  pivot_longer(
    everything(), names_to = "variable", values_to = "value"
  ) %>%
  group_by(variable, value) %>%
  summarise(
    count    = n(),
    pct      = n() / sum(n()) * 100,
    .groups  = "drop_last"
  ) %>%
  arrange(variable, value)


#### Visualization ####

## of continuous vars ##

hist(tab$YYY)
plot(density(tab$YYY))
boxplot(tab$YYY)

# boxplot with mean added as points
means <- aggregate(rendement ~ variete, data = pommes, FUN = mean)
boxplot(
  formula = rendement ~ variete, data = pommes,
  main = "Boite à moustaches"
)
points(1:3, means$rendement, col = "red")

## normal qqplot
qqnorm(pommes$rendement)
qqline(pommes$rendement)


## of factorial vars ##

barplot(table(tab3b$cyl))


## pair plot, 2 x 2
pairs(trees)
 

## save figures ##

# automatic loop on all variable numeric :
num_vars <- names(tab)[sapply(tab, is.numeric)]
pdf("rplot.pdf") # open the pdf
for (v in num_vars) {
  cat(">>> " , v)
  par(mfrow = c(1, 3)) # 3 graphiques côte à côte
  hist(tab[[v]], main = paste("Histogram of", v), xlab = v)
  plot(density(tab[[v]], na.rm = TRUE), main = paste("Density of", v), xlab = v)
  boxplot(tab[[v]], main = paste("Boxplot of", v), horizontal = TRUE)
  cat(" end ")
}
dev.off() # close the pdf
par(mfrow = c(1, 1)) # set back the config

# open a jpeg
jpeg(file.path(path_folder, "density_basic.jpg"))
plot(
  density(df_withoutna$bill_length_mm),
  main = "Density bill_length_mm"
)
dev.off() # close the jpeg


# open a png
png(file.path(path_folder, "density_basic.png"))
plot(
  density(df_withoutna$bill_length_mm),
  main = "Density bill_length_mm"
)
dev.off() # close



#### Statistical test ####

## normality ##

shapiro.test() 
# shapiro.test(dfiris2$Petal.Length[dfiris2$Species %in% "setosa"])
# shapiro.test(dfiris2$Petal.Length[dfiris2$Species %in% "virginica"])
# shapiro.test(df_withoutna$bill_length_mm[df_withoutna$sex %in% "female"])
# shapiro.test(df_withoutna$bill_length_mm[df_withoutna$sex %in% "male"])


## variance equal ? ##

bartlett.test(Petal.Length ~ Species, data = dfiris2)

library(car)
car::leveneTest()
# car::leveneTest(Petal.Length ~ Species, data = dfiris2)
# p = 8.871e-09 ***
# p > 0.05 → variances égales → Student's t-test
# p < 0.05 → variances inégales → Welch's t-test

## t-test ##

t.test(
  Petal.Length ~ Species,
  data = dfiris2, var.equal = FALSE
)

## wilcoxon test ##
wilcox.test(
  bill_length_mm ~ sex,
  data = df_withoutna
)

# mytest <- wilcox.test(bill_length_mm ~ sex, data = df_withoutna)
# mytab <-  broom::tidy(mytest)

#### gtsummary pkg ####

# describ and test, all in once. 

# install.packages("gtsummary")
library(gtsummary)

gtsummary::tbl_summary(
  data = penguins[, c("sex", "bill_length_mm")]
)

gtsummary::tbl_summary(
  data = penguins[, 
       c("species", "sex", "bill_length_mm", "bill_depth_mm")
  ],
  statistic = list(
    all_continuous() ~
      "Med = {median} (Q1 = {p25}, Q3 = {p75}) ; Mean = {mean}", 
    all_categorical() ~ 
      "n = {n} ({p} %)"
  )
)


gtsummary::tbl_summary(
  data = penguins[, c("sex", "bill_length_mm")],
  by = sex
)


table_res <- penguins %>% 
  tbl_summary(
    include = c(bill_length_mm), # select
    by = sex, # split table by group
    missing = "no" # don't list missing data separately, p.e.
  ) %>% 
  add_n() %>%  # add column with total number of non-missing observations
  add_p() %>%  # test for a difference between groups
  modify_header(label = "**My Variable**") %>%  # update the column header
  bold_labels()

# show it :
table_res


#### LM Linear Model ####

?stats::aov
?stats::lm
?stats::anova

##### ANOVA #####

## anova test, special case of LM ##
my_anova <- aov(formula = rendement ~ variete, data = pommes)
my_anova
summary(my_anova)
# diag plots 
par(mfrow=c(2, 2));
plot(my_anova)
par(mfrow=c(1,1))

## Test post-hoc, after anova
?TukeyHSD
TukeyHSD(x = my_anova)

pairwise.t.test(x = pommes$rendement, g = pommes$variete, p.adj = "bonf")

pairwise.t.test(x = pommes$rendement, g = pommes$variete, p.adj = "holm")

##### LM #####

# Volume = β_0 + β_1 × Girth + β_2 × Height + ε
reg <- lm(
  formula = Volume ~ Girth + Height, 
  data = trees
)
summary(reg)

# get coef = estimates
coef(reg) # returns coefficents
coef(reg)[3]  # returns directly Beta2 
# or
res <- summary(reg) # and  
res$coefficients

## get confidence interval of coefs
confint(reg, level = 0.95)
# or more precisly : 
confint(reg, level = 0.95)[3, ]

## R2 and adj R2
res <- summary(reg)
res$r.squared # and 
res$adj.r.squared

## Global F test < 0.05
# H0 : "all betas = 0" vs H1 : “At least one coefficient not null”.
# here we reject H0, it means the model is relevant 

## prediction 
predict(reg, data.frame(Girth = 8.3, Height = 70))
predict(
  reg, 
  data.frame(Girth = 8.3, Height = 70),
  interval = "confidence"
)

## Graph diag
par(mfrow = c(2, 2))
plot(reg)

## To focus on the 1st graph, if it seems to have a problem : 
par(mfrow = c(1,1));
plot(reg, 1)

# to check normality : QQplot 
par(mfrow = c(1,1)); plot(reg, 2)

# ## To go further about independence of errors : 
# par(mfrow = c(1, 2))
# acf(e)
# pacf(e)
# par(mfrow = c(1, 1)) ## get back params graphic 


## multi colinearity :the Klein's rule
cor(x = trees$Girth, y = trees$Height)^2
# to compare with
summary(reg)$r.squared

cor(trees) # compute all correlation 2 x 2


## Extrem values
## cook d
cook <- cooks.distance(reg)
cook[cook>1]

# to check cook d
par(mfrow = c(1,1)); plot(reg, 4)

# hat info
summary(influence.measures(reg))


##### Compare 2 LM models #####

reg1 = lm(Volume ~ Girth + Height + X3 + X4, data = mydata)
reg2 = lm(Volume ~ Girth + Height, data = mydata)
anova(reg1, reg2)

## AIC BIC

AIC(reg1) ; BIC(reg1) ;
AIC(reg2) ; BIC(reg2) ;


#### GLM ####

##### logistic regression #####

# binary information CC explain by weight
# boxplot(T2Ddata$weight ~ T2Ddata$CC) # better to get the labels

## glm binomial
reg <- glm(
  formula = CC ~ weight,
  data = T2Ddata,
  family = binomial(link = "logit")
)
summary(reg)

## 
OR <- exp(coef(reg)) # store it 
OR # show it
OR > 1 # check it > 1

## Diag
plot(reg)
plot(reg, 1) # horizontal line, good
plot(reg, 2)
plot(reg, 3)
plot(reg, 4)

