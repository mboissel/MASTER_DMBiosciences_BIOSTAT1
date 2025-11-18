# Practical 1 : 2025/10/24 

#### Envir ####

R.version
.libPaths()
getwd()
# setwd("your/folder")

#### calculus ####

2 + 3 # addition
8 - 12 # subtraction
5 ^ 2        # power (5 squared)
10 / 3       # division
14 * 25      # multiply 
sqrt(16)     # square root
# usual operator priority : 
exp((7 * 3 + 12 / 2 - 7^2 + sqrt(4) - abs(log(0.05))) / (2^3 + 1))


# store
a <- 5                       # assign numeric value
b <- "hello world"           # character string
c <- TRUE                    # logical (TRUE/FALSE)

# show 
a
print(a)                     # show content of object

# class = type
class(a); class(b); class(c) # check type of each object

# use the object
a^2
sqrt(a)

# modify the obj 
a <- 1
a <- a + 5
a

# copy the obj
y <- a
y
a <- 101
cat("y =", y, "and a =", a)

## Maths function available
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


# example of error  
5a + 5
# example of warning
log(-1)
# NaN = Not a Number
2 * (2 + log(-1))


#### dataset iris ####

data(iris)           # load dataset (optional)
?iris                # help
head(iris)           # show first rows
str(iris)            # structure of dataset
summary(iris)        # summary statistics
class(iris)          # type of object (data.frame)


# Select the first 10 rows
iris[1:10, ]
# Select sepal length column
iris$Sepal.Length
# Filter rows where Species = “setosa”
iris_setosa <- subset(iris, Species == "setosa")
# Or using indexing
iris[iris$Species == "setosa", ]
iris[iris$Species %in% "setosa", ]
mydf <- iris[iris$Species %in% "setosa", ]


#### save object ####

# Save an object to a file
saveRDS(mydf, file = "01-mydf.rds")
rm(mydf) # remove or hard clean : Session > Restart R
# Restore the object
mydf <- readRDS(file = "01-mydf.rds")

# Save multiple objects
save(mydf, a, file = "data.RData")
# Load all objects again
load("data.RData")
ls()

# Save the entire workspace
save.image(file = "my_work_space.RData")
# Restore your workspace
load("my_work_space.RData")
ls()





#### Les fonctions prédéfinies ####

# https://bookdown.org/ael/rexplor/chap1bis.html#chap1.7
# usage

# valeur par default 

# aide 
# help(<nom.fonct>)
# ou
# ?<nom.fonct>