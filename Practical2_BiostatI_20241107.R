# Practical 2 : 2025/11/07 

#### envir ####
getwd()
# setwd("...")
list.files()

#### Read ####

# tab1a <- read.table("yourfile.txt", header = TRUE)
tab1a <- read.table("iris_head80.txt", header = TRUE)

str(tab1a)
tab1a$Species <- as.factor(tab1a$Species)
# tab1a$Species <- factor(tab1a$Species, levels, labels...)
tab1a$Sepal.Length <- as.numeric(tab1a$Sepal.Length)
tab1a$Sepal.Width <- as.numeric(tab1a$Sepal.Width)
tab1a$Petal.Length <- as.numeric(tab1a$Petal.Length)
tab1a$Petal.Width <- as.numeric(tab1a$Petal.Width)
summary(tab1a)

# tab2a <- read.csv("english_format.csv")
# tab2b <- read.csv2("french_format.csv")
tab2a <- read.csv("iris_tail80_en.csv")
tab2b <- read.csv2("iris_tail80_fr.csv")

head(tab2a)
head(tab2b)

# install.packages("waldo")
library(waldo)
# compare(c("a", "b", "c"), c("a", "b"))
compare(tab2a, tab2b)

# install.packages("readxl")
library(readxl)
# tab3a <- read_excel("dataset.xlsx")
# tab3b <- read_excel("dataset.xlsx", sheet = 2)
# tab3b <- read_excel("dataset.xlsx", sheet = "put the name") # safer

tab3a <- read_excel("mycars.xlsx")
str(tab3a)
tab3b <- read_excel("mycars.xlsx", sheet = 2)
readxl::excel_sheets("mycars.xlsx")
tab3b <- read_excel("mycars.xlsx", sheet = "Sheet2")
str(tab3b)
compare(tab3a, tab3b)

## mycars.xlsx actually come from mtcars
?datasets::mtcars

str(tab3b)
table(tab3b$am)
tab3b$cyl <- as.factor(tab3b$cyl)
tab3b$vs <- as.factor(tab3b$vs)
tab3b$am <- as.factor(tab3b$am)

# tab4a <- read_excel("experiment.xlsx", range = "F20:I25")
tab4a <- read_excel("Luciferase_exp.xlsx", range = "G8:K11")
tab4a

tab4b <- reshape2::melt(
  tab4a,
  id.vars = c("condi"), ## what to keep as it is
  measure.vars = grep("variant", names(tab4a), value = TRUE), # what to put longer
  value.name = "luciferase_value"
)
tab4b

dt_4b <- data.table::melt(
  data = data.table::setDT(tab4a),
  id.vars = c("condi"), ## what to keep as it is
  measure.vars = grep("variant", names(tab4a), value = TRUE), # what to put longer
  value.name = "luciferase_value"
)
dt_4b

compare(tab4b, dt_4b)


#### Describe ####

dim(tab3b)

# summary(tab$YYY)
# mean(tab$YYY)
# median(tab$YYY)
# sd(tab$YYY)
# quantile(tab$YYY)
# sum(is.na(tab$YYY))
# length(tab$YYY)   # n
summary(tab3b$disp)
mean(tab3b$disp)
median(tab3b$disp)
sd(tab3b$disp)
quantile(tab3b$disp)
sum(is.na(tab3b$disp))
length(tab3b$disp)   # n

# install.packages("tidyverse")
library(tidyverse)

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


# hist(tab$YYY)
# plot(density(tab$YYY))
# boxplot(tab$YYY)
hist(tab3b$disp)
plot(density(tab3b$disp))
boxplot(tab3b$disp)

# tab <- tab3b
num_vars <- names(tab)[sapply(tab, is.numeric)]

pdf("rplot.pdf") 
for (v in num_vars) {
  cat(">>> " , v)
  par(mfrow = c(1, 3)) # 3 graphiques côte à côte
  hist(tab[[v]], main = paste("Histogram of", v), xlab = v)
  plot(density(tab[[v]], na.rm = TRUE), main = paste("Density of", v), xlab = v)
  boxplot(tab[[v]], main = paste("Boxplot of", v), horizontal = TRUE)
  cat(" end ")
}
dev.off()

par(mfrow = c(1, 1)) # set back the config


# table(tab$XXX)
# prop.table(table(tab$XXX)) * 100
# sum(is.na(tab$XXX))
# length(tab$YYY)   # n
# table(tab$XXX)
table(tab3b$cyl)
prop.table(table(tab3b$cyl)) * 100
sum(is.na(tab3b$cyl))
length(tab3b$cyl)   # n

barplot(table(tab3b$cyl))

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


#### Report ####

install.packages("writexl")
library(writexl)

df_quant <- tab3b %>%
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

df_qual <- tab3b %>%
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



write_xlsx(list(
  Quantitative = df_quant,
  Qualitative = df_qual
), "descriptive_report.xlsx")



#### Extra ####

install.packages("data.table")
install.packages("ggplot2")
# Fast import with fread
library(data.table)
tabDT <- fread("dataset.csv")
str(tabDT)
# Nice graphics with ggplot2
library(ggplot2)
ggplot(tabDT, aes(x = YYY)) + geom_histogram()
ggplot(tabDT, aes(x = XXX, y = YYY)) + geom_boxplot()
ggplot(tabDT, aes(x = XXX)) + geom_bar()

