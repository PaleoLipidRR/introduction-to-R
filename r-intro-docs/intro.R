# Howdy! Welcome to introduction to R

3 + 4

2349078/54

x <- 6

# new R object named "juice"
juice <- 4+4

#overwrite x to a series of numbers
x <- c(1,7,9)

apples <- 5
oranges <- 9
fruit <- apples + oranges

# View(x) ---> View() is to display x
# head(x) ---> is to show head of the data array

#  V = 4/3 π r³ volume equation
d <- 20
Vol <- ((4/3)*pi*((d/2)^3))

# to nest our parenthesis (adding parenthesis on both sides for the current text), 
# highlight the text and press "("














# Data classes and data frames

# ?class() # to check data class
# ?str() # to compactly display the internal structure of an R object,

Vol
juice <- 4+4
class(juice)
class(Vol)
class(4) # numeric is float (numbers with decimals)
class(8.3)
class(4L) # "L" to define that the values are integers (float in R is called "numeric" --- the data with decimals)
str(4L)

x <- c("1","7","9") # character in R is "str" in python
class(x)

numeric_vector <- c(1, 10, 49, 5)
character_vector <- c("a", "b", "c", "d")
str(numeric_vector)
str(character_vector)

names(numeric_vector) <- character_vector # named numeric array is the same as dict in python

# [rows, columns]. Using [] to index the data frame
numeric_vector[4]
numeric_vector
class(numeric_vector)

# Question --> is names a column names or just each element name/ how to display only column names

# Lists
list_example <- list(2^4, "cabbage", TRUE, 1+5i)
list_example2 <- list(2E2, "lettuce", VALUES = 1:25, VALUES_BINARY = FALSE)

str(list_example)
str(list_example2)

list_example2$VALUES # $ is used to select variables inside the list (equivalent to "." in python and MATLAB)
typeof(list_example2$VALUES_BINARY)

# Matrices
matrix(1:9, byrow = FALSE, nrow = 3)

q <- c(460, 314)
r <- c(290, 24)
w <- c(390, 165)

c(q, r, w)
# Put the data into a matrix
mydata <- matrix(c(q, r, w), byrow = TRUE, nrow = 3)


region <- c("one", "two")
category <- c("A", "B", "C")
rownames(mydata) <- category
colnames(mydata) <- region

totals <-rowSums(mydata)

mydata_better <- cbind(mydata, totals)
# View(mydata_better)


# Data frame
greetings <- c("hey", "hi", "howdy", "hello", "morning")
n <- c(99,15,324,54,12)
df <- data.frame(greetings, n)

df[3,] # to select the third row we need to put "," in the []
df[2] # by default, R will select column if we put in one index number

addition <- c("afternoon", 18)
?rbind() # add a new row
df <- rbind(df, addition)
head(df)

names(df)
colnames(df)[1:2] <- c("Greeting","Observed") # this is to change col names
rownames(df)[1:2] <- c("heys","hi2")

# Activity 4.4

df$Observed <- as.numeric(df$Observed)
class(df$Observed)

# Install R packages

# install.packages('vegan')
# install.packages("tidyverse")

library(vegan)
library(tidyverse)

# data import
greetings

# Path and working directories
getwd() # to get current working directory
setwd("/Users/ratta/Documents/GitHub/introduction-to-R/r-intro-docs/") # to set working directory

#Import Pizza
pizza <- read.csv("pizza.csv")
data("starwars")
# View(starwars)

asv_table <- read_table(file = "asv-table.txt")
asv_table2 <- read.table(file= "asv-table.txt",header = TRUE)

metadata <- read_delim(file = "metadata.txt")

#Export data
head(pizza)
# View(pizza)

?write_delim()
write_delim(pizza, file="pizza-plain.txt")
write_delim(pizza, file="pizza-tab.txt", delim = "\t")
write_delim(pizza, file="pizza-nocols.txt",quote = "all", col_names = FALSE)

?write_csv() # from tidyverse
?write.csv() # base R
write.csv(pizza, file = "pizza.csv", row.names = TRUE)
write_csv(pizza, file = "pizza2.csv")

# Tidyverse
head(pizza)


pizza %>% select(type)
types_of_pizza <- pizza %>% select(type)

types_of_pizza_char <- as.character(pizza %>% select(type))

# Filter()
# %>% this is called "pipe" which is used to string several commands that we wanted to perform on the R object
pizza %>% 
  select(style = type, everything(), -X) %>% 
  filter(num_toppings >= 2) %>% 
  filter(red_sauce != TRUE)

pizza %>% 
  select(style = type, everything(), -X) %>% ## the first position is to rename "type" to "style"; we tell R that we want to select "style" column (which does not exist yeat in our "pizza" data frame) that is the same thing as "type"
  filter(num_toppings >= 2 & red_sauce == TRUE)

pizza %>% 
  select(style = type, everything(), -X) %>% 
  filter(num_toppings >= 2 | red_sauce != TRUE)

# mutate() # Add a column
pizza$cost

totalcost <- pizza %>% 
  select(style = type, everything(), -X) %>% 
  mutate(TOTAL_COST = (cost * num_order))

sum(totalcost$TOTAL_COST)
write.csv(totalcost, file = "totalcost.csv")

# wide vs. long 
df_wide <- data.frame(Meal = c('Morning', 'Mid-morning', 'Lunch', 'Mid-afternoon', 'Dinner'),
                      Weight = c(9, 0, 12, 10, 6),
                      Dishes = c(3, 0, 1, 1, 4))
df_wide

df_long <- df_wide %>% pivot_longer(cols = c('Weight','Dishes'),names_to = "Attribute", values_to = "Value")

df_long

## Example from larger data table

names(asv_table2)

asv_long <- asv_table %>% pivot_longer(cols = 3:20)
asv_long <- asv_table %>% pivot_longer(starts_with("GordaRidge"))

head(asv_long)

# 5.1 Filtering data
## Subset ASV table so that ASV has greater than or equal to 100 sequences.
asv_long_100plus <- asv_long %>% 
  filter(value >=100)

dim(asv_long)
dim(asv_long_100plus)

# Star wars
head(starwars)
names(starwars)

## select data that for "Tatooine"
tatooine <- starwars %>% 
  filter(homeworld == "Tatooine") # select only those from Tatooine
head(tatooine)

unique(starwars$homeworld)


# Create a table of droids that are equal to or greater than 96 inches in height.
droids_tall <- starwars %>%
  filter(species == "Droid" & height >= 96)
droids_tall

# Summarizing data
# mutate() == adds a column
# summarize() == subsets, applys function

## What are the average masses of humans versus ewoks?
# group_by() summaries

starwars %>% 
  filter(species == "Human" | species == "Ewok") %>% 
  filter(!is.na(mass)) %>% 
  group_by(species) %>% 
    summarise(MEAN_MASS = mean(mass))

starwars %>% 
  filter(species == "Human" | species == "Ewok") %>% 
  filter(!is.na(mass)) %>% 
  group_by(species, homeworld) %>% 
    summarise(MEAN_MASS = mean(mass),
              MIN = min(mass),
              MAX = max(mass),
              )

# 3.0.1 
## Add a column that classifies tall vs. short based on the height of each species
height_class <- starwars %>% 
  mutate(HEIGHT_BIN = case_when(
    height > 100 ~ "tall", # "~" here means assign the rows where height >100 as "tall"
    height <= 100 ~ "short"
  ))

hist(starwars$height)

## use median_heights of each species

tmp <- starwars %>% 
  filter(!is.na(height)) %>% 
  group_by(species) %>% 
  mutate(MEDIAN_HEIGHT_SPECIES = median(height)) %>% 
  mutate(HEIGHT_BIN = case_when(
    height > MEDIAN_HEIGHT_SPECIES ~ "tall",
    height <= MEDIAN_HEIGHT_SPECIES ~ "short"
  ))


# Section 4
## Isolate non-droids on Alderaan, Naboo, Endor, Kamino, and Coruscant. Based on their reported sex, is there a relationship between their height, planet, and species?

sort(unique(starwars$homeworld)) #unique listed the elements in an order that represents which data come first in the original dataframe; adding sort() will sort the unique list alphabetically

planets <- c("Alderaan", "Naboo", "Endor", "Kamino", "Coruscant")

nondroids <- starwars %>% 
  filter(homeworld %in% planets) %>% 
  filter(species != "Droid")

ggplot(nondroids, aes(x = species, y = height, fill = sex)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(homeworld ~ .)
# try with "fill = sex" and without

# stat - "count" works when you don't have a y value. it counts occurrences

# Microbiology example
head(asv_long)
head(metadata)


?left_join()
# Combine with metadata
asv_withmetadata <- asv_long %>% 
  filter(value > 0) %>%
  left_join(metadata %>% 
              select(name = SAMPLE, everything())
            )

?separate()
taxon_classes = c("DOMAIN", "SUPERGROUP", "PHYLUM", "CLASS", "ELSE")
# parse taxon names
asv_withmetadata_tax <- asv_withmetadata %>% 
  separate(Taxon, into = taxon_classes, sep = ";")



# barplot to supergroup
asv_withmetadata_tax %>% 
  filter(SUPERGROUP != "Bacteria_X") %>% 
  filter(!is.na(SUPERGROUP)) %>% 
  group_by(name, SUPERGROUP) %>% 
    summarize(SUM_SUPERGROUP = sum(value)) %>% 
      ggplot(aes(x = name, y = SUM_SUPERGROUP, fill = SUPERGROUP)) +
      geom_bar(stat = "identity", position = "fill", color = "black") +
      theme(axis.text.x = element_text(angle = 90))

# PCA plot
head(asv_long)

# need matrix
asv_matrix <- asv_withmetadata_tax %>% 
  select(name, value, SAMPLETYPE, FeatureID) %>% 
  unite(SAMPLENAME, name, SAMPLETYPE, sep = ";") %>% 
  pivot_wider(names_from = FeatureID, values_fill = 0) %>% 
  column_to_rownames(var = "SAMPLENAME") %>% 
  as.matrix()

# Install packages for PCA analysis
# install.packages("compositions")
library(compositions)

?clr()
?prcomp()
asv_clr <- clr(asv_matrix)
View(t(asv_clr))
pca_output <- prcomp(asv_clr)
head(pca_output)

data.frame(pca_output$x, SAMPLE = rownames(pca_output$x)) %>%
  separate(SAMPLE, into = c("samplenames","sampletype"), sep = ";") %>%
    ggplot(aes(x = PC1, y = PC2, shape = sampletype, color)) +
    geom_point()
