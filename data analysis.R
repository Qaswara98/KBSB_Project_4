# Load necessary libraries
library(dplyr)
library(tidyverse)
library(rmcfs)
library(R.ROSETTA)
#library(tm)
#library(SnowballC)
#library(wordcloud)
#library(RColorBrewer)
#library(VisuNet)
#library(clusterProfiler)
#library(org.Hs.eg.db)
#library(tibble) 

# Set working directory 
#setwd("The Knowledge-based Systems for Bioinformatics/Proejcts/Project 4/")

# Read the dataset
project_4 <- read.csv("Project4.csv", header = TRUE, sep = "\t")

# Check dataset dimensions
dim(project_4)

# Inspect a subset of the dataset
project_4[1:5, 1:5]  # Display first 5 rows

# Display first n rows along with the last column
n <- 10  
project_4[1:n, c(1:n, dim(project_4)[2])]

# Check the structure and summary of the dataset
str(project_4)  
summary(project_4)  
class(project_4)

# Display frequency table for the "ethnicity" column
table(project_4$ethnicity)

# --------------------------------------
# Data Cleaning: Removing Unnecessary Columns
# --------------------------------------

# Remove specified columns
project_4 <- project_4 %>% 
  select(-c(age, height.cm., weight.kg., bmi, systolic, diastolic, sex))

# Verify the changes
head(project_4)

# Set "object_id" as row names
project_4 <- project_4 %>% column_to_rownames("object_id")


# Verify the final structure
head(project_4)
dim(project_4)

# Check for missing values
any(is.na(project_4))  

#------------------------------------------------------
# Monte Carlo Feature Selection (MCFS)
#------------------------------------------------------

# Check mcfs documentation
?mcfs

# Perform Monte Carlo feature selection (random reducts)
result <- mcfs(ethnicity ~ ., data = project_4, projections = 3500, projectionSize = 0.1, 
               splits = 5, splitSetSize = 365, balance = 2, cutoffPermutations = 8, 
               threadsNumber = 8)

print(result)
head(result$RI) 

# Render results
plot(result, type = "distances")

# Fetch significant results
result2 <- result$RI[1:result$cutoff_value, ]

# Check structure of results
dim(result2)  # Check number of rows and columns
summary(result2)  # See if there are NA values

# Print selected attributes
print(result2$attribute)

# Compute minimum significant RI
min_significant <- min(result2$RI, na.rm = TRUE)  

# Compute maximum insignificant RI safely
max_insignificant <- max(result$RI[result$RI < min_significant], na.rm = TRUE)

difference <- min_significant - max_insignificant

# Prepare ID graph plot
gid <- build.idgraph(result, size = 20)
plot.idgraph(gid, label_dist = 0.3)
