data <- read.table("Project4.csv", header = TRUE)

data <- subset(data, select = - c(age, height.cm., weight.kg., bmi, systolic, diastolic, sex)) 

library(tidyverse)
data <- data %>% column_to_rownames(.,'object_id')
               
