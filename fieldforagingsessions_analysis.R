# This file contains the script to perform an extra analysis  
# to compare field observations to a lab result in the article
# "Exploration-exploitation learning in foraging hummingbirds is driven by cumulative experience and energetic costs"
# by Courtney Donkersteeg and Roslyn Dakin

# For questions, please contact Courtney Donkersteeg (courtneydonkersteeg@cmail.carleton.ca)

library(stringr)
library(lubridate)
library(dplyr)

# First, get the data...
data <- read.csv('fieldforagingsessions_data.csv', stringsAsFactors=F)
head(data)
str(data)
nrow(data) 

# Filter for the 8x4 array with the random condition
# This indicates the size of the flower array used and the nectar configuration
# The array was 1/4 rewarded with positions changing each day
data <- data %>% filter(array_size == '8x4')
data <- data %>% filter(distribution == 'random')

# Distribution of flower IDs for comparison with lab data
# which is the distribution of visits across the array space
#### FIGURE A3
data$flowerID <- as.numeric(data$flowerID)
plot(density(data$flowerID), main = '', ylab = NA, xlab = "Probe distribution for field array (8x4)", yaxt = "n", zero.line = FALSE, frame.plot = FALSE); polygon(density(data$flowerID), col=rgb(0.8,0.8,0.8))


















