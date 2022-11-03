# This file contains the script to perform the null data permutations  
# used in the main analyses in the article
# "Exploration-exploitation learning in foraging hummingbirds is driven by cumulative experience and energetic costs"
# by Courtney Donkersteeg and Roslyn Dakin

# For questions, please contact Courtney Donkersteeg (courtneydonkersteeg@cmail.carleton.ca)

library(tidyverse)
library(lubridate)
library(stringr)
library(vegan)
library(lme4)

# Get the data:
data <- read_csv('labforagingsessions_data.csv')
data
nrow(data) # 3464 probes, 5 individuals in 7 trials

# First, set up the starttimes for each trial:
starttimes <- unique(data[,c('bird_rep','start','end')])
head(starttimes)
# start times are entered as 24-h clock time, so we convert to minutes since midnight:
starttimes <- starttimes %>%
  mutate(start_h = ifelse(nchar(start) > 3, substr(start, 1, 2), substr(start, 1, 1)),
         end_h = substr(end, 1, 2)) %>%
  mutate(start_m = ifelse(nchar(start) > 3, substr(start, 3, 4), substr(start, 2, 3)),
         end_m = substr(end, 3, 4)) %>% 
  mutate(start_val = as.numeric(start_h) * 60 + as.numeric(start_m),
         end_val = as.numeric(end_h) * 60 + as.numeric(end_m)) %>% 
  mutate(duration = end_val - start_val)
starttimes <- starttimes[!duplicated(starttimes$bird_rep),]
summary(starttimes$duration)

# Omit the session that was spread over two days from further analysis
data <- subset(data, bird_rep != 'bird_181-1') 
(meta <- unique(data[,c('bird_rep','birdID','rep','mass_g','sex')]))

#### Null data permutation analyses
# This permutation (perm 3) randomizes the flower IDs WITHIN a given bird-rep AND generates times from a uniform distribution.
# Here, we (generate 1,000 null permutated datasets, recalculating the experience column, and regenerating the response variables for each.
# In the null perm, we shuffle the IDs of flowers probed within a trial.
# And we regenerate times from a uniform distribution.
# This preserves each individual's distribution of flowers visited within a trial, but breaks (removes) any temporal associations within trials

data_perm <- vector('list', length = 1000) # where we will store 1,000 permuted datasets

for(j in 1:1000){
  pdata <- data
  for(i in unique(pdata$bird_rep)){
    temp <- pdata[pdata$bird_rep == i, ] 
    index <- sample(1:nrow(temp))
    temp$flowerID <- temp$flowerID[index]
    temp$rewarded <- temp$rewarded[index] # shuffle flowers
    strt <- as.numeric(starttimes[starttimes$bird_rep == i, 'start_val'])
    endt <- as.numeric(starttimes[starttimes$bird_rep == i, 'end_val'])
    temp$time_h <- runif(n = nrow(temp), min = 0, max = 240)/60 
    temp <- temp[with(temp, order(time_h)),]
    for(k in 2:nrow(temp)){
      temp$different[k] <- !(temp$flowerID[k] %in% temp$flowerID[(k - 1)])
    }
    pdata[pdata$bird_rep == i, ] <- temp 
    # print(i)
  }
  data_perm[[j]] <- pdata
  print(j)
}

mean(data$rewarded)
mean(data_perm[[sample(1:100, 1)]]$rewarded) # check: this should be stable (this aspect is maintained)

# Now we (re)determine the experience categories for each of the null perm data sets. This takes a few minutes.
for(k in 1:1000){
  
  pdata <- data_perm[[k]]
  pdata$explore <- F
  
  for(i in unique(pdata$birdID)){
    
    temp <- pdata[pdata$birdID==i, ]
    
    for(j in 2:nrow(temp)){
      fid <- temp$flowerID[j]
      curr_rep <- temp$rep[j]
      step1 <- temp[1:(j-1),]
      step2 <- subset(step1, rep == curr_rep)
      step3 <- subset(step2, flowerID == fid & rewarded == T)
      if(nrow(step3) > 0){
        temp$explore[j] <- F
      } else {
        temp$explore[j] <- T
      }
    }
    pdata[pdata$birdID==i, ] <- temp
    # print(i)
  }
  data_perm[[k]] <- pdata
  print(k)
  # Save the 1000 null datasets:
  save(data_perm, file = 'null_permuted_3_datasets.RData')
  message <- paste('You are on iteration', k, sep = ' ')
  save(message, file = 'null_data_progress_indicator.RData')
}

load('null_permuted_3_datasets.RData')

# Then, generate rdat for each:
t_int <- 20 # choose  a time interval in minutes
(times <- seq(0, 240, by = t_int)) 
rdat_base <- expand.grid(time = times, bird_rep = meta$bird_rep, stringsAsFactors = F) 
rdat_base$flr <- NA # Number of probes made
rdat_base$div <- NA # Diversity of flowers probed

rdat_perm <- vector('list', length = 1000)

for(k in 1:1000){
  rdat <- rdat_base
  for(i in 1:nrow(rdat)){
    temp_start <- starttimes[starttimes$bird_rep == rdat$bird_rep[i],]$start_val
    # get start time in min since midnight
    tm <- rdat[i,]$time + temp_start
    temp <- subset(data_perm[[k]], bird_rep == rdat$bird_rep[i])
    temp_prev <- subset(temp, min_midnight <= tm)
    temp_interval <- subset(temp_prev, (min_midnight >= (tm - t_int)) & (min_midnight <= tm))
    rdat$flr[i] <- nrow(temp_interval)
    if(nrow(temp_interval) > 0){
      rdat$div[i] <- diversity(table(temp_interval$flowerID), index = 'shannon')
    }
  }
  
  rdat <- merge(rdat, meta, by = 'bird_rep')
  rdat$time_h <- rdat$time/60
  rdat_perm[[k]] <- rdat
  print(k)
}

# Save the 1000 null rdat2 datasets:
save(rdat_perm, data_perm, file = 'null_permuted_3_rdat.RData')
load('null_permuted_3_rdat.RData')

# Fit and save 1000 models for the analyses with sig. time effects:
# Remember that there is no null permutation for choice diversity because the observed effects were not large or significant
null_mod_A <- vector('list', length = 1000) #non-reinforced choices
null_mod_B <- vector('list', length = 1000) #choice shift
null_mod_D <- vector('list', length = 1000) #foraging performance

# 1000 permutations for all three model predictions
# This will take a while
for(i in 1:1000){
  xxx <- data_perm[[i]]
  mod_A <- glmer(explore ~ rep * time_h + mass_g + (1|birdID), data = xxx, family = 'binomial')
  null_mod_A[[i]] <- mod_A
  mod_B <- glmer(different ~ rep * time_h + mass_g + (1|birdID), data = xxx, family = 'binomial')
  null_mod_B[[i]] <- mod_B
  mod_D <- glmer(rewarded ~ rep * time_h + mass_g + (1|birdID), data = xxx, family = 'binomial')
  null_mod_D[[i]] <- mod_D
  # null_predictions[[i]] <- visreg(mod, 'time_h', by = 'rep', overlay = T, scale = 'response', ylim = c(0, 1))
  print(i)
}

save(null_mod_A, null_mod_B, null_mod_D, file = 'null_permuted_3_predictions.RData')
load('null_permuted_3_predictions.RData')

mod_A_rate_1 <- rep(NA, 1000)
mod_A_rate_4 <- rep(NA, 1000)
mod_A_rate_7 <- rep(NA, 1000)
mod_B_rate_1 <- rep(NA, 1000)
mod_B_rate_4 <- rep(NA, 1000)
mod_B_rate_7 <- rep(NA, 1000)
mod_D_rate_1 <- rep(NA, 1000)
mod_D_rate_4 <- rep(NA, 1000)
mod_D_rate_7 <- rep(NA, 1000)
# 
for(i in 1:1000){
  xxx <- subset(data_perm[[i]], bird_rep != 'bird_181-1')
  
  pdata_null <- expand.grid(time_h = c(0, 4), rep = c(1, 4, 7), mass_g = 3.8)
  pdata_null$p_reward <- predict(null_mod_A[[i]], newdata = pdata_null, re.form = NA, type = 'response')
  pdata_null$p_explore <- predict(null_mod_B[[i]], newdata = pdata_null, re.form = NA, type = 'response')
  pdata_null$p_different <- predict(null_mod_D[[i]], newdata = pdata_null, re.form = NA, type = 'response')
  pdata_null_change <- pdata_null %>% group_by(rep) %>% summarize(sA = diff(p_reward), sB = diff(p_explore), sD = diff(p_different))
  pdata_null_change <- data.frame(pdata_null_change)
  mod_A_rate_1[i] <- pdata_null_change[1, 'sA']
  mod_A_rate_4[i] <- pdata_null_change[2, 'sA']
  mod_A_rate_7[i] <- pdata_null_change[3, 'sA']
  mod_B_rate_1[i] <- pdata_null_change[1, 'sB']
  mod_B_rate_4[i] <- pdata_null_change[2, 'sB']
  mod_B_rate_7[i] <- pdata_null_change[3, 'sB']
  mod_D_rate_1[i] <- pdata_null_change[1, 'sD']
  mod_D_rate_4[i] <- pdata_null_change[2, 'sD']
  mod_D_rate_7[i] <- pdata_null_change[3, 'sD']
  
  print(i)
}

save(mod_A_rate_1, mod_A_rate_4, mod_A_rate_7, mod_B_rate_1, mod_B_rate_4, mod_B_rate_7, mod_D_rate_1, mod_D_rate_4, mod_D_rate_7, file = 'null_permuted_3_change_estimates.RData')

#### Object 'null_permuted_3_change_estimates.RData' is loaded in main analysis and used for comparison with observed model predictions

