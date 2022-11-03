# This file contains the script to perform the main analyses   
# and generate the figures and tables in the article
# "Exploration-exploitation learning in foraging hummingbirds is driven by cumulative experience and energetic costs"
# by Courtney Donkersteeg and Roslyn Dakin

# For questions, please contact Courtney Donkersteeg (courtneydonkersteeg@cmail.carleton.ca)

library(tidyverse)
library(lubridate)
library(stringr)
library(vegan)
library(lme4)
library(lmerTest)
library(visreg)
library(rptR)

# Get the data:
data <- read_csv('labforagingsessions_data.csv')
data
nrow(data) # 3464 probes, 5 individuals in 7 trials

mass_vals <- unique(data[,c('birdID','bird_rep','sex','mass_g')])
mass_vals # range of bird mass
# model for sex difference in body mass
#### FIGURE A1
t.test(mass_g ~ sex, data = mass_vals)
boxplot(mass_vals$mass_g ~ mass_vals$birdID, xlab = "", ylab = "Body mass (g)", col = c("red", "blue", "blue", "blue", "red"), )
legend(x = "bottomright", legend = c("female", "male"), fill = c("red", "blue"))

(meta <- unique(data[,c('bird_rep','birdID','rep','mass_g','sex')]))
nrow(meta) # 35 trials

data %>% group_by(bird_rep) %>% summarize(nprobes_total = n()) %>% summarize(mean(nprobes_total), median(nprobes_total), sd(nprobes_total))
# median of 105, mean of 99 (SD 44.2) probes in 4-hour session, or 8-9 probes per 20 min bin
data$bird_rep_hr <- NA
for(i in 1:nrow(data)){
  data$bird_rep_hr[i] <- paste0(data$bird_rep[i], "-", data$hour[i])
}
data %>% group_by(bird_rep_hr) %>% summarize(nprobes_total = n()) %>% summarize(mean(nprobes_total), median(nprobes_total), sd(nprobes_total))
# median 19, mean 22.3 (SD 16.6) probes in 1 hour of a session

# Choice diversity metric requires temporal binning

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

t_int <- 20 # we choose a 20 min time interval (in minutes) for binning
#### can re-run this with 15 min or 25 min bins for sensitivity analysis
#### which would produce TABLE S1 & TABLE S2
(times <- seq(0, 240, by = t_int)) 
rdat <- expand.grid(time = times, bird_rep = meta$bird_rep, stringsAsFactors = F) # data frame where metrics of behaviour will be stored; one row for each bird-trial-time interval
# We will compute the following for each bin:
rdat$any_rewards <- NA
rdat$flr <- NA # Number of flowers probed in that time bin
rdat$div <- NA # Diversity of flowers probed in that time bin
rdat$uni <- NA # Number of unique flowers probed
rdat$all_prev <- NA # Total number of unique flowers probed up until that time

# Loop that calculates the metrics above for each bin:
for(i in 1:nrow(rdat)){
  temp_start <- starttimes[starttimes$bird_rep == rdat$bird_rep[i],]$start_val
  # get start time in min since midnight
  tm <- rdat[i,]$time + temp_start
  temp <- subset(data, bird_rep == rdat$bird_rep[i])
  temp_prev <- subset(temp, min_midnight <= tm)
  rdat$any_rewards[i] <- any(temp_prev$rewarded)
  rdat$all_prev[i] <- length(unique(temp_prev$flowerID))
  temp_interval <- subset(temp_prev, (min_midnight >= (tm - t_int)) & (min_midnight <= tm))
  rdat$flr[i] <- nrow(temp_interval)
  if(nrow(temp_interval) > 0){
    rdat$div[i] <- diversity(table(temp_interval$flowerID), index = 'shannon')
    # this is the Shannon-Weaver diversity index as a measure of diversity of flowers probed
    rdat$uni[i] <- nrow(temp_interval[!duplicated(temp_interval$flowerID), ])
  }
  print(i)
}

# Illustrate that SDI is dependent on the number of samples (probes)
#### FIGURE A2
plot(rdat$div ~ rdat$uni, xlab = 'Unique flowers probed', ylab = 'SDI of flowers probed', cex.lab = 1.4, pch = 16)

rdat <- merge(rdat, meta, by = 'bird_rep')
rdat$time_h <- rdat$time/60
# rdat2 <- subset(rdat, rep > 1 & any_rewards == T)
# nrow(rdat2)

# Distribution of flower IDs for comparison with field data
# which is the distribution of visits across the array space
# see script 'fieldforagingsessions_analysis.R' for other half of comparison
#### FIGURE A3
plot(density(data$flowerID), main = '', ylab = NA, xlab = "Probe distribution for lab array (10x10)", yaxt = "n", zero.line = FALSE, frame.plot = FALSE); polygon(density(data$flowerID), col=rgb(0.8,0.8,0.8))

# Omit the session that was spread over two days from further analysis
data <- subset(data, bird_rep != 'bird_181-1') 
rdat <- subset(rdat, bird_rep != 'bird_181-1') 

# Here we check repeatability (extent of among-individual differences) 
# for four behavioural metrics (three explorative, one exploitative) and body mass:
#### TABLE A4
set.seed(101) # NOTE we determined that setting seed may produce slightly different results depending on operating system

rptBinary(explore ~ 1 + (1|birdID), 'birdID', data = data, nboot = 50) # prob. of non-reinforced choices

rptBinary(different ~ 1 + (1|birdID), 'birdID', data = data, nboot = 50) # prob. of choice shift

rptGaussian(div ~ log(flr) + (1|birdID), 'birdID', data = rdat, nboot = 100) # choice diversity

rptBinary(rewarded ~ 1 + (1|birdID), 'birdID', data = data, nboot = 50) # foraging performance (rewarded probes)

rptGaussian(mass_g ~ 1 + (1|birdID), 'birdID', data = data, nboot = 50) # bird body mass


#### Analysis of relationship between metrics of behaviour and predictors
#### Predictors: 
#### rep = time across sessions 
#### time_h = time within sessions
#### mass_g = bird body mass 
#### rep:time_h = interaction between time across and within sessions


#### Non-reinforced choice (probability of probing unexperienced or unrewarded flower):
#### TABLE 1 & FIGURE 2
summary(modA <- glmer(c(explore) ~ rep * time_h + mass_g + (1|birdID), data = data, family = 'binomial'))
visreg(modA, 'mass_g', scale = 'response', ylim = c(0,1), xlab = 'Body Mass (g)', ylab = NA, line.par = list(col = 'black'))
par(mar = c(5,6,5,2), xpd = TRUE)
visreg(modA, 'time_h', by = 'rep', overlay = T, scale = 'response', ylim = c(0,1), breaks = 3, xlab = 'Time in Session (hr)', ylab = 'Non-reinforced choice\n(probability of unrewarded\nor unexperienced flower)', legend = FALSE)
legend("topright", inset = c(0.15, -0.15), legend = c("session 1", "session 4", "session 7"), lty = 1, lwd = 2, col = c("orangered2", "green3", "steelblue"), box.col = "white", horiz = TRUE)

AIC(modA) # check sex in place of body mass, because there is a sig. body mass effect:
AIC(modA_sex <- glmer(explore ~ rep * time_h + sex + (1|birdID), data = data, family = 'binomial')) 
# no significant delta-AIC

# 95% CI based on results of summary
#rep
print(0.11271 - (1.96*0.04472))
print(0.11271 + (1.96*0.04472))
#time_h
print(0.07006 - (1.96*0.08408))
print(0.07006 + (1.96*0.08408))
#mass_g
print(-0.35324 - (1.96*0.15155))
print(-0.35324 + (1.96*0.15155))
#rep:time_h
print(-0.08591 - (1.96*0.01816))
print(-0.08591 + (1.96*0.01816))


#### Choice shift (probability of switching to a different flower):
#### TABLE 2 & FIGURE 3
summary(modB <- glmer(c(different) ~ rep * time_h + mass_g + (1|birdID), data = data, family = 'binomial'))
visreg(modB, 'mass_g', scale = 'response', ylim = c(0,1), xlab = 'Body Mass (g)', ylab = NA, line.par = list(col = 'black'))
par(mar = c(5,5,6,2), xpd = TRUE)
visreg(modB, 'time_h', by = 'rep', overlay = T, scale = 'response', ylim = c(0,1), breaks = 3, xlab = 'Time in Session (hr)', ylab = 'Choice shift\n(probability of different flower)', legend = FALSE)
legend("topright", inset = c(0.15, -0.15), legend = c("session 1", "session 4", "session 7"), lty = 1, lwd = 2, col = c("orangered2", "green3", "steelblue"), box.col = "white", horiz = TRUE)

AIC(modB) # check sex in place of body mass, because there is a sig. body mass effect:
AIC(modB_sex <- glmer(different ~ rep * time_h + sex + (1|birdID), data = data, family = 'binomial'))
# no significant delta-AIC

# 95% CI based on results of summary
#rep
print(0.29333 - (1.96*0.05766))
print(0.29333 + (1.96*0.05766))
#time_h
print(0.70132 - (1.96*0.10909))
print(0.70132 + (1.96*0.10909))
#mass_g
print(-0.55066 - (1.96*0.31037))
print(-0.55066 + (1.96*0.31037))
#rep:time_h
print(-0.14614 - (1.96*0.02321))
print(-0.14614 + (1.96*0.02321))


#### Choice diversity (evenness of probes in 20-min time bin):
#### TABLE 3 & FIGURE 4
summary(modC <- lmer(div ~ log(flr) + rep * time_h + mass_g + (1|birdID), data = rdat))
visreg(modC, 'mass_g', xlab = 'Body Mass (g)', ylab = NA, line.par = list(col = 'black'))
par(mar = c(5,5,6,2), xpd = TRUE)
visreg(modC, 'time_h', by = 'rep', overlay = T, breaks = 3, xlab = 'Time in Session (hr)', ylab = 'Choice diversity', legend = FALSE)
legend("topright", inset = c(0.15, -0.15), legend = c("session 1", "session 4", "session 7"), lty = 1, lwd = 2, col = c("orangered2", "green3", "steelblue"), box.col = "white", horiz = TRUE)

AIC(modC) # check sex in place of body mass, because there is a sig. body mass effect:
AIC(modC_sex <- lmer(div ~ log(flr) + rep * time_h + sex + (1|birdID), data = rdat))
# no significant delta-AIC

# 95% CI based on results of summary
#rep
print(0.02731 - (1.96*0.02522))
print(0.02731 + (1.96*0.02522))
#time_h
print(0.06066 - (1.96*0.04777))
print(0.06066 + (1.96*0.04777))
#mass_g
print(-0.24862 - (1.96*0.14367))
print(-0.24862 + (1.96*0.14367))
#rep:time_h
print(-0.01838 - (1.96*0.01011))
print(-0.01838 + (1.96*0.01011))


#### Note: extra metric, foraging performance (probability of probing rewarded flower):
#### TABLE A3 & FIGURE A4
summary(modD <- glmer(rewarded ~ rep * time_h + mass_g + (1|birdID), data = data, family = 'binomial'))
visreg(modD, 'mass_g', scale = 'response', ylim = c(0,1), xlab = 'Body Mass (g)', ylab = NA, line.par = list(col = 'black'))
par(mar = c(5,5,6,2), xpd = TRUE)
visreg(modD, 'time_h', by = 'rep', overlay = T, scale = 'response', ylim = c(0,1), breaks = 3, xlab = 'Time in Session (hr)', ylab = 'Foraging performance\n(probability of rewarded flower)', legend = FALSE)
legend("topright", inset = c(0.15, -0.15), legend = c("session 1", "session 4", "session 7"), lty = 1, lwd = 2, col = c("orangered2", "green3", "steelblue"), box.col = "white", horiz = TRUE)

AIC(modD) # check sex in place of body mass, because there is a borderline sig. body mass effect:
AIC(modD_sex <- glmer(rewarded ~ rep * time_h + sex + (1|birdID), data = data, family = 'binomial'))
# no significant delta-AIC

#(95% CI)
#rep
print(-0.10925 - (1.96*0.04343))
print(-0.10925 + (1.96*0.04343))
#time_h
print(-0.15122 - (1.96*0.08259))
print(-0.15122 + (1.96*0.08259))
#mass_g
print(0.24930 - (1.96*0.12963))
print(0.24930 + (1.96*0.12963))
#rep:time_h
print(0.07045 - (1.96*0.01770))
print(0.07045 + (1.96*0.01770))

# Trend lines of each metric (modA, B, and D) for sessions 1, 4, and 7
# for comparison with the null permutations
# Choice diversity (modC) excluded from null comparison due to non-significant results (see Methods)

pdata <- expand.grid(time_h = c(0, 4), rep = c(1, 4, 7), mass_g = 3.8)
pdata$p_A <- predict(modA, newdata = pdata, re.form = NA, type = 'response')
pdata$p_B <- predict(modB, newdata = pdata, re.form = NA, type = 'response')
pdata$p_D <- predict(modD, newdata = pdata, re.form = NA, type = 'response')
pdata_change <- pdata %>% group_by(rep) %>% summarize(sA = diff(p_A), sB = diff(p_B), sD = diff(p_D))
pdata_change

# Note these trend lines for the null permutation comparison

#### Load null data permutations and model estimates
#### These are provided but can be generated again with script labforagingsessions_null.R
#### Note that every new permutation will be a little different each time
load('null_permuted_3_change_estimates.RData')

#ablines for observed data, see object pdata_change

par(mfrow = c(3,1), las = 1, mar = c(5,5,0.75,0.75), mgp = c(2, 0.5, 0))

#### FIGURE 2 NULL
plot(density(mod_A_rate_1), xlim = c(-0.5,0.1), main = '', ylab = NA, xlab = NA, yaxt = "n", zero.line = FALSE, frame.plot = FALSE); polygon(density(mod_A_rate_1), col=rgb(0.8,0.8,0.8)); abline(v = -0.0152, col = 'orangered2', lwd = 2)
plot(density(mod_A_rate_4), xlim = c(-0.5,0.1), main = '', ylab = NA, xlab = NA, yaxt = "n", zero.line = FALSE, frame.plot = FALSE); polygon(density(mod_A_rate_4), col=rgb(0.8,0.8,0.8)); abline(v = -0.264, col = 'green3', lwd = 2)
plot(density(mod_A_rate_7), xlim = c(-0.5,0.1), main = '', ylab = NA, xlab = NA, yaxt = "n", zero.line = FALSE, frame.plot = FALSE); polygon(density(mod_A_rate_7), col=rgb(0.8,0.8,0.8)); abline(v = -0.486, col = 'steelblue', lwd = 2)

#### FIGURE 3 NULL
plot(density(mod_B_rate_1), xlim = c(-0.3,0.4), main = '', ylab = NA, xlab = NA, yaxt = "n", zero.line = FALSE, frame.plot = FALSE); polygon(density(mod_B_rate_1), col=rgb(0.8,0.8,0.8)); abline(v = 0.351, col = 'orangered2', lwd = 2)
plot(density(mod_B_rate_4), xlim = c(-0.3,0.4), main = '', ylab = NA, xlab = NA, yaxt = "n", zero.line = FALSE, frame.plot = FALSE); polygon(density(mod_B_rate_4), col=rgb(0.8,0.8,0.8)); abline(v = 0.0735, col = 'green3', lwd = 2)
plot(density(mod_B_rate_7), xlim = c(-0.3,0.4), main = '', ylab = NA, xlab = NA, yaxt = "n", zero.line = FALSE, frame.plot = FALSE); polygon(density(mod_B_rate_7), col=rgb(0.8,0.8,0.8)); abline(v = -0.203, col = 'steelblue', lwd = 2)

#### FIGURE A4 NULL
plot(density(mod_D_rate_1), xlim = c(-0.2,0.4), main = '', ylab = NA, xlab = NA, yaxt = "n", zero.line = FALSE, frame.plot = FALSE); polygon(density(mod_D_rate_1), col=rgb(0.8,0.8,0.8)); abline(v = -0.0806, col = 'orangered2', lwd = 2)
plot(density(mod_D_rate_4), xlim = c(-0.2,0.4), main = '', ylab = NA, xlab = NA, yaxt = "n", zero.line = FALSE, frame.plot = FALSE); polygon(density(mod_D_rate_4), col=rgb(0.8,0.8,0.8)); abline(v = 0.129, col = 'green3', lwd = 2)
plot(density(mod_D_rate_7), xlim = c(-0.2,0.4), main = '', ylab = NA, xlab = NA, yaxt = "n", zero.line = FALSE, frame.plot = FALSE); polygon(density(mod_D_rate_7), col=rgb(0.8,0.8,0.8)); abline(v = 0.326, col = 'steelblue', lwd = 2)




