# Project: 2022 Squash Elicitor Field Experiment
# Final Project For ENT532 
# Nina Devine, December 2024 
# Data Cleaning and Wrangling ----
# importing and cleaning data ----
# GitHub: https://github.com/ninadevine/ent535
# install and load packages
library(tidyverse)
library(lubridate)

# import the fruit biomass data 
fruitbiomass <- read.csv("https://raw.githubusercontent.com/ninadevine/ent535/refs/heads/main/fruit_biomass_squash2022.csv")
# coerce date to be a date instead of a character
fruitbiomass$Date <- mdy(fruitbiomass$Date) 
# change column names to lower case
names(fruitbiomass) <- tolower(names(fruitbiomass))
# change name of treatment column from "type" to "trt"
fruitbiomass <- rename(fruitbiomass, trt = type)
# change the treatment column from a character to a factor 
fruitbiomass$trt <- as.factor(fruitbiomass$trt)
str(fruitbiomass) #looks good! 

#import the guild data 
guild.data <-  read.csv("https://raw.githubusercontent.com/ninadevine/ent535/refs/heads/main/guild_data_squash2022_clean.csv")
guild.data$Date <- mdy(guild.data$Date)
str(guild.data)
#do the same steps as before 
names(guild.data) <- tolower(names(guild.data)) 
guild.data <- rename(guild.data, trt = type)
str(guild.data)
guild.data$trt <- as.factor(guild.data$trt)
# remove all "unidentified" columns because they are not relevant to analysis
# but keep the total column which includes data from the unidentified columns 
guild.data <- guild.data %>% select(!(c(unidentified_black_insect,
                                        unidentified_red_abdomen,
                                        unidentified_eggs,
                                        unidentified_eggmass,
                                        unidentified_pupae)))

#import the herbivory data
herbivory.data <- read.csv("https://raw.githubusercontent.com/ninadevine/ent535/refs/heads/main/herbivory_squash2022.csv")
str(herbivory.data)
#lowercase names and coerce date column to date format 
names(herbivory.data) <- tolower(names(herbivory.data))
herbivory.data$date <- mdy(herbivory.data$date)
#rename treatment column 
herbivory.data <- rename(herbivory.data, trt = type)
herbivory.data$trt <- as.factor(herbivory.data$trt)
str(herbivory.data)
# create a new column expressing herbivory as a proportion rather than a percentage 
herbivory.data <- herbivory.data %>% mutate(prop.herb =  herbivory/100)
# write a CSV saving the cleaned data 
write_csv(herbivory.data, "herbivory_data_2024Nov06.csv", col_names = TRUE)
# import plant biomass and clean
plantbiomass <- read.csv("https://raw.githubusercontent.com/ninadevine/ent535/refs/heads/main/plant_biomass_squash2022_clean.csv")
names(plantbiomass) <- tolower(names(plantbiomass))                     
plantbiomass <- rename(plantbiomass, trt = type)
plantbiomass$trt <- as.factor(plantbiomass$trt)
str(plantbiomass)
write_csv(plantbiomass, "plant_biomass_data_2024Nov06.csv", col_names = TRUE)
#import taxon-level morphospecies data set
# and follow steps as above 
taxon.data <- read.csv("https://raw.githubusercontent.com/ninadevine/ent535/refs/heads/main/taxon_data_squash2022clean.csv")
names(taxon.data) <- tolower(names(taxon.data))                     
taxon.data <- rename(taxon.data, trt = type)
taxon.data$date <- mdy(taxon.data$date)
taxon.data$trt <- as.factor(taxon.data$trt)
str(taxon.data)

# i want to note that theoretically the above modifications could be changed during the import phase 
# but i've had issues with doing things that way, especially with dates 
# so I'm just doing everything manually to make sure it works 

#rename columns in the taxon data for easier coding
taxon.data <- rename(taxon.data, c(flea.beetle = chrysomelidae_flea_beetle, 
                                     SCB = cucumber_striped_beetle,
                                     SPCB = cucumber_spotted_beetle, 
                                     coccinellidae.adult = coccinelidae_ladybug))
#now merge the guild.data with taxon.data
taxon.guild.data <- full_join(guild.data, taxon.data)

# creating a summary data frame ----
# want to create an "averaged" datasheet that contains total abundances of each taxon and guild
# plus average herbivory 
# average fruit biomass per plant 
# total number of fruits per plant 
# and total plant biomass per plant
# ignoring date for now 
# creating a separate dataframe for the total number of fruits per plant
# create columns for total fruit weight and total fruit number per plant 
summed.data <- fruitbiomass %>% group_by(p_id, trt, block) %>%
  summarize(total_fruit_n = sum(fruits_n, na.rm = TRUE), 
            total_fruit_wt_g = sum(weight_g, na.rm = TRUE))
#create a column for average fruit weight 
summed.data <- summed.data %>% mutate(avg_fruit_wt = total_fruit_wt_g/total_fruit_n)
summed.data$avg_fruit_wt <- na_if(summed.data$avg_fruit_wt, NaN) #change NaNs to NAs

# need to deal with the fact that the total fruit wt and total fruit n columns 
# have zeroes for plants that died instead of NAs 
#change zeroes to NAs for when total fruit number = 0 and total fruit weight = 0 because plants died 
summed.data$total_fruit_n <- na_if(summed.data$total_fruit_n,0)
summed.data$total_fruit_wt_g <- na_if(summed.data$total_fruit_wt_g,0.0)

# add the summed taxon data
test.df <- taxon.guild.data %>% group_by(p_id, trt, block) %>% 
  summarise(across(2:99, sum, na.rm=TRUE))
summed.data <- full_join(summed.data,test.df) #join to summed dataframe
squash <- summed.data  #rename the object for easier coding 
squash$trt <- as.factor(squash$trt) # change treatment and block to factors 
squash$block <- as.factor(squash$block)
# create a data frame for summed herbivory data 
herb.2 <- herbivory.data %>% group_by(p_id) %>% 
  summarize(avg.herbivory = mean(prop.herb, na.rm = TRUE), 
            n.leaves = sum(leaves_n, na.rm = TRUE), 
            tot.herbivory = sum(herbivory/100, na.rm = TRUE))
herb.2$avg.herbivory <- na_if(herb.2$avg.herbivory, NaN) #change NaNs to NAs
# change zeroes to NAs because plants died
herb.2$tot.herbivory <- na_if(herb.2$tot.herbivory,0)
herb.2$n.leaves <- na_if(herb.2$n.leaves, 0)

# create the full dataset with all variables of interest
squash <- full_join(squash,herb.2)

# move some columns around for easier viewing 
squash <- squash %>% 
  relocate(avg.herbivory, .after = 6) %>% 
  relocate(tot.herbivory, .after = 7) %>% 
  relocate(n.leaves, .after = 8)

# add a column for all natural enemies
squash <- squash %>% mutate(natural_enemies = (predators+parasitoid))
# add a column for squash bug nymphs and adults combined because life history stages have the same feeding behavior 
squash <- squash %>% mutate(squash_bugs = sum(squash_bug_adult + squash_bug_nymph, na.rm = TRUE))
# move the column 
squash <- squash %>% relocate(squash_bugs, .after = 20)

# Data Exploration ----
# exploration of distributions ----
squash$trt <- as.factor(squash$trt)
squash$block <- as.factor(squash$block)
#using a loop function to plot multiple histograms at once
for(i in names(squash[5:22])){
  hist(squash[[i]], main = i, xlab = "Value")}
# everything is pretty zero-skewed
# except for average fruit weight, sq bug eggmass, dry weight, and SCB (sort of) 

# simple correlation testing ----
cor.test(squash$avg_fruit_wt, squash$suckers) #NS, but close 
plot(squash$avg_fruit_wt, squash$suckers)
cor.test(squash$avg_fruit_wt, squash$SCB) #NS, p = 0.1327
plot(squash$avg_fruit_wt, squash$SCB) #seems slightly positive? #probably related to overall plant size 
cor.test(squash$predators, squash$suckers) #S
# ^ increase in predators is correlated with increase in sucking insects 
plot(squash$predators, squash$suckers) #lots of dispersion here 
cor.test(squash$predators, squash$SCB) #S
plot(squash$predators, squash$SCB)
# more predators = more striped cucumber beetles 
cor.test(squash$natural_enemies, squash$squash_bugs) # S 
plot(squash$squash_bugs, squash$natural_enemies) 
# seems there is a positive relationship between predator and pest abundances 
# ANOVA: fruit number vs treatment ----
shapiro.test(squash$total_fruit_n) # data not normally distributed but we're gonna go with it 
fligner.test(total_fruit_n~trt, data = squash) # but variances are homogenous
fruit.trt <- aov(total_fruit_n~trt, data = squash, na.action = na.omit)
summary(fruit.trt)
TukeyHSD(fruit.trt)
shapiro.test(resid(fruit.trt)) #residuals normal 
# calculate number of fruits in each treatment group 
squash %>% group_by(trt) %>% summarize((mean(total_fruit_n, na.rm=TRUE)))
# 1 C                                      4.12
# 2 JA                                     4.61
# 3 SA                                     2.47
# calculate percent decrease 
((4.12-2.47)/4.12)*100
((4.61-2.47)/4.61)*100
# treatment with SA decreases fruit set by 40% compared to controls 
# and by 46% compared to JA 

# non-parametric tests---- 
# analyzing the effect of treatment on the abundance of 
# squash bugs and striped cucumber beetles
# and important guilds: pollinators and natural enemies 
# visualize distribution of variables by treatment ----
# to ensure within-treatment groups have similar distributions 
# natural enemies 
ggplot(squash, aes(x = natural_enemies)) + 
  geom_histogram(binwidth = 2) + 
  facet_grid(cols = vars(trt))
# pollinators 
ggplot(squash, aes(x = pollinators)) + 
  geom_histogram(binwidth = 2) + 
  facet_grid(cols = vars(trt))
# squash bugs 
ggplot(squash, aes(x = squash_bugs)) + 
  geom_histogram(binwidth = 30) + 
  facet_grid(cols = vars(trt))
#SCB 
ggplot(squash, aes(x = SCB)) + 
  geom_histogram(binwidth = 5) + 
  facet_grid(cols = vars(trt))
# all look pretty similar to me 

# run wilcoxon rank-sum tests ----
library(rstatix) # i like this package because it uses a pipe-friendly syntax 
# ran into trouble with the fact that my data is "grouped" by plant ID and by treatment 
# so I need to ungroup the data in order to analyze it 
squash <- ungroup(squash)
# natural enemies 
nat_en.trt <- squash %>% ungroup() %>%  kruskal_test(natural_enemies~trt)
nat_en.trt
dunn_test(squash, natural_enemies~trt, p.adjust.method = "BH")
# NS for all groups 

# squash bugs 
sqbug.trt <- squash %>% kruskal_test(squash_bugs~trt)
sqbug.trt
dunn_test(squash, squash_bugs~trt, p.adjust.method = "BH")
# NS for all groups 

# pollinators 
poll.trt <- squash %>% kruskal_test(pollinators~trt)
poll.trt
dunn_test(squash, pollinators~trt, p.adjust.method = "BH")
# C x JA significant 

# striped cucumber beetle 
scb.trt <- squash %>% kruskal_test(SCB~trt)
scb.trt
dunn_test(squash, SCB~trt, p.adjust.method = "BH")
# NS for all groups 

# herbivory 
herb.trt <- squash %>% kruskal_test(avg.herbivory~trt)
herb.trt
dunn_test(squash, avg.herbivory~trt, p.adjust.method = "BH")
# not significant 




# Making plots ----
library(patchwork)
library(ggpubr)
library(scales)
# don't need to load data - data already in R from running above code 
# correlation between natural enemies and squash bugs ----
cor.test(squash$squash_bugs, squash$natural_enemies) # p = 0.00689
ggplot(squash, aes(x = squash_bugs , y = natural_enemies)) + 
  geom_point() + 
  geom_smooth(method = "lm", alpha = 0.4, color = "black", fill = "darkmagenta") + 
  theme_light() + 
  stat_regline_equation(label.x = 15, label.y = 26) +
  stat_cor(aes(label= ..rr.label..),
           label.x= 10, label.y = 24)+
  stat_cor(aes(label = ..p.label..), label.x = 15, label.y = 28) +
  scale_y_continuous(breaks = seq(0,30, by = 5)) + 
  labs(x = "Squash bugs", 
       y = "Natural enemies", 
       title = "Natural Enemies Increase With Squash Bugs")

# correlation between natural enemies and squash bugs separated by treatment ----
ggplot(squash, aes(x = squash_bugs , y = natural_enemies, fill = trt)) + 
  scale_fill_manual(values = c("lightblue", 
                               "indianred3", 
                               "tan1")) + 
  geom_point() + 
  geom_smooth(method = "lm", alpha = 0.5, color = "black") + 
  theme_light() + 
  scale_y_continuous(breaks = seq(0,30, by = 5)) + 
  scale_x_continuous(breaks = seq(0,700, by = 100)) + 
  labs(x = "Squash bugs", 
       y = "Natural enemies", 
       fill = "Treatment",
       title = "Correlation Between Natural Enemies and Squash Bugs by Treatment") + 
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14), 
        title = element_text(size = 10))

# plotting SCB abundance vs. treatment ----
SCB_plot <- ggplot(squash, aes(x = trt, y = SCB, fill = trt)) + geom_boxplot() + 
  scale_fill_manual(values = c("lightblue", 
                               "indianred3", 
                               "tan1")) + 
  theme_classic() + 
  labs(title = "Striped cucumber beetles") + 
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.x = element_text(color = "black", size = 14), 
        axis.text.y = element_text(color = "black", size = 14),
        title = element_text(size = 14)) + 
  scale_y_continuous(breaks = seq(0,100, by = 20)) + 
  annotate("text", label = "ns", x = 0.75, y = 45, size = 5)+ 
  annotate("text", label = "ns", x = 1.75, y = 55, size = 5)+ 
  annotate("text", label = "ns", x = 2.75, y = 45, size = 5)
SCB_plot
# plotting natural enemies vs. trt ----
nat_enem_plot <- ggplot(squash, aes(x = trt, y = natural_enemies, fill = trt))+ 
  geom_boxplot() + 
  scale_fill_manual(values = c("lightblue", 
                               "indianred3", 
                               "tan1")) + 
  theme_classic() + 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 14), 
        axis.text.y = element_text(color = "black", size = 14), 
        title = element_text(size = 14)) + 
  labs(y = "Abundance", title = "Natural enemies")+ 
  annotate("text", x = 0.75, y = 14, label = "ns", size = 5) +
  annotate("text", x = 1.75, y = 13.5, label = "ns", size = 5)+
  annotate("text", x = 2.75, y = 10.5, label = "ns", size = 5)
nat_enem_plot

# pollinators vs. trt ----
poll_plot <- ggplot(squash, aes(x = trt, y = pollinators, fill = trt))+ 
  geom_boxplot() + 
  scale_fill_manual(values = c("lightblue", 
                               "indianred3", 
                               "tan1")) + 
  theme_classic() + 
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "black", size = 14), 
        axis.text.y = element_text(color = "black", size = 14), 
        title = element_text(size = 14)) + 
  labs(title = "Pollinators") + 
  scale_y_continuous(breaks = seq(0,20, by = 5)) + 
  annotate("text", label = "a", x = 0.75, y = 8.5, size = 5) + 
  annotate("text", label = "b", x = 1.75, y = 16, size = 5)+ 
  annotate("text", label = "ab", x = 2.75, y = 12, size = 5)
poll_plot

# suckers (squash bugs) vs. trt ----
sqbug_plot <- ggplot(squash, aes(x = trt, y = squash_bugs, fill = trt))+ 
  geom_boxplot() + 
  scale_fill_manual(values = c("lightblue", 
                               "indianred3", 
                               "tan1")) + 
  theme_classic() + 
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", size = 14),
        axis.text.y = element_text(colour = "black",  size = 14), 
        title = element_text(size = 14)) + 
  labs(y = "Abundance", title = "Squash bugs") + 
  scale_y_continuous(breaks = seq(0,500, by = 100)) + 
  annotate("text", label = "ns", x = 0.75, y = 350, size = 5)+ 
  annotate("text", label = "ns", x = 1.75, y = 430, size = 5)+ 
  annotate("text", label = "ns", x = 2.75, y = 290, size = 5)
sqbug_plot
# arrange the plots together ----
ggarrange(nat_enem_plot, poll_plot, sqbug_plot, SCB_plot, label.y = "Abundance")
ggsave("guild_plot.png", plot = last_plot())
# fruit number vs treatment plot ----
fruits_plot <- ggplot(squash, aes(x = trt, y = total_fruit_n, fill = trt)) + 
  geom_boxplot(color="black")+
  scale_fill_manual(values = c("lightblue", 
                               "indianred3", 
                               "tan1")) + 
  labs(x = NULL, 
       y = "Number of Fruits", 
       title = "Salicylic acid reduces fruit set") + 
  theme_classic() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(color = "black", size = 14.5),
        axis.text.y = element_text(colour = "black", size = 14.5), 
        axis.title.y = element_text(size = 14.5), 
        title = element_text(size = 15)) + 
  annotate("text", label = "ab", x = .75, y= 6.5, size = 5) + 
  annotate("text", label = "a", x = 1.75, y = 6.5, size = 5) + 
  annotate("text", label = "b", x = 2.75, y = 3.5, size = 5)
fruits_plot
ggsave("fruits_plot.png", plot = fruits_plot)
# herbivory by treatment ----
# plotting as a percentage for easier visualization 
herb.plot <- ggplot(squash, aes(x = trt, y = avg.herbivory*100, fill = trt)) + 
  geom_boxplot(color="black")+
  scale_fill_manual(values = c("lightblue", 
                               "indianred3", 
                               "tan1")) + 
  labs(x = NULL, 
       y = "Average Herbivory (%)") + 
  theme_classic() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(color = "black", size = 18),
        axis.text.y = element_text(colour = "black", size = 15), 
        axis.title.y = element_text(size = 16)) + 
  scale_y_continuous(limits = c(0,50), breaks = seq(0,50, by = 10))
herb.plot
ggsave("herbivory_vs_trt.png", plot = herb.plot)

# plot SCB abundance over time by treatment ----
scb_time <- ggplot(taxon.guild.data, aes(x = date, y = SCB, fill = trt)) + 
  geom_jitter(shape = 21, alpha = 0.8, show.legend = FALSE) + 
  facet_grid(.~trt) + 
  scale_fill_manual(values = c("lightblue", 
                               "indianred3", 
                               "tan1")) + 
  theme_linedraw() + 
  labs(y = "Striped Cucumber Beetle Abundance", x = "Date") + 
  ylim(0,50) + 
  theme( panel.grid.major = element_blank(), 
         axis.text.x = element_text(angle = 70, hjust = 1, size = 10), 
         strip.text.x = element_text(size = 14), 
         axis.title = element_text(size = 14))+
  scale_x_date(date_labels = ("%m/%d"),
               breaks = as_date(c("2022-06-15", 
                                  #"2022-06-17",
                                  "2022-06-22",
                                  #"2022-06-24",
                                  "2022-06-28",
                                  "2022-07-08",
                                  "2022-07-14",
                                  "2022-07-21",
                                  "2022-07-26",
                                  "2022-08-02")))
scb_time
ggsave("scb_time.png", plot = scb_time)


# PERMANOVA with taxon data ----
library(vegan)
# create a dataframe with all the counts for each taxa summed per plant throughout the season
perm.taxa <- taxon.data %>% group_by(p_id, trt) %>% select(-date) %>% 
  summarize(across(1:89, sum, na.rm=TRUE)) 

# make a table of dependent variables only 
bugs <- perm.taxa[ ,3:91]

# make a distance matrix using bray-curtis distances  
dist_matrix <- vegdist(bugs, method = "bray")

# do the PERMANOVA ! 
elicitor_permanova <- adonis2(dist_matrix~trt, permutations = 999, data = perm.taxa)
elicitor_permanova # not significant :(

# make a PCA plot 
library(ggfortify)
library(RColorBrewer)
pca_bug <- prcomp(bugs, scale. = FALSE)
autoplot(pca_bug, data = perm.taxa, color = "trt")+
  scale_color_brewer(palette = "Dark2") +
  theme_classic() + labs(color = "Treatment")

# make a PCA plot using ggplot instead
# build a data frame
pca <- as.data.frame(pca_bug$x[, 1:2]) # extract first two PCs
pca <- cbind(pca, perm.taxa$trt) # add trt to df 
colnames(pca) <- c("PC1", "PC2", "Species") # change column names
summary(pca_bug) #PC1 = 95.9%, PC2 = 2.51%
# plot
ggplot(pca, aes(PC1, PC2, color = perm.taxa$trt)) +
  geom_point(size = 2) + 
  scale_color_brewer(palette = "Dark2") + 
    theme_classic() + labs(x = "PC1: 95.9%", y = "PC2: 2.51%",  color = "Treatment") + 
  theme(axis.text = element_blank()) + 
  stat_ellipse(geom = "polygon", alpha = .1)
# lol that looks hilarious 

# survival analysis with plant mortality data ----
# read in the data frame and convert dates 
plant_mortality <- read.csv("https://raw.githubusercontent.com/ninadevine/ent535/refs/heads/main/plant_mortality.csv")
plant_mortality$plant.date <- mdy(plant_mortality$plant.date)
plant_mortality$mort.date <- mdy(plant_mortality$mort.date)
str(plant_mortality)
# add a new column for "days alive" 
plant_mortality <- mutate(plant_mortality, days.alive = as.numeric(mort.date-plant.date))
# try a survival analysis 
library(survival)
library(survminer)
# make the kaplan-meier estimate 
#function nested within function
# kaplan <- survfit(Surv(time, status)~response variable + covariate, data = data)
kaplan <- survfit(Surv(days.alive, status)~trt, data = plant_mortality)
summary(kaplan)

# log-rank test for significance differences in survival between groups 
survdiff(Surv(days.alive, status)~trt, data = plant_mortality)    
#pairwise comparison between groups 
pairwise_survdiff(Surv(days.alive, status)~trt, data = plant_mortality)
#not significant but a pattern between JA-SA and C-SA - seems like SA is affecting plant mortality
# make a graph
ggsurvplot(kaplan, data = plant_mortality,
           color = "strata",
           palette = c("#1B9E77", "#D95F02", "#7570B3"),
           pval = TRUE,
           linetype = "solid",
           linewidth = 0.6) + labs(x = "Days")

