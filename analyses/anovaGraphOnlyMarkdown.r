```{r, message=FALSE, warning = FALSE, echo = FALSE}
## Setup

#Load Packages and Read in Subject CSV files

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)
library(grid)
library(gridExtra)

# Get list of region CSV files
filelist = list.files(path = dataPath[1], recursive = TRUE, pattern = "table.csv", full.names = TRUE)

# Creates list of variables, comma separated with a header  
datalist = lapply(filelist, function(x)read.table(x, header=T, sep = ",")) 

# Aggregates into single variable
classAccuracies = do.call("rbind", datalist) 

## One-Way ANOVA

#Read accuracies into new table and compute ANOVA against chance (50%)

# Change variable types from integer to text
classAccuracies$subjectid = as.character(classAccuracies$subjectid)
classAccuracies$TrialTypeCombo = as.character(classAccuracies$TrialTypeCombo)

# Flag for bad subj ID - project specific
classAccuracies$subjectid[is.na(classAccuracies$subjectid)] <- "6666"

options(scipen=999) #Gers rid of scientific notation

classAccuraciesCov<-classAccuracies %>%                    #Create table named classAccuraciesCov
  group_by(roiid) %>%                                      #Within each ROI,
  summarise(pval = list(tidy(t.test(accuracy, mu=.5))),    #Perform a 1s t-test against chance, 
            sdAcc = sd(accuracy)) %>%                      #Then take sd of classifier acc,
  unnest() %>%                                             #Turn list to column,
  mutate_if(is.numeric, ~round(., 4))                      #Round variables to the fourth decimal, and

# Plot Table
classAccuraciesPlot = classAccuraciesCov[,c(1:5,7)]
colnames(classAccuraciesPlot) <- c("regionID", "avgClassAcc", "t-stat","pVal","N","maxAcc")
grid.table(classAccuraciesPlot)
```
Statistics Table

```{r, message=FALSE, warning = FALSE, fig.align = 'center', echo = FALSE}
grid.table(classAccuraciesPlot)
```
Graph of Classification Accuracy

```{r, message=FALSE, warning = FALSE, fig.align = 'center', echo = FALSE}
##Figure 
#Plot Group Accuracies with SEM
classAccuracies %>%
  ggplot(aes(x = roiid, y = accuracy, fill = roiid)) +                        #Create plot
  stat_summary(geom = "bar", fun.y = mean, position = "dodge") +              #Average classifier accuracy to make bars
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge") +   #Get SEM and display as error bars
  coord_cartesian(ylim = c(min(classAccuraciesCov$conf.low), max(classAccuraciesCov$conf.high))) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),                    #Here down is all formatting
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle(paste("Classification Accuracy For ", classAccuracies$TrialTypeCombo[1], sep="")) + 
  ylab("Classification Accuracy") + 
  geom_hline(yintercept=c(.5), linetype="dotted")
```
