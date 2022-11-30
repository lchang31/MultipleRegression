library("readxl")
library(tidyverse)
library("RColorBrewer")
install.packages("ggsignif")
library(ggsignif)
install.packages("purrr")
library(purrr)
            
#total
file <- 'Pub_graphs.xlsx'
sheets <- excel_sheets(file)
total <- map_df(sheets, ~ read_excel(file, sheet = .x))
total_long <- gather(total, subject, betavalues, Subject1:Subject18, factor_key=TRUE)
total_long<-group_by(total_long, Stimulus) 
total$features <- factor(total$features, levels = c("self", "others", "thing",
  "social_nonsocial", "mentalization", "visual", "auditory", "face", "action", "touch"))

# Organize ROI to display all ROI beta figures
#total$ROI <-factor(total$ROI, levels = c("STG", "MTG", "TPJ","TP", "Precu", "aMPFC", "pMPFC",
#                                         "IFG", "IPS", "Auditory", "Visual"))


# comment out one of the ones below to produce figures of specific network
# only show the social ROI
#total$ROI <-factor(total$ROI, levels = c("STG", "MTG"))

# only show the ToM ROI
#total$ROI <-factor(total$ROI, levels = c("TPJ","TP", "Precu", "aMPFC", "pMPFC"))

# only show the action observation ROI
#total$ROI <-factor(total$ROI, levels = c("IFG", "IPS"))

# only show the auditory visual ROI
#total$ROI <-factor(total$ROI, levels = c("Auditory", "Visual"))


#Remove rid of NA
total <- total[complete.cases(total[ , 1]),]      

# Group Stimulus
total$Stimulus <-factor(total$Stimulus, levels = c("Sherlock", "Summer"))

#if p value is greater than 0.05 then remove color
total$alpha = ifelse((is.na(total$"p.signif")), 0, 1)
p <- ggplot(total, aes(x = features, y = beta, fill = Stimulus, color = Stimulus)) +
 geom_bar(aes(features, fill = Stimulus), stat = "identity", position = position_dodge(0.8), 
           width = 0.8, alpha= total$alpha)+
  scale_fill_manual(values = c('#999999','#E69F00'))+
  scale_color_manual(values = c('#999999','#E69F00'))+
  
  # display standard deviation
  geom_errorbar(aes(ymin=beta-sd, ymax=beta+sd), color = "black", width=.2,
                position=position_dodge(.8))  + 
  # draw reliability lines
  geom_hline(data= total, aes(yintercept=reliability, color= Stimulus))  

# display final figure
p + theme_classic() + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
  facet_wrap(~ROI, nrow = 6,
             ncol = 2, scales = "free_y")




