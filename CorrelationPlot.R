library("readxl")
library(tidyverse)
library("RColorBrewer")

library(corrplot)
install.packages("purr")
corr <- as.matrix(read_excel("Summer_action.xlsx", sheet = "Neural Correlation"))

# Create the correlation of the summer graph
rownames(corr) <- c("STG", "MTG", "TPJ", 'TP', 'Precu',
                    'aMPFC','pMPFC', 'IFG oper','IPS','Auditory','Visual')
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA", 
                          "#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(corr, type = "upper", method="color",col = col(200),
         cl.lim = c(0, 1), addCoef.col = "black",tl.col="black",
         diag = FALSE)

predict_corr <- as.matrix(read_excel("Summer_action.xlsx", sheet = "Feature Correlation"))
rownames(predict_corr) <- c("self", "others", "things", "social", 
                            "mentalization", "touch","amplitude", "visual", "face", "action")
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(predict_corr,type="upper", 
         method="color", col=col(200), addCoef.col = "black",
         tl.col="black", diag = FALSE)



corr <- as.matrix(read_excel("Sherlock_with action.xlsx", sheet = "Neural Correlation"))
rownames(corr) <- c("STG", "MTG", "TPJ", 'TP', 'Precu',
                    'aMPFC','pMPFC', 'IFG oper','IPS','Auditory','Visual')
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA", 
                          "#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(corr, type = "upper", method="color",col = col(200),
         cl.lim = c(0, 1), addCoef.col = "black",tl.col="black",
         diag = FALSE)

predict_corr <- as.matrix(read_excel("Sherlock_with action.xlsx", sheet = "feature correlation"))
rownames(predict_corr) <- c("self", "others", "things", "social", 
                            "mentalization", "visual", "amplitude", "face", "action")
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(predict_corr, type="upper", method="color", col=col(200), 
         addCoef.col = "black",tl.col="black", diag = FALSE)

