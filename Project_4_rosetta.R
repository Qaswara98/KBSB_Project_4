set.seed(1)
library("usethis")
library(devtools)
install_github("komorowskilab/R.ROSETTA")
library(R.ROSETTA)


sig.features <- result2$attribute

sig.features <- c(sig.features, "ethnicity")

project_4_sig <- project_4[, sig.features]

# run rosetta
ros.P4 <- rosetta(project_4_sig, roc = TRUE, clroc = "Cau")
rules <- ros.P4$main
qual <- ros.P4$quality

#check quality
qual

# significant rules in the model
tabS <- table(rules[rules$pValue < 0.05,]$decision)
# fraction of significant rules, in [%]
tabS[1]/sum(rules$decision=="Afro")*100
tabS[2]/sum(rules$decision=="Asian")*100
tabS[3]/sum(rules$decision=="Cau")*100

# ROC curve
plotMeanROC(ros.P4)

# analyze the rule-based model
gf <- getFeatures(rules, filter = T, filterType = "pvalue", thr = 0.05)
genesAfro <- gf$features$Afro
genesAsian <- gf$features$Asian
genesCau <- gf$features$Cau

rec.ros.P4 <- recalculateRules(ros.P4, rules)

