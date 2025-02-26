library(devtools)
install_github("komorowskilab/VisuNet")
library(VisuNet)

# GO annotation
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")
library(clusterProfiler)
library(org.Hs.eg.db)
#word cloud
install.packages("tm")
install.packages("SnowballC")
install.packages("wordcloud")
install.packages("RColorBrewer")

#visunet
install.packages("DT")
vis <- visunet(ros.P4$main)


#visuarc
visuArc(vis,'Asian',feature='UTS2_Th17_48')


#UTS2_Th17_48
