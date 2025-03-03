# -------------------------------
# Load Required Libraries
# -------------------------------
library(dplyr)
library(tidyverse)
library(rmcfs)
library(R.ROSETTA)
library(VisuNet)
library(arcdiagram)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(tm)
library(SnowballC)
library(wordcloud)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(DT)
# -------------------------------
# Set Working Directory & Read Data
# -------------------------------
#setwd("The Knowledge-based Systems for Bioinformatics/Projects/Project 4/")
project_4 <- read.csv("Project4.csv", header = TRUE, sep = "\t")

# ------------------------------------
# Project inspection and Data Cleaning
# -------------------------------------

# Inspect a subset of the dataset
project_4[1:5, 1:5]  # First 5 rows and first 5 columns

# Display first n rows along with the last column
n <- 10  
project_4[1:n, c(1:n, dim(project_4)[2])]

# Check the structure and summary of the dataset
str(project_4)  
summary(project_4)  
class(project_4)

# Display frequency table for the "ethnicity" column
table(project_4$ethnicity)

# Remove unnecessary metadata columns
project_4 <- project_4 %>% 
  dplyr::select(-c(age, height.cm., weight.kg., bmi, systolic, diastolic, sex))
# Set "object_id" as row names
project_4 <- project_4 %>% column_to_rownames("object_id")
# Verify the final structure
head(project_4)
dim(project_4)
# Returns TRUE if there are any NAs, otherwise FALSE
any(is.na(project_4))  

# -------------------------------
# Feature Selection using Monte Carlo Feature Selection (MCFS)
# -------------------------------
result <- mcfs(ethnicity ~ ., data = project_4, 
               projections = 4000, projectionSize = 0.1, 
               splits = 5, splitSetSize = 365, balance = 2, 
               cutoffPermutations = 8, threadsNumber = 8)
print(result)
head(result$RI) 
# render results
plot(result, type = "distances")

# Extract significant features
result2 <- result$RI[1:result$cutoff_value,] 
dim(result2)  # Check the number of rows and columns
print(result2$attribute)

# prepare id graph plot
gid <- build.idgraph(result, size = 20)
plot.idgraph(gid, label_dist = 0.3)

sig.features <- c(result2$attribute, "ethnicity") #Include ethnicity for modeling


# Subset the dataset with significant features
project_4_sig <- project_4[, sig.features]


# -------------------------------
# Rule-Based Modeling with Rosetta
# -------------------------------

# Train a rule-based model using R.ROSETTA
ros.P4 <- rosetta(project_4_sig, roc = TRUE, clroc = "Afro")  
rules <- ros.P4$main  # Extract main rules
qual <- ros.P4$quality  # Extract model quality metrics
summary(ros.P4$quality) #Check performance of the Rosetta model



# ----------------------------------------
# Compute Significant Rules and Performance
# ----------------------------------------

# Count the number of significant rules per decision class
tabS <- table(rules[rules$pValue < 0.05,]$decision)

# Compute the percentage of significant rules per class
percent_afro <- (tabS["Afro"] / sum(tabS)) * 100
percent_asian <- (tabS["Asian"] / sum(tabS)) * 100
percent_caucasian <- (tabS["Cau"] / sum(tabS)) * 100

# Print the percentage of significant rules
print("Percentage of Significant Rules per Class:")
print(paste("Afro:", percent_afro, "%"))
print(paste("Asian:", percent_asian, "%"))
print(paste("Caucasian:", percent_caucasian, "%"))

# Plot the ROC curve for the rule-based model
plotMeanROC(ros.P4)

# -------------------------------
# Extract Features from Significant Rules
# -------------------------------

# Extract significant features based on p-value threshold
gf <- getFeatures(rules, filter = TRUE, filterType = "pvalue", thr = 0.05)

# Store feature lists for each class
genesAfro <- gf$features$Afro
genesAsian <- gf$features$Asian
genesCau <- gf$features$Cau

# Print extracted features
print("Significant Features for Afro:")
print(genesAfro)

print("Significant Features for Asian:")
print(genesAsian)

print("Significant Features for Caucasian:")
print(genesCau)

# ----------------------------------------
# Recalculate Rules for Better Interpretation
# ----------------------------------------

# Recalculate rules for alignment with the dataset
rec.ros.P4 <- recalculateRules(project_4_sig, rules)


# Display rule stability
table(rules$pValue < 0.05)  # Count significant rules
summary(rules$supportLHS)   # Check rule support
summary(rules$accuracyRHS)  # Rule accuracy


# Extract rules that classify as "Caucasian"
subset_data <- rec.ros.P4[rec.ros.P4$decision == 'Cau',]

# -------------------------------
# Extract Top 5 Rules Per Class
# -------------------------------

# Extract only significant rules
top_rules <- rec.ros.P4 %>% 
  filter(pValue < 0.05) %>%  # Select only significant rules
  arrange(pValue)  # Sort by significance (ascending p-value)

# Extract top 5 rules for each class
top_rules_afro <- head(top_rules[top_rules$decision == "Afro", ], 5)
top_rules_asian <- head(top_rules[top_rules$decision == "Asian", ], 5)
top_rules_caucasian <- head(top_rules[top_rules$decision == "Cau", ], 5)

# Print extracted rules
print("Top 5 Rules for Afro:")
print(top_rules_afro)

print("Top 5 Rules for Asian:")
print(top_rules_asian)

print("Top 5 Rules for Caucasian:")
print(top_rules_caucasian)

# -------------------------------
# Extract Features Used in Top Rules
# -------------------------------

# Function to extract unique features from the top rules
extract_features <- function(top_rules_class) {
  if (nrow(top_rules_class) == 0) {
    return(character(0))  # Return empty vector if no rules exist
  }
  return(unique(unlist(strsplit(paste(top_rules_class$features, collapse=","), ","))))
}

# Extract features used in the top 5 rules per class
features_top_afro <- extract_features(top_rules_afro)
features_top_asian <- extract_features(top_rules_asian)
features_top_caucasian <- extract_features(top_rules_caucasian)

# Print extracted features
print("Features from Top 5 Rules for Afro:")
print(features_top_afro)

print("Features from Top 5 Rules for Asian:")
print(features_top_asian)

print("Features from Top 5 Rules for Caucasian:")
print(features_top_caucasian)

# -------------------------------
# Identify Distinctive & Common Features
# -------------------------------

# Find features that are unique to each class
unique_top_afro <- setdiff(features_top_afro, c(features_top_asian, features_top_caucasian))
unique_top_asian <- setdiff(features_top_asian, c(features_top_afro, features_top_caucasian))
unique_top_caucasian <- setdiff(features_top_caucasian, c(features_top_afro, features_top_asian))

# Identify common features across all three classes
common_top_features <- intersect(features_top_afro, intersect(features_top_asian, features_top_caucasian))

# Print distinctive and common features
print("Unique features for Afro:")
print(unique_top_afro)

print("Unique features for Asian:")
print(unique_top_asian)

print("Unique features for Caucasian:")
print(unique_top_caucasian)

print("Common features across all classes:")
print(common_top_features)

# -------------------------------
# Feature Ranking Based on Occurrence
# -------------------------------

# Get all unique features across classes
all_features <- unique(c(features_top_afro, features_top_asian, features_top_caucasian))

# Function to count occurrences of each feature in a class
count_features <- function(features, all_features) {
  if (length(features) == 0) {
    return(rep(0, length(all_features)))  # Return zeroes if no features exist
  }
  return(as.numeric(table(factor(features, levels = all_features))))
}

# Create a feature ranking data frame
feature_ranking <- data.frame(
  Feature = all_features,
  Afro = count_features(features_top_afro, all_features),
  Asian = count_features(features_top_asian, all_features),
  Caucasian = count_features(features_top_caucasian, all_features)
)

# Ensure feature ranking is sorted correctly
feature_ranking <- feature_ranking %>% arrange(desc(rowSums(select(., -Feature))))

# Print feature ranking table
print("Feature Ranking Based on Occurrence:")
print(feature_ranking)

# -------------------------------
# Rule Visualization (Heatmap & Boxplot)
# -------------------------------

topRuleInd <- 1  # Define the index of the rule to visualize

# Generate heatmap visualization of rule-based model
plotRule(project_4_sig, rec.ros.P4, type = "heatmap", discrete = FALSE, ind = topRuleInd)

# Generate boxplot visualization of rule-based model
plotRule(project_4_sig, rec.ros.P4, type = "boxplot", discrete = FALSE, ind = topRuleInd)

#-----------------------------------
# Cluster Model Rules using Heatmap
# ----------------------------------

cluster_rules<-function(training_df,recal,support=7,fontsize=7,show_colnames=FALSE,show_rownames=FALSE){
  library(pheatmap)
  #process data for clustering, i.e. making a matrix with binary values 0 and 1
  dataM <- data.frame(matrix(ncol = length(recal$features[1:nrow(recal)]), nrow = length(row.names(training_df))))
  rownames(dataM)<-rownames(training_df)
  colnames(dataM)<-recal$features
  
  #assigning 1 to objects satisfying rules
  for(i in 1:length(recal$features[1:100])){
    list_of_features<-unlist(strsplit(recal$supportSetLHS[i],split=','))
    for(j in 1:length(list_of_features)){
      if(list_of_features[j] %in% rownames(dataM))(dataM[list_of_features[j],i]<-1)
    }
  }
  
  #replacing not satisfying rules as 0
  dataM[is.na(dataM)] <- 0
  
  #decision variable of the dataset
  decision_var<-names(training_df)[ncol(training_df)]
  ann <- data.frame( eval(parse(text=paste("training_df$", decision_var, sep = ""))))
  colnames(ann) <- 'Decision'
  rownames(ann)<- rownames(dataM)
  ann[] <- lapply( ann, factor)
  newdf<-dataM[rowSums(dataM[])>1,colSums(dataM[])>support]
  a <- filter(ann, rownames(ann) %in% rownames(newdf))
  annoCol <- list(category = unique(ann$Decision))
  setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
  print(pheatmap(as.matrix(newdf), annotation_row=a ,main='Clustering of Model Rules',fontsize = fontsize, border_color = 'white',annotation_colors = annoCol,cluster_cols = TRUE,show_rownames = show_rownames, show_colnames = show_colnames, cluster_rows = TRUE,color = c('grey88','gray39'),legend_breaks = c(0,1)))
  setHook("grid.newpage", NULL, "replace")
  grid.text("Rules", y=-0.07, gp=gpar(fontsize=15))
  grid.text("Visits", x=-0.07, rot=90, gp=gpar(fontsize=15))
}

cluster_rules(project_4_sig, rec.ros.P4 ,support=20)
# -------------------------------
# Visualizing the Model with VisuNet
# -------------------------------
vis <- visunet(ros.P4$main)
head(vis$Asian$nodes[order(-vis$Asian$nodes$NodeConnection), ])
# Extract unique gene symbols from vis$all$nodes$label
gene_symbols <- unique(as.character(vis$all$nodes$label))
# Check unique labels in ros.P4$main
print("Checking if weight.kg. exists in ros.P4$main:")
print(any(ros.P4$main$label == "weight.kg."))
# Remove "weight.kg." from vis$all$nodes
vis$all$nodes <- vis$all$nodes %>% filter(label != "weight.kg.")

# Verify that it is gone
print(any(vis$all$nodes$label == "weight.kg."))  # Should return FALSE

# -------------------------------
# Data Cleaning for Gene Analysis
# -------------------------------
# Extract unique gene symbols from VisuNet nodes
gene_symbols <- unique(as.character(vis$all$nodes$label))

# Remove experimental conditions from gene names
cleaned_gene_symbols <- gsub("_.*", "", gene_symbols)

# Check which genes match official symbols
valid_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")
matched_cleaned_genes <- cleaned_gene_symbols[cleaned_gene_symbols %in% valid_symbols]

# -------------------------------
# Convert Gene Symbols to Entrez IDs
# -------------------------------
set.seed(1000)
gene_entrez <- bitr(matched_cleaned_genes,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)

# -------------------------------
# Gene Ontology (GO) Enrichment Analysis
# -------------------------------
genes_GO <- groupGO(gene = unique(gene_entrez$ENTREZID),
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",   # Biological Process
                    level = 3,    # Less specific than level 5
                    readable = TRUE)
genes_GO
# -------------------------------
# Visualizing GO Enrichment
# -------------------------------
genes_GO_df <- as.data.frame(genes_GO)

# Filter GO terms where Count > 1
x1 <- genes_GO_df[genes_GO_df$Count > 1, c('ID', 'Count')]

# Check if x1 is empty before plotting
if (nrow(x1) > 0) {
  x <- as.numeric(x1$Count)
  names(x) <- x1$ID
  barplot(x, col = rainbow(20), las = 2)
} else {
  print("No GO terms have a count greater than 1. Barplot skipped.")
}

# -------------------------------
# Reactome Pathway Enrichment Analysis
# -------------------------------
# Perform Reactome pathway enrichment
genes_reactome <- enrichPathway(gene = unique(gene_entrez$ENTREZID),
                                organism = "human", 
                                pvalueCutoff = 0.05, 
                                readable = TRUE)

# Print Reactome results
print("Reactome Pathway Enrichment Results:")
print(genes_reactome)

# -------------------------------
# Reactome Pathway Visualization
# -------------------------------
barplot(genes_reactome, showCategory = 10, title = "Top 10 Enriched Reactome Pathways")
dotplot(genes_reactome, showCategory = 10, title = "Dotplot of Enriched Reactome Pathways")

########################################3
# -------------------------------
# KEGG Pathway Enrichment Analysis
# -------------------------------
genes_KEGG <- enrichKEGG(gene = unique(gene_entrez$ENTREZID),
                         organism = "hsa", # Human genes
                         pvalueCutoff = 0.05)

# Print KEGG results
print("KEGG Pathway Enrichment Results:")
print(genes_KEGG)

# -------------------------------
# Visualizing KEGG Pathway Enrichment
# -------------------------------
barplot(genes_KEGG, showCategory = 10, title = "KEGG Pathway Enrichment")
dotplot(genes_KEGG, showCategory = 10, title = "KEGG Dotplot")


# -------------------------------
# Disease Associations from DisGeNET
# -------------------------------
genes_disease <- enrichDO(
  gene = unique(gene_entrez$ENTREZID),  # Use the unique Entrez IDs from KEGG analysis
  ont = "DO",  # Disease Ontology
  pvalueCutoff = 0.05,  # Significance threshold
  readable = TRUE  # Convert gene IDs to gene symbols
)

# View results
print("Disease Ontology Enrichment Results:")
print(genes_disease)

# -------------------------------
# Visualization of Disease Enrichment
# -------------------------------
dotplot(genes_disease, showCategory = 10, title = "Disease Ontology Enrichment")
barplot(genes_disease, showCategory = 10, title = "Disease Ontology Enrichment")

# -------------------------------
# Highlighting Most Enriched Genes in VisuNet
# -------------------------------
genes_top_GO <- strsplit(genes_GO_df[which.max(genes_GO_df$Count), "geneID"], "/")[[1]]

# Update node visualization in VisuNet
nodes_RNO <- vis$all$nodes
nodes_RNO$shape <- rep("dot", length(nodes_RNO$label))
nodes_RNO$shape[which(as.character(nodes_RNO$label) %in% genes_top_GO)] <- "star"

# Create the node object list and rerun VisuNet
nodesL <- list(nodes = nodes_RNO, CustCol = c("shape"))
vis_out2 <- visunet(ros.P4$main, CustObjectNodes = nodesL)

visuArc <- function(net, decision, feature=NULL){
  #load libraries if you only copy the function
  library(VisuNet)
  if (!require(arcdiagram)) install.packages('arcdiagram')
  library(arcdiagram)
  
  #make df:  a data frame containing network connection in columns
  # df<- eval(parse(text=paste("net$", decision,'$nodes' ,sep = "")))
  df<-net[[decision]][['nodes']]
  cols_select<- c('label','DiscState','NodeConnection')
  df<- subset(df,select = cols_select)
  rownames(df)  <-make.names(df[,'label'], unique = TRUE)
  df <- df[,-1]
  
  #scale Node connection fro 0 to 1
  scale_values <- function(x){(x-min(x))/(max(x)-min(x))}
  df$NodeConnection <- scale_values(df$NodeConnection)
  
  
  #check for feature
  if(is.null(feature)){
    feature <- rownames(df[which.max(df[,'NodeConnection']),])
  }
  
  discLev<-df[feature,'DiscState']
  feature2 <- gsub("\\s*\\([^\\)]+\\)","",as.character(feature))
  edges <- net[[decision]]$edges[,1:4]
  edges_from <- as.character(edges$from)
  edges_to <- as.character(edges$to) ## gsub("\\=.*","",as.character(net$short$edges$to))
  indx <- unique(c(which(edges_from == paste0(feature2,"=",discLev)),which(edges_to == paste0(feature2,"=",discLev))))
  edges2 <- edges[indx,]
  M <- as.matrix(edges2[,1:2])
  
  #colors
  valsCols <- round(edges2$connNorm, digits = 1)*10
  edgesCon <- as.numeric(as.factor(as.character(valsCols)))
  colors <- colorRampPalette(c("gainsboro","lavender","darkorchid3"))(length(unique(valsCols)))[edgesCon]
  labelsNodes <- unique(c(t(M)))
  colNodes <- c("#56B4E9","#999999","#E69F00")[as.numeric(sub(".*=","",labelsNodes))] #EDIT THIS WITH YOUR COLOR CODING
  
  #nodes values
  nodesDlev<- df[,'DiscState']
  nodesNams <- gsub("\\s*\\([^\\)]+\\)","",as.character(rownames(df)))
  nodes2 <- paste0(nodesNams,"=",nodesDlev)
  nodeSize <- round(df[match(labelsNodes, nodes2),2],digits = 1)*10
  
  ordV <- c(1,order(edgesCon, decreasing = T)+1)
  
  colsLabs <- rep("black",length(labelsNodes))
  colsLabs[which(labelsNodes == paste0(feature2,"=",discLev))] <- "red"
  
  arcplot(M, lwd.arcs=edgesCon, col.arcs = colors, col.nodes = colNodes, labels=sub("=.*","",labelsNodes),
          ordering=ordV, col.labels=colsLabs, cex.labels=0.7, font=1, lwd.nodes = nodeSize, horizontal = FALSE,main = decision,adj=1)
}


#visuarc
visuArc(vis,'Cau',feature='IFIT2_IFNb_4')

#############################################
# Check if top rule features appear in VisuArc
visuArc(vis, 'Cau', feature='IFIT2_IFNb_4')
visuArc(vis, 'Asian', feature='UTS2_Activated_48')
visuArc(vis, 'Afro', feature='CYBB_Activated_48')

