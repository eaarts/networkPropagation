# network propagation for ciliopathies and mouse phenotypes

# Set working directory ----

setwd('W:/GROUP/Users/Ellen/NetworkPropagation/')

# Libraries ---- 

library(igraph)
library(pROC)
library(doParallel)
library(foreach)
library(org.Hs.eg.db)
library(tidyverse)
library(clusterProfiler)
library(ComplexHeatmap)

source('Code/networkPropagation.R')

'%notin%' = Negate('%in%')

# Interaction network full ----

#load open targets interaction network (IntAct, Reactome, SIGNOR, STRING)
intAll <- read.csv('./Datasets/interaction/interactionAll.csv')

#set threshold for STRING
intString = grep('string', intAll$sourceDatabase)
lowString = intString[intAll$scoring[intString] < 0.4]
intHigh = intAll[-lowString,]

intHigh = intHigh[grep('ENSG0', intHigh$targetA),]
intHigh = intHigh[grep('ENSG0', intHigh$targetB),]
intHigh = intHigh[!duplicated(intHigh[,c("targetA", "targetB")]),]

#create graph
intGraph = graph_from_data_frame(intHigh[,c("targetA", "targetB")], directed = F)
intGraphClean = igraph::simplify(intGraph, remove.loops = T, remove.multiple = T, edge.attr.comb = c(weight = 'max', 'ignore'))

#remove big files to save memory
rm(intAll, intGraph, intString, lowString)
gc()

# Run pageRank ----
variantsWithHPOandMP <- read.csv()

diseasesAll = unique(variantsWithHPOandMP$diseaseId)
diseasesAll = diseasesAll[grep('EFO|MONDO|HP|MP|Orphanet', diseasesAll)] #select IDs from DOID, EFO, HP, MONDO, Orphanet, MP, GO

pageRankScoresAllTraits <- lapply(diseasesAll, FUN = function(x){getPageRank(x, intGraphClean, variantsWithHPOandMP)})
names(pageRankScoresAllTraits) <- diseasesAll
pageRankScoresAllTraits = pageRankScoresAllTraits[!is.na(pageRankScoresAllTraits)]

pageRankDFAllTraits <- data.frame(t(sapply(pageRankScoresAllTraits, FUN = function(x){c(x$pageRank)})))

colnames(pageRankDFAllTraits) = V(intGraphClean)$name

pageRankDFAllTraits = pageRankDFAllTraits[!(rownames(pageRankDFAllTraits) %in% cilioEFO$ID & rownames(pageRankDFAllTraits) %notin% diseasesCilia),] #exclude extra ciliopathies

#cluster disorders based on pageRank scores
distTraits = dist(pageRankDFAllTraits, method = 'manhattan')