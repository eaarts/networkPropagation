# rank genes on top 10 related mouse phenotypes

'%notin%' = Negate('%in%')

library(igraph)
library(tidyverse)
library(pROC)
library(foreach)
library(doParallel)
library(ComplexHeatmap)
library(ggpubr)

source('Code/networkPropagation.R')

# load files ----

distTraitsMatrix = readRDS('data/distTraits.rds')

traitAnnotation = read.csv('data/traitOverview.csv')

variantsWithHPOandMP = read.csv('data/variantsCiliopathyMP.csv')

pageRankScores = readRDS('data/pagerankScores.rds')


# define ciliopathies ----

diseasesCilia = traitAnnotation$Var1[traitAnnotation$ciliopathy == T & !is.na(traitAnnotation$ciliopathy)]

# load network ----

#load open targets interaction network (IntAct, Reactome, SIGNOR, STRING)
intAll <- read.csv('./Datasets/interaction/interactionAll.csv') #data from open targets (https://ftp.ebi.ac.uk/pub/databases/IntAct/various/ot_graphdb/2022-07-22/)
#intAll <- read.csv('../../../Datasets/interaction/interactionAll.csv')

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

# define related mouse phenotypes per ciliopathy ----

distTraitsMatrix = distTraitsMatrix[rownames(distTraitsMatrix) %notin% c('MP_0002169','MP_0003171','MP_0003175','MP_0003176'),] #exclude normal phenotypes
distTraitsMatrix = distTraitsMatrix[,colnames(distTraitsMatrix) %notin% c('MP_0002169','MP_0003171','MP_0003175','MP_0003176')] #exclude normal phenotypes

#only focus on mouse phenotypes with at least 10 genes associated to the phenotype
variantsWithHPOandMP_10 = variantsWithHPOandMP[variantsWithHPOandMP$diseaseId %in% names(table(variantsWithHPOandMP$diseaseId)[table(variantsWithHPOandMP$diseaseId) >= 10]),]

distTraitsMatrix = distTraitsMatrix[rownames(distTraitsMatrix) %in% c(unique(variantsWithHPOandMP_10[grep('MP_', variantsWithHPOandMP_10$diseaseId),"diseaseId"]),diseasesCilia),] 
distTraitsMatrix = distTraitsMatrix[,colnames(distTraitsMatrix) %in% c(unique(variantsWithHPOandMP_10[grep('MP_', variantsWithHPOandMP_10$diseaseId),"diseaseId"]),diseasesCilia)]

relatedDiseasesPerCiliopathy <- list()
for (i in 1:length(diseasesCilia)) {
  distCiliopathy <-
    distTraitsMatrix[grep('MP_', rownames(distTraitsMatrix)), grep(diseasesCilia[i], colnames(distTraitsMatrix))]
  distCiliopathy <-
    distCiliopathy[names(distCiliopathy) %notin% diseasesCilia]
  if (length(distCiliopathy != 0)) {
    relatedDiseasesPerCiliopathy[[i]] <-
      names(distCiliopathy[order(distCiliopathy)][1:10])
  }
}

relatedDiseasesPerCiliopathy = lapply(relatedDiseasesPerCiliopathy, unique)
names(relatedDiseasesPerCiliopathy) = diseasesCilia
relatedDiseasesPerCiliopathy = relatedDiseasesPerCiliopathy[sapply(relatedDiseasesPerCiliopathy, length) > 0]

relatedMouse = unique(unlist(as.data.frame(relatedDiseasesPerCiliopathy)))

#select propagation scores of mouse phenotypes and ciliopathies
pageRankScoresSelected = pageRankScores[rownames(pageRankScores) %in% c(diseasesCilia, relatedMouse),]
pageRankScoresScaled = scale(t(pageRankScoresSelected))
pageRankScoresSelectedT = t(pageRankScoresSelected)

#check overlapping genes between ciliopathies and mouse phenotypes ----

variantsCiliopathy = read.csv('data/variantsCiliopathies.csv') 

variantsRelatedMP = variantsWithHPOandMP[variantsWithHPOandMP$diseaseId %in% relatedMouse,]

variantsRelatedMP$ciliopathy = variantsRelatedMP$targetId %in% variantsCiliopathy$targetId

fractionMPCiliopathy = sum(unique(variantsRelatedMP[,c("targetId","ciliopathy")])[,2])/length(unique(variantsRelatedMP[,c("targetId","ciliopathy")])[,2])*100

MPSeedCiliopathy = as.data.frame(table(variantsRelatedMP[,c("diseaseId","ciliopathy")]))

# run ranking ----
rankDFAll = data.frame()
for (i in 1:length(relatedDiseasesPerCiliopathy)) {
  pageRankScoresRelatedMP = pageRankScoresSelectedT[, colnames(pageRankScoresSelectedT) %in% relatedDiseasesPerCiliopathy[[i]]]
  
  pageRankScoresRank = apply(
    pageRankScoresRelatedMP,
    2,
    FUN = function(x) {
      rank(desc(x))
    }
  )
  
  pageRankScoresRank = as.data.frame(pageRankScoresRank)
  pageRankScoresRank$totalRank = unname(apply(
    pageRankScoresRank,
    1,
    FUN = function(x) {
      10 ^ (sum(log10(x)) / ncol(pageRankScoresRank)) #to prevent inf values, use sum of log10 instead of product
    }
  ))
  
  pageRankScoresRank$quantile1.5percent = pageRankScoresRank$totalRank < quantile(pageRankScoresRank$totalRank, 0.015)
  pageRankScoresRank$quantile5percent = pageRankScoresRank$totalRank < quantile(pageRankScoresRank$totalRank, 0.05)
  
  rankDF = data.frame(
    'gene' = rownames(pageRankScoresRank),
    'rank' = pageRankScoresRank$totalRank,
    'quantile1.5percent' = pageRankScoresRank$quantile1.5percent,
    'quantile5percent' = pageRankScoresRank$quantile5percent,
    'ciliopathy' = rep(names(relatedDiseasesPerCiliopathy[i]), nrow(pageRankScoresRank))
  )
  
  rankDFAll = rbind(rankDFAll, rankDF)
  
}

write.csv(rankDFAll, 'data/rankMP.csv')

rankDFAll_top100 = rankDFAll[rankDFAll$rank <= 100,]
rankDFAll_top100$ciliaryGene = rankDFAll_top100$gene %in% variantsCiliopathy$targetId

#how many genes are associated to the mouse phenotypes?
donutData = data.frame(category = c('yes', 'no'), count = c(sum(unique(rankDFAll_top100$gene) %in% variantsRelatedMP$targetId), sum(unique(rankDFAll_top100$gene) %notin% variantsRelatedMP$targetId)))
donutData$fraction = donutData$count/sum(donutData$count)
donutData$ymax = cumsum(donutData$fraction)
donutData$ymin = c(0, head(donutData$ymax, n=-1))
donutData$labelPosition <- (donutData$ymax + donutData$ymin) / 2
donutData$label <- paste0(donutData$category, "\n value: ", donutData$count)

ggplot(donutData, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  scale_fill_manual(values = c('#3A53A4', '#BE202E')) +
  coord_polar(theta="y") +
  xlim(c(1, 4)) +
  theme_void() +
  theme(legend.position = "none")

#how many genes are ciliopathy genes?
donutData = data.frame(category = c('yes', 'no'), count = c(sum(unique(rankDFAll_top100$gene) %in% variantsCiliopathy$targetId), sum(unique(rankDFAll_top100$gene) %notin% variantsCiliopathy$targetId)))
donutData$fraction = donutData$count/sum(donutData$count)
donutData$ymax = cumsum(donutData$fraction)
donutData$ymin = c(0, head(donutData$ymax, n=-1))
donutData$labelPosition <- (donutData$ymax + donutData$ymin) / 2
donutData$label <- paste0(donutData$category, "\n value: ", donutData$count)

ggplot(donutData, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  scale_fill_manual(values = c('#3A53A4', '#BE202E')) +
  coord_polar(theta="y") +
  xlim(c(1, 4)) +
  theme_void() +
  theme(legend.position = "none")

