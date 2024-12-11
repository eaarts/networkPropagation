# predict genes from related traits with different numbers of traits

'%notin%' = Negate('%in%')

library(igraph)
library(tidyverse)
library(pROC)
library(foreach)
library(doParallel)
library(ComplexHeatmap)
library(ggpubr)

source('code/0.networkPropagation.R')

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

distTraitsMatrix = distTraitsMatrix[rownames(distTraitsMatrix) %notin% c('MP_0002169','MP_0003171','MP_0003175','MP_0003176'),] #exclude normal phenotypes
distTraitsMatrix = distTraitsMatrix[,colnames(distTraitsMatrix) %notin% c('MP_0002169','MP_0003171','MP_0003175','MP_0003176')] #exclude normal phenotypes

#only focus on mouse phenotypes with at least 10 genes associated to the phenotype
variantsWithHPOandMP_10 = variantsWithHPOandMP[variantsWithHPOandMP$diseaseId %in% names(table(variantsWithHPOandMP$diseaseId)[table(variantsWithHPOandMP$diseaseId) >= 10]),]

distTraitsMatrix = distTraitsMatrix[rownames(distTraitsMatrix) %in% c(unique(variantsWithHPOandMP_10[grep('MP_', variantsWithHPOandMP_10$diseaseId),"diseaseId"]),diseasesCilia),] 
distTraitsMatrix = distTraitsMatrix[,colnames(distTraitsMatrix) %in% c(unique(variantsWithHPOandMP_10[grep('MP_', variantsWithHPOandMP_10$diseaseId),"diseaseId"]),diseasesCilia)]

# ROC with varying number of traits ----

rocResAll_traitNumber = list()
nTraits_traitNumber = list()
for (numberTraits in seq(10,100, 10)){
  
  relatedDiseasesPerCiliopathy <- list()
  for (i in 1:length(diseasesCilia)) {
    distCiliopathy <-
      distTraitsMatrix[grep('MP_', rownames(distTraitsMatrix)), grep(diseasesCilia[i], colnames(distTraitsMatrix))]
    distCiliopathy <-
      distCiliopathy[names(distCiliopathy) %notin% diseasesCilia]
    if (length(distCiliopathy != 0)) {
      relatedDiseasesPerCiliopathy[[i]] <-
        names(distCiliopathy[order(distCiliopathy)][1:numberTraits])
    }
  }
  
  relatedDiseasesPerCiliopathy = lapply(relatedDiseasesPerCiliopathy, unique)
  names(relatedDiseasesPerCiliopathy) = diseasesCilia
  relatedDiseasesPerCiliopathy = relatedDiseasesPerCiliopathy[sapply(relatedDiseasesPerCiliopathy, length) > 0]
  
  nTraits = c()
  rocResAll = c()
  for (i in 1:length(relatedDiseasesPerCiliopathy)){
    relatedTraits = relatedDiseasesPerCiliopathy[[i]]
    
    variantsMP = unique(variantsWithHPOandMP[variantsWithHPOandMP$diseaseId %in% relatedTraits, ])
    variantsMPNoCiliopathy = variantsMP[variantsMP$targetId %notin% variantsWithHPOandMP[variantsWithHPOandMP$diseaseId %in% names(relatedDiseasesPerCiliopathy), "targetId"],]
    variantsMPNoCiliopathy = unique(variantsMPNoCiliopathy)
    
    pageRankScoresRelatedMP <- lapply(relatedTraits, FUN = function(x){getPageRank(x, intGraphClean, variantsMPNoCiliopathy)})
    names(pageRankScoresRelatedMP) = relatedTraits
    pageRankScoresRelatedMP = pageRankScoresRelatedMP[!is.na(pageRankScoresRelatedMP)]
    
    pageRankDFRelatedMP <- as.data.frame(matrix(ncol = length(pageRankScoresRelatedMP), nrow = vcount(intGraphClean)))
    for (k in 1:length(pageRankScoresRelatedMP)){
      pageRankDFRelatedMP[,k] = pageRankScoresRelatedMP[[k]]$pageRank
    }
    
    rownames(pageRankDFRelatedMP) = V(intGraphClean)$name
    colnames(pageRankDFRelatedMP) = names(pageRankScoresRelatedMP)
    
    pageRankDFRelatedMPRank = apply(
      pageRankDFRelatedMP,
      2,
      FUN = function(x) {
        rank(dplyr::desc(x))
      }
    )
    
    pageRankDFRelatedMPRank = as.data.frame(pageRankDFRelatedMPRank)
    pageRankDFRelatedMPRank$totalRank = unname(apply(
      pageRankDFRelatedMPRank,
      1,
      FUN = function(x) {
        10^(sum(log10(x))/ncol(pageRankDFRelatedMPRank)) #to prevent inf values, use sum of log10 instead of product
      }
    ))
    
    pageRankDFRelatedMPRank = pageRankDFRelatedMPRank[order(pageRankDFRelatedMPRank$totalRank, decreasing = F),]
    
    geneList = pageRankDFRelatedMPRank$totalRank
    names(geneList) = rownames(pageRankDFRelatedMPRank)
    
    rocRes = roc(names(geneList) %in% variantsWithHPOandMP[variantsWithHPOandMP$diseaseId == names(relatedDiseasesPerCiliopathy)[[i]],"targetId"],geneList)$auc
    
    rocResAll = c(rocResAll, rocRes)
    
    nTraits = c(nTraits, length(pageRankDFRelatedMP))
  }
  
  rocResAll_traitNumber[[as.character(numberTraits)]] = rocResAll
  
  nTraits_traitNumber[[as.character(numberTraits)]] = nTraits
  
}

rocResAll_traitNumber_DF = as.data.frame(rocResAll_traitNumber)
rocResAll_traitNumber_DF = gather(rocResAll_traitNumber_DF, key = 'nTrait', value = 'ROC')
rocResAll_traitNumber_DF$nTrait = factor(rocResAll_traitNumber_DF$nTrait, levels = c('X10','X20', 'X30', 'X40', 'X50', 'X60', 'X70', 'X80','X90', 'X100'))
ggplot(rocResAll_traitNumber_DF, aes(x=nTrait, y=ROC)) + geom_boxplot() + theme_classic()

# define related mouse phenotypes per ciliopathy ----

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

# calculate ROC values for predicting ciliopathy genes based on rank ----

nTraits = c()
rocResAll = c()
for (i in 1:length(relatedDiseasesPerCiliopathy)){
  relatedTraits = relatedDiseasesPerCiliopathy[[i]]
  
  variantsMP = unique(variantsWithHPOandMP[variantsWithHPOandMP$diseaseId %in% relatedTraits, ])
  variantsMPNoCiliopathy = variantsMP[variantsMP$targetId %notin% variantsWithHPOandMP[variantsWithHPOandMP$diseaseId %in% names(relatedDiseasesPerCiliopathy), "targetId"],]
  variantsMPNoCiliopathy = unique(variantsMPNoCiliopathy)
  
  pageRankScoresRelatedMP <- lapply(relatedTraits, FUN = function(x){getPageRank(x, intGraphClean, variantsMPNoCiliopathy)})
  names(pageRankScoresRelatedMP) = relatedTraits
  pageRankScoresRelatedMP = pageRankScoresRelatedMP[!is.na(pageRankScoresRelatedMP)]
  
  pageRankDFRelatedMP <- as.data.frame(matrix(ncol = length(pageRankScoresRelatedMP), nrow = vcount(intGraphClean)))
  for (k in 1:length(pageRankScoresRelatedMP)){
    pageRankDFRelatedMP[,k] = pageRankScoresRelatedMP[[k]]$pageRank
  }
  
  rownames(pageRankDFRelatedMP) = V(intGraphClean)$name
  colnames(pageRankDFRelatedMP) = names(pageRankScoresRelatedMP)
  
  pageRankDFRelatedMPRank = apply(
    pageRankDFRelatedMP,
    2,
    FUN = function(x) {
      rank(dplyr::desc(x))
    }
  )
  
  pageRankDFRelatedMPRank = as.data.frame(pageRankDFRelatedMPRank)
  pageRankDFRelatedMPRank$totalRank = unname(apply(
    pageRankDFRelatedMPRank,
    1,
    FUN = function(x) {
      10^(sum(log10(x))/ncol(pageRankDFRelatedMPRank)) #to prevent inf values, use sum of log10 instead of product
    }
  ))
  
  pageRankDFRelatedMPRank = pageRankDFRelatedMPRank[order(pageRankDFRelatedMPRank$totalRank, decreasing = F),]
  
  geneList = pageRankDFRelatedMPRank$totalRank
  names(geneList) = rownames(pageRankDFRelatedMPRank)
  
  rocRes = roc(names(geneList) %in% variantsWithHPOandMP[variantsWithHPOandMP$diseaseId == names(relatedDiseasesPerCiliopathy)[[i]],"targetId"],geneList)$auc
  
  rocResAll = c(rocResAll, rocRes)
  
  nTraits = c(nTraits, length(pageRankDFRelatedMP))
}

rocResAllTrue = rocResAll

names(rocResAllTrue) = traitAnnotation$trait_label[match(names(relatedDiseasesPerCiliopathy), traitAnnotation$Var1)]
rocResAllTrueDF = as.data.frame(rocResAllTrue) %>% rownames_to_column(var = 'ciliopathy')
ggplot(rocResAllTrueDF, aes(y=reorder(ciliopathy, rocResAllTrue), x=rocResAllTrue)) + 
  geom_bar(stat = 'identity', color = 'black', fill = '#3A53A4') + theme_classic() + 
  geom_vline(xintercept = median(rocRandomvsTrue[rocRandomvsTrue$round == 'random', "roc"]), color = 'red', linewidth = 1)


#with random neighbor traits ----

mousePhenotypes = rownames(distTraitsMatrix)[grep('MP_', rownames(distTraitsMatrix))]

rocResAllDF = as.data.frame(matrix(nrow = 100, ncol = length(relatedDiseasesPerCiliopathy)))

for (j in 1:100){
  randomTraits = sample(mousePhenotypes, 10)
  
  rocResAll = c()
  for (i in names(relatedDiseasesPerCiliopathy)){
    variantsMP = unique(variantsWithHPOandMP[variantsWithHPOandMP$diseaseId %in% randomTraits, ])
    variantsMPNoCiliopathy = variantsMP[variantsMP$targetId %notin% variantsWithHPOandMP[variantsWithHPOandMP$diseaseId %in% names(relatedDiseasesPerCiliopathy), "targetId"],]
    variantsMPNoCiliopathy = unique(variantsMPNoCiliopathy)
    
    pageRankScoresRelatedMP <- lapply(randomTraits, FUN = function(x){getPageRank(x, intGraphClean, variantsMPNoCiliopathy)})
    names(pageRankScoresRelatedMP) = randomTraits
    pageRankScoresRelatedMP = pageRankScoresRelatedMP[!is.na(pageRankScoresRelatedMP)]
    
    pageRankDFRelatedMP <- as.data.frame(matrix(ncol = length(pageRankScoresRelatedMP), nrow = vcount(intGraphClean)))
    for (k in 1:length(pageRankScoresRelatedMP)){
      pageRankDFRelatedMP[,k] = pageRankScoresRelatedMP[[k]]$pageRank
    }
    
    rownames(pageRankDFRelatedMP) = V(intGraphClean)$name
    colnames(pageRankDFRelatedMP) = names(pageRankScoresRelatedMP)
    
    pageRankDFRelatedMPRank = apply(
      pageRankDFRelatedMP,
      2,
      FUN = function(x) {
        rank(dplyr::desc(x))
      }
    )
    
    pageRankDFRelatedMPRank = as.data.frame(pageRankDFRelatedMPRank)
    pageRankDFRelatedMPRank$totalRank = unname(apply(
      pageRankDFRelatedMPRank,
      1,
      FUN = function(x) {
        10^(sum(log10(x))/ncol(pageRankDFRelatedMPRank)) #to prevent inf values, use sum of log10 instead of product
      }
    ))
    
    #pageRankDFRelatedMPRank$totalRank[is.infinite(pageRankDFRelatedMPRank$totalRank)] = max(pageRankDFRelatedMPRank$totalRank[!is.infinite(pageRankDFRelatedMPRank$totalRank)])
    
    pageRankDFRelatedMPRank = pageRankDFRelatedMPRank[order(pageRankDFRelatedMPRank$totalRank, decreasing = F),]
    
    geneList = pageRankDFRelatedMPRank$totalRank
    names(geneList) = rownames(pageRankDFRelatedMPRank)
    
    rocRes = roc(names(geneList) %in% variantsWithHPOandMP[variantsWithHPOandMP$diseaseId == i,"targetId"],geneList)$auc
    
    rocResAll = c(rocResAll, rocRes)
  }
  rocResAllDF[j,] = rocResAll
}

colnames(rocResAllDF) = names(relatedDiseasesPerCiliopathy)

rocResAllDF_median = apply(rocResAllDF, 2, median)

rocRandomvsTrue = data.frame('roc' = c(rocResAllTrue, rocResAllDF_median), 'round' = c(rep('true', 21), rep('random', 21)))
ggplot(rocRandomvsTrue, aes(x = round, y=roc)) + geom_boxplot() + theme_classic() + stat_compare_means(method = 't.test')

