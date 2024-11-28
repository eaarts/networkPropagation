# candidate gene identification

setwd('W:/GROUP/Users/Ellen/NetworkPropagation/')

library(igraph)
library(doParallel)
library(foreach)
library(tidyverse)

traitAnnotation = read.csv('20231206_traitOverview_JBTSnoOverlap.csv')
variantsWithHPOandMP = read.csv('20240105_variantsWithHPOandMP_updatedJBTS.csv')
pageRankScores = readRDS('20231206_pageRankScores_JBTSnoOverlap.rds')

'%notin%' = Negate('%in%')

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

# select ciliopathies ----

disease = traitAnnotation$Var1[traitAnnotation$ciliopathy == T]

pageRankDF = pageRankScores[rownames(pageRankScores) %in% disease,]
pageRankDF = pageRankDF[match(disease, rownames(pageRankDF)),]

# run 1000x permutation ----

pageRankResRandom = as.data.frame(matrix(nrow = 1000, ncol = length(V(intGraphClean))))
colnames(pageRankResRandom) = V(intGraphClean)$name

randomPageRank <- function(disease, geneVariant, intGraphClean) {
  pageRankResRand = as.data.frame(matrix(nrow = 1000, ncol = length(V(intGraphClean))))
  colnames(pageRankResRand) = V(intGraphClean)$name
  for (i in 1:1000) {
    V(intGraphClean)$targetGene = 0
    variantList <-
      unique(geneVariant[geneVariant$diseaseId == paste0(disease), 'targetId'])
    randomTargets = sample(length(V(intGraphClean)), length(variantList))
    V(intGraphClean)$targetGene[randomTargets] = 1
    
    pageRank = page_rank(intGraphClean, personalized = V(intGraphClean)$targetGene)
    
    pageRankResRand[i, ] = pageRank$vector
  }
  return(pageRankResRand)
}

cores <- detectCores()
cl <- makeCluster(cores - 1)
registerDoParallel(cl)

pageRankResRandom <- foreach(i=iter(disease, by = 'cell'), .packages = c('igraph', 'dplyr')) %dopar% {
  source('Code/networkPropagation.R', local = T)
  randomPageRank(disease = i, geneVariant = variantsWithHPOandMP, intGraphClean = intGraphClean)
}

stopCluster(cl)

names(pageRankResRandom) = disease

pageRankDFPvalue = as.data.frame(matrix(nrow = nrow(pageRankDF), ncol = ncol(pageRankDF)))
for (i in 1:ncol(pageRankDF)){
  for (j in 1:nrow(pageRankDF)){
    pageRankDFPvalue[j,i] = sum(pageRankDF[j,i] > pageRankResRandom[[j]][,i]) / 1000
  }
}

colnames(pageRankDFPvalue) = colnames(pageRankDF)
rownames(pageRankDFPvalue) = rownames(pageRankDF)

pageRankDFPvalue = t(pageRankDFPvalue)

# select candidate genes ----

pageRankDFWoSeeds = pageRankDFPvalue[rownames(pageRankDFPvalue) %notin% variantsCiliopathy$targetId,]
candidateDFWoSeedsLong = pageRankDFWoSeeds %>% as.data.frame %>% rownames_to_column(var = 'geneID') %>% gather(., key = 'diseaseID', value = 'pvalue', -geneID)
candidateDFWoSeedsLong = candidateDFWoSeedsLong[candidateDFWoSeedsLong$pvalue >= 0.995, ]

candidateDFWoSeedsLong$symbol = AnnotationDbi::mapIds(org.Hs.eg.db, keys = candidateDFWoSeedsLong$gene, column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first')
candidateDFWoSeedsLong$diseaseName = traitAnnotation$trait_label[match(candidateDFWoSeedsLong$diseaseID, traitAnnotation$Var1)]

write.csv(candidateDFWoSeedsLong, '20240108_candidateGenes_pvalue0995_updatedJBTS.csv')

