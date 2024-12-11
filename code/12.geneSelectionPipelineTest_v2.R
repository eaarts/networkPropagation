# recover ciliary genes with pipeline for candidate gene selection

library(igraph)
library(doParallel)
library(foreach)
library(tidyverse)

traitAnnotation = read.csv('data/traitOverview.csv')
variantsWithHPOandMP = read.csv('data/variantsCiliopathyMP.csv')
pageRankScores = readRDS('data/pagerankScores.rds')
variantsCiliopathy = read.csv('data/variantsCiliopathies.csv') 

'%notin%' = Negate('%in%')

pageRankScoresMP = pageRankScores[grep('MP_', rownames(pageRankScores)),]

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

# select train and test genes ----

randomPageRank <- function(disease, geneVariant, intGraphClean) {
  pageRankResRand = sapply(c(1:1000), FUN = function(x){
    V(intGraphClean)$targetGene = 0
    randomTargets = sample(length(V(intGraphClean)), 10)
    V(intGraphClean)$targetGene[randomTargets] = 1
    
    pageRank = page_rank(intGraphClean, personalized = V(intGraphClean)$targetGene)
    
    return(pageRank$vector)
  })
}

pageRankResRandom <- randomPageRank(disease = disease, geneVariant = variantsCiliopathy, intGraphClean = intGraphClean)

rnaSingleCell = read_tsv('data/rna_single_cell_type.tsv')

rnaSingleCell_Wide = spread(rnaSingleCell[,c("Gene", "Cell type","nTPM")], key = 'Cell type', value = 'nTPM')
rnaSingleCell_Wide = column_to_rownames(rnaSingleCell_Wide, var = 'Gene')
rnaSingleCell_Wide = rnaSingleCell_Wide[rowSums(rnaSingleCell_Wide) != 0,]

rnaSingleCell_WideScaled = t(scale(t(rnaSingleCell_Wide)))
rnaSingleCell_WideScaled = as.data.frame(rnaSingleCell_WideScaled)

variantsWithHPOandMP_10 = variantsWithHPOandMP[variantsWithHPOandMP$diseaseId %in% names(table(variantsWithHPOandMP$diseaseId)[table(variantsWithHPOandMP$diseaseId) >= 10]),]

allScoresAll_allTraits = data.frame()

allScoresAllPosNeg_allTraits = data.frame()

allScoresAllPosNeg_allTraits_scaled = data.frame()

allScoresAll_allTraits_scaled = data.frame()

for (disease in c('MONDO_custom',
                  'MONDO_0018772',
                  'MONDO_0015229',
                  'MONDO_0018998')) {
  genes = unique(variantsCiliopathy[variantsCiliopathy$diseaseId == disease, "targetId"])
  
  set.seed(12345)
  negGenes = sample(colnames(pageRankScores)[colnames(pageRankScores) %notin% genes], 150)
  
  nSets = floor(length(genes) / 5)
  
  usedTrainSets = c()
  
  allScoresAll = data.frame()
  
  allScoresAll_scaled = data.frame()
  
  for (i in 1:nSets) {
    genes_train = sample(genes[genes %notin% usedTrainSets], 5)
    
    usedTrainSets = c(usedTrainSets, genes_train)
    
    genes_test = genes[genes %notin% genes_train]
    
    # run page rank ----
    
    V(intGraphClean)$targetGene = 0
    V(intGraphClean)$targetGene[V(intGraphClean)$name %in% genes_train] = 1
    
    pageRankTrue = page_rank(intGraphClean, personalized = V(intGraphClean)$targetGene)
    
    # run 1000x permutation ----
    
    pageRankDFPvalue = c()
    for (j in 1:length(pageRankTrue$vector)) {
      pageRankDFPvalue[j] = sum(pageRankTrue$vector[j] > pageRankResRandom[j, ]) / 1000
    }
    
    names(pageRankDFPvalue) = names(pageRankTrue$vector)
    
    # use mouse phenotypes for ranking ----
    pageRankTrueScaled = scale(pageRankTrue$vector)
    pageRankScoresMPScaled = scale(t(pageRankScoresMP))
    
    pageRankScoresMPScaled = pageRankScoresMPScaled[, colnames(pageRankScoresMPScaled) %notin% c('MP_0002169', 'MP_0003171', 'MP_0003175', 'MP_0003176')]
    pageRankScoresMPScaled = pageRankScoresMPScaled[, colnames(pageRankScoresMPScaled) %in% unique(variantsWithHPOandMP_10[grep('MP_', variantsWithHPOandMP_10$diseaseId),"diseaseId"])]
    
    distMP = apply(pageRankScoresMPScaled, 2, function(x)
      sqrt(sum((
        x - pageRankTrueScaled
      ) ^ 2)))
    
    relatedDiseasesPerCiliopathy <- names(distMP[order(distMP)][1:10])
    
    pageRankScoresRelated = pageRankScoresMP[rownames(pageRankScoresMP) %in% relatedDiseasesPerCiliopathy, ]
    pageRankScoresRelated = t(pageRankScoresRelated)
    
    pageRankScoresRank = apply(
      pageRankScoresRelated,
      2,
      FUN = function(x) {
        rank(dplyr::desc(x))
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
    
    # gene expression ----
    
    #expression of all ciliopathy genes
    
    rnaSingleCell_WideScaled$ciliopathyGene = rownames(rnaSingleCell_WideScaled) %in% variantsCiliopathy$targetId
    rnaSingleCell_WideScaled_train = rnaSingleCell_WideScaled[rownames(rnaSingleCell_WideScaled) %notin% genes_test, ]
    
    anovaRes = aov(ciliopathyGene ~ . , rnaSingleCell_WideScaled_train)
    tableAnova = summary(anovaRes)[[1]]
    cellsForModel = rownames(tableAnova)[tableAnova$`Pr(>F)` < 0.01]
    cellsForModel = gsub('`', '', str_trim(cellsForModel))
    
    rnaSingleCell_WideScaled_selected = rnaSingleCell_WideScaled_train[, colnames(rnaSingleCell_WideScaled_train) %in% cellsForModel]
    rnaSingleCell_WideScaled_selected$ciliopathyGene = rownames(rnaSingleCell_WideScaled_selected) %in% variantsCiliopathy$targetId
    
    model = glm(ciliopathyGene ~ .,
                family = binomial(link = 'logit'),
                data = rnaSingleCell_WideScaled_selected)
    
    expressionScore = predict(model, newdata = rnaSingleCell_WideScaled, type = 'response')
    expressionScore = expressionScore[match(names(pageRankDFPvalue), names(expressionScore))]
    expressionScore = expressionScore[!is.na(expressionScore)]
    
    # save scores for every round ----
    
    allGenes = Reduce(intersect, list(
      names(pageRankDFPvalue),
      names(expressionScore),
      rownames(pageRankScoresRank)
    ))
    
    allGenes = allGenes[allGenes %notin% genes_train]
    
    allScores = data.frame(
      geneExpression = expressionScore[match(allGenes, names(expressionScore))],
      pvalue = pageRankDFPvalue[match(allGenes, names(pageRankDFPvalue))],
      rank = pageRankScoresRank$totalRank[match(allGenes, rownames(pageRankScoresRank))],
      gene = allGenes,
      label = allGenes %in% genes_test
    )
    
    allScores_scaled = allScores
    allScores_scaled[, 1:3] = apply(allScores_scaled[, 1:3], 2, function(x) (x - min(x)) / (max(x) - min(x)))
    
    allScoresAll = rbind(allScoresAll, allScores)
    allScoresAll_scaled = rbind(allScoresAll_scaled, allScores_scaled)
    
  }
  
  allScoresAllPosNeg = allScoresAll[allScoresAll$gene %in% c(genes, negGenes), ]
  allScoresAllPosNeg_allTraits = rbind(allScoresAllPosNeg_allTraits, allScoresAllPosNeg)
  
  allScoresAllPosNeg_scaled = allScoresAll_scaled[allScoresAll_scaled$gene %in% c(genes, negGenes), ]
  allScoresAllPosNeg_allTraits_scaled = rbind(allScoresAllPosNeg_allTraits_scaled,
                                              allScoresAllPosNeg_scaled)
  
  allScoresAll_allTraits = rbind(allScoresAll_allTraits, allScoresAll)
  
  allScoresAll_allTraits_scaled = rbind(allScoresAll_allTraits_scaled, allScoresAll_scaled)
  
}

# overall selection of genes ----

allScoresAllPosNeg_allTraits_scaled$geneExpressionLog = log(allScoresAllPosNeg_allTraits_scaled$geneExpression)
allScoresAllPosNeg_allTraits_scaled$geneExpressionLog[is.infinite(allScoresAllPosNeg_allTraits_scaled$geneExpressionLog)] = min(allScoresAllPosNeg_allTraits_scaled$geneExpressionLog[!is.infinite(allScoresAllPosNeg_allTraits_scaled$geneExpressionLog)])

modelOverall = glm(label ~ geneExpressionLog + pvalue + rank, family = binomial(link = 'logit'), data = allScoresAllPosNeg_allTraits_scaled) 

modelOverallPvalue = glm(label ~ pvalue, family = binomial(link = 'logit'), data = allScoresAllPosNeg_allTraits_scaled) 
modelOverallExpr = glm(label ~ geneExpressionLog, family = binomial(link = 'logit'), data = allScoresAllPosNeg_allTraits_scaled) 
modelOverallRank = glm(label ~ rank, family = binomial(link = 'logit'), data = allScoresAllPosNeg_allTraits_scaled) 
modelOverallRankPvalue = glm(label ~ pvalue + rank, family = binomial(link = 'logit'), data = allScoresAllPosNeg_allTraits_scaled) 

#write.csv(allScoresAllPosNeg_allTraits_scaled, '20240517_trainScoresPosNeg_scaled_10MP_5train.csv')
#write.csv(allScoresAll_allTraits_scaled, '20240517_trainScoresAllGenes_scaled_10MP_5train.csv')

write_rds(modelOverallRankPvalue, 'data/candidateModel.rds')

# test model with other ciliopathies ----

testCiliopathies = list()
allScoresScaledTest = list()
allScoresTest = list()

for (disease in c(
  'MONDO_0019005',
  'MONDO_0018921',
  'MONDO_0015375',
  #'EFO_0008620',
  'MONDO_0015993',
  'MONDO_0015461',
  'MONDO_0017842',
  'MONDO_custom',
  'MONDO_0018772',
  'MONDO_0015229',
  'MONDO_0018998'
)) {
  genes = unique(variantsCiliopathy[variantsCiliopathy$diseaseId == disease, "targetId"])
  
  prAUCAll = data.frame()
  partialAUCAll = data.frame()
  AUCAll = data.frame()
  
  for (rundidx in 1:10) {
    genes_train = sample(genes, 5)
    
    genes_test = genes[genes %notin% genes_train]
    
    negGenes = sample(colnames(pageRankScores)[colnames(pageRankScores) %notin% genes], 1000)
    
    # run page rank ----
    
    V(intGraphClean)$targetGene = 0
    V(intGraphClean)$targetGene[V(intGraphClean)$name %in% genes_train] = 1
    
    pageRankTrue = page_rank(intGraphClean, personalized = V(intGraphClean)$targetGene)
    
    # run 1000x permutation ----
    
    pageRankDFPvalue = c()
    for (i in 1:length(pageRankTrue$vector)) {
      pageRankDFPvalue[i] = sum(pageRankTrue$vector[i] > pageRankResRandom[i, ]) / 1000
    }
    
    names(pageRankDFPvalue) = names(pageRankTrue$vector)
    
    # use mouse phenotypes for ranking ----
    pageRankTrueScaled = scale(pageRankTrue$vector)
    pageRankScoresMPScaled = scale(t(pageRankScoresMP))
    
    pageRankScoresMPScaled = pageRankScoresMPScaled[, colnames(pageRankScoresMPScaled) %notin% c('MP_0002169', 'MP_0003171', 'MP_0003175', 'MP_0003176')]
    pageRankScoresMPScaled = pageRankScoresMPScaled[, colnames(pageRankScoresMPScaled) %in% unique(variantsWithHPOandMP_10[grep('MP_', variantsWithHPOandMP_10$diseaseId),"diseaseId"])]
    
    distMP = apply(pageRankScoresMPScaled, 2, function(x)
      sqrt(sum((
        x - pageRankTrueScaled
      ) ^ 2)))
    
    relatedDiseasesPerCiliopathy <- names(distMP[order(distMP)][1:10])
    
    pageRankScoresRelated = pageRankScoresMP[rownames(pageRankScoresMP) %in% relatedDiseasesPerCiliopathy, ]
    pageRankScoresRelated = t(pageRankScoresRelated)
    
    pageRankScoresRank = apply(
      pageRankScoresRelated,
      2,
      FUN = function(x) {
        rank(dplyr::desc(x))
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
    
    # gene expression ----
    
    #expression of all ciliopathy genes
    
    rnaSingleCell_WideScaled$ciliopathyGene = rownames(rnaSingleCell_WideScaled) %in% variantsCiliopathy$targetId
    rnaSingleCell_WideScaled_train = rnaSingleCell_WideScaled[rownames(rnaSingleCell_WideScaled) %notin% genes_test, ]
    
    anovaRes = aov(ciliopathyGene ~ . , rnaSingleCell_WideScaled_train)
    tableAnova = summary(anovaRes)[[1]]
    cellsForModel = rownames(tableAnova)[tableAnova$`Pr(>F)` < 0.01]
    cellsForModel = gsub('`', '', str_trim(cellsForModel))
    
    rnaSingleCell_WideScaled_selected = rnaSingleCell_WideScaled_train[, colnames(rnaSingleCell_WideScaled_train) %in% cellsForModel]
    rnaSingleCell_WideScaled_selected$ciliopathyGene = rownames(rnaSingleCell_WideScaled_selected) %in% variantsCiliopathy$targetId
    
    model = glm(ciliopathyGene ~ .,
                family = binomial(link = 'logit'),
                data = rnaSingleCell_WideScaled_selected)
    
    expressionScore = predict(model, newdata = rnaSingleCell_WideScaled, type = 'response')
    expressionScore = expressionScore[match(names(pageRankDFPvalue), names(expressionScore))]
    expressionScore = expressionScore[!is.na(expressionScore)]
    
    # overall score
    
    allGenes = Reduce(intersect, list(
      names(pageRankDFPvalue),
      names(expressionScore),
      rownames(pageRankScoresRank)
    ))
    
    allGenes = allGenes[allGenes %notin% genes_train]
    
    allScores = data.frame(
      geneExpression = expressionScore[match(allGenes, names(expressionScore))],
      pvalue = pageRankDFPvalue[match(allGenes, names(pageRankDFPvalue))],
      rank = pageRankScoresRank$totalRank[match(allGenes, rownames(pageRankScoresRank))],
      gene = allGenes,
      label = allGenes %in% genes_test
    )
    
    allScores_scaled = allScores
    allScores_scaled[, 1:3] = apply(allScores_scaled[, 1:3], 2, function(x) (x - min(x)) / (max(x) - min(x)))
    
    allScores_scaled$geneExpressionLog = log(allScores_scaled$geneExpression)
    
    # allScores$pvalue = -log10(1 - allScores$pvalue + 0.001)
    # allScores$geneExpression = log(allScores$geneExpression)
    # allScores$rank = (allScores$rank-min(allScores$rank))/(min(allScores$rank)-max(allScores$rank))
    
    allScores$overallScore = predict(modelOverall, newdata = allScores_scaled)
    allScores$overallScorePvalue = predict(modelOverallPvalue, newdata = allScores_scaled)
    allScores$overallScoreExpr = predict(modelOverallExpr, newdata = allScores_scaled)
    allScores$overallScoreRank = predict(modelOverallRank, newdata = allScores_scaled)
    allScores$overallScoreRankPvalue = predict(modelOverallRankPvalue, newdata = allScores_scaled)
    
    allScores$overallScore[is.infinite(allScores$overallScore)] = min(allScores$overallScore[!is.infinite(allScores$overallScore)])
    allScores$overallScoreExpr[is.infinite(allScores$overallScoreExpr)] = min(allScores$overallScoreExpr[!is.infinite(allScores$overallScoreExpr)])
    
    allScoresPosNeg = allScores[rownames(allScores) %in% c(genes, negGenes),]
    
    prAUCAll[rundidx, 'overallScore'] = PRROC::pr.curve(scores.class0 = allScoresPosNeg$overallScore, weights.class0 =  allScoresPosNeg$label)$auc.integral
    prAUCAll[rundidx, 'scorePvalue'] = PRROC::pr.curve(scores.class0 = allScoresPosNeg$overallScorePvalue, weights.class0 =  allScoresPosNeg$label)$auc.integral
    prAUCAll[rundidx, 'scoreExpr'] = PRROC::pr.curve(scores.class0 = allScoresPosNeg$overallScoreExpr, weights.class0 =  allScoresPosNeg$label)$auc.integral
    prAUCAll[rundidx, 'scoreRank'] = PRROC::pr.curve(scores.class0 = allScoresPosNeg$overallScoreRank, weights.class0 =  allScoresPosNeg$label)$auc.integral
    prAUCAll[rundidx, 'scoreRankPvalue'] = PRROC::pr.curve(scores.class0 = allScoresPosNeg$overallScoreRankPvalue, weights.class0 =  allScoresPosNeg$label)$auc.integral
    
    partialAUCAll[rundidx, 'overallScore'] = pROC::roc(allScoresPosNeg$label, allScoresPosNeg$overallScore, partial.auc = c(1, 0.95), partial.auc.correct = T)$auc[1]
    partialAUCAll[rundidx, 'scorePvalue'] = pROC::roc(allScoresPosNeg$label,allScoresPosNeg$overallScorePvalue,partial.auc = c(1, 0.95), partial.auc.correct = T)$auc[1]
    partialAUCAll[rundidx, 'scoreExpr'] = pROC::roc(allScoresPosNeg$label, allScoresPosNeg$overallScoreExpr, partial.auc = c(1, 0.95), partial.auc.correct = T)$auc[1]
    partialAUCAll[rundidx, 'scoreRank'] = pROC::roc(allScoresPosNeg$label,allScoresPosNeg$overallScoreRank,partial.auc = c(1, 0.95), partial.auc.correct = T)$auc[1]
    partialAUCAll[rundidx, 'scoreRankPvalue'] = pROC::roc(allScoresPosNeg$label,allScoresPosNeg$overallScoreRankPvalue,partial.auc = c(1, 0.95), partial.auc.correct = T)$auc[1]
    
    AUCAll[rundidx, 'overallScore'] = pROC::roc(allScoresPosNeg$label, allScoresPosNeg$overallScore)$auc[1]
    AUCAll[rundidx, 'scorePvalue'] = pROC::roc(allScoresPosNeg$label, allScoresPosNeg$overallScorePvalue)$auc[1]
    AUCAll[rundidx, 'scoreExpr'] = pROC::roc(allScoresPosNeg$label, allScoresPosNeg$overallScoreExpr)$auc[1]
    AUCAll[rundidx, 'scoreRank'] = pROC::roc(allScoresPosNeg$label, allScoresPosNeg$overallScoreRank)$auc[1]
    AUCAll[rundidx, 'scoreRankPvalue'] = pROC::roc(allScoresPosNeg$label, allScoresPosNeg$overallScoreRankPvalue)$auc[1]
    
  }
  
  testCiliopathies[[disease]] = list(prAUCAll, partialAUCAll, AUCAll)
  
  allScoresTest[[disease]] = allScores
  allScoresScaledTest[[disease]] = allScores_scaled
}

library(ggpubr)
#my_comparisons <- list( c('overallScore', 'scoreExpr'), c('overallScore', 'scorePvalue'), c('overallScore', 'scoreRank'), c('overallScore', 'scoreRankPvalue') )
my_comparisons <- list(c('scoreRankPvalue', 'scorePvalue'), c('scoreRankPvalue', 'scoreRank'), c('scoreExpr', 'scoreRankPvalue') , c('overallScore', 'scoreRankPvalue'))

boxplotPRAUC = lapply(testCiliopathies[1:6], function(x) x[[1]] %>% 
                        gather(key = 'scoretype', value = 'scorevalue') %>% 
                        mutate(scoretype = fct_relevel(scoretype, 'scoreExpr', 'scoreRank', 'scorePvalue', 'scoreRankPvalue', 'overallScore')) %>%
                        ggplot(., aes(x=scoretype, y=scorevalue)) + geom_boxplot() + theme_classic() + 
                        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
                        scale_y_continuous(limits = c(0,1.5)) +
                        stat_compare_means(comparisons = my_comparisons))
cowplot::plot_grid(plotlist = boxplotPRAUC, ncol = 3)

boxplotpartialAUC = lapply(testCiliopathies[1:6], function(x) x[[2]] %>% 
                             gather(key = 'scoretype', value = 'scorevalue') %>% 
                             mutate(scoretype = fct_relevel(scoretype, 'scoreExpr', 'scoreRank', 'scorePvalue', 'scoreRankPvalue', 'overallScore')) %>%
                             ggplot(., aes(x=scoretype, y=scorevalue)) + geom_boxplot() + theme_classic() + 
                             theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
                             scale_y_continuous(limits = c(0.5,1.25)) + 
                             stat_compare_means(comparisons = my_comparisons))
cowplot::plot_grid(plotlist = boxplotpartialAUC, ncol = 4)

boxplotAUC = lapply(testCiliopathies[1:6], function(x) x[[3]] %>% 
                      gather(key = 'scoretype', value = 'scorevalue') %>% 
                      mutate(scoretype = fct_relevel(scoretype, 'scoreExpr', 'scoreRank', 'scorePvalue', 'scoreRankPvalue', 'overallScore')) %>%
                      ggplot(., aes(x=scoretype, y=scorevalue)) + geom_boxplot() + theme_classic() + 
                      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
                      scale_y_continuous(limits = c(0.5,1.25)) + 
                      stat_compare_means(comparisons = my_comparisons))
cowplot::plot_grid(plotlist = boxplotAUC, ncol = 4)

modelTest = list('accuracies' = testCiliopathies, 'allscores' = allScoresTest, 'allscoresscaled' = allScoresScaledTest)
#write_rds(modelTest, '20240723_accuraciesTestCiliopathies_train5genes_min10Genes.rds')

#plot for all ciliopathies in one plot
bind_rows(lapply(testCiliopathies[1:6], function(x) x[[3]])) %>%
  gather(key = 'scoretype', value = 'scorevalue') %>% 
  mutate(scoretype = fct_relevel(scoretype, 'scoreExpr', 'scoreRank', 'scorePvalue', 'scoreRankPvalue', 'overallScore')) %>%
  ggplot(., aes(x=scoretype, y=scorevalue)) + geom_boxplot() + theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  stat_compare_means(comparisons = my_comparisons)

