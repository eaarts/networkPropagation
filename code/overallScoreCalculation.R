# calculate overall scores ciliopathies

setwd('W:/GROUP/Users/Ellen/NetworkPropagation/')

library(tidyverse)
library(org.Hs.eg.db)

'%notin%' = Negate('%in%')

pageRankPvalue = read.csv('20240105_pagerankDFPvalue_updatedJBTS.csv', row.names = 1)

rankMP = read.csv('20240527_rankMP.csv')

expressionScore = read.csv('20240329_HPAExpressionLocalization.csv', row.names = 1)

#allGenes = Reduce(intersect, list(rownames(pageRankPvalue),expressionScore$gene,rankMP$gene))
allGenes = Reduce(intersect, list(rownames(pageRankPvalue),rankMP$gene))

allScores = data.frame()
for (disease in colnames(pageRankPvalue)){
  rankCiliopathy = rankMP[rankMP$ciliopathy == disease, ]
  allSCoresCiliopathy = data.frame(pvalue = pageRankPvalue[match(allGenes, rownames(pageRankPvalue)), disease], 
                                   rank = rankCiliopathy[match(allGenes, rankCiliopathy$gene), "rank"],
                                   #expressionScore = expressionScore[match(allGenes, expressionScore$gene),"expression"],
                                   gene = allGenes,
                                   ciliopathy = disease)
  allSCoresCiliopathy[,1:2] = apply(allSCoresCiliopathy[,1:2], 2, function(x) (x-min(x))/(max(x)-min(x)))
  allScores = rbind(allScores, allSCoresCiliopathy)
}

#allScores$overallScore = 0.4741494 * allScores$pvalue - 0.3603909 * allScores$rank + 0.1654598 * allScores$expressionScore #weights based on PR AUC value of separate scores 

model = read_rds('20240529_modelPvalueRank.rds')

allScores$overallScore = predict(model, newdata = allScores)

#allScores$localization = expressionScore$localization[match(allScores$gene, expressionScore$gene)]

# ciliary gene? ----

ciliaryGenes = readxl::read_xlsx('ciliaGenes.xlsx')

allScores$ciliary = allScores$gene %in% ciliaryGenes$Ensembl.Gene.ID

# some uncertain evidence? ----

uncertainGenes = read.csv('20240105_uncertainGenesCiliopathies.csv')

allScores$uncertainEvidence = allScores$gene %in% uncertainGenes$targetId

# ciliopathy gene? ----

variantsCiliopathy = read.csv('20231206_variantsCiliopathy_JBTSnoOverlap.csv')

allScores$ciliopathyGene = allScores$gene %in% variantsCiliopathy$targetId

# ciliopathy gene, specific ----

for (i in unique(allScores$ciliopathy)){
  allScores[allScores$ciliopathy == i, 'ciliopathyGeneSpecific'] = allScores[allScores$ciliopathy == i, 'gene'] %in% variantsCiliopathy[variantsCiliopathy$diseaseId == i, 'targetId']
} 

# disease names ----

traitAnnotation = read.csv('20231206_traitOverview_JBTSnoOverlap.csv')

allScores$diseaseName = traitAnnotation$trait_label[match(allScores$ciliopathy, traitAnnotation$Var1)]

# gene names ----

allScores$geneName = AnnotationDbi::mapIds(org.Hs.eg.db, keys = allScores$gene, column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first')

# rank scores per ciliopathy ----

for (i in unique(allScores$ciliopathy)){
  allScores[allScores$ciliopathy == i, 'overallRank'] = 19307 - rank(allScores[allScores$ciliopathy == i, "overallScore"], ties.method = 'max')
}

# subcellular localization data ----

subcellular = read_tsv('Datasets/humanProteinAtlas/subcellular_location.tsv/subcellular_location.tsv')

allScores$localization = subcellular$`Main location`[match(allScores$gene, subcellular$Gene)]

# expression score ----

allScores$expression = expressionScore$expression[match(allScores$gene, expressionScore$gene)]


  