# calculate overall scores ciliopathies

library(tidyverse)
library(org.Hs.eg.db)

'%notin%' = Negate('%in%')

pageRankPvalue = read.csv('data/pvaluesPropagation.csv', row.names = 1)
rankMP = read.csv('data/rankMP.csv')
expressionScore = read.csv('data/HPAExpressionLocalization.csv', row.names = 1)

allGenes = Reduce(intersect, list(rownames(pageRankPvalue),rankMP$gene))

allScores = data.frame()
for (disease in colnames(pageRankPvalue)){
  rankCiliopathy = rankMP[rankMP$ciliopathy == disease, ]
  allSCoresCiliopathy = data.frame(pvalue = pageRankPvalue[match(allGenes, rownames(pageRankPvalue)), disease], 
                                   rank = rankCiliopathy[match(allGenes, rankCiliopathy$gene), "rank"],
                                   gene = allGenes,
                                   ciliopathy = disease)
  allSCoresCiliopathy[,1:2] = apply(allSCoresCiliopathy[,1:2], 2, function(x) (x-min(x))/(max(x)-min(x)))
  allScores = rbind(allScores, allSCoresCiliopathy)
}


model = read_rds('data/candidateModel.rds')

allScores$overallScore = predict(model, newdata = allScores)

# ciliary gene? ----

ciliaryGenes = readxl::read_xlsx('data/ciliaGenes.xlsx')

allScores$ciliary = allScores$gene %in% ciliaryGenes$Ensembl.Gene.ID

# some uncertain evidence? ----

uncertainGenes = read.csv('data/uncertainGenesCiliopathies.csv')

allScores$uncertainEvidence = allScores$gene %in% uncertainGenes$targetId

# ciliopathy gene? ----

variantsCiliopathy = read.csv('data/variantsCiliopathies.csv')

allScores$ciliopathyGene = allScores$gene %in% variantsCiliopathy$targetId

# ciliopathy gene, specific ----

for (i in unique(allScores$ciliopathy)){
  allScores[allScores$ciliopathy == i, 'ciliopathyGeneSpecific'] = allScores[allScores$ciliopathy == i, 'gene'] %in% variantsCiliopathy[variantsCiliopathy$diseaseId == i, 'targetId']
} 

# disease names ----

traitAnnotation = read.csv('data/traitOverview.csv')

allScores$diseaseName = traitAnnotation$trait_label[match(allScores$ciliopathy, traitAnnotation$Var1)]

# gene names ----

allScores$geneName = AnnotationDbi::mapIds(org.Hs.eg.db, keys = allScores$gene, column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first')

# rank scores per ciliopathy ----

for (i in unique(allScores$ciliopathy)){
  allScores[allScores$ciliopathy == i, 'overallRank'] = 19307 - rank(allScores[allScores$ciliopathy == i, "overallScore"], ties.method = 'max')
}

# subcellular localization data ----

subcellular = read_tsv('data/subcellular_location.tsv')

allScores$localization = subcellular$`Main location`[match(allScores$gene, subcellular$Gene)]

# expression score ----

allScores$expression = expressionScore$expression[match(allScores$gene, expressionScore$gene)]

write.csv(allScores, 'data/finalScores.csv')