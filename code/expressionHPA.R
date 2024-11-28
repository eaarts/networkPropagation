# data from human protein atlas
setwd('W:/GROUP/Users/Ellen/NetworkPropagation/')

library(ComplexHeatmap)
library(ggpubr)
library(tidyverse)

candidateGenes = read.csv('20240108_candidateGenes_pvalue0995_updatedJBTS.csv')
variantsCiliopathy = read.csv('20231206_variantsCiliopathy_JBTSnoOverlap.csv')
traitAnnotation = read.csv('20231206_traitOverview_JBTSnoOverlap.csv')

'%notin%' = Negate('%in%')


# load consensus tissue data ----
# 
# rnaTissue = read_tsv('Datasets/humanProteinAtlas/rna_tissue_consensus.tsv/rna_tissue_consensus.tsv')
# 
# rnaTissue_ciliopathies = rnaTissue[rnaTissue$Gene %in% variantsCiliopathy$targetId,]
# 
# rnaTissue_ciliopathiesWide = spread(rnaTissue_ciliopathies[,c("Gene", "Tissue","nTPM")], key = 'Tissue', value = 'nTPM')
# rnaTissue_ciliopathiesWide = column_to_rownames(rnaTissue_ciliopathiesWide, var = 'Gene')
# rnaTissue_ciliopathiesWide = rnaTissue_ciliopathiesWide[rowSums(rnaTissue_ciliopathiesWide, na.rm = T) != 0,]
# rnaTissue_ciliopathiesWide = rnaTissue_ciliopathiesWide[colSums(apply(rnaTissue_ciliopathiesWide,1,is.na)) == 0,]
# 
# rnaTissue_ciliopathiesWide_selectedGenes = rnaTissue_ciliopathiesWide[apply(rnaTissue_ciliopathiesWide, 1, max) > 10,]
# 
# ciliopathies = traitAnnotation$Var1[traitAnnotation$ciliopathy == T]
# expressionCiliopathies = as.data.frame(matrix(nrow = length(ciliopathies), ncol = ncol(rnaTissue_ciliopathiesWide)))
# for (i in 1:length(ciliopathies)){
#   if(sum(unique(variantsCiliopathy$targetId[variantsCiliopathy$diseaseId == ciliopathies[i]]) %in% rownames(rnaTissue_ciliopathiesWide_selectedGenes)) > 1){
#   expressionCiliopathies[i,] = apply(rnaTissue_ciliopathiesWide[rownames(rnaTissue_ciliopathiesWide) %in% variantsCiliopathy$targetId[variantsCiliopathy$diseaseId == ciliopathies[i]],], 2, mean)
#   } else {
#     expressionCiliopathies[i,] = NA
#   }
# }
# 
# rownames(expressionCiliopathies) = traitAnnotation$trait_label[match(ciliopathies, traitAnnotation$Var1)]
# colnames(expressionCiliopathies) = colnames(rnaTissue_ciliopathiesWide)
# 
# expressionCiliopathies = expressionCiliopathies[!is.na(expressionCiliopathies[,1]),]
# 
# expressionCiliopathies_selected = expressionCiliopathies[,apply(expressionCiliopathies, 2, max) > 1]
# 
# Heatmap(scale(t(expressionCiliopathies)))
# 
# rownames(rnaTissue_ciliopathiesWide) = variantsCiliopathy$geneName[match(rownames(rnaTissue_ciliopathiesWide), variantsCiliopathy$targetId)]
# 
# Heatmap(t(scale(t(rnaTissue_ciliopathiesWide))))

# load single cell consensus data ----

rnaSingleCell = read_tsv('Datasets/humanProteinAtlas/rna_single_cell_type.tsv/rna_single_cell_type.tsv')

rnaSingleCell_Wide = spread(rnaSingleCell[,c("Gene", "Cell type","nTPM")], key = 'Cell type', value = 'nTPM')
rnaSingleCell_Wide = column_to_rownames(rnaSingleCell_Wide, var = 'Gene')
rnaSingleCell_Wide = rnaSingleCell_Wide[rowSums(rnaSingleCell_Wide) != 0,]

rnaSingleCell_WideScaled = t(scale(t(rnaSingleCell_Wide)))
rnaSingleCell_WideScaled = as.data.frame(rnaSingleCell_WideScaled)



# use supervised learning model to define weights

ciliopathyGenes = unique(variantsCiliopathy$targetId)  
nonCiliopathyGenes = rownames(rnaSingleCell_WideScaled)[rownames(rnaSingleCell_WideScaled) %notin% variantsCiliopathy$targetId]  

rnaSingleCell_WideScaled$ciliopathyGene = rownames(rnaSingleCell_WideScaled) %in% ciliopathyGenes

set.seed(12345)
testGenes = c(sample(ciliopathyGenes, 0.1*length(ciliopathyGenes)), sample(nonCiliopathyGenes, 0.1*length(nonCiliopathyGenes)))
trainGenes = rownames(rnaSingleCell_WideScaled)[rownames(rnaSingleCell_WideScaled) %notin% testGenes]

trainData = rnaSingleCell_WideScaled[rownames(rnaSingleCell_WideScaled) %in% trainGenes,]
testData = rnaSingleCell_WideScaled[rownames(rnaSingleCell_WideScaled) %in% testGenes,]

anovaRes = aov(ciliopathyGene ~ . , trainData)
tableAnova = summary(anovaRes)[[1]]
cellsForModel = rownames(tableAnova)[tableAnova$`Pr(>F)` < 0.01]
cellsForModel = gsub('`', '', str_trim(cellsForModel))

trainData_selected = trainData[,colnames(trainData) %in% cellsForModel]
trainData_selected$ciliopathyGene = rownames(trainData_selected) %in% ciliopathyGenes

model = glm(ciliopathyGene ~ ., family = binomial(link = 'logit'), data = trainData_selected)

testPrediction = predict(model, newdata = testData, type = 'response')

plot(0.79)

car::Anova(model, type = 2)

# calculate final expression score ----

anovaRes = aov(ciliopathyGene ~ . , rnaSingleCell_WideScaled)
tableAnova = summary(anovaRes)[[1]]
cellsForModel = rownames(tableAnova)[tableAnova$`Pr(>F)` < 0.01]
cellsForModel = gsub('`', '', str_trim(cellsForModel))

rnaSingleCell_WideScaled_selected = rnaSingleCell_WideScaled[, colnames(rnaSingleCell_WideScaled) %in% cellsForModel]
rnaSingleCell_WideScaled_selected$ciliopathyGene = rownames(rnaSingleCell_WideScaled_selected) %in% ciliopathyGenes

model = glm(ciliopathyGene ~ ., family = binomial(link = 'logit'), data = rnaSingleCell_WideScaled_selected)

expressionScore = predict(model, type = 'response')

HPAdata = data.frame(expression = expressionScore)

# expression per ciliopathy ----

ciliopathies = traitAnnotation$Var1[traitAnnotation$ciliopathy == T]

expressionCiliopathies = as.data.frame(matrix(nrow = length(ciliopathies), ncol = ncol(rnaSingleCell_Wide)))
for (i in 1:length(ciliopathies)){
  expressionCiliopathies[i,] = apply(rnaSingleCell_Wide[rownames(rnaSingleCell_Wide) %in% variantsCiliopathy$targetId[variantsCiliopathy$diseaseId == ciliopathies[i]],], 2, mean)

}

rownames(expressionCiliopathies) = traitAnnotation$trait_label[match(ciliopathies, traitAnnotation$Var1)]
colnames(expressionCiliopathies) = colnames(rnaSingleCell_Wide)

expressionCiliopathies = scale(t(expressionCiliopathies))

#expressionCiliopathies_selected = expressionCiliopathies[apply(expressionCiliopathies, 1, max) > 1,]
expressionCiliopathies_selected = expressionCiliopathies[rownames(expressionCiliopathies) %in% cellsForModel,]

Heatmap(expressionCiliopathies_selected)

rnaSingleCell_Wide_ciliopathies = rnaSingleCell_Wide[rownames(rnaSingleCell_Wide) %in% variantsCiliopathy$targetId, ]
rownames(rnaSingleCell_Wide_ciliopathies) = variantsCiliopathy$geneName[match(rownames(rnaSingleCell_Wide_ciliopathies), variantsCiliopathy$targetId)]

rnaSingleCell_Wide_ciliopathies = t(scale(t(rnaSingleCell_Wide_ciliopathies)))

rnaSingleCell_Wide_ciliopathies = rnaSingleCell_Wide_ciliopathies[,colnames(rnaSingleCell_Wide_ciliopathies) %in% cellsForModel]

Heatmap(rnaSingleCell_Wide_ciliopathies)

# load subcellular localization data ----

subcellular = read_tsv('Datasets/humanProteinAtlas/subcellular_location.tsv/subcellular_location.tsv')

HPAdata$localization = subcellular$`Main location`[match(HPAdata$gene, subcellular$Gene)]

write.csv(HPAdata, '20240329_HPAExpressionLocalization.csv')

