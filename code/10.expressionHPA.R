# data from human protein atlas

library(ComplexHeatmap)
library(ggpubr)
library(tidyverse)

variantsCiliopathy = read.csv('data/variantsCiliopathies.csv') 
traitAnnotation = read.csv('data/traitOverview.csv')

'%notin%' = Negate('%in%')

# load single cell consensus data ----

rnaSingleCell = read_tsv('data/rna_single_cell_type.tsv') #HPA dataset

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

pROC::roc(testData$ciliopathyGene,testPrediction)$auc

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

ciliopathies = traitAnnotation$Var1[traitAnnotation$ciliopathy == T & !is.na(traitAnnotation$ciliopathy)]

expressionCiliopathies = as.data.frame(matrix(nrow = length(ciliopathies), ncol = ncol(rnaSingleCell_Wide)))
for (i in 1:length(ciliopathies)){
  expressionCiliopathies[i,] = apply(rnaSingleCell_Wide[rownames(rnaSingleCell_Wide) %in% variantsCiliopathy$targetId[variantsCiliopathy$diseaseId == ciliopathies[i]],], 2, mean)

}

rownames(expressionCiliopathies) = traitAnnotation$trait_label[match(ciliopathies, traitAnnotation$Var1)]
colnames(expressionCiliopathies) = colnames(rnaSingleCell_Wide)

expressionCiliopathies = scale(t(expressionCiliopathies))

expressionCiliopathies_selected = expressionCiliopathies[rownames(expressionCiliopathies) %in% cellsForModel,]

Heatmap(expressionCiliopathies_selected)

# load subcellular localization data ----

subcellular = read_tsv('data/subcellular_location.tsv') #from HPA dataset

HPAdata$localization = subcellular$`Main location`[match(rownames(HPAdata), subcellular$Gene)]

HPAdata = rownames_to_column(HPAdata, var = 'gene')

write.csv(HPAdata, 'HPAExpressionLocalization.csv')
