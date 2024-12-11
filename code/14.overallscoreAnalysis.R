# analyze overall score

traitAnnotation = read.csv('data/traitOverview.csv')
PPIClusters = read.csv('data/PPIFullNetworkClusters.csv')
variantsCiliopathy = read.csv('data/variantsCiliopathies.csv') 

'%notin%' = Negate('%in%')

allScores = read.csv('data/finalScores.csv', row.names = 1)

allScoresNoKnown = allScores[allScores$ciliopathyGene == F,]

# correlate gene ranking ----
lowGenes = table(allScoresNoKnown[allScoresNoKnown$overallRank > 1000, "gene"])
lowGenes = names(lowGenes[lowGenes==21])

allScoresNoKnownhighGenes = allScoresNoKnown[allScoresNoKnown$gene %notin% lowGenes,]

corRes = matrix(nrow = length(unique(allScores$ciliopathy)), ncol = length(unique(allScores$ciliopathy)))
colnames(corRes) = unique(allScores$ciliopathy)
rownames(corRes) = unique(allScores$ciliopathy)

for (i in unique(allScores$ciliopathy)){
  for (j in unique(allScores$ciliopathy)){
    corRes[i,j] = cor(allScoresNoKnownhighGenes[allScoresNoKnownhighGenes$ciliopathy == i, "overallScore"], 
                      allScoresNoKnownhighGenes[allScoresNoKnownhighGenes$ciliopathy == j, "overallScore"], method = 'pearson')
  }}

colnames(corRes) = traitAnnotation$trait_label[match(colnames(corRes), traitAnnotation$Var1)]
rownames(corRes) = traitAnnotation$trait_label[match(rownames(corRes), traitAnnotation$Var1)]

library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(0.5, 1), c("white", "red"))
Heatmap(corRes, col = col_fun)

# enrichment of several things ----
gsearesAll = list()
for (ciliopathy in unique(allScoresNoKnown$ciliopathy)) {
  scoresCiliopathy = allScoresNoKnown[allScoresNoKnown$ciliopathy == ciliopathy, ]
  
  ciliaryGenes = unique(scoresCiliopathy[scoresCiliopathy$ciliary == T, "gene"])
  uncertainGenes = unique(scoresCiliopathy[scoresCiliopathy$uncertainEvidence == T, "gene"])
  localization = unique(scoresCiliopathy[, c("localization", "gene")])
  colnames(localization) = c('set', 'gene')
  
  doubleLoc = grep(';', localization$set)
  for (i in doubleLoc) {
    newLoc = data.frame('set' = str_split(localization$set[doubleLoc[i]], ';')[[1]],
                        'gene' = localization$gene[doubleLoc[i]])
    localization = rbind(localization, newLoc)
  }
  localization = localization[-doubleLoc, ]
  localization = unique(localization)
  localization = localization[!is.na(localization$set), ]
  
  geneSets = data.frame('set' = c(
    rep('ciliary', each = length(ciliaryGenes)),
    rep('uncertain', each = length(uncertainGenes))
  ),
  'gene' = c(ciliaryGenes, uncertainGenes))
  geneSets = rbind(geneSets, localization)
  
  scoreVector = scale(scoresCiliopathy$overallScore)
  names(scoreVector) = scoresCiliopathy$gene
  scoreVector = sort(scoreVector, decreasing = T)
  
  gseares = clusterProfiler::GSEA(geneList = scoreVector, TERM2GENE = geneSets, pvalueCutoff = 0.05, maxGSSize = 2000)
  
  gsearesAll[[ciliopathy]] = gseares@result
}


table(unlist(lapply(gsearesAll, function(x) x[[1]]))) %>% as.data.frame() %>% ggplot(., aes(x=Freq, y=reorder(Var1, Freq))) + geom_bar(stat = 'identity') + theme_classic()

# module enrichment - fraction ----

#modules
allScoresNoKnown$module = PPIClusters$walktrapCluster[match(allScoresNoKnown$gene, PPIClusters$gene)]

#candidate
candidateGenes = allScoresNoKnown[allScoresNoKnown$overallRank <= 100,]

candidateGenes_singleGenes = candidateGenes[!duplicated(candidateGenes$gene), ]
candidateGenes_singleGenes_withModule = candidateGenes_singleGenes[!is.na(candidateGenes_singleGenes$module),]

allGenesModule = PPIClusters[,c("gene", "walktrapCluster")]
unknownGenesModule = allGenesModule[allGenesModule$gene %notin% variantsCiliopathy$targetId,]

tableAllGenesModule = as.data.frame(table(allGenesModule$walktrapCluster))

tableAllGenesModule$candidateGenes = sapply(tableAllGenesModule$Var1, function(x) sum(candidateGenes_singleGenes_withModule$module %in% x))

tableAllGenesModule$unknownGenes = sapply(tableAllGenesModule$Var1, function(x) sum(unknownGenesModule$walktrapCluster %in% x))

for( i in 1:nrow(tableAllGenesModule)){
  tableAllGenesModule$fractionCandidate[i] = tableAllGenesModule$candidateGenes[i]/tableAllGenesModule$unknownGenes[i]
}

for (i in 1:nrow(tableAllGenesModule)){
  tableAllGenesModule$NciliopathyVariant[i] = sum(PPIClusters[PPIClusters$walktrapCluster == as.character(tableAllGenesModule$Var1[i]), "gene"] %in% variantsCiliopathy$targetId)
}

for( i in 1:nrow(tableAllGenesModule)){
  tableAllGenesModule$fractionKnown[i] = tableAllGenesModule$NciliopathyVariant[i]/tableAllGenesModule$Freq[i]
}

allGenesModule$ciliary = allScores$ciliary[match(allGenesModule$gene, allScores$gene)]
ciliaryGenesModule = allGenesModule[allGenesModule$ciliary == T,]
tableAllGenesModule$ciliaryGenes = sapply(tableAllGenesModule$Var1, function(x) sum(ciliaryGenesModule$walktrapCluster %in% x))

for( i in 1:nrow(tableAllGenesModule)){
  tableAllGenesModule$fractionCiliary[i] = tableAllGenesModule$ciliaryGenes[i]/tableAllGenesModule$Freq[i]
}

tableAllGenesModule = tableAllGenesModule[rowSums(tableAllGenesModule[,c("fractionCandidate", "fractionKnown", "fractionCiliary")], na.rm = T) != 0 ,]


tableAllGenesModuleHighFrac = tableAllGenesModule[apply(tableAllGenesModule, 1, FUN = function(x){max(x[c("fractionCandidate", "fractionKnown", "fractionCiliary")], na.rm = T)}) >= 0.6,]

for (i in 1:nrow(tableAllGenesModuleHighFrac)){
  tableAllGenesModuleHighFrac$GO[i] = clusterProfiler::enrichGO(PPIClusters[PPIClusters$walktrapCluster == tableAllGenesModuleHighFrac$Var1[i], "gene"], org.Hs.eg.db, keyType = 'ENSEMBL', minGSSize = 5)@result$Description[1]
}

tableAllGenesModuleHighFrac$clusterName = paste0(tableAllGenesModuleHighFrac$Var1, '_', tableAllGenesModuleHighFrac$GO)

tableAllGenesModuleHighFrac_long = gather(tableAllGenesModuleHighFrac[,c("clusterName", "fractionCandidate", "fractionKnown", "fractionCiliary")], key = 'category', value = 'fraction', -clusterName)

ggplot(tableAllGenesModuleHighFrac_long, aes(x = category, y=clusterName, size = fraction, color = fraction)) + 
  geom_point() + 
  hrbrthemes::theme_ipsum(base_family = 'sans') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_gradient2(low = "white", high = "red", limit = c(0, max(tableAllGenesModuleHighFrac_long$fraction))) + 
  scale_size(range = c(0, 5))


# ciliary genes, low scores ----

ggplot(allScoresNoKnown[allScoresNoKnown$ciliary == T, ], aes(x=overallRank)) + geom_histogram(fill = 'lightgrey', color = 'black') + theme_classic()

hist(table(allScoresNoKnown[allScoresNoKnown$overallRank <= 100, "geneName"]), breaks = 20)


