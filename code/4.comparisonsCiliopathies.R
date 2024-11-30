# similarity analysis of ciliopathies using propagation scores

# Set working directory ----

setwd('W:/GROUP/Users/Ellen/NetworkPropagation/')

# Libraries ---- 

library(tidyverse)
library(ComplexHeatmap)

# Load necessary files ----

pageRankScores = readRDS('20231206_pageRankScores_JBTSnoOverlap.rds') #calculated network propagation scores for ciliopathies
traitAnnotation = read.csv('20231206_traitOverview_JBTSnoOverlap.csv') #labels for ciliopathies
variantsCiliopathy = read.csv('20231206_variantsCiliopathy_JBTSnoOverlap.csv') #curated list of ciliopathy genes

# Categorize ciliopathies with their phenotypes ----

phenotype = data.frame('ciliopathy' = unique(traitAnnotation[traitAnnotation$ciliopathy == T, "trait_label"]), 'renal' = F, 'retina' = F, 'CNS' = F, 'polydactyly' = F, 'skeletal' = F, 'hearing' = F)
phenotype$renal[grep('kidney|meckel|renal|nephro|senior|orofacio|bardet', phenotype$ciliopathy, ignore.case = T)] = T
phenotype$retina[grep('ocular|bardet|senior|usher|leber|retinitis|cone', phenotype$ciliopathy, ignore.case = T)] = T
phenotype$CNS[grep('joubert|meckel|acrocallosal|hydrolethalus|orofacio|bardet', phenotype$ciliopathy, ignore.case = T)] = T
phenotype$polydactyly[grep('meckel|orofacio|acrocallosal|bardet|hydrolethalus|ellis|short rib', phenotype$ciliopathy, ignore.case = T)] = T
phenotype$skeletal[grep('Jeune|short rib|ellis|cranioe', phenotype$ciliopathy, ignore.case = T)] = T
phenotype$hearing[grep('usher', phenotype$ciliopathy, ignore.case = T)] = T

# Calculate distance between ciliopathies based on propagation scores ----

pageRankDFCilia = pageRankScores[rownames(pageRankScores) %in% traitAnnotation[traitAnnotation$ciliopathy == T, "Var1"],]
pageRankDFCilia_scaled = t(scale(t(pageRankDFCilia)))
rownames(pageRankDFCilia_scaled) = traitAnnotation$trait_label[match(rownames(pageRankDFCilia_scaled), traitAnnotation$Var1)]
propagationDist = as.matrix(dist(pageRankDFCilia_scaled, method = 'euclidean')) #best clustering with euclidean distance and scaled propagation scores

# Calculate distance between ciliopathies based on overlapping genes ----

variantsTable = variantsCiliopathy[,c("targetId", "diseaseId")]
variantsTable = unique(variantsTable[,1:2])
variantsTable$mutated = 1
variantsTable = spread(variantsTable,key = diseaseId, value = mutated)
variantsTable = column_to_rownames(variantsTable, var = 'targetId' )
variantsTable[is.na(variantsTable)] = 0
variantsTable = t(variantsTable)
rownames(variantsTable) = traitAnnotation$trait_label[match(rownames(variantsTable), traitAnnotation$Var1)]
seedDist = as.matrix(dist(variantsTable, method = 'binary'))
  
# Plot clustering of ciliopathies ----

# seed gene based
hclustTraits = hclust(as.dist(seedDist))
pdf('Figures/Fig1/dendrogramSeedGenes.pdf', width = 8, height = 7)
plot(hclustTraits)
dev.off()

phenotype_df = phenotype
rownames(phenotype_df) = NULL
phenotype_df = column_to_rownames(phenotype_df, var = 'ciliopathy')
phenotype_df[phenotype_df == T] = 1

phenotype_df = phenotype_df[match(hclustTraits$labels, rownames(phenotype_df)),]
pdf('Figures/Fig1/dendrogramSeedGenes_phenotypeAnno.pdf', width = 7, height = 4.5)
Heatmap(t(phenotype_df), col = c('white', 'black'), rect_gp = gpar(col = "white", lwd = 2), cluster_columns = hclustTraits)
dev.off()

# propagation score based
hclustTraits = hclust(as.dist(propagationDist))
pdf('Figures/Fig1/dendrogramPropagationScores.pdf', width = 8, height = 7)
plot(hclustTraits)
dev.off()

phenotype_df = phenotype
rownames(phenotype_df) = NULL
phenotype_df = column_to_rownames(phenotype_df, var = 'ciliopathy')
phenotype_df[phenotype_df == T] = 1

phenotype_df = phenotype_df[match(hclustTraits$labels, rownames(phenotype_df)),]
pdf('Figures/Fig1/dendrogramPropagationScores_phenotypeAnno.pdf', width = 7, height = 4.5)
Heatmap(t(phenotype_df), col = c('white', 'black'), rect_gp = gpar(col = "white", lwd = 2), cluster_columns = hclustTraits)
dev.off()

# Calculate AUC values of ciliopathies with similar phenotypes ----

# create dataframe with ciliopathy pairs, indicating if they have a similar phenotype
phenotypeDist = combn(phenotype$ciliopathy, 2) #create ciliopathy pairs
phenotypeDist = as.data.frame(t(phenotypeDist))
phenotypeDist$renal = F
phenotypeDist$CNS = F
phenotypeDist$hearing = F
phenotypeDist$polydactyly = F
phenotypeDist$retina = F
phenotypeDist$skeletal = F

for (i in 1:nrow(phenotypeDist)) { #annotate if pairs have a similar phenotype
  if (phenotype[phenotype$ciliopathy == phenotypeDist[i, 1], "renal"] == T &
      phenotype[phenotype$ciliopathy == phenotypeDist[i, 2], "renal"] == T) {
    phenotypeDist$renal[i] = T
  } 
  if (phenotype[phenotype$ciliopathy == phenotypeDist[i, 1], "CNS"] == T &
      phenotype[phenotype$ciliopathy == phenotypeDist[i, 2], "CNS"] == T) {
    phenotypeDist$CNS[i] = T
  } 
  if (phenotype[phenotype$ciliopathy == phenotypeDist[i, 1], "hearing"] == T &
      phenotype[phenotype$ciliopathy == phenotypeDist[i, 2], "hearing"] == T) {
    phenotypeDist$hearing[i] = T
  } 
  if (phenotype[phenotype$ciliopathy == phenotypeDist[i, 1], "polydactyly"] == T &
      phenotype[phenotype$ciliopathy == phenotypeDist[i, 2], "polydactyly"] == T) {
    phenotypeDist$polydactyly[i] = T
  } 
  if (phenotype[phenotype$ciliopathy == phenotypeDist[i, 1], "retina"] == T &
      phenotype[phenotype$ciliopathy == phenotypeDist[i, 2], "retina"] == T) {
    phenotypeDist$retina[i] = T
  } 
  if (phenotype[phenotype$ciliopathy == phenotypeDist[i, 1], "skeletal"] == T &
      phenotype[phenotype$ciliopathy == phenotypeDist[i, 2], "skeletal"] == T) {
    phenotypeDist$skeletal[i] = T
  }
}

# calculate distance between pairs based on seed gene overlap and propagation scores
for (i in 1:nrow(phenotypeDist)){phenotypeDist$propDist[i] = propagationDist[phenotypeDist[i,1], phenotypeDist[i,2]]}
for (i in 1:nrow(phenotypeDist)){phenotypeDist$seedDist[i] = seedDist[phenotypeDist[i,1], phenotypeDist[i,2]]}

# split dataframe up per phenotype so distance between ciliopathies with the same phenotype can be compared to the distance between these ciliopathies and ciliopathies with another phenotype
phenotypeDistRenal = phenotypeDist[(phenotypeDist$V1 %in% phenotype[phenotype$renal == T, "ciliopathy"] | phenotypeDist$V2 %in% phenotype[phenotype$renal == T, "ciliopathy"]), ]
phenotypeDistCNS = phenotypeDist[(phenotypeDist$V1 %in% phenotype[phenotype$CNS == T, "ciliopathy"] | phenotypeDist$V2 %in% phenotype[phenotype$CNS == T, "ciliopathy"]), ]
phenotypeDisthearing = phenotypeDist[(phenotypeDist$V1 %in% phenotype[phenotype$hearing == T, "ciliopathy"] | phenotypeDist$V2 %in% phenotype[phenotype$hearing == T, "ciliopathy"]), ]
phenotypeDistretina = phenotypeDist[(phenotypeDist$V1 %in% phenotype[phenotype$retina == T, "ciliopathy"] | phenotypeDist$V2 %in% phenotype[phenotype$retina == T, "ciliopathy"]), ]
phenotypeDistpolydactyly = phenotypeDist[(phenotypeDist$V1 %in% phenotype[phenotype$polydactyly == T, "ciliopathy"] | phenotypeDist$V2 %in% phenotype[phenotype$polydactyly == T, "ciliopathy"]), ]
phenotypeDistskeletal = phenotypeDist[(phenotypeDist$V1 %in% phenotype[phenotype$skeletal == T, "ciliopathy"] | phenotypeDist$V2 %in% phenotype[phenotype$skeletal == T, "ciliopathy"]), ]

# calculate AUC values per phenotype and based on either the seed gene distance or propagation score distance
rocRes = data.frame('seedDist' = NA, 'propDist' = NA, phenotype = c('renal', 'retina', 'CNS', 'hearing', 'polydactyly', 'skeletal'))
rocRes[rocRes$phenotype == 'CNS', 'seedDist'] = pROC::roc(phenotypeDistCNS$CNS, phenotypeDistCNS$seedDist)$auc
rocRes[rocRes$phenotype == 'CNS', "propDist"] = pROC::roc(phenotypeDistCNS$CNS, phenotypeDistCNS$propDist)$auc
rocRes[rocRes$phenotype == 'renal', "propDist"] = pROC::roc(phenotypeDistRenal$renal, phenotypeDistRenal$propDist)$auc
rocRes[rocRes$phenotype == 'renal', "seedDist"] = pROC::roc(phenotypeDistRenal$renal, phenotypeDistRenal$seedDist)$auc
rocRes[rocRes$phenotype == 'retina', "seedDist"] = pROC::roc(phenotypeDistretina$retina, phenotypeDistretina$seedDist)$auc
rocRes[rocRes$phenotype == 'retina', "propDist"] = pROC::roc(phenotypeDistretina$retina, phenotypeDistretina$propDist)$auc
rocRes[rocRes$phenotype == 'hearing', "propDist"] = pROC::roc(phenotypeDisthearing$hearing, phenotypeDisthearing$propDist)$auc
rocRes[rocRes$phenotype == 'hearing', "seedDist"] = pROC::roc(phenotypeDisthearing$hearing, phenotypeDisthearing$seedDist)$auc
rocRes[rocRes$phenotype == 'polydactyly', "seedDist"] = pROC::roc(phenotypeDistpolydactyly$polydactyly, phenotypeDistpolydactyly$seedDist)$auc
rocRes[rocRes$phenotype == 'polydactyly', "propDist"] = pROC::roc(phenotypeDistpolydactyly$polydactyly, phenotypeDistpolydactyly$propDist)$auc
rocRes[rocRes$phenotype == 'skeletal', "propDist"] = pROC::roc(phenotypeDistskeletal$skeletal, phenotypeDistskeletal$propDist)$auc
rocRes[rocRes$phenotype == 'skeletal', "seedDist"] = pROC::roc(phenotypeDistskeletal$skeletal, phenotypeDistskeletal$seedDist)$auc

# create a long dataframe for plotting
rocRes = gather(rocRes, key = 'distType', value = 'dist', -phenotype)

# plot the AUC values grouped by phenotype and split by distance metrix
pltAUC = ggplot(rocRes, aes(x=phenotype, y=dist, fill = distType)) + 
  geom_bar(stat = 'identity', position = 'dodge', color = 'black') +
  theme_classic(base_size = 15) +
  scale_fill_manual(values=c("#D2A8AC", "#BE1E2D"), name = 'Distance metric', labels = c('Propagation scores', 'Seed genes')) +
  xlab('Phenotype') + ylab('AUROC') + 
  scale_x_discrete(labels=c('CNS disease', 'Hearing loss', 'Polydactyly', 'Renal disease', 'Eye disease', 'Skeletal dysplasia')) + 
  coord_flip()

ggsave('Figures/Fig1/ROCAUCperPhenotype.pdf', plot = pltAUC, width = 17, height = 10, units = 'cm')

# Calculate AUC values for all phenotypes combined ----

for (i in 1:nrow(phenotypeDist)){
  phenotypeDist$comb[i] = sum(phenotypeDist[i,c("renal","CNS", "hearing","polydactyly","skeletal", "retina")]) > 0
} 

overallROCseed = pROC::roc(phenotypeDist$comb, phenotypeDist$seedDist, direction = '>')
overallROCprop = pROC::roc(phenotypeDist$comb, phenotypeDist$propDist, direction = '>')

pdf('Figures/Fig1/ROCAUCOverall.pdf', width = 5.5, height = 5)
plot(overallROCseed, col = 'red')
plot(overallROCprop, add = TRUE, col = 'blue')
legend(0.75, 0.15, legend=c(paste0("Seed genes", ', ROCAUC = ', round(overallROCseed$auc, 2)), paste0("Propagation scores", ', ROCAUC = ', round(overallROCprop$auc, 2))),
       col=c("red", "blue"),lty = 1, cex=0.8)
dev.off()


