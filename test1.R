# PageRank for ciliopathies with genes from Open Targets Direct Genetic Associations

# Set working directory ----

setwd('W:/GROUP/Users/Ellen/NetworkPropagation/')

# Libraries ---- 

library(igraph)
library(pROC)
library(doParallel)
library(foreach)
library(org.Hs.eg.db)
library(tidyverse)
library(clusterProfiler)
library(ComplexHeatmap)

source('Code/networkPropagation.R')

'%notin%' = Negate('%in%')

# Interaction network full ----

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

#remove big files to save memory
rm(intAll, intGraph, intString, lowString)
gc()

# load variants from Open Targets ----

#get EFO IDs for ciliopathies
cilioEFO = read.csv('./EFOIdsCiliopathies.csv')
cilioEFO = cilioEFO[-grep('Orphanet', cilioEFO$ID),] #kick out orphanet IDs -> mainly obsolete traits or double traits with MONDO Ids

#variants = read.csv('Datasets/associationByDatatypeDirect/associationAll.csv')
variants = read.csv('Datasets/associationByDatasourceDirect/associationAll.csv')

# variants[variants$diseaseId %in% cilioEFO$ID[grep('Bardet', cilioEFO$description)],"diseaseId"] = 'MONDO_0015229' #all BBS subtypes together
# variants[variants$diseaseId %in% cilioEFO$ID[grep('^orofaciodigital', cilioEFO$description, ignore.case = T)],"diseaseId"] = 'MONDO_0015375' #all orofaciodigital subtypes together
# variants[variants$diseaseId %in% cilioEFO$ID[grep('nephronophthisis', cilioEFO$description)],"diseaseId"] = 'MONDO_0019005' #all nephronophthisis subtypes together
# variants[variants$diseaseId %in% cilioEFO$ID[grep('Usher syndrome type 1', cilioEFO$description)],"diseaseId"] = 'MONDO_0010168' 
# variants[variants$diseaseId %in% cilioEFO$ID[grep('Usher syndrome type 2', cilioEFO$description)],"diseaseId"] = 'MONDO_0016484' 
# variants[variants$diseaseId %in% cilioEFO$ID[grep('short', cilioEFO$description, ignore.case = T)],"diseaseId"] = 'MONDO_0015461' 
# variants[variants$diseaseId %in% cilioEFO$ID[grep('retinitis', cilioEFO$description)],"diseaseId"] = 'MONDO_0019200' 
# variants[variants$diseaseId %in% cilioEFO$ID[grep('^asphy|Jeune syndrome', cilioEFO$description)],"diseaseId"] = 'MONDO_0018770'
# variants[variants$diseaseId %in% cilioEFO$ID[grep('polycystic', cilioEFO$description, ignore.case = T)],"diseaseId"] = 'EFO_0008620'
# variants[variants$diseaseId %in% cilioEFO$ID[grep('Ellis', cilioEFO$description, ignore.case = T)],"diseaseId"] = 'MONDO_0009162'
# variants[variants$diseaseId %in% cilioEFO$ID[grep('hydrolet', cilioEFO$description, ignore.case = T)],"diseaseId"] = 'MONDO_0006037' 

variantCertain = variants[(variants$datasourceId == 'uniprot_literature' & variants$score > 0.5) | 
                            (variants$datasourceId == 'gene2phenotype' & variants$score > 0.5) | 
                            (variants$datasourceId == 'clingen' & variants$score > 0.5) | 
                            (variants$datasourceId == 'uniprot_variants' & variants$score > 0.5) | 
                            (variants$datasourceId == 'genomics_england' & variants$score > 0.5) | 
                            (variants$datasourceId == 'orphanet' & variants$score > 0.5) | 
                            (variants$datasourceId == 'gene_burden' & variants$score > 0.5) | 
                            (variants$datasourceId == 'ot_genetics_portal' & variants$score > 0.5) | 
                            (variants$datasourceId == 'eva' & variants$score > 0.8),]

variantCertainCilia = read.csv('variantsCilia.csv')
diseasesCilia = unique(variantCertainCilia$diseaseId)
variantCertain = variantCertain[variantCertain$diseaseId %notin% diseasesCilia,]
variantCertain = variantCertain[,-1]
variantCertain = rbind(variantCertain, variantCertainCilia[,1:6])

cilioEFO[nrow(cilioEFO) + 1,'ID'] = 'MONDO_custom'
cilioEFO[nrow(cilioEFO),"description"] = 'cilia-associated RP' 

variantUncertain = variants[variants$X %notin% variantCertain$X,]
variantUncertainCilia = variantUncertain[variantUncertain$diseaseId %in% diseasesCilia,]
variantUncertainCilia = variantUncertainCilia[variantUncertainCilia$diseaseId %in% variantCertainCilia$diseaseId,]
variantUncertainCilia$diseaseName = cilioEFO$description[match(variantUncertainCilia$diseaseId, cilioEFO$ID)]

rm(variants)
gc()

#plot number of seed genes per ciliopathy
variantCertainCilia$diseaseName = cilioEFO$description[match(variantCertainCilia$diseaseId, cilioEFO$ID)]
par(mar = c(25,4,2,1))
seedNumber = table(unique(variantCertainCilia[,c("targetId", "diseaseName")])[,2])
seedNumber = seedNumber[seedNumber > 1]
seedNumber = seedNumber[order(seedNumber, decreasing = T)]
barplot(seedNumber, las =2, ylab = 'Number of seed genes')
#plot data sources for ciliopathy genes
pie(table(variantCertainCilia$datasourceId), col = RColorBrewer::brewer.pal(8, 'Set3'))

rm(variants)
gc()

# Add human phenotype ontology 

#HPO = read.table('Datasets/HPO/phenotype.hpoa', sep = '\t', comment.char = '', skip = 4, header = T)
HPOGenes = read.table('Datasets/HPO/phenotype_to_genes.txt', sep = '\t')
colnames(HPOGenes) = c('diseaseId', 'HPOLabel', 'EntrezId', 'GeneSymbol', 'SourceInfo', 'Source', 'DiseaseId')

HPOGenes$HPOLabel = gsub('$', '_HP',HPOGenes$HPOLabel)
HPOGenes$diseaseId = gsub(':', '_', HPOGenes$diseaseId)

HPOGenes$EntrezId = as.character(HPOGenes$EntrezId)
HPOGenes$targetId = mapIds(x = org.Hs.eg.db, keys = HPOGenes$EntrezId, column = 'ENSEMBL', keytype = 'ENTREZID', multiVals = 'first')
# 
# variantsWithHPO = variantCertain[,c("diseaseId", "targetId")]
# variantsWithHPO = rbind(variantsWithHPO, HPOGenes[,c("diseaseId", "targetId")])

## other HPO file (more HPO terms included)
HPOGenes2 = read.table('Datasets/HPO/genes_to_phenotype.txt', sep = '\t', header = T)
HPOGenes2$hpo_id = gsub(':', '_', HPOGenes2$hpo_id)

HPOGenes2$ncbi_gene_id = as.character(HPOGenes2$ncbi_gene_id)
HPOGenes2$targetId = mapIds(x = org.Hs.eg.db, keys = HPOGenes2$ncbi_gene_id, column = 'ENSEMBL', keytype = 'ENTREZID', multiVals = 'first')

HPOGenes2 = HPOGenes2[,c("hpo_id", "targetId", "hpo_name")]
HPOGenes2$hpo_name = gsub('$', '_HP',HPOGenes2$hpo_name)
colnames(HPOGenes2)[1] = "diseaseId"

# variantsWithHPO = rbind(variantsWithHPO, HPOGenes2[,c("diseaseId", "targetId")])

# Add mouse phenotypes

mouse = read.csv('Datasets/mousePhenotypes/phenotypeAll.csv')
mouse$modelPhenotypeLabel = gsub('$', '_MP',mouse$modelPhenotypeLabel)
mouse$modelPhenotypeId = gsub(':', '_', mouse$modelPhenotypeId)

colnames(mouse)[c(2,4)] = c("diseaseId", "targetId")

variantsWithHPOandMP = variantCertain[,c("diseaseId", "targetId")]
variantsWithHPOandMP = rbind(variantsWithHPOandMP, HPOGenes2[,c("diseaseId", "targetId")], mouse[,c("diseaseId", "targetId")],HPOGenes[,c("diseaseId", "targetId")])
variantsWithHPOandMP = unique(variantsWithHPOandMP)

# Run pageRank ----
pageRankScores <- lapply(diseasesCilia, FUN = function(x){getPageRank(x, intGraphClean, variantCertain)})
names(pageRankScores) <- diseasesCilia
pageRankScores = pageRankScores[!is.na(pageRankScores)]

pageRankDF <- as.data.frame(matrix(nrow = length(pageRankScores), ncol = vcount(intGraphClean)))
for (i in 1:length(pageRankScores)){
  pageRankDF[i,] = pageRankScores[[i]]$pageRank
}

colnames(pageRankDF) = V(intGraphClean)$name
rownames(pageRankDF) = cilioEFO$description[match(names(pageRankScores), cilioEFO$ID)]

pageRankDF = pageRankDF[-grep('spondy', rownames(pageRankDF)),] #exclude spondylometaphyseal dysplasia-cone-rod dystrophy syndrome

# Run random page rank for every ciliopathy

pageRankResRandom = as.data.frame(matrix(nrow = 100, ncol = length(V(intGraphClean))))
colnames(pageRankResRandom) = V(intGraphClean)$name

randomPageRank <- function(disease, geneVariant) {
  pageRankResRand = as.data.frame(matrix(nrow = 100, ncol = length(V(intGraphClean))))
  colnames(pageRankResRand) = V(intGraphClean)$name
  for (i in 1:100) {
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

pageRankResRandom <- lapply(names(pageRankScores), FUN = function(x){randomPageRank(x, variantCertain)})
names(pageRankResRandom) = rownames(pageRankDF)

pageRankDFPvalue = as.data.frame(matrix(nrow = nrow(pageRankDF), ncol = ncol(pageRankDF)))
for (i in 1:ncol(pageRankDF)){
  for (j in 1:nrow(pageRankDF)){
    pageRankDFPvalue[j,i] = t.test(mu = pageRankDF[j,i], pageRankResRandom[[j]][,i], alternative = 'less')$p.value #wilcox.test better as not perse normal distributions
  }
}

pageRankDFPvalue = -log10(pageRankDFPvalue)

colnames(pageRankDFPvalue) = colnames(pageRankDF)
rownames(pageRankDFPvalue) = rownames(pageRankDF)

pageRankDFPvalueAdjusted = apply(pageRankDFPvalue, 1, function(x) p.adjust(10^-x))
pageRankDFPvalueAdjusted  = -log10(pageRankDFPvalueAdjusted)
pageRankDFPvalueAdjusted = t(pageRankDFPvalueAdjusted)

#cluster disorders based on pageRank scores
distTraits = dist(pageRankDFCilia, method = 'manhattan')
hclustTraits = hclust(distTraits)
plot(hclustTraits)
clustTraits = cutree(hclustTraits, h = 0.7)

#also cluster based on seed gene similarity
nVariants = table(unique(variantCertainCilia[,c("diseaseId", "targetId")])[,1])
variantCertainCilia = variantCertainCilia[variantCertainCilia$diseaseId %in% names(nVariants[nVariants > 1]),]

variantsTable = variantCertainCilia[,c("targetId", "diseaseId")]
variantsTable = unique(variantsTable[,1:2])
variantsTable$mutated = 1
variantsTable = spread(variantsTable,key = diseaseId, value = mutated)
variantsTable = column_to_rownames(variantsTable, var = 'targetId' )
variantsTable[is.na(variantsTable)] = 0
variantsTable = t(variantsTable)
rownames(variantsTable) = cilioEFO$description[match(rownames(variantsTable), cilioEFO$ID)]
plot(hclust(dist(variantsTable, method = 'binary')))

#also cluster based on phenotypes
phenotype = read.csv('Datasets/diseaseToPhenotype/phenotypeAll.csv')
phenotype = phenotype[,-1]
phenotypeCilia = phenotype[phenotype$disease %in% diseases,]

phenotypeTable = phenotypeCilia
phenotypeTable$present = 1
phenotypeTable = spread(phenotypeTable,key = disease, value = present)
phenotypeTable = column_to_rownames(phenotypeTable, var = 'phenotype' )
phenotypeTable[is.na(phenotypeTable)] = 0
phenotypeTable = t(phenotypeTable)
rownames(phenotypeTable) = cilioEFO$description[match(rownames(phenotypeTable), cilioEFO$ID)]
plot(hclust(dist(phenotypeTable, method = 'binary')))

# PCA ----
pcaPageRank = prcomp(t(scale(t(pageRankDF))))
ggplot(as.data.frame(pcaPageRank$x), aes(x=PC1, y=PC2, label=rownames(pcaPageRank$x))) + geom_point() + geom_label()
pc1Contrib = factoextra::fviz_contrib(pcaPageRank,choice = 'var', top = 10, axes = 1)
pc2Contrib = factoextra::fviz_contrib(pcaPageRank,choice = 'var', top = 10, axes = 2)

pc1Contrib = pc1Contrib$data[order(pc1Contrib$data$contrib, decreasing = T),]
pc2Contrib = pc2Contrib$data[order(pc2Contrib$data$contrib, decreasing = T),]

pc1Contrib$name = AnnotationDbi::mapIds(org.Hs.eg.db, keys = rownames(pc1Contrib), column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first')
pc2Contrib$name = AnnotationDbi::mapIds(org.Hs.eg.db, keys = rownames(pc2Contrib), column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first')

pc1Contrib$variant = 0
pc1Contrib[rownames(pc1Contrib) %in% variantsGeneticAssocCilia$targetId, "variant"] = 1
pc2Contrib$variant = 0
pc2Contrib[rownames(pc2Contrib) %in% variantsGeneticAssocCilia$targetId, "variant"] = 1

ggplot(pc1Contrib[1:100,], aes(y=contrib,x=reorder(name, -contrib), fill=variant)) + 
  geom_bar(stat = 'identity') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(pc2Contrib[1:100,], aes(y=contrib,x=reorder(name, -contrib), fill=variant)) + 
  geom_bar(stat = 'identity') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# UMAP/TSNE ----

pageRankDF = pageRankDF[,apply(pageRankDF,2,sum) != 0]
umapDisease = umap::umap(scale(pageRankDF), min_dist = 0.01, spread = 0.1)
colnames(umapDisease$layout) = c('UMAP1', 'UMAP2')
ggplot(as.data.frame(umapDisease$layout), aes(x=UMAP1, y=UMAP2, label = rownames(umapDisease$layout))) + geom_point() + geom_label()

tsneDisease = Rtsne::Rtsne(scale(pageRankDF), perplexity = 15, check_duplicates = F)
colnames(tsneDisease$Y) = c('TSNE1', 'TSNE2')
ggplot(as.data.frame(tsneDisease$Y), aes(x=TSNE1, y=TSNE2, label = rownames(pageRankDF))) + geom_point() + geom_label()

ggplot(as.data.frame(umapDisease$layout), aes(x=UMAP1, y=UMAP2, label = rownames(umapDisease$layout), color = as.character(kmeansDisease$cluster))) + geom_point() + geom_label()
ggplot(as.data.frame(tsneDisease$Y), aes(x=TSNE1, y=TSNE2, label = rownames(pageRankDF), color = as.character(kmeansDisease$cluster))) + geom_point() + geom_label()


# New gene discovery ----

#genes with highest pageRank scores which were not seed genes
pageRankDFWoSeeds = pageRankDF[,colnames(pageRankDF) %notin% variantCertainCilia$targetId]
candidateDFWoSeedsLong = pageRankDFWoSeeds %>% as.data.frame %>% rownames_to_column(var = 'disease') %>% gather(., key = 'gene', value = 'pageRank', -disease)
candidateDFWoSeedsLong$symbol = AnnotationDbi::mapIds(org.Hs.eg.db, keys = candidateDFWoSeedsLong$gene, column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first')

ggplot(candidateDFWoSeedsLong, aes(group = disease, y = pageRank, x= disease)) + 
  geom_boxplot() + 
  geom_label(data=candidateDFWoSeedsLong[candidateDFWoSeedsLong$pageRank > 0.001,], aes(label = symbol), label.padding = unit(0.05, units = 'cm'), size = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#most significant genes which were not seeds 
pageRankDFWoSeeds = pageRankDFPvalueAdjusted[,colnames(pageRankDFPvalueAdjusted) %notin% variantCertainCilia$targetId]
candidateDFWoSeedsLong = pageRankDFWoSeeds %>% as.data.frame %>% rownames_to_column(var = 'disease') %>% gather(., key = 'gene', value = 'pvalue', -disease)
candidateDFWoSeedsLong$symbol = AnnotationDbi::mapIds(org.Hs.eg.db, keys = candidateDFWoSeedsLong$gene, column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first')

ggplot(candidateDFWoSeedsLong, aes(group = disease, y = pvalue, x= disease)) + 
  geom_boxplot() + 
  geom_label(data=candidateDFWoSeedsLong[candidateDFWoSeedsLong$pvalue > 150,], aes(label = symbol), label.padding = unit(0.05, units = 'cm'), size = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#not including ciliopathy subtypes
candidateDFWoSeedsLong = candidateDFWoSeedsLong[-grep('syndrome 1|dystrophy 3|related disorders|sis 1|syndrome I|^Polycy|type 2A',candidateDFWoSeedsLong$disease),]

ggplot(candidateDFWoSeedsLong, aes(group = disease, y = pvalue, x= disease)) + 
  geom_boxplot() + 
  geom_label(data=candidateDFWoSeedsLong[candidateDFWoSeedsLong$pvalue > 200,], aes(label = symbol), label.padding = unit(0.05, units = 'cm'), size = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#check if genes are cilia-related
ciliaGenes <- readxl::read_xlsx('../Literature/SCGSv2.xlsx', skip = 1)
candidateDFWoSeedsLong$cilia = 0
candidateDFWoSeedsLong$cilia[candidateDFWoSeedsLong$gene %in% ciliaGenes$`Ensembl ID`] = 1
ggplot(candidateDFWoSeedsLong, aes(group = cilia, y= pvalue)) + geom_boxplot()
ggplot(candidateDFWoSeedsLong, aes(group = disease, y = pvalue, x= disease)) + 
  geom_boxplot() + 
  geom_label(data=candidateDFWoSeedsLong[candidateDFWoSeedsLong$pvalue > 150,], aes(label = symbol, color = cilia), label.padding = unit(0.05, units = 'cm'), size = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#plot median and SD of gene pageRank scores to identify genes with high scores for all disorders or only for specific disorders 
geneMedian = apply(pageRankDFPvalueAdjusted, 2, median)
geneSD =  apply(pageRankDFPvalueAdjusted, 2, sd)
geneSignificant = apply(pageRankDFPvalueAdjusted, 2, FUN = function(x) sum(x > 3))
geneMax = apply(pageRankDFPvalueAdjusted, 2, max)
geneOverview = data.frame(median = geneMedian, SD = geneSD, name = names(geneMedian), significant = geneSignificant, max = geneMax)
geneOverview$symbol = AnnotationDbi::mapIds(org.Hs.eg.db, keys = geneOverview$name, column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first')
geneOverview$seed = 0 
geneOverview$seed[geneOverview$name %in% variantCertainCilia$targetId] = 1
ggplot(geneOverview, aes(x=max, y=significant, color = seed)) + geom_point() + ggrepel::geom_text_repel(data = geneOverview[geneOverview$max > 300 | geneOverview$significant > 25,], aes(label = symbol))

ggplot(geneOverview[geneOverview$seed == 0,], aes(x=max, y=significant)) + 
  geom_point() + 
  ggrepel::geom_text_repel(data = geneOverview[(geneOverview$max > 250 | geneOverview$significant > 22) & geneOverview$seed == 0,], aes(label = symbol))

#check scores for other genes in open targets related to ciliopathies: either other source or low score for genetic associations
variantUncertainCilia = variantUncertainCilia[variantUncertainCilia$targetId %notin% variantCertainCilia$targetId,]
variantUncertainCilia = variantUncertainCilia[variantUncertainCilia$diseaseName %in% rownames(pageRankDFPvalueAdjusted),]
for (i in 1:nrow(variantUncertainCilia)){
  if (variantUncertainCilia$targetId[i] %in% colnames(pageRankDFPvalueAdjusted)){
    variantUncertainCilia$pvalue[i] <- pageRankDFPvalueAdjusted[variantUncertainCilia$diseaseName[i], variantUncertainCilia$targetId[i]]
  } else {
    variantUncertainCilia$pvalue[i] <- 0
  }
}

variantUncertainCilia$geneSymbol = AnnotationDbi::mapIds(org.Hs.eg.db, keys = variantUncertainCilia$targetId, column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first')

ggplot(variantUncertainCilia, aes(group = diseaseName, y = pvalue, x= diseaseName)) + 
  geom_boxplot() + 
  geom_label(data=variantUncertainCilia[variantUncertainCilia$pvalue > 3,], aes(label = geneSymbol), label.padding = unit(0.05, units = 'cm'), size = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#color by cross validation significance 
for (i in 1:nrow(variantUncertainCilia)){
  if (variantUncertainCilia$diseaseId[i] %in% names(pageRankScoresCrossVal)){
    variantUncertainCilia$crossVal[i] <- pageRankScoresCrossVal[[variantUncertainCilia$diseaseId[i]]][variantUncertainCilia$targetId[i],26]
  } else {
    variantUncertainCilia$crossVal[i] <- FALSE
  }
}

ggplot(variantUncertainCilia, aes(group = diseaseName, y = pvalue, x= diseaseName)) + 
  geom_boxplot() + 
  geom_label(data=variantUncertainCilia[variantUncertainCilia$pvalue > 3,], aes(label = geneSymbol, color = crossVal), label.padding = unit(0.05, units = 'cm'), size = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

##check if variants are in modules
variantNew = variantUncertainCilia[variantUncertainCilia$pvalue > 3 & variantUncertainCilia$crossVal == T,]

for (i in 1:nrow(variantNew)){
  variantNew$inModules[i] = sum(modulesAllComb$diseaseId == variantNew$diseaseId[i] & modulesAllComb$gene == variantNew$targetId[i])
}

for (i in 1:nrow(variantNew)){
  if (variantNew$inModules[i] == 1){
    variantNew$module[i] = unlist(modulesAllComb[modulesAllComb$diseaseId == variantNew$diseaseId[i] & modulesAllComb$gene == variantNew$targetId[i], 'moduleId'])
  } else {
    variantNew$module[i] = NA
  }}

for (i in 1:nrow(variantNew)){
  if (!is.na(variantNew$module[i]) & unlist(variantNew$module[i]) %in% overlapModules$moduleName){
    variantNew$overlapModule[i] = overlapModules[overlapModules$moduleName %in% unlist(variantNew$module[i]),"moduleNo"]
    variantNew$overlapModuleGO[i] = overlapModules[overlapModules$moduleName %in% unlist(variantNew$module[i]),"GO"]
  } else {
    variantNew$overlapModule[i] =  NA
    variantNew$overlapModuleGO[i] = NA
  }
}

#match DOID names to OT names
DOIDtoOT = data.frame('DOID' = rownames(pageRankDFPvalueAdjustedDOID))
for (i in 1:nrow(DOIDtoOT)){
  if (DOIDtoOT$DOID[i] %in% rownames(pageRankDFPvalueAdjusted)){
    DOIDtoOT$OT[i] = DOIDtoOT$DOID[i]
  } else {
    DOIDtoOT$OT[i] = rownames(pageRankDFPvalueAdjusted)[grep(paste0(DOIDtoOT$DOID[i], '$'), rownames(pageRankDFPvalueAdjusted), ignore.case = T)]
  }
}

for (i in 1:nrow(variantUncertainCilia)){
  if ((variantUncertainCilia$targetId[i] %in% colnames(pageRankDFPvalueAdjustedDOID)) & (sum(DOIDtoOT$OT == variantUncertainCilia$diseaseName[i]) > 0)){
    variantUncertainCilia$pvalueDOID[i] <- pageRankDFPvalueAdjustedDOID[DOIDtoOT$DOID[DOIDtoOT$OT == variantUncertainCilia$diseaseName[i]], variantUncertainCilia$targetId[i]]
  } else {
    variantUncertainCilia$pvalueDOID[i] <- 0
  }
}


variantUncertainCilia$DOIDSignificant = 'no'
variantUncertainCilia$DOIDSignificant[variantUncertainCilia$pvalueDOID > 3] = 'yes'

ggplot(variantUncertainCilia, aes(group = diseaseName, y = pvalue, x= diseaseName)) + 
  geom_boxplot() + 
  geom_label(data=variantUncertainCilia[variantUncertainCilia$pvalue > 3,], aes(label = geneSymbol, color = DOIDSignificant), label.padding = unit(0.05, units = 'cm'), size = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

plot(variantUncertainCilia$score, variantUncertainCilia$pvalue)

# Ciliopathies with other phenotypes ----

variantsWithHPOandMP <- variantsWithHPOandMP[variantsWithHPOandMP$targetId %in% V(intGraphClean)$name,] #TODO check why some genes are not in network

diseasesAll = unique(variantsWithHPOandMP$diseaseId)
diseasesAll = diseasesAll[grep('EFO|MONDO|HP|MP|Orphanet', diseasesAll)] #select IDs from DOID, EFO, HP, MONDO, Orphanet, MP, GO

pageRankScoresAllTraits <- lapply(diseasesAll, FUN = function(x){getPageRank(x, intGraphClean, variantsWithHPOandMP)})
names(pageRankScoresAllTraits) <- diseasesAll
pageRankScoresAllTraits = pageRankScoresAllTraits[!is.na(pageRankScoresAllTraits)]

pageRankDFAllTraits <- data.frame(t(sapply(pageRankScoresAllTraits, FUN = function(x){c(x$pageRank)})))

colnames(pageRankDFAllTraits) = V(intGraphClean)$name

pageRankDFAllTraits = pageRankDFAllTraits[!(rownames(pageRankDFAllTraits) %in% cilioEFO$ID & rownames(pageRankDFAllTraits) %notin% diseasesCilia),] #exclude extra ciliopathies

#cluster disorders based on pageRank scores
distTraits = dist(pageRankDFAllTraits, method = 'manhattan')

distTraitsMatrix <- as.matrix(distTraits)

relatedDiseases <- c()
for (i in 1:length(diseasesCilia)){
  distCiliopathy <- distTraitsMatrix[,grep(diseasesCilia[i], colnames(distTraitsMatrix))]
  distCiliopathy <- distCiliopathy[names(distCiliopathy) %notin% diseasesCilia]
  if (length(distCiliopathy != 0)){
    relatedDiseases <- c(relatedDiseases, names(distCiliopathy[order(distCiliopathy)][1:50]))
  }
}

relatedDiseases <- c(relatedDiseases, diseasesCilia)
relatedDiseases <- unique(relatedDiseases)

relatedDiseasesNoHP <- c()
for (i in 1:length(diseasesCilia)){
  distCiliopathy <- distTraitsMatrix[-grep('HP_', rownames(distTraitsMatrix)),grep(diseasesCilia[i], colnames(distTraitsMatrix))]
  distCiliopathy <- distCiliopathy[names(distCiliopathy) %notin% diseasesCilia]
  if (length(distCiliopathy != 0)){
    relatedDiseasesNoHP <- c(relatedDiseasesNoHP, names(distCiliopathy[order(distCiliopathy)][1:50]))
  }
}

relatedDiseasesNoHP <- c(relatedDiseasesNoHP, diseasesCilia)
relatedDiseasesNoHP <- unique(relatedDiseasesNoHP)


relatedDiseasesNoHPMP <- c()
for (i in 1:length(diseasesCilia)){
  distCiliopathy <- distTraitsMatrix[-grep('HP_|MP_', rownames(distTraitsMatrix)),grep(diseasesCilia[i], colnames(distTraitsMatrix))]
  distCiliopathy <- distCiliopathy[names(distCiliopathy) %notin% diseasesCilia]
  if (length(distCiliopathy != 0)){
    relatedDiseasesNoHPMP <- c(relatedDiseasesNoHPMP, names(distCiliopathy[order(distCiliopathy)][1:50]))
  }
}

relatedDiseasesNoHPMP <- c(relatedDiseasesNoHPMP, diseasesCilia)
relatedDiseasesNoHPMP <- unique(relatedDiseasesNoHPMP)


relatedDiseasesAll = unique(c(relatedDiseases,relatedDiseasesNoHP,relatedDiseasesNoHPMP))

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

jaccardRelated = as.data.frame(matrix(nrow = length(relatedDiseasesAll), ncol = length(diseasesCilia)))
for (i in 1:(length(relatedDiseasesAll))){
  for (j in 1:length(diseasesCilia)){
    jaccardRelated[i,j] = jaccard(unique(variantsWithHPOandMP[variantsWithHPOandMP$diseaseId == relatedDiseasesAll[i],"targetId"]), 
                                  unique(variantsWithHPOandMP[variantsWithHPOandMP$diseaseId == diseasesCilia[j],"targetId"]))
  }}

relatedDiseasesAll = relatedDiseasesAll[-which(jaccardRelated > 0.7, arr.ind = T)[,1]] #exclude traits with perfect seed gene overlap
relatedDiseasesAll = unique(c(relatedDiseasesAll, diseasesCilia))

distTraitsMatrixRelatedDiseases <- distTraitsMatrix[colnames(distTraitsMatrix) %in% relatedDiseasesAll, rownames(distTraitsMatrix) %in% relatedDiseasesAll]
# distTraitsMatrixRelatedDiseases <- distTraitsMatrix[colnames(distTraitsMatrix) %in% unique(c(relatedDiseasesNoHPMP)), rownames(distTraitsMatrix) %in% unique(c(relatedDiseasesNoHPMP))]
# distTraitsMatrixRelatedDiseases <- distTraitsMatrix[colnames(distTraitsMatrix) %in% unique(c(relatedDiseasesNoHPMP, relatedDiseases, relatedDiseasesNoHP)), rownames(distTraitsMatrix) %in% unique(c(relatedDiseasesNoHPMP, relatedDiseases, relatedDiseasesNoHP))]

hclustTraits = hclust(as.dist(distTraitsMatrixRelatedDiseases))
if (exists("EFOAnno") == F) {
  EFOAnno = ontologyIndex::get_ontology('./Datasets/annotation/efo.obo')
}
EFOAnno$id = gsub(':', '_', EFOAnno$id)

isCiliopathy = rep(0, length(hclustTraits$labels))
isCiliopathy[hclustTraits$labels %notin% diseasesCilia] = 1
isCiliopathy[hclustTraits$labels %in% diseasesCilia] = 2

hclustTraits$labels[hclustTraits$labels %in% EFOAnno$id] = unname(EFOAnno$name[match(hclustTraits$labels[hclustTraits$labels %in% EFOAnno$id], EFOAnno$id)])
hclustTraits$labels[hclustTraits$labels %in% HPOGenes$diseaseId] = unname(HPOGenes$HPOLabel[match(hclustTraits$labels[hclustTraits$labels %in% HPOGenes$diseaseId], HPOGenes$diseaseId)])
hclustTraits$labels[hclustTraits$labels %in% HPOGenes2$diseaseId] = unname(HPOGenes2$hpo_name[match(hclustTraits$labels[hclustTraits$labels %in% HPOGenes2$diseaseId], HPOGenes2$diseaseId)])
hclustTraits$labels[hclustTraits$labels %in% mouse$diseaseId] = unname(mouse$modelPhenotypeLabel[match(hclustTraits$labels[hclustTraits$labels %in% mouse$diseaseId], mouse$diseaseId)])
hclustTraits$labels[hclustTraits$labels == 'MONDO_custom'] = 'cilia-associated retinitis pigmentosa'

clustTraits = cutree(hclustTraits, h = 0.65)

library(ape)
plot(as.phylo(hclustTraits), tip.color = c('black','red')[isCiliopathy], cex = 0.8, label.offset = 0.01, direction = 'downwards', edge.width = 2)


# Identify gene modules for ciliopathies and related diseases----

cores <- detectCores()
cl <- makeCluster(cores - 1)
registerDoParallel(cl)

modulesAll <- foreach(i=iter(relatedDiseases, by = 'cell'), .packages = c('igraph', 'dplyr')) %dopar% {
  pageRankModules(disease = i, intGraph = intGraphClean, geneVariant = variantsWithHPOandMP)
}

stopCluster(cl)

names(modulesAll) <- relatedDiseases

modulesAll = modulesAll[unlist(lapply(modulesAll, FUN = function(x){sum(is.na(x)) < 2}))]

#test significance of modules

ksRes <- vector(mode = "list", length = length(modulesAll))
for(i in 1:length(modulesAll)){
  for (j in names(table(modulesAll[[i]][[2]]$walktrapCluster))){
    ksRes[[i]] <- c(ksRes[[i]], ks.test(modulesAll[[i]][[1]][modulesAll[[i]][[1]]$walktrapCluster == paste0(j), "pageRank"], 
                                        modulesAll[[i]][[1]]$pageRank,
                                        alternative = 'less')$p.value)
  }}

ksRes <- lapply(ksRes, p.adjust)

for (i in 1:length(modulesAll)){
  modulesAll[[i]][[3]] <- modulesAll[[i]][[2]][modulesAll[[i]][[2]]$walktrapCluster %in% names(table(modulesAll[[i]][[2]]$walktrapCluster))[which(ksRes[[i]] < 0.05)],]
}



jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

for (i in 1:length(modulesAll)){
  modulesAll[[i]][[3]]$moduleId = paste(modulesAll[[i]][[3]]$walktrapCluster, names(modulesAll)[i], sep = '_')
}

for (i in 1:length(modulesAll)){
  modulesAll[[i]][[1]]$moduleId = paste(modulesAll[[i]][[1]]$walktrapCluster, names(modulesAll)[i], sep = '_')
}

modulesAllComb <- NULL
for (i in 1:length(modulesAll)){
  modulesAllComb <- rbind(modulesAllComb, modulesAll[[i]][[3]])
}

modulesAllCombInclNoSeed <- NULL
for (i in 1:length(modulesAll)){
  modulesAllCombInclNoSeed <- rbind(modulesAllCombInclNoSeed, modulesAll[[i]][[1]])
}

#to add extra modules to old matrix
jaccardResMat = matrix(nrow = length(unique(modulesAllComb$moduleId)), ncol = length(unique(modulesAllComb$moduleId)))
rownames(jaccardResMat) = unique(modulesAllComb$moduleId)
colnames(jaccardResMat) = unique(modulesAllComb$moduleId)

for(i in 1:ncol(jaccardResMatOld)){
  jaccardResMat[2735:4590,colnames(jaccardResMatOld)[i]] = jaccardResMatOld[,colnames(jaccardResMatOld)[i]]
  jaccardResMat[colnames(jaccardResMatOld)[i],2735:4590] = jaccardResMatOld[colnames(jaccardResMatOld)[i],]
}

jaccardResMat[lower.tri(jaccardResMat)] = Inf

jaccardResMatLong = jaccardResMat %>% as.data.frame() %>% rownames_to_column(., var = 'moduleId2') %>% gather(., value = 'jaccard', key = 'moduleId1', -'moduleId2')
jaccardResMatLong = jaccardResMatLong[jaccardResMatLong$moduleId1 != jaccardResMatLong$moduleId2 & is.na(jaccardResMatLong$jaccard),]

jaccardResMatLong$jaccard = apply(jaccardResMatLong, 1, FUN = function(x){jaccard(unlist(modulesAllComb[modulesAllComb$moduleId == x['moduleId1'],"gene"]), unlist(modulesAllComb[modulesAllComb$moduleId == x['moduleId2'], 'gene']))})

for (i in 1:nrow(jaccardResMatLong)){
  jaccardResMat[jaccardResMatLong[i,"moduleId2"], jaccardResMatLong[i,"moduleId1"]] = jaccardResMatLong[i,"jaccard"]
}

diag(jaccardResMat) = 1
jaccardResMat[lower.tri(jaccardResMat)] = t(jaccardResMat)[lower.tri(jaccardResMat)]



jaccardResMat <- matrix(nrow = length(unique(modulesAllComb$moduleId)), ncol = length(unique(modulesAllComb$moduleId)))
rownames(jaccardResMat) = unique(modulesAllComb$moduleId)
colnames(jaccardResMat) = unique(modulesAllComb$moduleId)
for (i in rownames(jaccardResMat)){
  jaccardResMat[paste0(i),] <- sapply(colnames(jaccardResMat),FUN=function(x){jaccard(unlist(modulesAllComb[modulesAllComb$moduleId == x,"gene"]), unlist(modulesAllComb[modulesAllComb$moduleId == paste0(i), 'gene']))})
}

ComplexHeatmap::Heatmap(jaccardResMat)

set.seed(12345)
jaccardUMAP = umap::umap(jaccardResMat)
jaccardkMeans = kmeans(jaccardResMat, centers = 25)
ggplot(as.data.frame(jaccardUMAP$layout), aes(x=jaccardUMAP$layout[,1], y= jaccardUMAP$layout[,2], color = as.character(jaccardkMeans$cluster))) + geom_point()


combinedModules = rep(list(sapply(c('genes', 'modules', 'diseaseId', 'diseases'), function(x) NULL)), max(jaccardkMeans$cluster))
for (i in 1:max(jaccardkMeans$cluster)){
  modules = names(jaccardkMeans$cluster[jaccardkMeans$cluster == i])
  genes =  unique(modulesAllComb[modulesAllComb$moduleId %in% modules, 'gene'])
  diseaseId = gsub('^.*\\d_','', modules)
  diseases = EFOAnno$name[match(diseaseId, EFOAnno$id)]
  combinedModules[[i]]$genes = genes
  combinedModules[[i]]$modules = modules
  combinedModules[[i]]$diseaseId = diseaseId
  combinedModules[[i]]$diseases = unname(diseases)
}

moduleGOBP <- NULL
moduleGOBP <- lapply(combinedModules, FUN = function(x){
  enrichGO(gene = x$genes, 
           OrgDb = 'org.Hs.eg.db', 
           keyType = 'ENSEMBL', 
           pvalueCutoff = 0.01, 
           pAdjustMethod = 'BH', 
           ont = 'BP')@result
})

moduleGOBP3 <- lapply(moduleGOBP, FUN = function(x){x$Description[3]})

moduleMatrix = as.data.frame(matrix(nrow = length(modulesAll), ncol = length(combinedModules)))
rownames(moduleMatrix) = EFOAnno$name[match(names(modulesAll), EFOAnno$id)]
colnames(moduleMatrix) = unlist(moduleGOBP1)

for ( i in 1:ncol(moduleMatrix)){
  for (j in rownames(moduleMatrix)){
    if (paste0(j) %in% combinedModules[[i]]$diseases){
      moduleMatrix[paste0(j),i] = 1
    } else
      moduleMatrix[paste0(j),i] = 0
  }
}

Heatmap(t(moduleMatrix), column_names_gp = gpar(fontsize = 9), row_names_gp = gpar(fontsize = 10), col = c('lightgrey', 'darkblue'))

# moduleGOBP <- list()
# for (i in unique(modulesAllComb$moduleId)) { #hard to parallelize: 'database is locked' error
#   moduleGOBP[[i]] <- enrichGO(gene = unname(unlist(modulesAllComb[modulesAllComb$moduleId == paste0(i), 'gene'])), 
#                                OrgDb = 'org.Hs.eg.db', 
#                                keyType = 'ENSEMBL', 
#                                pvalueCutoff = 0.01, 
#                                pAdjustMethod = 'BH', 
#                                ont = 'BP')@result
# }
# 
# 
# 


#Other way to identify overlapping modules
moduleInteractions = jaccardResMat %>% as.data.frame() %>% rownames_to_column(var = 'B') %>% gather(key = 'A', value = 'jaccard', -B)
moduleInteractions = moduleInteractions[moduleInteractions$jaccard >= 0.7,]
moduleGraph = graph_from_data_frame(moduleInteractions[,c('A','B')], directed = F)
moduleGraph = igraph::simplify(moduleGraph, remove.loops = T)
overlapModules = components(moduleGraph)
overlapModules = overlapModules$membership[overlapModules$membership %in% which(overlapModules$csize > 2)]
overlapModules = data.frame(moduleName = names(overlapModules), moduleNo = unname(overlapModules))

overlapModulesGenes = list()
for (i in unique(overlapModules$moduleNo)){
  overlapModulesGenes[[paste0(i)]] = unique(modulesAllComb$gene[modulesAllComb$moduleId %in% overlapModules$moduleName[overlapModules$moduleNo == paste0(i)]])
}


overlapModulesGOBP <- NULL
overlapModulesGOBP <- lapply(overlapModulesGenes, FUN = function(x){
  enrichGO(gene = x, 
           OrgDb = 'org.Hs.eg.db', 
           keyType = 'ENSEMBL', 
           pvalueCutoff = 0.01, 
           pAdjustMethod = 'BH', 
           ont = 'BP')@result
})

overlapModules$disease = gsub('.*\\d_', '', overlapModules$moduleName)

if (exists("EFOAnno") == F) {
  EFOAnno = ontologyIndex::get_ontology('./Datasets/annotation/efo.obo')
}
EFOAnno$id = gsub(':', '_', EFOAnno$id)
EFOAnno$name[grep('Orpha',EFOAnno$id)] = gsub('$', '_Orpha',EFOAnno$name[grep('Orpha',EFOAnno$id)])

overlapModules$diseaseName[overlapModules$disease %in% EFOAnno$id] = unname(EFOAnno$name[match(overlapModules$disease[overlapModules$disease %in% EFOAnno$id], EFOAnno$id)])
overlapModules$diseaseName[overlapModules$disease %in% HPOGenes$diseaseId] = unname(HPOGenes$HPOLabel[match(overlapModules$disease[overlapModules$disease %in% HPOGenes$diseaseId], HPOGenes$diseaseId)])
overlapModules$diseaseName[overlapModules$disease %in% HPOGenes2$diseaseId] = unname(HPOGenes2$hpo_name[match(overlapModules$disease[overlapModules$disease %in% HPOGenes2$diseaseId], HPOGenes2$diseaseId)])
overlapModules$diseaseName[overlapModules$disease %in% mouse$diseaseId] = unname(mouse$modelPhenotypeLabel[match(overlapModules$disease[overlapModules$disease %in% mouse$diseaseId], mouse$diseaseId)])
overlapModules$diseaseName[overlapModules$disease == 'MONDO_custom'] = 'cilia-associated retinitis pigmentosa'
overlapModules = overlapModules[-grep('\t', overlapModules$diseaseName),] #exclude trait with weird name

for (i in 1:nrow(overlapModules)){
  overlapModules[i,'GO'] = overlapModulesGOBP[[paste0(overlapModules$moduleNo[i])]]$Description[[1]]
}

overlapModules$moduleGO = paste0(overlapModules$moduleNo, '_', overlapModules$GO)

for (i in 1:nrow(overlapModules)){
  overlapModules[i,'GO_ID'] = overlapModulesGOBP[[paste0(overlapModules$moduleNo[i])]]$ID[[1]]
}

for (i in 1:nrow(overlapModules)){
  overlapModules[i,'GO2'] = overlapModulesGOBP[[paste0(overlapModules$moduleNo[i])]]$Description[[2]]
}

for (i in 1:nrow(overlapModules)){
  overlapModules[i,'GO3'] = overlapModulesGOBP[[paste0(overlapModules$moduleNo[i])]]$Description[[3]]
}

for (i in 1:nrow(overlapModules)){
  overlapModules[i,'GO4'] = overlapModulesGOBP[[paste0(overlapModules$moduleNo[i])]]$Description[[4]]
}

for (i in 1:nrow(overlapModules)){
  seed = unique(variantsWithHPOandMP[variantsWithHPOandMP$diseaseId == overlapModules[i,"disease"],"targetId"])
  overlapModules$seed[i] = list(seed[seed %in% overlapModulesGenes[[paste0(overlapModules[i,"moduleNo"])]]])
}



overlapModules$selected = 1

moduleMatrix = spread(overlapModules[,c('diseaseName', 'moduleGO', 'selected')], key = diseaseName, value = selected) %>% column_to_rownames(var = 'moduleGO')
moduleMatrix[is.na(moduleMatrix)] = 0

Heatmap((moduleMatrix), column_names_gp = gpar(fontsize = 9), row_names_gp = gpar(fontsize = 10), col = c('lightgrey', 'darkblue'))

Heatmap((moduleMatrix[grep('Usher|deaf|blind|leber|ear', colnames(moduleMatrix), ignore.case = T)]), column_names_gp = gpar(fontsize = 9), row_names_gp = gpar(fontsize = 10), col = c('lightgrey', 'darkblue'))

Heatmap((moduleMatrix[grep('joubert|meckel|bardet', colnames(moduleMatrix), ignore.case = T)]), column_names_gp = gpar(fontsize = 9), row_names_gp = gpar(fontsize = 10), col = c('lightgrey', 'darkblue'))


exocytosis = intHigh[intHigh$targetA %in% overlapModulesGenes[names(overlapModulesGenes) == '1'][[1]] & intHigh$targetB %in% overlapModulesGenes[names(overlapModulesGenes) == '1'][[1]],]
exocytosis = exocytosis[(exocytosis$sourceDatabase == 'string' & exocytosis$scoring >= 0.9) | (exocytosis$sourceDatabase == 'intact' & exocytosis$scoring >= 0.4),]
exocytosis$targetA = AnnotationDbi::mapIds(org.Hs.eg.db, keys =exocytosis$targetA, column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first')
exocytosis$targetB = AnnotationDbi::mapIds(org.Hs.eg.db, keys =exocytosis$targetB, column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first')
write.table(exocytosis, 'module1Exocytosis.txt')

ciliumAssembly = intHigh[intHigh$targetA %in% overlapModulesGenes[names(overlapModulesGenes) == '16'][[1]] & intHigh$targetB %in% overlapModulesGenes[names(overlapModulesGenes) == '16'][[1]],]
ciliumAssembly = ciliumAssembly[(ciliumAssembly$sourceDatabase == 'string' & ciliumAssembly$scoring >= 0.9) | (ciliumAssembly$sourceDatabase == 'intact' & ciliumAssembly$scoring >= 0.4),]
ciliumAssembly$targetA = AnnotationDbi::mapIds(org.Hs.eg.db, keys =ciliumAssembly$targetA, column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first')
ciliumAssembly$targetB = AnnotationDbi::mapIds(org.Hs.eg.db, keys =ciliumAssembly$targetB, column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first')
write.table(ciliumAssembly, 'module16ciliumAssembly.csv')

intraCiliaryTransport = intHigh[intHigh$targetA %in% overlapModulesGenes[names(overlapModulesGenes) == '2'][[1]] & intHigh$targetB %in% overlapModulesGenes[names(overlapModulesGenes) == '2'][[1]],]
intraCiliaryTransport = intraCiliaryTransport[(intraCiliaryTransport$sourceDatabase == 'string' & intraCiliaryTransport$scoring >= 0.9) | (intraCiliaryTransport$sourceDatabase == 'intact' & intraCiliaryTransport$scoring >= 0.4),]
intraCiliaryTransport$targetA = AnnotationDbi::mapIds(org.Hs.eg.db, keys =intraCiliaryTransport$targetA, column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first')
intraCiliaryTransport$targetB = AnnotationDbi::mapIds(org.Hs.eg.db, keys =intraCiliaryTransport$targetB, column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first')
write.table(intraCiliaryTransport, 'module2IntraciliaryTransport.txt')

ciliumOrganization = intHigh[intHigh$targetA %in% overlapModulesGenes[names(overlapModulesGenes) == '25'][[1]] & intHigh$targetB %in% overlapModulesGenes[names(overlapModulesGenes) == '25'][[1]],]
ciliumOrganization = ciliumOrganization[(ciliumOrganization$sourceDatabase == 'string' & ciliumOrganization$scoring >= 0.9) | (ciliumOrganization$sourceDatabase == 'intact' & ciliumOrganization$scoring >= 0.4),]
ciliumOrganization$targetA = AnnotationDbi::mapIds(org.Hs.eg.db, keys =ciliumOrganization$targetA, column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first')
ciliumOrganization$targetB = AnnotationDbi::mapIds(org.Hs.eg.db, keys =ciliumOrganization$targetB, column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first')
write.table(ciliumOrganization, 'module25ciliumOrganization.txt')

ciliumorganization25 = read.table('module25ciliumOrganization.txt')
ciliumorganization25Seed = data.frame('symbol' = unique(c(ciliumorganization25$targetA, ciliumorganization25$targetB)))
ciliumorganization25Seed$ensembl = AnnotationDbi::mapIds(org.Hs.eg.db, keys = ciliumorganization25Seed$symbol, column = 'ENSEMBL', keytype = 'SYMBOL', multiVals = 'first')
ciliumorganization25Seed = cbind(ciliumorganization25Seed, as.data.frame(matrix(nrow = nrow(ciliumorganization25Seed), ncol = length(overlapModules[overlapModules$moduleNo == '25',"diseaseName"]))))
colnames(ciliumorganization25Seed)[3:ncol(ciliumorganization25Seed)] =  overlapModules[overlapModules$moduleNo == '25',"diseaseName"]
for (i in 3:ncol(ciliumorganization25Seed)){
  ciliumorganization25Seed[,i] = ciliumorganization25Seed$ensembl %in% unlist(overlapModules[overlapModules$moduleNo == '25' & overlapModules$diseaseName == colnames(ciliumorganization25Seed)[i], "seed"])
}
write.table(ciliumorganization25Seed, 'module25_ciliumorganization_seedGenes.txt', quote = F, row.names = F, sep = '\t')


ciliumOrganization = intHigh[intHigh$targetA %in% overlapModulesGenes[names(overlapModulesGenes) == '41'][[1]] & intHigh$targetB %in% overlapModulesGenes[names(overlapModulesGenes) == '41'][[1]],]
ciliumOrganization = ciliumOrganization[(ciliumOrganization$sourceDatabase == 'string' & ciliumOrganization$scoring >= 0.9) | (ciliumOrganization$sourceDatabase == 'intact' & ciliumOrganization$scoring >= 0.4),]
ciliumOrganization$targetA = AnnotationDbi::mapIds(org.Hs.eg.db, keys =ciliumOrganization$targetA, column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first')
ciliumOrganization$targetB = AnnotationDbi::mapIds(org.Hs.eg.db, keys =ciliumOrganization$targetB, column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first')
write.table(ciliumOrganization, 'module41ciliumOrganization.txt')

microtubuleOrganization = intHigh[intHigh$targetA %in% overlapModulesGenes[names(overlapModulesGenes) == '17'][[1]] & intHigh$targetB %in% overlapModulesGenes[names(overlapModulesGenes) == '17'][[1]],]
microtubuleOrganization = microtubuleOrganization[(microtubuleOrganization$sourceDatabase == 'string' & microtubuleOrganization$scoring >= 0.9) | (microtubuleOrganization$sourceDatabase == 'intact' & microtubuleOrganization$scoring >= 0.4),]
microtubuleOrganization$targetA = AnnotationDbi::mapIds(org.Hs.eg.db, keys =microtubuleOrganization$targetA, column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first')
microtubuleOrganization$targetB = AnnotationDbi::mapIds(org.Hs.eg.db, keys =microtubuleOrganization$targetB, column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first')
write.csv(microtubuleOrganization, 'module17microtubuleOrganization.csv')

polyubiquination = intHigh[intHigh$targetA %in% overlapModulesGenes[names(overlapModulesGenes) == '39'][[1]] & intHigh$targetB %in% overlapModulesGenes[names(overlapModulesGenes) == '39'][[1]],]
polyubiquination = polyubiquination[(polyubiquination$sourceDatabase == 'string' & polyubiquination$scoring >= 0.9) | (polyubiquination$sourceDatabase == 'intact' & polyubiquination$scoring >= 0.4),]
polyubiquination$targetA = AnnotationDbi::mapIds(org.Hs.eg.db, keys =polyubiquination$targetA, column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first')
polyubiquination$targetB = AnnotationDbi::mapIds(org.Hs.eg.db, keys =polyubiquination$targetB, column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first')
write.csv(polyubiquination, 'module39polyubiquination.csv')

phosphatidylinositol = intHigh[intHigh$targetA %in% overlapModulesGenes[names(overlapModulesGenes) == '53'][[1]] & intHigh$targetB %in% overlapModulesGenes[names(overlapModulesGenes) == '53'][[1]],]
phosphatidylinositol = phosphatidylinositol[(phosphatidylinositol$sourceDatabase == 'string' & phosphatidylinositol$scoring >= 0.9) | (phosphatidylinositol$sourceDatabase == 'intact' & phosphatidylinositol$scoring >= 0.4),]
phosphatidylinositol$targetA = AnnotationDbi::mapIds(org.Hs.eg.db, keys =phosphatidylinositol$targetA, column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first')
phosphatidylinositol$targetB = AnnotationDbi::mapIds(org.Hs.eg.db, keys =phosphatidylinositol$targetB, column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first')

proteinOxidation = intHigh[intHigh$targetA %in% overlapModulesGenes[names(overlapModulesGenes) == '319'][[1]] & intHigh$targetB %in% overlapModulesGenes[names(overlapModulesGenes) == '319'][[1]],]
proteinOxidation = proteinOxidation[(proteinOxidation$sourceDatabase == 'string' & proteinOxidation$scoring >= 0.9) | (proteinOxidation$sourceDatabase == 'intact' & proteinOxidation$scoring >= 0.4),]
proteinOxidation$targetA = AnnotationDbi::mapIds(org.Hs.eg.db, keys =proteinOxidation$targetA, column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first')
proteinOxidation$targetB = AnnotationDbi::mapIds(org.Hs.eg.db, keys =proteinOxidation$targetB, column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first')

moduleNetwork = overlapModulesSelectedInt[,c("diseaseName", "moduleGO")]
write.csv(moduleNetwork, 'moduleTraitInteractions_SelectedInt_top25Genes.csv')

for (i in 1:nrow(overlapModules)){
  overlapModules$GO_offspring[i] = length(as.list(GO.db::GOBPOFFSPRING[paste0(overlapModules$GO_ID[i])])[[1]])
}

for (i in 1:nrow(overlapModules)){
  overlapModules$GO_size[i] = as.numeric(gsub('/.*','',overlapModulesGOBP[names(overlapModulesGOBP) == overlapModules$moduleNo[i]][[1]]$BgRatio[[1]]))
}

# test for enrichment of overlapping modules per page rank cluster
hclustTraits = hclust(as.dist(distTraitsMatrixRelatedDiseases[rownames(distTraitsMatrixRelatedDiseases) %in% overlapModules$disease,rownames(distTraitsMatrixRelatedDiseases) %in% overlapModules$disease]), method = 'ward.D2')
clustTraits = cutree(hclustTraits, h = 1)

overlapModules$pagerankClust = unname(clustTraits)[match(overlapModules$disease, names(clustTraits))]
overlapModules = overlapModules[!is.na(overlapModules$pagerankClust),]
modulePerClust = unique(overlapModules[,c("moduleGO", "pagerankClust")])

for (i in 1:nrow(modulePerClust)){
  modulePerClust$enrichP[i] = phyper(nrow(overlapModules[overlapModules$pagerankClust == modulePerClust[i,"pagerankClust"] & overlapModules$moduleGO == modulePerClust[i,"moduleGO"],]) - 1,  
                                     nrow(overlapModules[overlapModules$moduleGO == modulePerClust[i,"moduleGO"],]), 
                                     nrow(overlapModules) -  nrow(overlapModules[overlapModules$moduleGO == modulePerClust[i,"moduleGO"],]), 
                                     nrow(overlapModules[overlapModules$pagerankClust == modulePerClust[i,"pagerankClust"],]), lower.tail = F)
}

modulePerClust$enrichP = -log10(modulePerClust$enrichP)

for (i in 1:nrow(overlapModules)){
  overlapModules[i, 'enrichP'] = modulePerClust[modulePerClust$moduleGO == overlapModules[i,"moduleGO"] & modulePerClust$pagerankClust == overlapModules[i, "pagerankClust"], "enrichP"]
}

overlapModulesWide = spread(overlapModules[,c("diseaseName", "moduleGO", "pagerankClust", "enrichP")], key = moduleGO, value = enrichP)
overlapModulesWide = column_to_rownames(overlapModulesWide, var = 'diseaseName')
overlapModulesWide[is.na(overlapModulesWide)] = 0

pagerankClustHeatmap = overlapModulesWide$pagerankClust
overlapModulesWide = overlapModulesWide[,-1]

overlapModulesWide[overlapModulesWide > 0 & overlapModulesWide < 1.3] = 0.5
overlapModulesWide[overlapModulesWide >= 1.3] = 1

hclustIds = hclustTraits$labels
hclustTraits$labels[hclustIds %in% EFOAnno$id] = unname(EFOAnno$name[match(hclustIds[hclustIds %in% EFOAnno$id], EFOAnno$id)])
hclustTraits$labels[hclustIds %in% HPOGenes$diseaseId] = unname(HPOGenes$HPOLabel[match(hclustIds[hclustIds %in% HPOGenes$diseaseId], HPOGenes$diseaseId)])
hclustTraits$labels[hclustIds %in% HPOGenes2$diseaseId] = unname(HPOGenes2$hpo_name[match(hclustIds[hclustIds %in% HPOGenes2$diseaseId], HPOGenes2$diseaseId)])
hclustTraits$labels[hclustIds %in% mouse$diseaseId] = unname(mouse$modelPhenotypeLabel[match(hclustIds[hclustIds %in% mouse$diseaseId], mouse$diseaseId)])
hclustTraits$labels[hclustIds == 'MONDO_custom'] = 'cilia-associated retinitis pigmentosa'

overlapModulesWide = overlapModulesWide[match(hclustTraits$labels, rownames(overlapModulesWide)),]

ComplexHeatmap::Heatmap(overlapModulesWide, split = max(clustTraits), col =structure(c('white', 'grey', 'red'), names = c("0", "0.5", "1")), gap = unit(0.1, 'cm'), cluster_rows = hclustTraits, row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6))


modulePerClustWide = spread(modulePerClust, key = 'moduleGO', value = enrichP)
modulePerClustWide = column_to_rownames(modulePerClustWide, var = 'pagerankClust')
modulePerClustWide[is.na(modulePerClustWide)] = 0

ComplexHeatmap::Heatmap(modulePerClustWide)

# test for enrichment of page rank scores in cilia modules

ciliaModule = read.csv('modulesCiliaRelatedSysciliaGenes.csv')
pageRankDFCilia = pageRankDFAllTraits[rownames(pageRankDFAllTraits) %in% diseasesCilia, ]
pageRankDFRelated = pageRankDFAllTraits[rownames(pageRankDFAllTraits) %in% relatedDiseases, ]

ciliaModuleEnrich = list()
for (i in 1:nrow(pageRankDFRelated)) {
  ciliaModuleEnrich[[i]] = clusterProfiler::GSEA(
    unlist(pageRankDFRelated[i, order(pageRankDFRelated[i, ], decreasing = T), drop = T]),
    TERM2GENE = ciliaModule[, c("V1", "V2")],
    minGSSize = 5,
    verbose = F
  )@result
}

names(ciliaModuleEnrich) = rownames(pageRankDFRelated)
moduleMatrixEnrich = matrix(nrow = nrow(pageRankDFRelated), ncol = length(unique(ciliaModule$V1)))
colnames(moduleMatrixEnrich) = unique(ciliaModule$V1)
for (i in 1:nrow(moduleMatrixEnrich)) {
  moduleMatrixEnrich[i, colnames(moduleMatrixEnrich) %in% ciliaModuleEnrich[[i]]$Description] = 1
}
moduleMatrixEnrich[is.na(moduleMatrixEnrich)] = 0
rownames(moduleMatrixEnrich) =  EFOAnno$name[match(names(ciliaModuleEnrich), EFOAnno$id)]

moduleMatrixEnrich = moduleMatrixEnrich[,colSums(moduleMatrixEnrich) != 0]
moduleMatrixEnrich = moduleMatrixEnrich[rowSums(moduleMatrixEnrich) != 0,]

ComplexHeatmap::Heatmap(moduleMatrixEnrich, row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 8))



hclustTraits = hclust(as.dist(distTraitsMatrixRelatedDiseases[rownames(distTraitsMatrixRelatedDiseases) %in% ciliaModuleEnrich$diseaseId,rownames(distTraitsMatrixRelatedDiseases) %in% ciliaModuleEnrich$diseaseId]), method = 'ward.D2')
clustTraits = cutree(hclustTraits, h = 1)

ciliaModuleEnrich$pagerankClust = unname(clustTraits)[match(ciliaModuleEnrich$diseaseId, names(clustTraits))]
ciliaModuleEnrich = ciliaModuleEnrich[!is.na(ciliaModuleEnrich$pagerankClust),]
modulePerClust = unique(ciliaModuleEnrich[,c("complex", "pagerankClust")])

for (i in 1:nrow(modulePerClust)){
  modulePerClust$enrichP[i] = phyper(nrow(ciliaModuleEnrich[ciliaModuleEnrich$pagerankClust == modulePerClust[i,"pagerankClust"] & ciliaModuleEnrich$complex == modulePerClust[i,"complex"],]) - 1,  
                                     nrow(ciliaModuleEnrich[ciliaModuleEnrich$complex == modulePerClust[i,"complex"],]), 
                                     nrow(ciliaModuleEnrich) -  nrow(ciliaModuleEnrich[ciliaModuleEnrich$complex == modulePerClust[i,"complex"],]), 
                                     nrow(ciliaModuleEnrich[ciliaModuleEnrich$pagerankClust == modulePerClust[i,"pagerankClust"],]), lower.tail = F)
}

modulePerClust$enrichP = -log10(modulePerClust$enrichP)

for (i in 1:nrow(ciliaModuleEnrich)){
  ciliaModuleEnrich[i, 'enrichP'] = modulePerClust[modulePerClust$complex == ciliaModuleEnrich[i,"complex"] & modulePerClust$pagerankClust == ciliaModuleEnrich[i, "pagerankClust"], "enrichP"]
}

hclustIds = hclustTraits$labels
hclustTraits$labels[hclustIds %in% EFOAnno$id] = unname(EFOAnno$name[match(hclustIds[hclustIds %in% EFOAnno$id], EFOAnno$id)])
hclustTraits$labels[hclustIds %in% HPOGenes$diseaseId] = unname(HPOGenes$HPOLabel[match(hclustIds[hclustIds %in% HPOGenes$diseaseId], HPOGenes$diseaseId)])
hclustTraits$labels[hclustIds %in% HPOGenes2$diseaseId] = unname(HPOGenes2$hpo_name[match(hclustIds[hclustIds %in% HPOGenes2$diseaseId], HPOGenes2$diseaseId)])
hclustTraits$labels[hclustIds %in% mouse$diseaseId] = unname(mouse$modelPhenotypeLabel[match(hclustIds[hclustIds %in% mouse$diseaseId], mouse$diseaseId)])
hclustTraits$labels[hclustIds == 'MONDO_custom'] = 'cilia-associated retinitis pigmentosa'

ciliaModuleEnrich$diseaseName = hclustTraits$labels[match(ciliaModuleEnrich$diseaseId, hclustIds)]

ciliaModuleEnrichWide = spread(ciliaModuleEnrich[,c("diseaseName", "complex", "pagerankClust", "enrichP")], key = complex, value = enrichP)
ciliaModuleEnrichWide = column_to_rownames(ciliaModuleEnrichWide, var = 'diseaseName')
ciliaModuleEnrichWide[is.na(ciliaModuleEnrichWide)] = 0

pagerankClustHeatmap = ciliaModuleEnrichWide$pagerankClust
ciliaModuleEnrichWide = ciliaModuleEnrichWide[,-1]

ciliaModuleEnrichWide[ciliaModuleEnrichWide > 0 & ciliaModuleEnrichWide < 1.3] = 0.5
ciliaModuleEnrichWide[ciliaModuleEnrichWide >= 1.3] = 1

ciliaModuleEnrichWide = ciliaModuleEnrichWide[match(hclustTraits$labels, rownames(ciliaModuleEnrichWide)),]

ComplexHeatmap::Heatmap(ciliaModuleEnrichWide, split = max(clustTraits), col =structure(c('white', 'grey', 'red'), names = c("0", "0.5", "1")), gap = unit(0.1, 'cm'), cluster_rows = hclustTraits, row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6))
