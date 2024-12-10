# look further into mouse phenotypes

setwd('W:/GROUP/USERS/Ellen/NetworkPropagation/')

# load libraries ----

library(igraph)
library(tidyverse)
library(pROC)
library(foreach)
library(doParallel)
library(ComplexHeatmap)
library(clusterProfiler)
library(RColorBrewer)
library(dendextend)

source('Code/networkPropagation.R')

'%notin%' = Negate('%in%')

# load files ----

distTraits = readRDS('20231215_distTraits_scaledPropagationScores_updatedJBTS.rds')
distTraitsMatrix = as.matrix(distTraits)

traitAnnotation = read.csv('20231206_traitOverview_JBTSnoOverlap.csv')

variantsWithHPOandMP = read.csv('20240105_variantsWithHPOandMP_updatedJBTS.csv')

pageRankScores = readRDS('20231206_pageRankScores_JBTSnoOverlap.rds')

# define ciliopathies ----

diseasesCilia = traitAnnotation$Var1[traitAnnotation$ciliopathy == T]

# load network ----

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

# define related mouse phenotypes per ciliopathy ----

distTraitsMatrix = distTraitsMatrix[rownames(distTraitsMatrix) %notin% c('MP_0002169','MP_0003171','MP_0003175','MP_0003176'),] #exclude normal phenotypes
distTraitsMatrix = distTraitsMatrix[,colnames(distTraitsMatrix) %notin% c('MP_0002169','MP_0003171','MP_0003175','MP_0003176')] #exclude normal phenotypes

#only focus on mouse phenotypes with at least 10 genes associated to the phenotype
variantsWithHPOandMP_10 = variantsWithHPOandMP[variantsWithHPOandMP$diseaseId %in% names(table(variantsWithHPOandMP$diseaseId)[table(variantsWithHPOandMP$diseaseId) >= 10]),]

distTraitsMatrix = distTraitsMatrix[rownames(distTraitsMatrix) %in% c(unique(variantsWithHPOandMP_10[grep('MP_', variantsWithHPOandMP_10$diseaseId),"diseaseId"]),diseasesCilia),] 
distTraitsMatrix = distTraitsMatrix[,colnames(distTraitsMatrix) %in% c(unique(variantsWithHPOandMP_10[grep('MP_', variantsWithHPOandMP_10$diseaseId),"diseaseId"]),diseasesCilia)]

relatedDiseasesPerCiliopathy <- list()
for (i in 1:length(diseasesCilia)) {
  distCiliopathy <-
    distTraitsMatrix[grep('MP_', rownames(distTraitsMatrix)), grep(diseasesCilia[i], colnames(distTraitsMatrix))]
  distCiliopathy <-
    distCiliopathy[names(distCiliopathy) %notin% diseasesCilia]
  if (length(distCiliopathy != 0)) {
    relatedDiseasesPerCiliopathy[[i]] <-
      names(distCiliopathy[order(distCiliopathy)][1:20])
  }
}

relatedDiseasesPerCiliopathy = lapply(relatedDiseasesPerCiliopathy, unique)
names(relatedDiseasesPerCiliopathy) = diseasesCilia
relatedDiseasesPerCiliopathy = relatedDiseasesPerCiliopathy[sapply(relatedDiseasesPerCiliopathy, length) > 0]

# create visualization ----

relatedMouse = unique(unlist(as.data.frame(relatedDiseasesPerCiliopathy)))

#select propagation scores of mouse phenotypes and ciliopathies
pageRankScoresSelected = pageRankScores[rownames(pageRankScores) %in% c(diseasesCilia, relatedMouse),]
pageRankScoresScaled = scale(t(pageRankScoresSelected))

#run umap on scaled scores
# set.seed(12345)
# umapRes = umap::umap(t(pageRankScoresScaled))

#obtain mouse phenotype ancestors for umap coloring
# MPOntology = ontologyIndex::get_ontology('Datasets/mousePhenotypes/MP_ontology.txt')
# MPOntology$id = gsub(':','_', MPOntology$id)
# MPOntology$ancestors = lapply(MPOntology$ancestors, FUN = function(x){gsub(':','_', x)})
# 
# MPontologyAncestors = unname(gsub(':','_',MPOntology$id[sapply(MPOntology$ancestors, length) == 2])) #main MP categories only have two ancestors
# 
# relatedMouseAncestors = MPOntology$ancestors[match(relatedMouse, MPOntology$id)]
# relatedMouseAncestors = lapply(relatedMouseAncestors, function(x) x[x %in% MPontologyAncestors])
# 
# orderMP = factor( #manual order of ancestors in case a phenotype has multiple ancestors (ordered by importance of phenotype to ciliopathies)
#   levels = c(
#     'MP_0005391',
#     'MP_0005377',
#     'MP_0005394',
#     'MP_0005390',
#     'MP_0005367',
#     'MP_0005379',
#     'MP_0005382',
#     'MP_0005371',
#     'MP_0005370',
#     'MP_0003631',
#     'MP_0005369',
#     'MP_0005385',
#     'MP_0005381',
#     'MP_0005388',
#     'MP_0005376',
#     'MP_0005387',
#     'MP_0005389',
#     'MP_0005397',
#     'MP_0010771',
#     'MP_0005378',
#     'MP_0005380',
#     'MP_0005384'
#   )
# )
# 
# for (i in 1:length(relatedMouseAncestors)){
#   if (length(relatedMouseAncestors[[i]]) == 1){
#     next
#   } else {
#     ancestorOrder = match(relatedMouseAncestors[[i]], levels(orderMP))
#     relatedMouseAncestors[[i]] = relatedMouseAncestors[[i]][ancestorOrder == min(ancestorOrder)]
#   }
# }
# 
# relatedMouseAncestors = unlist(relatedMouseAncestors)
# names(relatedMouseAncestors)  = relatedMouse

# #create dataframe for umap plot
# umapResDF = as.data.frame(umapRes$layout)
# umapResDF$ancestor = unname(relatedMouseAncestors)[match(rownames(umapResDF), names(relatedMouseAncestors))]
# umapResDF$ancestorName = MPOntology$name[match(umapResDF$ancestor, MPOntology$id)]
# umapResDF$traitLabel = traitAnnotation$trait_label[match(rownames(umapResDF), traitAnnotation$Var1)]
# umapResDF$ancestorName[umapResDF$ancestorName %in% names(table(umapResDF$ancestorName)[table(umapResDF$ancestorName) < 4])] = 'other'
# umapResDF$ancestorName[is.na(umapResDF$ancestorName)] = 'ciliopathy'
# umapResDF$ciliopathy = umapResDF$ancestorName == 'ciliopathy'
# 
# #define colors
# nb.cols <- length(unique(umapResDF$ancestorName))
# mycolors <- colorRampPalette(brewer.pal(8, "Accent"))(nb.cols)
# 
# #plot umap
# ggplot(umapResDF, aes(x=V1, y=V2, color = ancestorName, shape = ciliopathy)) + geom_point(size = 4) + 
#   ggrepel::geom_text_repel(data=subset(umapResDF,traitLabel %in% traitAnnotation$trait_label[traitAnnotation$ciliopathy == T]), 
#             aes(x=V1,y=V2,label=traitLabel), color = 'black') + theme_classic() + scale_color_manual(values = mycolors)
# 
# #try clustering
# distRelatedMP = dist(t(pageRankScoresScaled))
# 
# set.seed(1234)
# kmeansRes = kmeans(t(pageRankScoresScaled), centers = 10)
# umapResDF$cluster = kmeansRes$cluster[match(rownames(umapResDF), names(kmeansRes$cluster))]
# umapResDF$cluster = as.character(umapResDF$cluster)
# 
# nb.cols <- length(unique(umapResDF$cluster))
# mycolors <- colorRampPalette(brewer.pal(9, "Set1"))(nb.cols)
# 
# ggplot(umapResDF, aes(x=V1, y=V2, color = cluster, shape = ciliopathy)) + geom_point(size = 4) +
#   ggrepel::geom_text_repel(data=subset(umapResDF,traitLabel %in% traitAnnotation$trait_label[traitAnnotation$ciliopathy == T]), 
#                            aes(x=V1,y=V2,label=traitLabel), color = 'black') + theme_classic() + scale_color_manual(values = mycolors)

# hierarchical tree instead of umap ----
MPOntology = ontologyIndex::get_ontology('Datasets/mousePhenotypes/MP_ontology.txt')
MPOntology$id = gsub(':','_', MPOntology$id)
MPOntology$ancestors = lapply(MPOntology$ancestors, FUN = function(x){gsub(':','_', x)})

MPontologyAncestors = unname(gsub(':','_',MPOntology$id[sapply(MPOntology$ancestors, length) == 2])) #main MP categories only have two ancestors

relatedMouseAncestors = MPOntology$ancestors[match(relatedMouse, MPOntology$id)]
relatedMouseAncestors = lapply(relatedMouseAncestors, function(x) x[x %in% MPontologyAncestors])

relatedMouseAncestorsDF = stack(relatedMouseAncestors)
relatedMouseAncestorsDF$ind = gsub(':','_',relatedMouseAncestorsDF$ind)
colnames(relatedMouseAncestorsDF) = c('ancestor', 'MP')

distRelatedMP = dist(t(pageRankScoresScaled))

hclustRelatedMP = hclust(distRelatedMP, method = 'ward.D2')

hclustClusters = cutree(hclustRelatedMP, k = 11)

hclustDF = as.data.frame(hclustClusters)
hclustDF$traitLabel = traitAnnotation$trait_label[match(rownames(hclustDF), traitAnnotation$Var1)]

hclustAncestry = expand.grid(unique(hclustDF$hclustClusters), unique(relatedMouseAncestorsDF$ancestor))
for (i in 1:nrow(hclustAncestry)){
  hclustAncestry$fraction[i] = sum(relatedMouseAncestorsDF$MP %in% rownames(hclustDF)[hclustDF$hclustClusters == hclustAncestry$Var1[i]] & 
                                     relatedMouseAncestorsDF$ancestor == hclustAncestry$Var2[i]) / 
    sum(hclustDF$hclustClusters == hclustAncestry$Var1[i])
}

hclustAncestry$ancestryLabel = MPOntology$name[match(hclustAncestry$Var2, MPOntology$id)]

hclustRelatedMP$labels = traitAnnotation$trait_label[match(hclustRelatedMP$labels, traitAnnotation$Var1)]
hclustRelatedMPDend = as.dendrogram(hclustRelatedMP)
hclustRelatedMPDend = color_branches(hclustRelatedMPDend, k = 11)
plot(hclustRelatedMPDend)

for (i in 1:nrow(hclustDF)){
  orderAncestry = hclustAncestry[hclustAncestry$Var1 == hclustDF$hclustClusters[i], "Var2"][order(hclustAncestry[hclustAncestry$Var1 == hclustDF$hclustClusters[i], "fraction"], decreasing = T)]
  hclustDF$ancestor[i] = as.character(orderAncestry[min(match(relatedMouseAncestorsDF[relatedMouseAncestorsDF$MP == rownames(hclustDF)[i], "ancestor"], orderAncestry))])
}

hclustDF$ancestorName = MPOntology$name[match(hclustDF$ancestor, MPOntology$id)]
hclustDF$ancestorName[is.na(hclustDF$ancestorName)] = 'ciliopathy'

hclustDF = hclustDF[match(hclustRelatedMP$labels[hclustRelatedMP$order], hclustDF$traitLabel),]

hclustDF$ancestorName[hclustDF$ancestorName %in% names(table(hclustDF$ancestorName))[table(hclustDF$ancestorName) < 5]] = 'other'

nb.cols <- length(unique(hclustDF$ancestorName))
#mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
mycolors <- rcartocolor::carto_pal(12, "Safe")
mycolors = c(mycolors, '#000000')

Heatmap(hclustDF$ancestorName, row_split = factor(hclustDF$hclustClusters, levels = unique(hclustDF$hclustClusters)), col = mycolors)

hclustAncestry_wide = spread(hclustAncestry[,c("Var1", "Var2", "fraction")], key = Var2, value = fraction)
hclustAncestry_wide = hclustAncestry_wide[,apply(hclustAncestry_wide, 2, max) > 0.2]

hclustAncestry_wide = column_to_rownames(hclustAncestry_wide, var = 'Var1')

colnames(hclustAncestry_wide) = MPOntology$name[match(colnames(hclustAncestry_wide), MPOntology$id)]

Heatmap(hclustAncestry_wide, cluster_rows = F, col = c('white', 'red'))

# modules of related mouse phenotypes ----

PPIClusters = read.csv('20230929_PPIFullNetworkClusters_allEdges_max20Nodes.csv') #clusters

#wilcoxon test
cores <- detectCores()
cl <- makeCluster(cores - 1)
registerDoParallel(cl)

pageRankScoresSelectedT = t(pageRankScoresSelected)

wilcoxRes = foreach (i = 1:ncol(pageRankScoresSelectedT), .combine = cbind) %dopar% {
  sapply(
    unique(PPIClusters$walktrapCluster),
    FUN = function(x) {
      wilcox.test(pageRankScoresSelectedT[, i][rownames(pageRankScoresSelectedT) %in% PPIClusters[PPIClusters$walktrapCluster == x, "gene"]],
                  pageRankScoresSelectedT[, i], 
                  alternative = 'greater')$p.value
    }
  )
  
}

stopCluster(cl)

rownames(wilcoxRes) = unique(PPIClusters$walktrapCluster)
colnames(wilcoxRes) = colnames(pageRankScoresSelectedT)

wilcoxRes_adjusted = apply(wilcoxRes, 2, function(x){p.adjust(x, method = 'bonferroni')})
wilcoxRes_adjusted = -log10(wilcoxRes_adjusted)

wilcoxResDF = wilcoxRes_adjusted %>% as.data.frame(.) %>% rownames_to_column(., var = 'module') %>% gather(., key = 'trait', value = 'pvalue', -module)
wilcoxResDF = wilcoxResDF[wilcoxResDF$pvalue > 1.3,]

for (i in 1:nrow(wilcoxResDF)){
  wilcoxResDF$seed[i] = sum(PPIClusters[PPIClusters$walktrapCluster == wilcoxResDF[i,"module"], "gene"] %in% 
                            variantsWithHPOandMP[variantsWithHPOandMP$diseaseId == wilcoxResDF[i,"trait"], "targetId"])
}

wilcoxResDF = wilcoxResDF[wilcoxResDF$seed > 0,]

# #enrichment of modules in clusters?
# wilcoxResDF$cluster = umapResDF$cluster[match(wilcoxResDF$trait, rownames(umapResDF))]
# 
# modulePercluster = as.data.frame(table(wilcoxResDF[,c("module", "cluster")]))
# 
# for (i in unique(modulePercluster$cluster)){
#   modulePercluster[modulePercluster$cluster == i, 'freqByClusSize'] = modulePercluster[modulePercluster$cluster == i, "Freq"]/table(kmeansRes$cluster)[i]
# }
# 
# modulePresent = rownames(umapResDF) %in% wilcoxResDF[wilcoxResDF$module == '7;29;1;3;1', "trait"]
# 
# ggplot(umapResDF, aes(x=V1, y=V2, color = modulePresent, shape = ciliopathy)) + geom_point(size = 4) + 
#   ggrepel::geom_text_repel(data=subset(umapResDF,traitLabel %in% traitAnnotation$trait_label[traitAnnotation$ciliopathy == T]), 
#                            aes(x=V1,y=V2,label=traitLabel), color = 'black') + th eme_classic() + scale_color_manual(values = c('#3A53A4', '#BE202E'))

#for hclust
hclustModules = expand.grid(unique(hclustDF$hclustClusters), unique(wilcoxResDF$module))

for (i in 1:nrow(hclustModules)){
  hclustModules$fraction[i] = sum(wilcoxResDF$trait %in% rownames(hclustDF)[hclustDF$hclustClusters == hclustModules$Var1[i]] & 
                                    wilcoxResDF$module == hclustModules$Var2[i]) / 
    sum(hclustDF$hclustClusters == hclustModules$Var1[i])
}

for (i in 1:nrow(hclustDF)){
  orderModule = hclustModules[hclustModules$Var1 == hclustDF$hclustClusters[i], "Var2"][order(hclustModules[hclustModules$Var1 == hclustDF$hclustClusters[i], "fraction"], decreasing = T)]
  hclustDF$module[i] = as.character(orderModule[min(match(wilcoxResDF[wilcoxResDF$trait == rownames(hclustDF)[i], "module"], orderModule))])
}

hclustDF$module[hclustDF$module %in% names(table(hclustDF$module))[table(hclustDF$module) < 5]] = 'other'

nb.cols <- length(unique(hclustDF$module))
#mycolors <- colorRampPalette(brewer.pal(8, "Accent"))(nb.cols)
mycolors <- rcartocolor::carto_pal(nb.cols, "Bold")

Heatmap(hclustDF$module, split = factor(hclustDF$hclustClusters, levels = unique(hclustDF$hclustClusters)), col = mycolors)

hclustModules_wide = spread(hclustModules, key = Var2, value = fraction)
hclustModules_wide = hclustModules_wide[,apply(hclustModules_wide, 2, max) >= 0.5]

hclustModules_wide = column_to_rownames(hclustModules_wide, var = 'Var1')

Heatmap(hclustModules_wide, cluster_rows = F, col = c('white', 'red'))

# also run pca ----

# pcaRes = prcomp(t(pageRankScoresScaled))
# 
# pcaDF = as.data.frame(pcaRes$x)
# 
# pcaDF$ancestor = relatedMouseAncestorsDF$ancestor[match(rownames(pcaDF), relatedMouseAncestorsDF$MP)]
# pcaDF$ancestorName = MPOntology$name[match(pcaDF$ancestor, MPOntology$id)]
# pcaDF$traitLabel = traitAnnotation$trait_label[match(rownames(pcaDF), traitAnnotation$Var1)]
# pcaDF$ancestorName[pcaDF$ancestorName %in% names(table(pcaDF$ancestorName)[table(pcaDF$ancestorName) < 4])] = 'other'
# pcaDF$ancestorName[is.na(pcaDF$ancestorName)] = 'ciliopathy'
# pcaDF$ciliopathy = pcaDF$ancestorName == 'ciliopathy'
# 
# nb.cols <- length(unique(pcaDF$ancestorName))
# mycolors <- colorRampPalette(brewer.pal(8, "Accent"))(nb.cols)
# 
# ggplot(pcaDF, aes(x=PC1, y=PC2, color = ancestorName, shape = ciliopathy)) + geom_point(size = 4) + 
#   ggrepel::geom_text_repel(data=subset(pcaDF,traitLabel %in% traitAnnotation$trait_label[traitAnnotation$ciliopathy == T]), 
#                            aes(x=PC1,y=PC2,label=traitLabel), color = 'black') + theme_classic() + scale_color_manual(values = mycolors)
# 
# #can the genes that mainly split the phenotypes be identified?
# geneContrib = get_pca_var(pcaRes)$contrib
# indContrib = get_pca_ind(pcaRes)$contrib
# 
# geneMPcontrib = matrix(nrow = nrow(geneContrib), ncol = nrow(indContrib))
# for (i in 1:ncol(indContrib)){
#   geneMPcontrib[,i] = apply(geneContrib, 1, function(x) indContrib[i,]%*%x)
# }
# 
# rownames(geneMPcontrib) = rownames(geneContrib)
# colnames(geneMPcontrib) = rownames(indContrib)
# 
# #ciliopathyGenesHighContrib = unique(variantsCiliopathy[variantsCiliopathy$targetId %in% unique(names(which(geneMPcontrib > 50, arr.ind = T)[, 1])), "targetId"])
# ciliopathyGenesHighContrib = unique(names(which(geneMPcontrib > 300, arr.ind = T)[, 1]))
# 
# geneMPcontrib_ciliopathy = geneMPcontrib[rownames(geneMPcontrib) %in% ciliopathyGenesHighContrib,]
# 
# colnames(geneMPcontrib_ciliopathy) = traitAnnotation$trait_label[match(colnames(geneMPcontrib_ciliopathy), traitAnnotation$Var1)]
# rownames(geneMPcontrib_ciliopathy) = variantsCiliopathy$geneName[match(rownames(geneMPcontrib_ciliopathy), variantsCiliopathy$targetId)]
# 
# Heatmap(geneMPcontrib_ciliopathy)

#