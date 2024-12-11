# look further into mouse phenotypes

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

traitAnnotation = read.csv('data/traitOverview.csv')

variantsWithHPOandMP = read.csv('data/variantsCiliopathyMP.csv')

pageRankScores = readRDS('data/pagerankScores.rds')

# define ciliopathies ----

diseasesCilia = traitAnnotation$Var1[traitAnnotation$ciliopathy == T & !is.na(traitAnnotation$ciliopathy)]

# define related mouse phenotypes per ciliopathy ----

pageRankScores_scaled = scale(t(pageRankScores))
distTraitsMatrix = as.matrix(dist(t(pageRankScores_scaled), method = 'euclidean')) 

write_rds(distTraitsMatrix, 'data/distTraits.rds')

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

# hierarchical tree instead of umap ----
MPOntology = ontologyIndex::get_ontology('data/MP_ontology.txt')
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

PPIClusters = read.csv('data/PPIFullNetworkClusters.csv') #clusters

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
