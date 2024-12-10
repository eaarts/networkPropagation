# clusters related to ciliopathies

# Libraries ---- 

library(tidyverse)
library(ComplexHeatmap)
library(doParallel)
library(circlize)
library(igraph)
library(ggraph)
library(clusterProfiler)

'%notin%' = Negate('%in%')

# Load files ----

PPIClusters = read.csv('data/PPIFullNetworkClusters.csv') #clusters

pageRankScores = readRDS('data/pagerankScores.rds') #calculated network propagation scores for ciliopathies

traitAnnotation = read.csv('data/traitOverview.csv') #labels for ciliopathies

variantsCiliopathy = read.csv('data/variantsCiliopathies.csv') #curated list of ciliopathy genes

# cilia page rank scores ----

pageRankDFCilia = pageRankScores[rownames(pageRankScores) %in% traitAnnotation[traitAnnotation$ciliopathy == T, "Var1"],]
pageRankDFCilia = t(pageRankDFCilia)

pageRankDFCilia_scaled = t(scale(pageRankDFCilia))
rownames(pageRankDFCilia_scaled) = traitAnnotation$trait_label[match(rownames(pageRankDFCilia_scaled), traitAnnotation$Var1)]
propagationDist = as.matrix(dist(pageRankDFCilia_scaled, method = 'euclidean')) 

hclustTraits = hclust(as.dist(propagationDist))

# wilcoxon test for enrichment ----

cores <- detectCores()
cl <- makeCluster(cores - 1)
registerDoParallel(cl)

wilcoxRes = foreach (i = 1:ncol(pageRankDFCilia), .combine = cbind) %dopar% {
  sapply(
    unique(PPIClusters$walktrapCluster),
    FUN = function(x) {
      wilcox.test(pageRankDFCilia[, i][rownames(pageRankDFCilia) %in% PPIClusters[PPIClusters$walktrapCluster == x, "gene"]],
                  pageRankDFCilia[, i], #[rownames(pageRankDFCilia) %notin% PPIClusters[PPIClusters$walktrapCluster == x, "gene"]]
                  alternative = 'greater')$p.value
    }
  )
  
}

stopCluster(cl)

rownames(wilcoxRes) = unique(PPIClusters$walktrapCluster)
colnames(wilcoxRes) = colnames(pageRankDFCilia)

wilcoxRes_adjusted = apply(wilcoxRes, 2, function(x){p.adjust(x, method = 'bonferroni')})
wilcoxRes_adjusted = -log10(wilcoxRes_adjusted)

wilcoxRes_adjusted = wilcoxRes_adjusted[apply(wilcoxRes_adjusted,1,function(x){max(x, na.rm=T)}) > 1.3,]

wilcoxRes_variable = wilcoxRes_adjusted[apply(wilcoxRes_adjusted,1,FUN = function(x){sd(x) / mean(x) * 100}) > 200,]

wilcoxRes_binary = wilcoxRes_variable
for (i in 1:nrow(wilcoxRes_binary)) {
  for (j in 1:ncol(wilcoxRes_binary)) {
    moduleGenes =  PPIClusters[PPIClusters$walktrapCluster == rownames(wilcoxRes_variable)[i], "gene"]
    diseaseGenes = variantsCiliopathy[variantsCiliopathy$diseaseId == colnames(wilcoxRes_variable)[j], "targetId"]
    if (wilcoxRes_variable[i, j] > 1.3 &
        sum(diseaseGenes %in% moduleGenes) > 0) {
      wilcoxRes_binary[i, j] = 1
    } else {
      wilcoxRes_binary[i, j] = 0
    }
  }
}

colnames(wilcoxRes_variable) = traitAnnotation$trait_label[match(colnames(wilcoxRes_variable), traitAnnotation$Var1)]

# only take modules with known variant
wilcoxRes_binary = wilcoxRes_adjusted
for (i in 1:nrow(wilcoxRes_binary)) {
  for (j in 1:ncol(wilcoxRes_binary)) {
    moduleGenes =  PPIClusters[PPIClusters$walktrapCluster == rownames(wilcoxRes_adjusted)[i], "gene"]
    diseaseGenes = variantsCiliopathy[variantsCiliopathy$diseaseId == colnames(wilcoxRes_adjusted)[j], "targetId"]
    if (wilcoxRes_adjusted[i, j] > 1.3 &
        sum(diseaseGenes %in% moduleGenes) > 0) {
      wilcoxRes_binary[i, j] = 1
    } else {
      wilcoxRes_binary[i, j] = 0
    }
  }
}

wilcoxRes_binary = wilcoxRes_binary[rowSums(wilcoxRes_binary) > 0,]

colnames(wilcoxRes_binary) = traitAnnotation$trait_label[match(colnames(wilcoxRes_binary), traitAnnotation$Var1)]

Heatmap(wilcoxRes_binary, cluster_columns = hclustTraits, col = c('white', 'black'), rect_gp = gpar(col = "white", lwd = 2))

# cluster protein modules using re-clustering pattern ----

clusteringDF = PPIClusters[,"walktrapCluster", drop =F]

clusteringDF$clust1 = 1
clusteringDF$clust2 = gsub(';.*','',clusteringDF$walktrapCluster)
clusteringDF$clust3 = sapply(strsplit(clusteringDF$walktrapCluster,';'), FUN = function(x) paste(x[1:2], collapse = ';'))
clusteringDF$clust4 = sapply(strsplit(clusteringDF$walktrapCluster,';'), FUN = function(x) paste(x[1:3], collapse = ';'))
clusteringDF$clust5 = sapply(strsplit(clusteringDF$walktrapCluster,';'), FUN = function(x) paste(x[1:4], collapse = ';'))
clusteringDF$clust6 = sapply(strsplit(clusteringDF$walktrapCluster,';'), FUN = function(x) paste(x[1:5], collapse = ';'))

clusteringDF_selected = clusteringDF[clusteringDF$walktrapCluster %in% rownames(wilcoxRes_binary)[rowSums(wilcoxRes_binary) > 0],]

d1 = clusteringDF_selected[,c("clust1", "clust2")]
d1 = unique(d1)
colnames(d1) = c('from', 'to')
d2 = clusteringDF_selected[,c("clust2", "clust3")]
d2 = unique(d2)
colnames(d2) = c('from', 'to')
d3 = clusteringDF_selected[,c("clust3", "clust4")]
d3 = unique(d3)
colnames(d3) = c('from', 'to')
d4 = clusteringDF_selected[,c("clust4", "clust5")]
d4 = unique(d4)
colnames(d4) = c('from', 'to')
d5 = clusteringDF_selected[,c("clust5", "clust6")]
d5 = unique(d5)
colnames(d5) = c('from', 'to')

edges <- rbind(d1, d2, d3, d4, d5)
mygraph <- graph_from_data_frame(edges)

V(mygraph)$nGene = table(unlist(clusteringDF[,2:7]))[match(V(mygraph)$name, names(table(unlist(clusteringDF[,2:7]))))]
V(mygraph)$ciliopathyCluster = V(mygraph)$name %in% rownames(wilcoxRes_binary)[rowSums(wilcoxRes_binary) > 0]

ggraph(mygraph, layout = 'dendrogram', circular = FALSE) + 
  geom_edge_link(arrow = arrow(length = unit(2, 'mm')), 
                 end_cap = circle(3, 'mm')) +
  geom_node_point(aes(size = nGene), color = '#93C1E8') +
  theme_void() + geom_node_text(aes(label = name, color = ciliopathyCluster), size = 3) + scale_size_continuous(breaks = seq(1000,16000,2000))

# annotate protein modules with GO terms ----

PPIClustersCiliopathies = list()
for (i in unique(rownames(wilcoxRes_binary)[rowSums(wilcoxRes_binary) > 0])){
  PPIClustersCiliopathies[[i]] = enrichGO(PPIClusters[PPIClusters$walktrapCluster == i,"gene"], OrgDb = 'org.Hs.eg.db', keyType = 'ENSEMBL', ont = "BP")
}