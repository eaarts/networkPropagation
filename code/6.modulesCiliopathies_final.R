# clusters related to ciliopathies

# Set working directory ----

setwd('W:/GROUP/Users/Ellen/NetworkPropagation/')

# Libraries ---- 

library(tidyverse)
library(ComplexHeatmap)
library(doParallel)
library(circlize)

'%notin%' = Negate('%in%')

# Load files ----

PPIClusters = read.csv('20230929_PPIFullNetworkClusters_allEdges_max20Nodes.csv') #clusters

pageRankScores = readRDS('20231206_pageRankScores_JBTSnoOverlap.rds') #calculated network propagation scores for ciliopathies

traitAnnotation = read.csv('20231206_traitOverview_JBTSnoOverlap.csv') #labels for ciliopathies

variantsCiliopathy = read.csv('20231206_variantsCiliopathy_JBTSnoOverlap.csv') #curated list of ciliopathy genes

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


Heatmap(
  wilcoxRes_variable,
  col = colorRamp2(c(0, 3), c("white", "red")),
  rect_gp = gpar(col = "white", lwd = 2),
  cluster_columns = hclustTraits,
  cell_fun = function(j, i, x, y, w, h, fill) {
    if (wilcoxRes_variable[i, j] > 3) {
      grid.text("***", x, y)
    } else if (wilcoxRes_variable[i, j] > 2) {
      grid.text("**", x, y) 
    } else if (wilcoxRes_variable[i, j] > 1.3) {
      grid.text("*", x, y)
    }
    
  }
)

Heatmap(
  wilcoxRes_variable,
  col = colorRamp2(c(0, 3), c("white", "red")),
  rect_gp = gpar(col = "white", lwd = 2),
  cluster_columns = hclustTraits,
  cell_fun = function(j, i, x, y, w, h, fill) {
    if (wilcoxRes_binary[i, j] == 1) {
      grid.text("*", x, y)
    }
  }
)


Heatmap(
  t(wilcoxRes_variable),
  col = colorRamp2(c(0, 3), c("white", "red")),
  rect_gp = gpar(col = "white", lwd = 2),
  cluster_rows = hclustTraits,
  cell_fun = function(j, i, x, y, w, h, fill) {
    if (wilcoxRes_variable[j, i] > 3) {
      grid.text("***", x, y)
    } else if (wilcoxRes_variable[j, i] > 2) {
      grid.text("**", x, y) 
    } else if (wilcoxRes_variable[j, i] > 1.3) {
      grid.text("*", x, y)
    }
    
  }
)

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

# enrich per phenotype ----

phenotype = data.frame('ciliopathy' = unique(traitAnnotation[traitAnnotation$ciliopathy == T, "trait_label"]), 'renal' = F, 'retina' = F, 'CNS' = F, 'polydactyly' = F, 'skeletal' = F, 'hearing' = F)
phenotype$renal[grep('kidney|meckel|renal|nephro|senior|orofacio|bardet', phenotype$ciliopathy, ignore.case = T)] = T
phenotype$retina[grep('ocular|bardet|senior|usher|leber|retinitis|cone', phenotype$ciliopathy, ignore.case = T)] = T
phenotype$CNS[grep('joubert|meckel|acrocallosal|hydrolethalus|orofacio|bardet', phenotype$ciliopathy, ignore.case = T)] = T
phenotype$polydactyly[grep('meckel|orofacio|acrocallosal|bardet|hydrolethalus|ellis|short rib', phenotype$ciliopathy, ignore.case = T)] = T
phenotype$skeletal[grep('Jeune|short rib|ellis|cranioe', phenotype$ciliopathy, ignore.case = T)] = T
phenotype$hearing[grep('usher', phenotype$ciliopathy, ignore.case = T)] = T

wilcoxRes_adjusted_long = wilcoxRes_binary %>% as.data.frame() %>% rownames_to_column(var = 'complex') %>% gather(key = 'ciliopathy', value = 'sign', -complex)
wilcoxRes_adjusted_long = merge(wilcoxRes_adjusted_long, phenotype, by = 'ciliopathy')
wilcoxRes_adjusted_long_sig = wilcoxRes_adjusted_long[wilcoxRes_adjusted_long$sign == T,]

enrichComplex = unique(wilcoxRes_adjusted_long[, c("complex"), drop = F])
for (j in colnames(phenotype)[-1]) {
  for (i in 1:nrow(enrichComplex)) {
    enrichComplex[i, j] = phyper(
      nrow(wilcoxRes_adjusted_long[wilcoxRes_adjusted_long$complex == enrichComplex$complex[i] &
                                     wilcoxRes_adjusted_long[, j] == T & wilcoxRes_adjusted_long$sign == T,]) - 1,
      sum(phenotype[,j] == T),
      sum(phenotype[,j] == F),
      sum(wilcoxRes_adjusted_long$complex == enrichComplex$complex[i] & wilcoxRes_adjusted_long$sign == T),
      lower.tail = F
    )
  }
}

#enrichComplex[,-1] = apply(enrichComplex[,-1], 2, function(x)p.adjust(x,method = 'BH'))
enrichComplex[,-1] = -log10(enrichComplex[,-1])

enrichComplex_sig = enrichComplex[apply(enrichComplex[,-1],1,max) > 1.3,]

enrichComplex_sig = gather(enrichComplex_sig, key = 'phenotype', value = 'pvalue', -complex)

ggplot(enrichComplex_sig,
       aes(
         x = phenotype,
         y = complex,
         size = pvalue,
         color = pvalue
       )) + geom_point() + theme_classic() + 
  scale_size_continuous(range = c(-1, 10), breaks = seq(0.5, 4, 0.5)) + 
  scale_color_gradient(low = 'white', high = 'red')






