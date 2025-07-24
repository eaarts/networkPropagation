# prioritization of disease genes

# set paths to files ----
variant_file = 'test_genes.txt' # put file path for your seed genes here
# example file:
# ENSG00000169126
# ENSG00000185658
# ENSG00000167131
# ENSG00000105479
# ENSG00000198003
# ENSG00000157856 

# PageRank score calculation ---- 
getPageRank <-
  function(disease, intGraph, geneVariant) {
    variantList <-
      unique(geneVariant[geneVariant$diseaseId == paste0(disease), 'targetId'])
    
    if (length(variantList) < 2 || sum(V(intGraph)$name %in% variantList) < 2) {
      pageRankRes = NA
    } else {
      V(intGraph)$targetGene = 0
      V(intGraph)$targetGene[V(intGraph)$name %in% variantList] = 1
      
      pageRank = page_rank(intGraph, personalized = V(intGraph)$targetGene)
      
      pageRankRes = data.frame(
        gene = V(intGraph)$name,
        pageRank = pageRank$vector,
        target = V(intGraph)$targetGene
      )
      
    }
    return(pageRankRes)
  }

# load and/or install libraries ----
packages <- c("igraph", "stringr", "tidyverse", "dplyr", "AnnotationDbi", "org.Hs.eg.db")
installed <- packages %in% rownames(installed.packages())
if (any(!installed)) {
    install.packages(packages[!installed])
}
lapply(packages, library, character.only = TRUE)

'%notin%' = Negate('%in%')

# Interaction network full ----

#load open targets interaction network (IntAct, Reactome, SIGNOR, STRING)
intGraphClean = readRDS('intGraphClean.rds')

# Run pageRank ----
seed_genes = read.table(variant_file)

#strip from spaces
seed_genes = str_trim(seed_genes[,1])

V(intGraphClean)$targetGene = 0
V(intGraphClean)$targetGene[V(intGraphClean)$name %in% seed_genes] = 1

pageRank = page_rank(intGraphClean, personalized = V(intGraphClean)$targetGene)

pageRankRes = data.frame(
        gene = V(intGraphClean)$name,
        pageRank = pageRank$vector,
        target = V(intGraphClean)$targetGene)

pageRankRes$pageRank_scaled = scale(pageRankRes$pageRank)

# inlcude mouse phenotypes ----
variantsWithHPOandMP = read.csv('variantsCiliopathyMP.csv')

# only keep mouse phenotypes with at least 10 genes associated to the phenotype
variantsWithHPOandMP_10 = variantsWithHPOandMP[variantsWithHPOandMP$diseaseId %in% names(table(variantsWithHPOandMP$diseaseId)[table(variantsWithHPOandMP$diseaseId) >= 10]),]
variantsWithHPOandMP_10 = variantsWithHPOandMP_10[grep('MP_', variantsWithHPOandMP_10$diseaseId),]
variantsWithHPOandMP_10 = variantsWithHPOandMP_10[!variantsWithHPOandMP_10$diseaseId %in% c('MP_0002169','MP_0003171','MP_0003175','MP_0003176'),] # exclude normal phenotypes

mouse_phenotypes = unique(variantsWithHPOandMP_10$diseaseId)

# run network propagation for all mouse phenotypes
pageRankScoresAllTraits <- lapply(mouse_phenotypes, FUN = function(x){getPageRank(x, intGraphClean, variantsWithHPOandMP_10)})
names(pageRankScoresAllTraits) <- mouse_phenotypes
pageRankScoresAllTraits = pageRankScoresAllTraits[!is.na(pageRankScoresAllTraits)]

pageRankDFAllTraits <- data.frame(sapply(pageRankScoresAllTraits, FUN = function(x){c(x$pageRank)}))

rownames(pageRankDFAllTraits) = V(intGraphClean)$name

# scale scores
pageRankScores_scaled = scale(pageRankDFAllTraits)

# calculate distance to mouse phenotypes ----
pr_vec <- t(as.matrix(pageRankRes$pageRank_scaled))
pagerankscores_MP_mat <- t(as.matrix(pageRankScores_scaled))
distTraits = proxy::dist(pagerankscores_MP_mat, pr_vec, method = "Euclidean")

relatedTraits = rownames(distTraits)[order(distTraits)][1:10]

# mouse phenotype names
traitAnnotation = read.csv('traitOverview.csv')
traitAnnotation[traitAnnotation$Var1 %in% relatedTraits, 'trait_label']

distDF = data.frame('dist' = distTraits[,1], 'mouse_pheno' = rownames(distTraits))
distDF$mouse_pheno_label = traitAnnotation$trait_label[match(distDF$mouse_pheno, traitAnnotation$Var1)]

# rank genes based on related mouse phenotypes ----

pageRankScoresSelected = as.data.frame(pageRankScores_scaled[,colnames(pageRankScores_scaled) %in% relatedTraits])

pageRankScoresRank = apply(
    pageRankScoresSelected,
    2,
    FUN = function(x) {
        rank(dplyr::desc(x))
    }
)
  
pageRankScoresRank = as.data.frame(pageRankScoresRank)
pageRankScoresRank$totalRank = unname(apply(
    pageRankScoresRank,
    1,
    FUN = function(x) {
        10 ^ (sum(log10(x)) / ncol(pageRankScoresRank)) #to prevent inf values, use sum of log10 instead of product
    }
))

pageRankRes$rank = pageRankScoresRank$totalRank

# calculate p-value for network propagation scores ----

pageRankResRand = as.data.frame(matrix(nrow = 1000, ncol = length(V(intGraphClean))))
colnames(pageRankResRand) = V(intGraphClean)$name
set.seed(1234)
for (i in 1:1000) {
    V(intGraphClean)$targetGene = 0
    randomTargets = sample(length(V(intGraphClean)), sum(pageRankRes$target))
    V(intGraphClean)$targetGene[randomTargets] = 1

    pageRank = page_rank(intGraphClean, personalized = V(intGraphClean)$targetGene)

    pageRankResRand[i, ] = pageRank$vector
}

for (i in 1:nrow(pageRankRes)) {
    pageRankRes$pvalue[i] = sum(pageRankRes$pageRank[i] > pageRankResRand[, i]) / 1000
}

# calculate overall score based ciliopathy model ----
model = readRDS('candidateModel.rds')

pageRankRes[,c('pvalue', 'rank')] = apply(pageRankRes[,c('pvalue', 'rank')], 2, function(x) (x-min(x))/(max(x)-min(x)))

pageRankRes$overallScore = predict(model, newdata = pageRankRes)

# add gene names ----
pageRankRes$geneName = AnnotationDbi::mapIds(org.Hs.eg.db, keys = pageRankRes$gene, column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first')

# check if gene is mouse phenotype gene ----
for (i in 1:length(relatedTraits)){
    mouseGenes = variantsWithHPOandMP_10[variantsWithHPOandMP_10$diseaseId %in% relatedTraits[i], 'targetId']
    pageRankRes[, relatedTraits[i]] = FALSE
    pageRankRes[pageRankRes$gene %in% mouseGenes, relatedTraits[i]] = TRUE
}

colnames(pageRankRes)[9:18] = traitAnnotation$trait_label[match(colnames(pageRankRes)[9:18], traitAnnotation$Var1)]

# save results ----
write.csv(pageRankRes, 'gene_ranking.csv', row.names = FALSE)

write.csv(distDF, 'MP_dists.csv', row.names = FALSE)
