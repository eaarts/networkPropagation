# network propagation for ciliopathies and mouse phenotypes

# Libraries ---- 

library(igraph)

source('code/networkPropagation.R')

'%notin%' = Negate('%in%')

# Interaction network full ----

#load open targets interaction network (IntAct, Reactome, SIGNOR, STRING)
#intAll <- read.csv('./Datasets/interaction/interactionAll.csv') #data from open targets (https://ftp.ebi.ac.uk/pub/databases/IntAct/various/ot_graphdb/2022-07-22/)
intAll <- read.csv('../../../Datasets/interaction/interactionAll.csv')

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

# Run pageRank ----
variantsWithHPOandMP <- read.csv('data/variantsCiliopathyMP.csv')

diseasesAll = unique(variantsWithHPOandMP$diseaseId)
diseasesAll = diseasesAll[grep('EFO|MONDO|MP', diseasesAll)] #select IDs from EFO, MONDO, MP

pageRankScoresAllTraits <- lapply(diseasesAll, FUN = function(x){getPageRank(x, intGraphClean, variantsWithHPOandMP)})
names(pageRankScoresAllTraits) <- diseasesAll
pageRankScoresAllTraits = pageRankScoresAllTraits[!is.na(pageRankScoresAllTraits)]

pageRankDFAllTraits <- data.frame(t(sapply(pageRankScoresAllTraits, FUN = function(x){c(x$pageRank)})))

colnames(pageRankDFAllTraits) = V(intGraphClean)$name

variantCertainCilia = read.csv('./data/variantsCiliopathies.csv') #includes same variants as data from open targets but with some manual curation
diseasesCilia = unique(variantCertainCilia$diseaseId)

cilioEFO = read.csv('data/EFOIdsCiliopathies.csv')

pageRankDFAllTraits = pageRankDFAllTraits[!(rownames(pageRankDFAllTraits) %in% cilioEFO$ID & rownames(pageRankDFAllTraits) %notin% diseasesCilia),] #exclude extra ciliopathies

saveRDS(pageRankDFAllTraits, 'data/pagerankScores.rds')
