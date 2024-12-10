# cluster entire protein network

# Load libraries ----

library(igraph)

# Create PPI network ----

# load open targets interaction network (IntAct, Reactome, SIGNOR, STRING)
intAll <- read.csv('./Datasets/interaction/interactionAll.csv') #data from open targets (https://ftp.ebi.ac.uk/pub/databases/IntAct/various/ot_graphdb/2022-07-22/)

# set threshold for STRING
intString = grep('string', intAll$sourceDatabase)
lowString = intString[intAll$scoring[intString] < 0.4]
intHigh = intAll[-lowString,]

intHigh = intHigh[grep('ENSG0', intHigh$targetA),]
intHigh = intHigh[grep('ENSG0', intHigh$targetB),]
intHigh = intHigh[!duplicated(intHigh[,c("targetA", "targetB")]),]

# create graph
intGraph = graph_from_data_frame(intHigh[,c("targetA", "targetB")], directed = F)
intGraphClean = igraph::simplify(intGraph, remove.loops = T, remove.multiple = T, edge.attr.comb = c(weight = 'max', 'ignore'))

# remove big files to save memory
rm(intAll, intGraph, intString, lowString)
gc()

# Initial clustering of network ----
cwt=cluster_walktrap(intGraphClean, 
                     steps = 6,
                     merges = TRUE, 
                     modularity = TRUE, 
                     membership = TRUE)

clusterAll = data.frame('gene' = V(intGraphSub)$name)
clusterAll$walktrapCluster = cwt$membership

clust=as.matrix(as.data.frame(table(cwt$membership)))

# Recluster network to obtain smaller clusters ----

round = 1

while(sum(as.numeric(clust[,2])>=20)>0 & round != 5){ #recluster until every cluster has 20 genes or less, or until 5 rounds of reclustering have been performed
  
  temp=clust[as.numeric(clust[,2])>=20,,drop=F]
  
  for (i in 1:nrow(temp)){
    
    reClust=clusterAll[clusterAll[,"walktrapCluster"]==temp[i,1],]	
    
    intGraphSubRe = subgraph(intGraphSub, V(intGraphSub)$name %in% reClust$gene)
    
    cwtRe=cluster_walktrap(intGraphSubRe, 
                           #weights = E(intGraphSubRe)$weight, 
                           steps = 6,
                           merges = TRUE, 
                           modularity = TRUE, 
                           membership = TRUE)
    
    reClust[,"walktrapCluster"]=paste(reClust[,"walktrapCluster"],cwtRe$membership,sep=";")

    clusterAll[clusterAll[,"walktrapCluster"]%in%temp[i,1],"walktrapCluster"]=reClust[,"walktrapCluster"]

  }
  
  clust=as.matrix(as.data.frame(table(clusterAll[,"walktrapCluster"])))
  round = round + 1
  
}

write_rds(clusterAll, 'PPIFullNetworkClusters.csv')
