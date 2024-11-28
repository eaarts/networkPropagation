# PageRank scores ---- 
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

# Modules ----
# pageRankModules <-
#   function(disease = NULL, variantList = NULL, intGraph, geneVariant = NULL) {
#     if (is.null(variantList)){
#     variantList <-
#       unique(geneVariant[geneVariant$diseaseId == paste0(disease), 'targetId'])
#     } else
#       variantList = variantList
#     
#     if (length(variantList) < 2) {
#       highRank = NA
#       highRankBest = NA
#     } else {
#       V(intGraph)$targetGene = 0
#       V(intGraph)$targetGene[V(intGraph)$name %in% variantList] = 1
#       
#       pageRank = page_rank(intGraph, personalized = V(intGraph)$targetGene)
#       
#       pageRankRes = data.frame(
#         gene = V(intGraph)$name,
#         pageRank = pageRank$vector
#       )
#       
#       ##Network filter
#       
#       highRank = pageRankRes[as.numeric(pageRankRes[,"pageRank"]) > quantile(as.numeric(pageRankRes[,"pageRank"]))[4],]
#       
#       intGraphSub = subgraph(intGraph, V(intGraph)$name %in% highRank$gene)
#       
#       cwt=cluster_walktrap(intGraphSub, 
#                            #weights = E(intGraphSub)$weight, 
#                            steps = 6,
#                            merges = TRUE, 
#                            modularity = TRUE, 
#                            membership = TRUE)
#       
#       nodeDegree = igraph::degree(intGraphSub)
#       
#       highRank=cbind(highRank,nodeDegree,cwt$membership)
#       
#       colnames(highRank)[(ncol(highRank)-1):ncol(highRank)]=c("nodeDegree","walktrapCluster")
#       
#       clust=as.matrix(as.data.frame(table(cwt$membership)))
#       
#       ####Recluster
#       
#       if(sum(as.numeric(clust[,2])>=100)>0){
#         
#         #temp=rbind(c("0","0"),clust[as.numeric(clust[,2])>=100,])
#         temp=clust[as.numeric(clust[,2])>=100,,drop=F]
#         
#         for (i in 1:nrow(temp)){
#           
#           reClust=highRank[highRank[,"walktrapCluster"]==temp[i,1],]	
#           
#           intGraphSubRe = subgraph(intGraphSub, V(intGraphSub)$name %in% reClust$gene)
#           
#           cwtRe=cluster_walktrap(intGraphSubRe, 
#                                  #weights = E(intGraphSubRe)$weight, 
#                                  steps = 6,
#                                  merges = TRUE, 
#                                  modularity = TRUE, 
#                                  membership = TRUE)
#           
#           reClust[,"walktrapCluster"]=paste(reClust[,"walktrapCluster"],cwtRe$membership,sep=";")
#           #reClust[,"walktrapModularity"]=cwtRe$modularity
#           
#           highRank[highRank[,"walktrapCluster"]%in%temp[i,1],"walktrapCluster"]=reClust[,"walktrapCluster"]
#           #highRank[highRank[,"walktrapCluster"]%in%temp[i,1],"walktrapModularity"]=reClust[,"walktrapModularity"]
#           
#         }
#         
#       }
#       
#       ####Re-Recluster
#       
#       clust=as.matrix(as.data.frame(table(highRank[,"walktrapCluster"])))
#       
#       if(sum(as.numeric(clust[,2])>=100)>0){
#         
#         temp=clust[as.numeric(clust[,2])>=100,,drop=F]
#         
#         for (i in 1:nrow(temp)){
#           
#           reClust=highRank[highRank[,"walktrapCluster"]==temp[i,1],]	
#           
#           intGraphSubRe = subgraph(intGraphSub, V(intGraphSub)$name %in% reClust$gene)
#           
#           cwtRe=cluster_walktrap(intGraphSubRe, 
#                                  #weights = E(intGraphSubRe)$weight, 
#                                  steps = 6,
#                                  merges = TRUE, 
#                                  modularity = TRUE, 
#                                  membership = TRUE)
#           
#           reClust[,"walktrapCluster"]=paste(reClust[,"walktrapCluster"],cwtRe$membership,sep=";")
#           #reClust[,"walktrapModularity"]=cwtRe$modularity
#           
#           highRank[highRank[,"walktrapCluster"]%in%temp[i,1],"walktrapCluster"]=reClust[,"walktrapCluster"]
#           #highRank[highRank[,"walktrapCluster"]%in%temp[i,1],"walktrapModularity"]=reClust[,"walktrapModularity"]
#           
#         }
#         
#       }
#       
#       ####Re-Re-Recluster
#       
#       clust=as.matrix(as.data.frame(table(highRank[,"walktrapCluster"])))
#       
#       if(sum(as.numeric(clust[,2])>=100)>0){
#         
#         temp=clust[as.numeric(clust[,2])>=100,,drop=F]
#         
#         for (i in 1:nrow(temp)){
#           
#           reClust=highRank[highRank[,"walktrapCluster"]==temp[i,1],]	
#           
#           intGraphSubRe = subgraph(intGraphSub, V(intGraphSub)$name %in% reClust$gene)
#           
#           cwtRe=cluster_walktrap(intGraphSubRe, 
#                                  #weights = E(intGraphSubRe)$weight, 
#                                  steps = 6,
#                                  merges = TRUE, 
#                                  modularity = TRUE, 
#                                  membership = TRUE)
#           
#           reClust[,"walktrapCluster"]=paste(reClust[,"walktrapCluster"],cwtRe$membership,sep=";")
#           #reClust[,"walktrapModularity"]=cwtRe$modularity
#           
#           highRank[highRank[,"walktrapCluster"]%in%temp[i,1],"walktrapCluster"]=reClust[,"walktrapCluster"]
#           #highRank[highRank[,"walktrapCluster"]%in%temp[i,1],"walktrapModularity"]=reClust[,"walktrapModularity"]
#           
#         }
#         
#       }
#       
#       highRankBest = highRank[highRank$walktrapCluster %in% names(table(highRank$walktrapCluster)[table(highRank$walktrapCluster) > 9]),]
#       highRankBest = group_by(highRankBest, walktrapCluster) %>% filter(sum(gene %in% variantList) > 0)
#     }
#     return(list(highRank, highRankBest))
#   }
# 
# 




# Modules ----
pageRankModules <-
  function(disease = NULL, variantList = NULL, intGraph, geneVariant = NULL) {
    if (is.null(variantList)){
      variantList <-
        unique(geneVariant[geneVariant$diseaseId == paste0(disease), 'targetId'])
    } else
      variantList = variantList
    
    if (length(variantList) < 2) {
      highRank = NA
      highRankBest = NA
    } else {
      V(intGraph)$targetGene = 0
      V(intGraph)$targetGene[V(intGraph)$name %in% variantList] = 1
      
      pageRank = page_rank(intGraph, personalized = V(intGraph)$targetGene)
      
      pageRankRes = data.frame(
        gene = V(intGraph)$name,
        pageRank = pageRank$vector
      )
      
      ##Network filter
      
      highRank = pageRankRes[as.numeric(pageRankRes[,"pageRank"]) > quantile(as.numeric(pageRankRes[,"pageRank"]), 0.75),]
      # 
      intGraphSub = intHigh[intHigh$targetA %in% highRank$gene & intHigh$targetB %in% highRank$gene,]
      intGraphSub=intGraphSub[(intGraphSub$sourceDatabase == 'string' & intGraphSub$scoring >= 0.9) | (intGraphSub$sourceDatabase == 'intact' & intGraphSub$scoring >= 0.4),]
      intGraphSub = graph_from_data_frame(intGraphSub[,c("targetA", "targetB")], directed = F)
      intGraphSub  = igraph::simplify(intGraphSub , remove.loops = T, remove.multiple = T, edge.attr.comb = c(weight = 'max', 'ignore'))
# 
#       intGraphSub = subgraph(intGraph, V(intGraph)$name %in% highRank$gene)
      
      cwt=cluster_walktrap(intGraphSub, 
                           #weights = E(intGraphSub)$weight, MAYBE TRY WITH WEIGHTS
                           steps = 6,
                           merges = TRUE, 
                           modularity = TRUE, 
                           membership = TRUE)
      
      nodeDegree = igraph::degree(intGraphSub)
      
      highRank = highRank[match(V(intGraphSub)$name, rownames(highRank)), ]
      highRank=cbind(highRank,nodeDegree,cwt$membership)
      
      highRank = as.data.frame(highRank)
      
      colnames(highRank)[(ncol(highRank)-1):ncol(highRank)]=c("nodeDegree","walktrapCluster")
      #colnames(highRank)=c("gene", "nodeDegree","walktrapCluster")
      
      clust=as.matrix(as.data.frame(table(cwt$membership)))
      
      ####Recluster
      
      round = 1
      
      while(sum(as.numeric(clust[,2])>=100)>0 & round != 5){
        
        #temp=rbind(c("0","0"),clust[as.numeric(clust[,2])>=100,])
        temp=clust[as.numeric(clust[,2])>=100,,drop=F]
        
        for (i in 1:nrow(temp)){
          
          reClust=highRank[highRank[,"walktrapCluster"]==temp[i,1],]	
          
          intGraphSubRe = subgraph(intGraphSub, V(intGraphSub)$name %in% reClust$gene)
          
          cwtRe=cluster_walktrap(intGraphSubRe, 
                                 #weights = E(intGraphSubRe)$weight, 
                                 steps = 6,
                                 merges = TRUE, 
                                 modularity = TRUE, 
                                 membership = TRUE)
          
          reClust[,"walktrapCluster"]=paste(reClust[,"walktrapCluster"],cwtRe$membership,sep=";")
          #reClust[,"walktrapModularity"]=cwtRe$modularity
          
          highRank[highRank[,"walktrapCluster"]%in%temp[i,1],"walktrapCluster"]=reClust[,"walktrapCluster"]
          #highRank[highRank[,"walktrapCluster"]%in%temp[i,1],"walktrapModularity"]=reClust[,"walktrapModularity"]
          
        }
        
        clust=as.matrix(as.data.frame(table(highRank[,"walktrapCluster"])))
        round = round + 1
        
      }
      
      highRankBest = highRank[highRank$walktrapCluster %in% names(table(highRank$walktrapCluster)[table(highRank$walktrapCluster) > 4]),]
      highRankBest = group_by(highRankBest, walktrapCluster) %>% filter(sum(gene %in% variantList) > 0)
    }
    return(list(highRank, highRankBest))
  }

