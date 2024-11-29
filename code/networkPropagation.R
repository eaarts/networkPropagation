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