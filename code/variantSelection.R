# trait variants

'%notin%' = Negate('%in%')
library(AnnotationDbi)
library(org.Hs.eg.db)

# load variants from Open Targets ----

variants = read.csv('Datasets/associationByDatasourceDirect/associationAll.csv', row.names = 1) #data from open targets (ftp.ebi.ac.uk/pub/databases/opentargets/platform/24.09/output/etl/json/associationByOverallDirect)

variantCertain = variants[(variants$datasourceId == 'uniprot_literature' & variants$score > 0.5) | 
                            (variants$datasourceId == 'gene2phenotype' & variants$score > 0.5) | 
                            (variants$datasourceId == 'clingen' & variants$score > 0.5) | 
                            (variants$datasourceId == 'uniprot_variants' & variants$score > 0.5) | 
                            (variants$datasourceId == 'genomics_england' & variants$score > 0.5) | 
                            (variants$datasourceId == 'orphanet' & variants$score > 0.5) | 
                            (variants$datasourceId == 'gene_burden' & variants$score > 0.5) | 
                            (variants$datasourceId == 'ot_genetics_portal' & variants$score > 0.5) | 
                            (variants$datasourceId == 'eva' & variants$score > 0.8),]

variantCertainCilia = read.csv('./data/variantsCiliopathies.csv')
diseasesCilia = unique(variantCertainCilia$diseaseId)
variantCertain = variantCertain[variantCertain$diseaseId %notin% diseasesCilia,]
variantCertain = variantCertain[,-1]
variantCertain = rbind(variantCertain, variantCertainCilia[,1:6])

variantUncertain = variants[variants$X %notin% variantCertain$X,]
variantUncertainCilia = variantUncertain[variantUncertain$diseaseId %in% diseasesCilia,]

rm(variants)
gc()

#plot number of seed genes per ciliopathy ----
variantCertainCilia$diseaseName = traitAnnotation$trait_label[match(variantCertainCilia$diseaseId, traitAnnotation$Var1)]
par(mar = c(25,4,2,1))
seedNumber = table(unique(variantCertainCilia[,c("targetId", "diseaseName")])[,2])
seedNumber = seedNumber[seedNumber > 1]
seedNumber = seedNumber[order(seedNumber, decreasing = T)]
barplot(seedNumber, las =2, ylab = 'Number of seed genes')

# Add human phenotype ontology ----

#HPO = read.table('Datasets/HPO/phenotype.hpoa', sep = '\t', comment.char = '', skip = 4, header = T)
HPOGenes = read.table('Datasets/HPO/phenotype_to_genes.txt', sep = '\t')
colnames(HPOGenes) = c('diseaseId', 'HPOLabel', 'EntrezId', 'GeneSymbol', 'SourceInfo', 'Source', 'DiseaseId')

HPOGenes$HPOLabel = gsub('$', '_HP',HPOGenes$HPOLabel)
HPOGenes$diseaseId = gsub(':', '_', HPOGenes$diseaseId)

HPOGenes$EntrezId = as.character(HPOGenes$EntrezId)
HPOGenes$targetId = mapIds(x = org.Hs.eg.db, keys = HPOGenes$EntrezId, column = 'ENSEMBL', keytype = 'ENTREZID', multiVals = 'first')

## other HPO file (more HPO terms included)
HPOGenes2 = read.table('Datasets/HPO/genes_to_phenotype.txt', sep = '\t', header = T)
HPOGenes2$hpo_id = gsub(':', '_', HPOGenes2$hpo_id)

HPOGenes2$ncbi_gene_id = as.character(HPOGenes2$ncbi_gene_id)
HPOGenes2$targetId = mapIds(x = org.Hs.eg.db, keys = HPOGenes2$ncbi_gene_id, column = 'ENSEMBL', keytype = 'ENTREZID', multiVals = 'first')

HPOGenes2 = HPOGenes2[,c("hpo_id", "targetId", "hpo_name")]
HPOGenes2$hpo_name = gsub('$', '_HP',HPOGenes2$hpo_name)
colnames(HPOGenes2)[1] = "diseaseId"

# Add mouse phenotypes ----

mouse = read.csv('Datasets/mousePhenotypes/phenotypeAll.csv')
mouse$modelPhenotypeLabel = gsub('$', '_MP',mouse$modelPhenotypeLabel)
mouse$modelPhenotypeId = gsub(':', '_', mouse$modelPhenotypeId)

colnames(mouse)[c(2,4)] = c("diseaseId", "targetId")

variantsWithHPOandMP = variantCertain[,c("diseaseId", "targetId")]
variantsWithHPOandMP = rbind(variantsWithHPOandMP, HPOGenes2[,c("diseaseId", "targetId")], mouse[,c("diseaseId", "targetId")],HPOGenes[,c("diseaseId", "targetId")])
variantsWithHPOandMP = unique(variantsWithHPOandMP)

variantsWithHPOandMP = variantsWithHPOandMP[variantsWithHPOandMP$diseaseId %in% traitAnnotation$Var1,]

# Save files ----

write.csv(variantsWithHPOandMP, '20240105_variantsWithHPOandMP_updatedJBTS.csv')
write.csv(variantUncertainCilia, '20240105_uncertainGenesCiliopathies.csv')
