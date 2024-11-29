# trait variants

'%notin%' = Negate('%in%')
library(AnnotationDbi)
library(org.Hs.eg.db)

# load variants from Open Targets ----

variants = read.csv('Datasets/associationByDatasourceDirect/associationAll.csv', row.names = 1) #data from open targets (ftp.ebi.ac.uk/pub/databases/opentargets/platform/24.09/output/etl/json/associationByOverallDirect)
#variants = read.csv('../../../Datasets/associationByDatasourceDirect/associationAll.csv', row.names = 1) #data from open targets (ftp.ebi.ac.uk/pub/databases/opentargets/platform/24.09/output/etl/json/associationByOverallDirect)

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
variantCertain = rbind(variantCertain, variantCertainCilia[,1:6])

variantUncertain = variants[paste0(variants$diseaseId, variants$targetId) %notin% paste0(variantCertain$diseaseId, variantCertain$targetId),]
variantUncertainCilia = variantUncertain[variantUncertain$diseaseId %in% diseasesCilia,]

rm(variants)
gc()

#plot number of seed genes per ciliopathy ----
traitAnnotation = read.csv('./data/traitOverview.csv')

variantCertainCilia$diseaseName = traitAnnotation$trait_label[match(variantCertainCilia$diseaseId, traitAnnotation$Var1)]
par(mar = c(25,4,2,1))
seedNumber = table(unique(variantCertainCilia[,c("targetId", "diseaseName")])[,2])
seedNumber = seedNumber[seedNumber > 1]
seedNumber = seedNumber[order(seedNumber, decreasing = T)]
barplot(seedNumber, las =2, ylab = 'Number of seed genes')

# Add mouse phenotypes ----

mouse = read.csv('Datasets/mousePhenotypes/phenotypeAll.csv', row.names = 1) # data from open targets (ftp.ebi.ac.uk/pub/databases/opentargets/platform/24.09/output/etl/json/mousePhenotypes)
#mouse = read.csv('../../../Datasets/mousePhenotypes/phenotypeAll.csv')

mouse$modelPhenotypeLabel = gsub('$', '_MP',mouse$modelPhenotypeLabel)
mouse$modelPhenotypeId = gsub(':', '_', mouse$modelPhenotypeId)

colnames(mouse)[c(2,4)] = c("diseaseId", "targetId")

variantsWithHPOandMP = variantCertain[,c("diseaseId", "targetId")]
variantsWithHPOandMP = rbind(variantsWithHPOandMP,mouse[,c("diseaseId", "targetId")])
variantsWithHPOandMP = unique(variantsWithHPOandMP)

variantsWithHPOandMP = variantsWithHPOandMP[variantsWithHPOandMP$diseaseId %in% traitAnnotation$Var1,]

# Save files ----

write.csv(variantsWithHPOandMP, 'variantsCiliopathyMP.csv')
write.csv(variantUncertainCilia, 'uncertainGenesCiliopathies.csv')
