# ciliopathy IDs

# Set working directory ----

setwd('W:/GROUP/Users/Ellen/NetworkPropagation/')

# Libraries ---- 

library(tidyverse)

# Open Targets EFO annotation ----

EFOAnno = ontologyIndex::get_ontology('./Datasets/annotation/efo.obo')

## select ciliopathies

cilioEFO = EFOAnno$name[grep(
  'cone-rod|acrocallosal|alstrom|bardet|polycystic kidney|carpenter|coach|cranioectodermal|creveld|hydrolethalus|joubert|leber congenital|oculocerebrorenal|mckusick-kaufman|meckel syndrome|nphp|nephronophthisis|orofaciodigital|loken|short-rib|short rib|syndactyly-telecanthus|stromme|^usher|_usher|^retinitis pigmentosa \\d|^retinitis pigmentosa$|jeune|asphyxiating thoracic dystrophy',
  EFOAnno$name,
  ignore.case = T
)]

#remove cole-carpenter syndrome
cilioEFO = cilioEFO[-grep('cole', unname(cilioEFO), ignore.case = T)] 

cilioEFO = data.frame(ID = names(cilioEFO), description = unname(cilioEFO))
cilioEFO$ID = gsub(':', '_', cilioEFO$ID)

write.csv(cilioEFO, './EFOIdsCiliopathies.csv')

# DOID disease IDs ----

jensenText = read_tsv('Datasets/DISEASES_Jensen/human_disease_textmining_filtered.tsv', col_names = F)

cilioDOID = unique(jensenText[grep(
  'acrocallosal|alstrom|bardet|^cone-rod|polycystic kidney|carpenter|coach|cranioectodermal|creveld|hydrolethalus|joubert|leber congenital|oculocerebrorenal|mckusick-kaufman|meckel syndrome|nphp|nephronophthisis|orofaciodigital|loken|short-rib|syndactyly-telecanthus|stromme|^usher|_usher|retinitis pigmentosa|jeune|asphyxiating thoracic dystrophy',
  jensenText$X4,
  ignore.case = T
),3:4])

cilioDOID = cilioDOID[-grep('cole|miles', cilioDOID$X4, ignore.case = T),]
colnames(cilioDOID) = c('ID', 'description')
cilioDOID$ID = gsub(':', '_', cilioDOID$ID)

write.csv(cilioDOID, './DOIDIdsCiliopathies.csv')
