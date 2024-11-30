# ciliopathy IDs

# Libraries ---- 

library(tidyverse)

# Open Targets EFO annotation ----

EFOAnno = ontologyIndex::get_ontology('./Datasets/annotation/efo.obo')
#EFOAnno = ontologyIndex::get_ontology('../../../Datasets/annotation/efo.obo') #from https://www.ebi.ac.uk/efo/

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

write.csv(cilioEFO, 'data/EFOIdsCiliopathies.csv')
