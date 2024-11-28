# prepare files for network propagation

# Set working directory ----

setwd('W:/GROUP/Users/Ellen/NetworkPropagation/')

# Libraries ----

library(dplyr)
library(sparklyr)
library(sparklyr.nested)

# Open Targets interactions ----

interactionPath <- 'W:/GROUP/Users/Ellen/NetworkPropagation/Datasets/interaction/'

## establish connection
sc <- spark_connect(master = "local", log = "console",
                    config = list(sparklyr.verbose = TRUE))

## read interaction dataset
interaction <- spark_read_parquet(sc,
                                  path = interactionPath)

## define necessary columns
columns <- interaction %>%
  sdf_schema() %>%
  lapply(function(x) do.call(tibble, x)) %>% 
  rbind()

## read dataset and select columns (loop for memory issues)
fileTable = cbind("",list.files(path = interactionPath, pattern= "snappy"))
colnames(fileTable)=c("file_name","parquet")
fileTable[,"file_name"]=paste0(interactionPath,
                               "part",
                               c(1:nrow(fileTable)),
                               ".csv",sep="")
for (i in 1:nrow(fileTable)){
  print(i)
  SetPath <- fileTable[i,"parquet"]
  Set <- spark_read_parquet(sc, path = paste0(interactionPath, SetPath), overwrite = TRUE)
  df <- Set %>%
    select(sourceDatabase,
           targetA,
           targetB,
           scoring) %>%
    collect()
  write.csv(df,fileTable[i,"file_name"],row.names=F)
}

spark_disconnect(sc)

## create dataframe with all interactions
fileList = list.files(path = interactionPath, pattern = "part.*\\.csv")

intAll = as.data.frame(matrix(ncol = 4))
colnames(intAll) = colnames(read.csv(paste0(interactionPath, fileList[[1]])))
for (i in 1:length(fileList)){
  print(i)
  df <- read.csv(paste0(interactionPath, fileList[[i]]))
  intAll <- rbind(intAll,df)
}

intAll = intAll[!is.na(intAll$targetA) & !is.na(intAll$targetB),] #filter targets wo ID
intAll = intAll[lapply(intAll$targetA, FUN = nchar) == 15,]
intAll = intAll[lapply(intAll$targetB, FUN = nchar) == 15,]

write.csv(intAll, './Datasets/interaction/interactionAll.csv')

# Common variants - Open Targets Genetics Portal ----

variantPath <- 'W:/GROUP/Users/Ellen/NetworkPropagation/Datasets/sourceId=ot_genetics_portal/'

## establish connection
sc <- spark_connect(master = "local", log = "console",
                    config = list(sparklyr.verbose = TRUE))

## read variant dataset
variant <- spark_read_parquet(sc,
                              path = variantPath)

## define necessary columns
columns <- variant %>%
  sdf_schema() %>%
  lapply(function(x) do.call(tibble, x)) %>% 
  rbind()

## read dataset and select columns (loop for memory issues)
fileTable = cbind("",list.files(path = variantPath, pattern= "snappy"))
colnames(fileTable)=c("file_name","parquet")
fileTable[,"file_name"]=paste(variantPath,
                              "part",
                              c(1:nrow(fileTable)),
                              ".csv",sep="")
for (i in 1:nrow(fileTable)){
  print(i)
  SetPath <- fileTable[i,"parquet"]
  Set <- spark_read_parquet(sc, path = paste0(variantPath, SetPath), overwrite = TRUE)
  df <- Set %>%
    select(targetId,
           datatypeId,
           diseaseFromSource,
           diseaseFromSourceMappedId) %>%
    collect()
  write.csv(df,fileTable[i,"file_name"],row.names=F)
}

spark_disconnect(sc)

## create dataframe with all variants
fileList = list.files(path = variantPath, pattern = "part.*\\.csv")

variantAll = as.data.frame(matrix(ncol = 4))
colnames(variantAll) = colnames(read.csv(paste0(variantPath, fileList[[1]])))
for (i in 1:length(fileList)){
  print(i)
  df <- read.csv(paste0(variantPath, fileList[[i]]))
  variantAll <- rbind(variantAll,df)
}

variantAll = variantAll[-1,]

write.csv(variantAll, './Datasets/sourceId=ot_genetics_portal/variantAll.csv')

# Open Targets ChEMBL ----

chemblPath = 'W:/GROUP/Users/Ellen/NetworkPropagation/Datasets/chembl/'

## read chembl dataset
chembl <- spark_read_parquet(sc,
                             path = chemblPath)

## define necessary columns
columns <- chembl %>%
  sdf_schema() %>%
  lapply(function(x) do.call(tibble, x)) %>% 
  rbind()

## read dataset and select columns (loop for memory issues)
fileTable = cbind("",list.files(path = chemblPath, pattern= "snappy"))
colnames(fileTable)=c("file_name","parquet")
fileTable[,"file_name"]=paste0(chemblPath,
                               "part",
                               c(1:nrow(fileTable)),
                               ".csv",sep="")
for (i in 1:nrow(fileTable)){
  print(i)
  SetPath <- fileTable[i,"parquet"]
  Set <- spark_read_parquet(sc, path = paste0(chemblPath, SetPath), overwrite = TRUE)
  df <- Set %>%
    select(targetId,
           clinicalStatus,
           clinicalPhase,
           diseaseFromSource,
           diseaseFromSourceMappedId) %>%
    collect()
  write.csv(df,fileTable[i,"file_name"],row.names=F)
}

spark_disconnect(sc)

## create dataframe with all chembl targets
fileList = list.files(path = chemblPath, pattern= "part.*\\.csv")

chemblAll = as.data.frame(matrix(ncol = 5))
colnames(chemblAll) = colnames(read.csv(paste0(chemblPath, fileList[[1]])))
for (i in 1:length(fileList)){
  print(i)
  df <- read.csv(paste0(chemblPath, fileList[[i]]))
  chemblAll <- rbind(chemblAll,df)
}

chemblAll = chemblAll[-1,]

write.csv(chemblAll, './Datasets/chembl/chemblAll.csv')

# Open Targets Associations by datatype - direct interactions ----

variantPath <- 'W:/GROUP/Users/Ellen/NetworkPropagation/Datasets/associationByDatatypeDirect/'

## establish connection
sc <- spark_connect(master = "local", log = "console",
                    config = list(sparklyr.verbose = TRUE))

## read variant dataset
variant <- spark_read_parquet(sc,
                              path = variantPath)

## define necessary columns
columns <- variant %>%
  sdf_schema() %>%
  lapply(function(x) do.call(tibble, x)) %>% 
  rbind()

## read dataset and select columns (loop for memory issues)
fileTable = cbind("",list.files(path = variantPath, pattern= "snappy"))
colnames(fileTable)=c("file_name","parquet")
fileTable[,"file_name"]=paste(variantPath,
                              "part",
                              c(1:nrow(fileTable)),
                              ".csv",sep="")
for (i in 1:nrow(fileTable)){
  print(i)
  SetPath <- fileTable[i,"parquet"]
  Set <- spark_read_parquet(sc, path = paste0(variantPath, SetPath), overwrite = TRUE)
  df <- Set %>%
    select(diseaseId,
           targetId,
           datatypeId,
           score,
           evidenceCount) %>%
    collect()
  write.csv(df,fileTable[i,"file_name"],row.names=F)
}

spark_disconnect(sc)

## create dataframe with all variants
fileList = list.files(path = variantPath, pattern = "part.*\\.csv")

variantAll = as.data.frame(matrix(ncol = 5))
colnames(variantAll) = colnames(read.csv(paste0(variantPath, fileList[[1]])))
for (i in 1:length(fileList)){
  print(i)
  df <- read.csv(paste0(variantPath, fileList[[i]]))
  variantAll <- rbind(variantAll,df)
}

variantAll = variantAll[-1,]

write.csv(variantAll, './Datasets/associationByDatatypeDirect/associationAll.csv')

# Open Targets Associations by data source - direct interactions ----

variantPath <- 'W:/GROUP/Users/Ellen/NetworkPropagation/Datasets/associationByDatasourceDirect/'

## establish connection
sc <- spark_connect(master = "local", log = "console",
                    config = list(sparklyr.verbose = TRUE))

## read variant dataset
variant <- spark_read_parquet(sc,
                              path = variantPath)

## define necessary columns
columns <- variant %>%
  sdf_schema() %>%
  lapply(function(x) do.call(tibble, x)) %>% 
  rbind()

## read dataset and select columns (loop for memory issues)
fileTable = cbind("",list.files(path = variantPath, pattern= "snappy"))
colnames(fileTable)=c("file_name","parquet")
fileTable[,"file_name"]=paste(variantPath,
                              "part",
                              c(1:nrow(fileTable)),
                              ".csv",sep="")
for (i in 1:nrow(fileTable)){
  print(i)
  SetPath <- fileTable[i,"parquet"]
  Set <- spark_read_parquet(sc, path = paste0(variantPath, SetPath), overwrite = TRUE)
  df <- Set %>%
    select(diseaseId,
           targetId,
           datatypeId,
           datasourceId,
           score,
           evidenceCount) %>%
    collect()
  write.csv(df,fileTable[i,"file_name"],row.names=F)
}

spark_disconnect(sc)

## create dataframe with all variants
fileList = list.files(path = variantPath, pattern = "part.*\\.csv")

variantAll = as.data.frame(matrix(ncol = 6))
colnames(variantAll) = colnames(read.csv(paste0(variantPath, fileList[[1]])))
for (i in 1:length(fileList)){
  print(i)
  df <- read.csv(paste0(variantPath, fileList[[i]]))
  variantAll <- rbind(variantAll,df)
}

variantAll = variantAll[-1,]

write.csv(variantAll, './Datasets/associationByDatasourceDirect/associationAll.csv')


# Open Targets Disease Phenotypes ----

variantPath <- 'W:/GROUP/Users/Ellen/NetworkPropagation/Datasets/diseaseToPhenotype/'

## establish connection
sc <- spark_connect(master = "local", log = "console",
                    config = list(sparklyr.verbose = TRUE))

## read variant dataset
variant <- spark_read_parquet(sc,
                              path = variantPath)

## define necessary columns
columns <- variant %>%
  sdf_schema() %>%
  lapply(function(x) do.call(tibble, x)) %>% 
  rbind()

## read dataset and select columns (loop for memory issues)
fileTable = cbind("",list.files(path = variantPath, pattern= "snappy"))
colnames(fileTable)=c("file_name","parquet")
fileTable[,"file_name"]=paste(variantPath,
                              "part",
                              c(1:nrow(fileTable)),
                              ".csv",sep="")
for (i in 1:nrow(fileTable)){
  print(i)
  SetPath <- fileTable[i,"parquet"]
  Set <- spark_read_parquet(sc, path = paste0(variantPath, SetPath), overwrite = TRUE)
  df <- Set %>%
    select(disease,
           phenotype) %>%
    collect()
  write.csv(df,fileTable[i,"file_name"],row.names=F)
}

spark_disconnect(sc)

## create dataframe with all variants
fileList = list.files(path = variantPath, pattern = "part.*\\.csv")

variantAll = as.data.frame(matrix(ncol = 2))
colnames(variantAll) = colnames(read.csv(paste0(variantPath, fileList[[1]])))
for (i in 1:length(fileList)){
  print(i)
  df <- read.csv(paste0(variantPath, fileList[[i]]))
  variantAll <- rbind(variantAll,df)
}

variantAll = variantAll[-1,]

write.csv(variantAll, './Datasets/diseaseToPhenotype/phenotypeAll.csv')

# Open Targets Mouse Phenotypes ----

variantPath <- 'W:/GROUP/Users/Ellen/NetworkPropagation/Datasets/mousePhenotypes/'

## establish connection
sc <- spark_connect(master = "local", log = "console",
                    config = list(sparklyr.verbose = TRUE))

## read variant dataset
variant <- spark_read_parquet(sc,
                              path = variantPath)

## define necessary columns
columns <- variant %>%
  sdf_schema() %>%
  lapply(function(x) do.call(tibble, x)) %>% 
  rbind()

## read dataset and select columns (loop for memory issues)
fileTable = cbind("",list.files(path = variantPath, pattern= "snappy"))
colnames(fileTable)=c("file_name","parquet")
fileTable[,"file_name"]=paste(variantPath,
                              "part",
                              c(1:nrow(fileTable)),
                              ".csv",sep="")
for (i in 1:nrow(fileTable)){
  print(i)
  SetPath <- fileTable[i,"parquet"]
  Set <- spark_read_parquet(sc, path = paste0(variantPath, SetPath), overwrite = TRUE)
  df <- Set %>%
    select(modelPhenotypeId,
           modelPhenotypeLabel,
           targetFromSourceId,
           targetInModel,
           targetInModelEnsemblId,
           targetInModelMgiId) %>%
    collect()
  write.csv(df,fileTable[i,"file_name"],row.names=F)
}

spark_disconnect(sc)

## create dataframe with all variants
fileList = list.files(path = variantPath, pattern = "part.*\\.csv")

variantAll = as.data.frame(matrix(ncol = 6))
colnames(variantAll) = colnames(read.csv(paste0(variantPath, fileList[[1]])))
for (i in 1:length(fileList)){
  print(i)
  df <- read.csv(paste0(variantPath, fileList[[i]]))
  variantAll <- rbind(variantAll,df)
}

variantAll = variantAll[-1,]

write.csv(variantAll, './Datasets/mousePhenotypes/phenotypeAll.csv')
