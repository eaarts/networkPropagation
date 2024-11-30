# plot data sources for ciliopathy genes

library(ggplot2)

variantsCiliopathy = read.csv('data/variantsCiliopathies.csv') #curated list of ciliopathy genes

freqtable = as.data.frame(table(unique(variantsCiliopathy[,c("datasourceId", "diseaseId", "targetId")])[,1]))

ggplot(freqtable, aes(x=reorder(Var1, -Freq), y=Freq)) + geom_bar(stat = 'identity') + theme_classic()
