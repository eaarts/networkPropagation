setwd('W:/GROUP/Users/Ellen/NetworkPropagation/')

variantsCiliopathy = read.csv('20231206_variantsCiliopathy_JBTSnoOverlap.csv') #curated list of ciliopathy genes

freqtable = as.data.frame(table(unique(variantsCiliopathy[,c("datasourceId", "diseaseId", "targetId")])[,1]))

ggplot(freqtable, aes(x=reorder(Var1, -Freq), y=Freq)) + geom_bar(stat = 'identity') + theme_classic()



# donutData$fraction = donutData$Freq/sum(donutData$Freq)
# donutData$ymax = cumsum(donutData$fraction)
# donutData$ymin = c(0, head(donutData$ymax, n=-1))
# donutData$labelPosition <- (donutData$ymax + donutData$ymin) / 2
# donutData$label <- paste0(donutData$Var1, "\n value: ", donutData$Freq)
# 
# ggplot(donutData, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Var1)) +
#   geom_rect() +
#   geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
#   scale_fill_manual(values = rcartocolor::carto_pal(9,"Safe")) +
#   coord_polar(theta="y") +
#   xlim(c(1, 4)) +
#   theme_void() +
#   theme(legend.position = "none")
