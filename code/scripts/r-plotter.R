library(readxl)
library(ggplot2)

a<-read_excel(file.choose())
p<-ggplot(data=a, aes(x=Name, y=`Pubmed Count`)) + geom_bar(stat="identity", fill="steelblue")+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p
