Name<-list.files("D:/pa_microbiome/wholeGenome_project/exiguobacterium/PathwaysGenome/sumGapSeq")
ID<-read.delim2("D:/pa_microbiome/wholeGenome_project/exiguobacterium/PathwaysGenome/sumGapSeq/GCA_000019905-all-Pathways.tbl",comment.char="#")[,1]
Results<-data.frame(ID=ID)
  
for(i in Name){
  path<-paste("D:/pa_microbiome/wholeGenome_project/exiguobacterium/PathwaysGenome/sumGapSeq/",i,sep="")
  Df<-read.delim2(path,comment.char="#")[,c(1,3)]
  Df[Df$Prediction=="false",2]<-0
  Df[Df$Prediction=="true",2]<-1
  colnames(Df)[2]<-i
  Results<-merge(Results,Df,by="ID",all=TRUE)
}

write.csv(Results,"D:/pa_microbiome/wholeGenome_project/exiguobacterium/PathwaysGenome/sumGapSeq/Results.csv")
