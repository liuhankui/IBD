#single transcriptome matrix
scRNA<-'D:/tmp/tmp'
rdf<-read.table(scRNA,head=T,as.is=T,sep="\t",row.names = 1)
rdf<-rdf[rowSums(rdf,na.rm=T)>0,]
sdx<-apply(rdf,1,function(x){sd(x,na.rm=T)})
rn<-toupper(row.names(rdf))

geneList<-c('AARS','AARS2','AARSL','CARS','CARS2','DARS','DARS2','EARS2','EPRS','EPRS','FARS2','FARSA','FARSB','FARSL','FARSLA','FARSLB','GARS','HARS','HARS2','HARSL','IARS','IARS2','IFI53','KARS','LARS','LARS2','MARS','MARS2','NARS','NARS2','PARS2','QARS','QARS','QARS','QPRS','QPRS','RARS','RARS2','RARSL','SARS','SARS2','SARSM','TARS','TARS2','USH3B','VARS','VARS2','VARS2','VARS2L','VARSL','WARS','WARS2','YARS','YARS2')
geneList<-geneList[geneList %in% rn]

i<-10000
sdf<-data.frame(n=rep(length(geneList),i))
sdf$s<-apply(sdf,1,function(x){mean(sample(sdx,x[1]),na.rm=T)})

sd_value_distribution<-sd(sdf$s,na.rm=T)
mean_value_distribution<-mean(sdf$s,na.rm=T)
mean_value_geneList<-mean(sdx[rn %in% geneList],na.rm=T)
p_value<-1-pnorm(mean_value_geneList,mean=mean_value_distribution,sd=sd_value_distribution)
print(p_value)
