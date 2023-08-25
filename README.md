# IBD loci

# gut scRNA
```
wget --no-check-certificate https://cellgeni.cog.sanger.ac.uk/gutcellatlas/Full_obj_raw_counts_nosoupx_v2.h5ad
warning: 5.8G file size
```

# h5ad2txt
```
library(rhdf5)
library(Matrix)

h5<-H5Fopen("Full_obj_raw_counts_nosoupx.h5ad")
indptr<-h5$'X/indptr'
gdf<-data.frame(
     symbol=h5$'var/_index',
     id=h5$'var/gene_ids'
)
cdf<-data.frame(
      Diagnosis=as.character(factor(h5$'obs/Diagnosis',levels=0:(length(h5$'obs/__categories/Diagnosis')-1),labels=h5$'obs/__categories/Diagnosis'))
      Category=as.character(factor(h5$'obs/category',levels=0:(length(h5$'obs/__categories/category')-1),labels=h5$'obs/__categories/category'))
      Subtype=as.character(factor(h5$'obs/Integrated_05',levels=0:(length(h5$'obs/__categories/Integrated_05')-1),labels=h5$'obs/__categories/Integrated_05'))
)
cdf$sample<-as.character(factor(cdf$Diagnosis,levels=c('fetal','Healthy adult','Pediatric Crohn Disease','Pediatric healthy'),labels=c('HF','HA','PCD','HP')))
h5closeAll()

write.table(gdf,file='gene.txt',quote=F,row.names=F,col.names=F,sep=' ')
write.table(cdf,file='celltype.txt',quote=F,row.names=F,col.names=F,sep=' ')

# $ data   : num [1:760344941(1d)] 1 1 1 4 1 1 1 1 1 1 ...
# $ indices: int [1:760344941(1d)] 15 53 102 154 216 223 244 269 271 326 ...
# $ indptr : int [1:428470(1d)] 0 970 1672 2394 3141 4369 5545 6131 6906 7533 ...

cells<-diff(indptr)
marker<-seq(0,length(cells),by=1000)
index<-1
for(k in marker){
  start<-k+1
  end<-k+1000
  if(end>length(cells)){end<-length(cells)}
  rowN<-rep(1:(end-start+1),cells[start:end])
  num<-length(rowN)
  colN<-h5read("Full_obj_raw_counts_nosoupx.h5ad", "/X/indices", index=list(index:(index+num-1)))
  counts<-h5read("Full_obj_raw_counts_nosoupx.h5ad", "/X/data", index=list(index:(index+num-1)))
  index<-index+num
  df<-as.matrix(sparseMatrix(i = rowN, j = colN + 1,x = as.numeric(counts),dims=c(end-start+1,33538)))
  write.table(df,file="counts.txt",append=T,quote=F,row.names=F,col.names=F,sep=' ')
}
```

# subset data

```
cat celltype.txt|'{print $2,$3 >> $4".celltype.txt"}'
cat celltype.txt|cut -d ' ' -f 4|paste - counts.txt|awk -F '\t' '{print $2 >> $1".counts.txt"}'
```

# rda
```
library(EWCE)
library(ewceData)
library(sctransform)

df<-read.table('PCD.counts.txt')
df<-t(df)
names(df)<-paste0("I",seq(ncol(df)))
gf<-read.table('gene.txt')
row.names(df)<-gf[,1]

af<-read.table('PCD.celltype.txt')
annotLevels <- list(level1class = af[,2], level2class = af[,3])

ctd_file <- generate_celltype_data(
    exp=as.matrix(df),
    annotLevels=annotLevels,
    groupName='PCD',savePath='./'
)
```


# IBD loci
```
zcat dbsnp_138.b37.vcf.gz|awk '!/#/{print $1"_"$2,$3,$4,$5}' | sort -k1,1 -T ./tmp > pos2allele.txt
awk -F '\t' 'NR>1{print $1"_"$3,$2,$10,$11}' tableS2.txt|sort -k1,1 > ibd.pos
join ibd.pos pos2allele.txt|awk '{if($1!=s){print $1,$2,$6,$7,"D1="$3";D2="$4};s=$1}'|tr '_' ' '|awk '{print $1,$2,$3,$4,$5,". .",$6}'|tr ' ' '\t'  > ibd.vcf
vep --assembly GRCh37 --fork 4 -i ibd.vcf -o ibd.vep --vcf --no_stats --merged --force_overwrite --offline --use_given_ref --per_gene --symbol --canonical --protein --biotype --nearest symbol --fasta hg19.masked.fa.gz --dir_cache ./cache
```
# enrichment

```
source('./bin/bootstrap_enrichment_test.r')
source('./bin/cell_list_dist.r')
source('./bin/generate_controlled_bootstrap_geneset.r')
source('./bin/get_summed_proportions.r')

gdf<-read.table('data/ibd.gene')
x<-unique(gdf$V1)
load('./data/CellTypeData_PCD.rda')

bg<-attr(ctd[[2]]$specificity,'dimnames')[[1]]
hits<-x[x %in% bg]
set.seed(2023)
rdf<-bootstrap_enrichment_test(sct_data=ctd,
                               hits=hits,
                               bg=bg,
                               reps=10000,
                               annotLevel=2
                              )
rdf$results$celltype<-row.names(rdf$results)
rdf$results$FDR<-p.adjust(rdf$results$p,method='fdr')
write.table(rdf$results,file='enrichment.txt',sep='\t',quote=F,row.names=F,col.names=T)
```
# plot
```
df<-read.table('enrichment.txt',sep='\t',head=T)
df$sign<-NA
df$sign[df$FDR<0.05]<-'*'
df$sign[df$FDR<0.01]<-'**'
df$sign[df$FDR<0.005]<-'***'

df<-df[order(df$p),]
df<-df[order(df$V2),]
df$celltype<-factor(df$celltype,levels=unique(df$celltype),order=T)

Fig1A<-ggplot(df,aes(celltype,-log10(p)))+
  geom_histogram(stat='identity')+
  geom_text(aes(label=sign),hjust=0,vjust=0.75,size=5)+
  coord_flip()+
  xlab('')+ggtitle('A')+
  scale_y_continuous(expression(paste(-log[10],"P-value")),
                     limits = c(0,5),expand=c(0,0))+
  theme_classic()+
  theme(legend.position = c(0.85,0.9),
        plot.title = element_text(size = 15,colour="black"),
        axis.text = element_text(size=12,colour="black"),
        #axis.text.x = element_text(angle=-45,hjust=0),
        axis.title = element_text(size=15,colour="black"),
        strip.text = element_text(size=12,colour="black"))
```

cell counts
```
for i in HF HP HA PCD
do
cat $i.celltype.txt|cut -d ' ' -f 3|paste 0 

done
```
















