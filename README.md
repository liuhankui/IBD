# IBD analysis script 

## Pre-stpep 1: IBD loci [shell]
```
zcat dbsnp_138.b37.vcf.gz|awk '!/#/{print $1"_"$2,$3,$4,$5}' | sort -k1,1 -T ./tmp > pos2allele.txt
awk -F '\t' 'NR>1{print $1"_"$3,$2,$10,$11}' tableS2.txt|sort -k1,1 > ibd.pos
join ibd.pos pos2allele.txt|awk '{if($1!=s){print $1,$2,$6,$7,"D1="$3";D2="$4};s=$1}'|tr '_' ' '|awk '{print $1,$2,$3,$4,$5,". .",$6}'|tr ' ' '\t'  > ibd.vcf
vep --assembly GRCh37 --fork 4 -i ibd.vcf -o ibd.vep --vcf --no_stats --merged --force_overwrite --offline --use_given_ref --per_gene --symbol --canonical --protein --biotype --nearest symbol --fasta hg19.masked.fa.gz --dir_cache ./cache
```

## Pre-stpep 2: get gut scRNA data [shell]
```
wget --no-check-certificate https://cellgeni.cog.sanger.ac.uk/gutcellatlas/Full_obj_raw_counts_nosoupx_v2.h5ad
warning: 5.8G file size
```

## Pre-stpep 3: get counts from h5 file format [R]
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
      Diagnosis=as.character(factor(h5$'obs/Diagnosis',levels=0:(length(h5$'obs/__categories/Diagnosis')-1),labels=h5$'obs/__categories/Diagnosis')),
      Category=as.character(factor(h5$'obs/category',levels=0:(length(h5$'obs/__categories/category')-1),labels=h5$'obs/__categories/category')),
      Subtype=as.character(factor(h5$'obs/Integrated_05',levels=0:(length(h5$'obs/__categories/Integrated_05')-1),labels=h5$'obs/__categories/Integrated_05'))
)
cdf$sample<-as.character(factor(cdf$Diagnosis,levels=c('fetal','Healthy adult','Pediatric Crohn Disease','Pediatric healthy'),labels=c('HF','HA','IBD','HP')))
h5closeAll()

write.table(gdf,file='gene.txt',quote=F,row.names=F,col.names=F,sep='\t')
write.table(cdf,file='celltype.txt',quote=F,row.names=F,col.names=F,sep='\t')

# $ data   : num [1:760344941(1d)] 1 1 1 4 1 1 1 1 1 1 ...
# $ indices: int [1:760344941(1d)] 15 53 102 154 216 223 244 269 271 326 ...
# $ indptr : int [1:428470(1d)] 0 970 1672 2394 3141 4369 5545 6131 6906 7533 ...

unlink('counts.txt')
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

## Pre-stpep 4: split data [shell]
```
#rm -f *.celltype.txt *.counts.txt *.Monocytes.txt
cat celltype.txt|awk -F '\t' '{print $2"\t"$3 >> $4".celltype.txt"}'
cat celltype.txt|cut -f 4|paste - counts.txt|awk -F '\t' '{print $2 >> $1".counts.txt"}'

for i in HF HP HA IBD
do
cat $i.celltype.txt|paste - $i.counts.txt|awk -v i=$i -F '\t' '$2=="Monocytes"{print $3 >> i".Monocytes.txt"}'
for j in Monocytes ILC3 Th1
do
echo -n "$i $j "
cat $i.celltype.txt|awk -F '\t' -v c=$j -v a=0  -v b=0 '{if($2==c){a++};b++}END{print a,b}'
done
done > cell.counts.txt
```

## Pre-stpep 5: calculate cell-type exrpression specificty [R]
```
library(EWCE)

df<-read.table('IBD.counts.txt')
cdf<-read.table('IBD.celltype.txt',sep='\t')
annotLevels <- list(level1class = cdf[,2], level2class =cdf[,1])
df<-t(df[,-1])
#names(df)<-paste0("I",seq(ncol(df)))
gdf<-read.table('gene.txt')
row.names(df)<-gdf[,1]

ctd_file <- generate_celltype_data(
    exp=as.matrix(df),
    annotLevels=annotLevels,
    groupName='IBD',savePath='./'
)
```

## required R packages
```
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(plyr)
library(reshape2)
library(igraph)
library(scales)
library(Matrix)
library(huge)
library(RColorBrewer)
library(tidyverse)
library(dplyr)
library(GGally)
library(network)
library(ggnetwork)
library(geomtextpath)
library(intergraph)
library(dendextend)
library(polynom)
```

## cell-type expression enrichment

```
source('./bin/bootstrap_enrichment_test.r')
source('./bin/cell_list_dist.r')
source('./bin/generate_controlled_bootstrap_geneset.r')
source('./bin/get_summed_proportions.r')
load('./data/CellTypeData_IBD.rda')

gdf<-read.table('data/ibd.gene')
x<-unique(gdf$V1)
bg<-attr(ctd[[1]]$specificity,'dimnames')[[1]]
hits<-x[x %in% bg]
set.seed(2023)
rdf<-bootstrap_enrichment_test(sct_data=ctd,
                               hits=hits,
                               bg=bg,
                               reps=10000,
                               annotLevel=1
                              )
rdf$results$celltype<-row.names(rdf$results)
rdf$results$FDR<-p.adjust(rdf$results$p,method='fdr')
write.table(rdf$results,file='enrichment.txt',sep='\t',quote=F,row.names=F,col.names=T)
```



## Fig. 1A
```
df<-read.table('enrichment.txt',sep='\t',head=T)
df$sign<-NA
df$sign[df$FDR<0.05]<-'*'
df$sign[df$FDR<0.01]<-'**'
df$sign[df$FDR<0.005]<-'***'

df<-df[order(df$p),]
df<-df[order(df$V2),]
df$celltype<-factor(df$celltype,levels=unique(df$celltype),order=T)

ggplot(df,aes(celltype,-log10(p)))+
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


## Fig. 1B
```
df<-read.table('cell.counts.txt')
ggplot(df,aes(V1,V3/V4*100))+
  geom_histogram(stat='identity',aes(fill=V2),colour='black',width=0.8,position='dodge')+
  scale_fill_brewer('',palette = 'Set2')+
  theme_classic()+
  coord_flip()+
  xlab('')+
  ylab('Cells proportion (%)')+
  ggtitle('B')+
  theme(legend.position = c(0.7,0.4),
        plot.title = element_text(size = 15),
        axis.text = element_text(size=12,colour="black"),
        axis.title = element_text(size=15,colour="black"),
        strip.text = element_text(size=12,colour="black"))
```

## Fig. 1C
```
df<-read.table('IBD.Monocytes.txt')
ss<-colSums(df)>0
df<-df[,ss]
df<-log(df+1)
sdf<-read.table('gene.txt')
sdf<-sdf[ss,]
sdf$sd<-apply(df,2,function(x){1/sd(x/mean(x))})

IBD<-read.table('./data/ibd.gene')
IBD<-unique(IBD$V1)
HI<-c('TP73','DFFB','KCNAB2','CHD5','CAMTA1','PINK1','SAM68','KCNQ4','GLUT1','MYH','FOXE3','HUD','INK4C','NFIA','CCN1','ABCA4','WNT2B','ADAR','ATP1A2','MPZ','MYOC','HRPT2','LRH-1','IRF6','PROX1','TP53BP2','NLRP3','ID2','MYCN','GCKR','SPAST','MSH6','FSHR','SPR','PAX8','SMADIP1','RPRM','SCN1A','HOXD13','COL3A1','SLC40A1','SATB2','SUMO1','BMPR2','XRCC5','PAX3','STK25','CHL1','SRGAP3','VHL','GHRL','PPARG','SRG3','RASSF1A','TKT','MITF','FOXP1','ROBO1','DIRC2','ATP2C1','FOXL2','ATR','SI','TERC','SOX2','OPA1','TFRC','FGFR3','LETM1','SH3BP2','MSX1','RBPJ','PHOX2B','ENAM','MAPK10','PKD2','SNCA','RIEG','ANK2','MAD2L1','PLK4','FBXW7','TERT','SEMA5A','GDNF','FGF10','PIK3R1','APC','RAD50','SMAD5','EGR1','TCOF1','NPM1','NKX2-5','MSX2','NSD1','FOXC1','DSP','EEF1E1','TNXA','TNX','HMGA1','RUNX2','CD2AP','ELOVL4','NT5E','SIM1','COL10A1','PARK2','TWIST1','GLI3','GCK','FKBP6','ELN','LIMK1','RFC2','GTF3','GTF2I','NCF1','KRIT1','COL1A2','SHFM1','RELN','FOXP2','CAV1','ST7','BRAF','SHH','HLXB9','GATA4','NKX3-1','FGFR1','CHD7','CSN5','EYA1','TRPS1','DMRT1','DMRT2','MLLT3','ARF','CDKN2B','BAG1','PAX5','GCNT1','ROR2','PTCH1','NR5A1','LMX1B','ENG','TSC1','COL5A1','NOTCH1','EHMT1','KLF6','GATA3','ANX7','PTEN','PAX2','FGF8','BUB3','CDKN1C','NUP98','PAX6','WT1','EXT2','ALX4','FEN1','SF1','FGF3','FZD4','ATM','H2AX','FLI1','NFRKB','PHB2','ETV6','CDKN1B','COL2A1','KRT5','MYF6','IGF1','SERCA2','TBX5','TBX3','HNF1A','BRCA2','FKHR','RB1','ZIC2','LIG4','COCH','NPAS3','NKX2-1','PAX9','BMP4','GCH1','SIX6','RAD51B','BCL11B','SPRED1','BUBR1','DLL4','FBN1','ALDH1A2','TPM1','P450SCC','BLM','COUP-TFII','SOX8','TSC2','PKD1','CBP','SOCS1','PRM2','PRM1','ABCC6','ERAF','SALL1','CBFB','CTCF','WWOX','FOXF1','FOXC2','YWHAE','HIC1','LIS1','P53','PMP22','COPS3','RAI1','TOP3A','SHMT1','RNF135','NF1','SUZ12','MEL-18','KLHL10','STAT5B','STAT5A','BECN1','BRCA1','PGRN','MAPT','CSH1','POLG2','PRKAR1A','SOX9','NHERF1','FSCN2','DSG1','DSG2','TCF4','FECH','MC4R','GALR1','SALL3','LKB1','PNPLA6','RYR1','TGFB1','RPS19','DMPK','CRX','PRPF31','JAG1','PAX1','GDF5','HNF4A','SALL4','MC3R','RAE1','GNAS','EDN3','KCNQ2','SOX18','SLC5A3','RUNX1','DYRK1A','COL6A1','PRODH','DGCR2','HIRA','TBX1','COMT','RTN4R','PCQAP','LZTR1','INI1','MYH9','SOX10','FBLN1','PPARA','PROSAP2','SHOX','P2RY8','NLGN4X','TRAPPC2','RPS4X','CSF2RA')
LOFT<-c('ALDH1B1','ANKRD30A','ANKRD35','ANO5','AP1G2','APIP','APOBEC3A','APOBEC3B','ASB15','ASPSCR1','ATP10B','ATP11A','ATP12A','ATP2C2','BPHL','BPIFA3','BTN3A3','BTNL2','BTNL8','BTNL9','C10orf53','C11orf40','C1orf127','C2orf40','C3orf14','C4orf46','C4orf50','C6','CABYR','CAPN9','CARD6','CCDC121','CCDC13','CCDC60','CCDC66','CD180','CD36','CD96','CDH19','CDK11A','CDKL2','CELA1','CEP72','CES1','CES5A','CFHR1','CFHR2','CFHR3','CHD1L','CHIT1','CHPF2','CLCN1','CLYBL','CNKSR1','COL16A1','COL6A2','COL6A5','CPXM2','CROT','CRYGN','CRYZ','CSH1','CTSE','CYP2A6','CYP2C8','CYP2D6','CYP2F1','CYP3A5','CYP4B1','DCDC2B','DCHS2','DDX60','DEFB126','DHDH','DMBT1','DNAH7','DQX1','DUOX2','ECT2L','EFCAB13','EFCAB3','EFCAB5','EFCAB6','ENOSF1','ENPEP','EPPK1','EPX','ERAP1','ERV3-1','EXO5','FAM129C','FAM151A','FAM187B','FAM45A','FAM81B','FCGBP','FCGR2A','FCN3','FLG','FLG2','FMO2','FRG2B','FUT2','FUT6','GADL1','GBGT1','GBP3','GBP4','GCFC2','GH2','GJB4','GLB1L2','GMPR','GOLGA8S','GP6','GPATCH2L','GRIN3B','GRK7','GYPB','HELB','HK3','HLA-B','HLA-DPA1','HPSE','HRG','HRNR','IDI2','IFIH1','IFNK','IL17RC','IL3RA','IQCH','ITIH1','KIAA0753','KIAA1257','KIAA1586','KIR3DL1','KLK14','KLK3','KRT4','KRT77','KRT83','LMF2','LMO7','LPA','LRRC39','LRTM1','MANEA','MAP3K4','MAZ','MCF2L','MCOLN3','MFSD9','MGAM','MLANA','MMP10','MOGAT1','MOK','MOXD1','MS4A6A','MST1','MUC17','MUC6','MUTYH','MYBBP1A','MYH1','MYH13','MYH8','MYO1A','MYOC','MYOF','NAALAD2','NBPF14','NBPF15','NEIL1','NLRP13','NLRP9','NOP16','NUDT8','OARD1','OBSCN','OCEL1','OR8S1','PAPLN','PDE11A','PDIA2','PGPEP1L','PHRF1','PKD1L2','PKHD1L1','PLA2G2C','PLA2G4D','PLA2R1','PLEKHG7','PLIN4','PNLIPRP3','POLM','POTEH','PPEF2','PPL','PPP1R3A','PRAMEF2','PRB1','PRB2','PRB4','PSG1','PSG11','PSG4','PSG9','PTCHD3','PTGDR','PXDNL','PZP','RAI1','RERGL','RETSAT','RFPL1','RGPD4','RGS11','RHD','RNF32','ROPN1B','RP1L1','RPTN','RTKN2','RTP1','SAMD11','SEMG2','SERHL2','SERPINA10','SERPINA9','SERPINB3','SFI1','SIGLEC1','SIGLEC5','SLC17A9','SLC22A10','SLC22A14','SLC22A25','SLC26A10','SLC5A4','SLCO1B1','SLFN13','SPATA31A6','SPATA4','SPATC1','SPNS3','SULT1A2','SULT1C4','SYNM','SYTL2','TAF6','TCF3','TCHHL1','TEKT3','TGM4','THBS4','THEM5','TIGD6','TLR10','TLR5','TMC2','TMEM82','TMIE','TMPRSS7','TNN','TRIM22','TRIM45','TRIM48','TRIM59','TRMT10B','TRMT2A','TTC38','TTN','UGT2B10','UGT2B17','UGT2B28','UMODL1','UNC93A','UPB1','UPK3A','UPP2','USP45','USP6','VILL','VWA3B','VWA7','WDR27','WDR90','XIRP1','XRRA1','ZAN','ZNF223','ZNF229','ZNF257','ZNF30','ZNF343','ZNF396','ZNF417','ZNF486','ZNF528','ZNF544','ZNF587','ZNF599','ZNF611','ZNF790','ZNF83','ZNF831','ZNF844','ZNF846','ZNF860','ZNF878','ZNF92','ZRANB3')

gdf<-data.frame(gene=c(IBD,HI,LOFT),type=c(rep('IBD',length(IBD)),rep('HI',length(HI)),rep('LOFT',length(LOFT))))

udf<-merge(gdf,sdf,by.x='gene',by.y='V1')

ggplot(udf,aes(type,sd))+
  geom_boxplot(aes(fill=type),outlier.colour = NA)+
  xlab('Genes')+ggtitle('C')+
  scale_y_continuous('Strictness')+
  scale_fill_brewer('',palette = 'Set2')+
  theme_classic()+
  theme(legend.position = 'none',
        plot.title = element_text(size = 15),
        axis.text = element_text(size=12,colour="black"),
        axis.title = element_text(size=15,colour="black"),
        strip.text = element_text(size=12,colour="black"))+
  coord_cartesian(ylim=c(0,1.5))+
  geom_signif(comparisons =  list(c("IBD","HI"),
                                  c('HI','LoFT'),
                                  c('IBD','LoFT')),
              parse = T,y_position = c(0.8,1,1.2),tip_length=0.01,
              map_signif_level=function(w) {if(w<0.001){
                f=round(w/10^floor(log10(abs(w))),2)
                d=floor(log10(abs(w)))
                paste0(f,"*x*10^",d)}else{round(w,3)}},
              step_increase = 0)
```

## Fig. 1D
```
pdf<-c(1)
num_ibd<-length(IBD)
for(i in seq(500000)){
  pdf[i]<-mean(sample(x=sdf$sd,size=num_ibd))
}
pdf<-data.frame(mstrictness=pdf)

ibd_s<-mean(sdf$sd[sdf3$gene %in% IBD])
dd<-density(pdf$mstrictness)
mdf<-data.frame(mstrictness=dd$x,den=dd$y,source='Observed')
ddf<-data.frame(mstrictness=seq(0.15,0.35,length=1000))
ddf$den<-dnorm(ddf$mstrictness,mean=mean(pdf$mstrictness),sd=sd(pdf$mstrictness))
idf<-ddf[ddf$mstrictness>ibd_s,]
ddf$source<-'Fitted'
ddf<-rbind(ddf,mdf)

ggplot()+
  geom_textpath(data=ddf,aes(mstrictness,den,label=source,colour=source),
                linewidth=0.8,size=3,fontface=2,hjust=0.4,vjust=0.3)+
  geom_textvline(xintercept =ibd_s,colour=brewer.pal(3, 'Set2')[3],
                 linetype=2,size=3,linewidth=0.5,
                 label='# of IBD genes')+
  scale_colour_brewer(palette = "Set2")+
  xlab('Mean strictness')+ylab('Density')+ggtitle('D')+
  coord_cartesian(xlim=c(0.18,0.34))+
  theme_classic()+
  theme(legend.position = 'none',
        plot.title = element_text(size = 15),
        axis.text = element_text(size=12,colour="black"),
        axis.title = element_text(size=15,colour="black"),
        strip.text = element_text(size=12,colour="black"))

```

## Fig. 1E
```
zdf<-read.table('gene.zscore',head=T)

zdf<-merge(gdf,zdf[,c(1,4,6)],by='gene')

ggplot(zdf,aes(x=lof_z,group = type))+
  geom_textdensity(aes(colour=type,label=type),linewidth=0.8,size=3,fontface=2,hjust=0.25,vjust=0.3)+
  theme_classic()+
  scale_colour_brewer('Diseases',palette = "Set2")+
  coord_cartesian(xlim = c(-5,15))+
  xlab('LoF Z-score')+
  ylab('Density')+
  ggtitle('E')+
  theme(legend.position = 'none',
        plot.title = element_text(size = 15),
        axis.text = element_text(size=12,colour="black"),
        axis.title = element_text(size=15,colour="black"),
        strip.text = element_text(size=12,colour="black"))
```

## Fig. 2
```
gdf<-read.table('ibd.gene')
gene<-read.table('gene.txt')
gene<-gene$V1

#-----------network IBD----------------------
df<-read.table('IBD.Monocytes.txt')
names(df)<-gene
df<-df[,gene %in% gdf$V1]
df<-df[,colSums(df)>0]
df<-df[,colSums(df==0)/(nrow(df)-1)<0.9]
df<-log(df+1)
exp<-colSums(df)/nrow(df)
cor.matrix<- cor(df,use = "pairwise.complete.obs")
mbModel <- huge(as.matrix(df), method = "mb")
mbOptRIC = huge.select(mbModel,criterion="stars")
mbOptRICGraph = mbOptRIC$refit
net<-network(mbOptRICGraph)
gnet1<-as.undirected(set_vertex_attr(asIgraph(net),"name", value = names(df)))
network.vertex.names(net)<-names(df)

tmp<-as.matrix(sparseMatrix(i=mbOptRIC$refit@i+1,p=mbOptRIC$refit@p,x=mbOptRIC$refit@x,dims=c(ncol(df),ncol(df))),diag=T)
net %e% "weights"<-abs(cor.matrix[which(tmp==1 & lower.tri(tmp),arr.ind=TRUE)])
net %e% "type"<- as.character(factor(sign(cor.matrix[which(tmp==1 & lower.tri(tmp),arr.ind=TRUE)]),levels=c(-1,0,1),labels=c('N','P','P')))
gid<-names(df)
net %v% "size" = as.numeric(exp)
net %v% "colour" = as.character(factor(gid,levels=gdf$V1[gdf$V1 %in% gid],labels=gdf$V2[gdf$V1 %in% gid]))

ndf1<-ggnetwork(net)
ndf1$source<-'Paediatric IBD'

#-----------network HP----------------------
df<-read.table('HP.Monocytes.txt')
names(df)<-gene
df<-df[,gene %in% gdf$V1]
df<-df[,colSums(df)>0]
df<-df[,colSums(df==0)/(nrow(df)-1)<0.9]
df<-log(df+1)
exp<-colSums(df)/nrow(df)
cor.matrix<- cor(df,use = "pairwise.complete.obs")
mbModel <- huge(as.matrix(df), method = "mb")
mbOptRIC = huge.select(mbModel,criterion="stars")
mbOptRICGraph = mbOptRIC$refit
net<-network(mbOptRICGraph)
gnet2<-as.undirected(set_vertex_attr(asIgraph(net),"name", value = names(df)))
network.vertex.names(net)<-names(df)

tmp<-as.matrix(sparseMatrix(i=mbOptRIC$refit@i+1,p=mbOptRIC$refit@p,x=mbOptRIC$refit@x,dims=c(ncol(df),ncol(df))),diag=T)
net %e% "weights"<-abs(cor.matrix[which(tmp==1 & lower.tri(tmp),arr.ind=TRUE)])
net %e% "type"<- as.character(factor(sign(cor.matrix[which(tmp==1 & lower.tri(tmp),arr.ind=TRUE)]),levels=c(-1,0,1),labels=c('N','P','P')))
gid<-names(df)
net %v% "size" = as.numeric(exp)
net %v% "colour" = as.character(factor(gid,levels=gdf$V1[gdf$V1 %in% gid],labels=gdf$V2[gdf$V1 %in% gid]))

ndf2<-ggnetwork(net)
ndf2$source<-'Paediatric Healthy'

#-----------network HF----------------------
df<-read.table('HF.Monocytes.txt')
names(df)<-gene
df<-df[,gene %in% gdf$V1]
df<-df[,colSums(df)>0]
df<-df[,colSums(df==0)/(nrow(df)-1)<0.9]
df<-log(df+1)
exp<-colSums(df)/nrow(df)
cor.matrix<- cor(df,use = "pairwise.complete.obs")
mbModel <- huge(as.matrix(df), method = "mb")
mbOptRIC = huge.select(mbModel,criterion="stars")
mbOptRICGraph = mbOptRIC$refit
net<-network(mbOptRICGraph)
gnet3<-as.undirected(set_vertex_attr(asIgraph(net),"name", value = names(df)))
#net<-network(mbModel$path[[3]])
network.vertex.names(net)<-names(df)

tmp<-as.matrix(sparseMatrix(i=mbOptRIC$refit@i+1,p=mbOptRIC$refit@p,x=mbOptRIC$refit@x,dims=c(ncol(df),ncol(df))),diag=T)
net %e% "weights"<-abs(cor.matrix[which(tmp==1 & lower.tri(tmp),arr.ind=TRUE)])
net %e% "type"<- as.character(factor(sign(cor.matrix[which(tmp==1 & lower.tri(tmp),arr.ind=TRUE)]),levels=c(-1,0,1),labels=c('N','P','P')))
gid<-names(df)
net %v% "size" = as.numeric(exp)
net %v% "colour" = as.character(factor(gid,levels=gdf$V1[gdf$V1 %in% gid],labels=gdf$V2[gdf$V1 %in% gid]))

ndf3<-ggnetwork(net)
ndf3$source<-'Fetal Healthy'

#-----------network HA----------------------
df<-read.table('HA.Monocytes.txt')
names(df)<-gene
df<-df[,gene %in% gdf$V1]
df<-df[,colSums(df)>0]
df<-df[,colSums(df==0)/(nrow(df)-1)<0.9]
df<-log(df+1)
exp<-colSums(df)/nrow(df)
cor.matrix<- cor(df,use = "pairwise.complete.obs")
mbModel <- huge(as.matrix(df), method = "mb")
mbOptRIC = huge.select(mbModel,criterion="stars")
mbOptRICGraph = mbOptRIC$refit
net<-network(mbOptRICGraph)
gnet4<-as.undirected(set_vertex_attr(asIgraph(net),"name", value = names(df)))
network.vertex.names(net)<-names(df)

tmp<-as.matrix(sparseMatrix(i=mbOptRIC$refit@i+1,p=mbOptRIC$refit@p,x=mbOptRIC$refit@x,dims=c(ncol(df),ncol(df))),diag=T)
net %e% "weights"<-abs(cor.matrix[which(tmp==1 & lower.tri(tmp),arr.ind=TRUE)])
net %e% "type"<- as.character(factor(sign(cor.matrix[which(tmp==1 & lower.tri(tmp),arr.ind=TRUE)]),levels=c(-1,0,1),labels=c('N','P','P')))
gid<-names(df)
net %v% "size" = as.numeric(exp)
net %v% "colour" = as.character(factor(gid,levels=gdf$V1[gdf$V1 %in% gid],labels=gdf$V2[gdf$V1 %in% gid]))

ndf4<-ggnetwork(net)
ndf4$source<-'Adult Healthy'

#----network merge----------------
ndf<-rbind(ndf1,ndf2,ndf3,ndf4)
ndf$source<-factor(ndf$source,levels=unique(ndf$source),order=T)
ggplot(ndf,aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(aes(linewidth=weights),colour='gray') +
  geom_nodes(aes(colour=colour,size=size),alpha = 0.5) +
  geom_nodes(aes(colour=colour,size=size/5))+
  geom_nodetext_repel(aes(label=vertex.names),size=2,box.padding = 0.1)+
  scale_color_brewer('Diseases',palette = "Set2") +
  scale_linewidth_continuous('Correlations',range = c(0.1,1))+
  scale_size('Expression(log[C+1])')+
  facet_wrap(~source,scale='free')+
  theme_blank()
```


# Fig3
```
cl<-function(g){
  A<-as.matrix(get.adjacency(g))
  S<-A+t(A)
  deg<-igraph::degree(g,mode=c("total"))
  num<-diag(S %*% S %*% S)
  denom<-diag(A %*% A)
  denom<-2*(deg*(deg-1)-2*denom)
  cc<-mean(num[denom!=0]/denom[denom!=0])
  return(cc)
}

tmp1<-data.frame(V=length(V(gnet1)),E=length(E(gnet1)),D=graph.density(gnet1),T=transitivity(gnet1),C=cl(gnet1),S='IBD')
tmp2<-data.frame(V=length(V(gnet2)),E=length(E(gnet2)),D=graph.density(gnet2),T=transitivity(gnet2),C=cl(gnet2),S='PH')
tmp3<-data.frame(V=length(V(gnet3)),E=length(E(gnet3)),D=graph.density(gnet3),T=transitivity(gnet3),C=cl(gnet3),S='FH')
tmp4<-data.frame(V=length(V(gnet4)),E=length(E(gnet4)),D=graph.density(gnet4),T=transitivity(gnet4),C=cl(gnet4),S='AH')
tdf<-rbind(tmp1,tmp2,tmp3,tmp4)
tdf<-melt(tdf[,-4],id='S')

tdf$variable<-factor(tdf$variable,
                     levels=c('V','E','D','C'),
                     labels=c('Nodes','Edges','Density','Cluster'))

#fig3A
ggplot()+
  geom_histogram(data=tdf,aes(S,value,fill=variable),stat='identity')+
  scale_fill_brewer(palette = "Set2")+
  xlab('')+ylab('# summary')+ggtitle('A')+
  facet_wrap(~variable,scale='free',nrow=2)+
  theme_classic()+
  theme(legend.position = 'none',
        plot.title = element_text(size = 15),
        axis.text = element_text(size=12,colour="black"),
        axis.text.x = element_text(angle=-30,hjust=0),
        axis.title = element_text(size=15,colour="black"),
        strip.text = element_text(size=12,colour="black"))


tmp1<-data.frame(degree=igraph::degree(gnet1),center=betweenness(gnet1))
tmp1$gene<-row.names(tmp1)
tmp1$source<-'PCD'

tmp2<-data.frame(degree=igraph::degree(gnet2),center=betweenness(gnet2))
tmp2$gene<-row.names(tmp2)
tmp2$source<-'PH'

tmp3<-data.frame(degree=igraph::degree(gnet3),center=betweenness(gnet3))
tmp3$gene<-row.names(tmp3)
tmp3$source<-'FH'

tmp4<-data.frame(degree=igraph::degree(gnet4),center=betweenness(gnet4))
tmp4$gene<-row.names(tmp4)
tmp4$source<-'AH'

adf<-merge(tmp1[,c(1,3)],tmp2[,c(1,3)],by='gene')
bdf<-merge(tmp1[,c(1,3)],tmp4[,c(1,3)],by='gene')
wilcox.test(adf$degree.x,adf$degree.y,paired = T)$p.value
wilcox.test(bdf$degree.x,bdf$degree.y,paired = T)$p.value


ddf<-rbind(tmp1,tmp2,tmp3,tmp4)

#fig3B
  ggplot()+
  geom_histogram(data=ddf,aes(degree,fill=source),binwidth = 1)+
  scale_fill_brewer(palette = "Set2")+
  xlab('Node degree')+ylab('Counts')+ggtitle('B')+
  facet_wrap(~source,scale='free',nrow=2)+
  theme_classic()+
  theme(legend.position = 'none',
        plot.title = element_text(size = 15),
        axis.text = element_text(size=12,colour="black"),
        axis.title = element_text(size=15,colour="black"),
        strip.text = element_text(size=12,colour="black"))


ldf<-ddf[ddf$center>180 & ddf$degree>5,]
#fig3C
ggplot(ddf,aes(degree,center))+
  geom_point()+
  xlab('Degree')+ylab('Central score')+ggtitle('C')+
  facet_wrap(~source,scale='free',nrow=1)+
  theme_classic()+
  theme(legend.position = 'none',
        plot.title = element_text(size = 15),
        axis.text = element_text(size=12,colour="black"),
        axis.title = element_text(size=15,colour="black"),
        strip.text = element_text(size=12,colour="black"))+
  geom_text_repel(data=ldf,aes(label=gene),size=3,box.padding = 0.1,max.overlaps=20)
dev.off()
```

## Fig. 4
```
kc1<-fastgreedy.community(gnet1)
kc2<-fastgreedy.community(gnet2)
kc3<-fastgreedy.community(gnet3)
kc4<-fastgreedy.community(gnet4)

dnd1 <- as.dendrogram(kc1)
dnd2 <- as.dendrogram(kc2)
dnd3 <- as.dendrogram(kc3)
dnd4 <- as.dendrogram(kc4)
dnd1 <- ladder(dnd1)
dnd2 <- ladder(dnd2)
dnd3 <- ladder(dnd3)
dnd4 <- ladder(dnd4)

dndlist <- dendextend::dendlist(dnd2, dnd1)
dendextend::tanglegram(dndlist, fast = TRUE, margin_inner = 5)

tanglegram(rank_branches(dnd2), rank_branches(dnd1), edge.lwd = 2,
           margin_inner = 5, type = "t", center = TRUE,
           axes=F)
```



```
df<-read.table('IBD.counts.txt')
cdf<-read.table('IBD.celltype.txt',sep='\t')
names(cdf)[3]<-'CT'
df<-cbind(cdf[,3,drop=F],df)
df<-ddply(df,.(CT),function(x){colSums(x[-1])})
annotLevels <- list(level1class = df[,1], level2class =df[,1])
df<-t(df[,-1])
names(df)<-paste0("I",seq(ncol(df)))
gdf<-read.table('gene.txt')
row.names(df)<-gdf[,1]

ctd_file <- generate_celltype_data(
    exp=as.matrix(df),
    annotLevels=annotLevels,
    groupName='IBD',savePath='./'
)
```
