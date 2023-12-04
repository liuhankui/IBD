# IBD analysis script 

## 1. prepare IBD genes and gut scRNA data
## Pre-step 1: IBD loci [shell]
```
# Liu, J. Z. et al. Association analyses identify 38 susceptibility loci for inflammatory bowel disease and highlight shared genetic risk across populations. Nature genetics 47, 979–986 (2015)
wget --no-check-certificate https://static-content.springer.com/esm/art%3A10.1038%2Fng.3359/MediaObjects/41588_2015_BFng3359_MOESM22_ESM.xlsx

# transform 41588_2015_BFng3359_MOESM22_ESM.xlsx to tableS2.txt by yourself

wget --no-check-certificate https://data.broadinstitute.org/snowman/hg19/variant_calling/vqsr_resources/Exome/v2/dbsnp_138.b37.vcf.gz
# warning: 1.8G file size

zcat dbsnp_138.b37.vcf.gz|awk '!/#/{print $1"_"$2,$3,$4,$5}' | sort -k1,1 -T ./tmp > pos2allele.txt
awk -F '\t' 'NR>1{print $1"_"$3,$2,$10,$11}' tableS2.txt|sort -k1,1 > ibd.pos
join ibd.pos pos2allele.txt|awk '{if($1!=s){print $1,$2,$6,$7,"D1="$3";D2="$4};s=$1}'|tr '_' ' '|awk '{print $1,$2,$3,$4,$5,". .",$6}'|tr ' ' '\t'  > ibd.vcf

# VEP, https://github.com/Ensembl/ensembl-vep
# warning: install vep and download cache

vep --assembly GRCh37 -i ibd.vcf -o ibd.vep --vcf --merged --offline --use_given_ref --per_gene --symbol --canonical --protein --biotype --nearest symbol --fasta hg19.fa.gz --dir_cache ./cache
```

## Pre-step 2: get gene z-score from gnomAD [shell]
```
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
zcat gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz|cut -f '1,2,5,21,33,34' > ./data/gene.zscore
```

## Pre-step 3: get discovery scRNA data [shell]
```
# Elmentaite, R. et al. Cells of the human intestinal tract mapped across space and time. Nature 597, 250–255 (2021)
wget --no-check-certificate https://cellgeni.cog.sanger.ac.uk/gutcellatlas/Full_obj_raw_counts_nosoupx_v2.h5ad -O discovery.h5ad
# warning: 5.8G file size
```

## Pre-step 4: get replication scRNA data [shell]
```
# colon immune
wget https://datasets.cellxgene.cziscience.com/18ddb5f7-d6fe-4fb4-9621-f231fe19088b.h5ad -O replication.immune.h5ad
# colon epithelial
wget https://datasets.cellxgene.cziscience.com/6afb8034-4a11-4a08-8ec2-048968865db9.h5ad -O replication.epithelial.h5ad
# colon stromal
wget https://datasets.cellxgene.cziscience.com/9531d09f-bcd9-48b6-a545-3bdcc5c7cfe1.h5ad -O replication.stromal.h5ad
```

## Pre-step 5: get counts from discovery h5ad file format [R]
```
library(rhdf5)
library(Matrix)

h5<-H5Fopen("discovery.h5ad")
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
cdf$sample<-as.character(factor(cdf$Diagnosis,levels=c('fetal','Healthy adult','Pediatric Crohn Disease','Pediatric healthy'),labels=c('HF','HA','pCD','HP')))
h5closeAll()

write.table(gdf,file='discovery.gene.txt',quote=F,row.names=F,col.names=F,sep='\t')
write.table(cdf,file='discovery.celltype.txt',quote=F,row.names=F,col.names=F,sep='\t')

# $ data   : num [1:760344941(1d)] 1 1 1 4 1 1 1 1 1 1 ...
# $ indices: int [1:760344941(1d)] 15 53 102 154 216 223 244 269 271 326 ...
# $ indptr : int [1:428470(1d)] 0 970 1672 2394 3141 4369 5545 6131 6906 7533 ...

unlink('discovery.counts.txt')
cells<-diff(indptr)
marker<-seq(0,length(cells),by=1000)
index<-1
for(k in marker){
  start<-k+1
  end<-k+1000
  if(end>length(cells)){end<-length(cells)}
  rowN<-rep(1:(end-start+1),cells[start:end])
  num<-length(rowN)
  colN<-h5read("discovery.h5ad", "/X/indices", index=list(index:(index+num-1)))
  counts<-h5read("discovery.h5ad", "/X/data", index=list(index:(index+num-1)))
  index<-index+num
  df<-as.matrix(sparseMatrix(i = rowN, j = colN + 1,x = as.numeric(counts),dims=c(end-start+1,nrow(gdf))))
  write.table(df,file="discovery.counts.txt",append=T,quote=F,row.names=F,col.names=F,sep=' ')
}
```

## Pre-step 6: split discovery data [shell]
```
cat discovery.celltype.txt|awk -F '\t' '{print $2"\t"$3 >> $4".celltype.txt"}'
cat discovery.celltype.txt|cut -f 4|paste - discovery.counts.txt|awk -F '\t' '{print $2 >> $1".counts.txt"}'

for i in HF HP HA pCD
do
cat $i.celltype.txt|paste - $i.counts.txt|awk -v i=$i -F '\t' '$2=="Monocytes"{print $3 >> i".Monocytes.txt"}'
cat $i.Monocytes.txt|gzip -f > ./data/$i.Monocytes.txt.gz
for j in Monocytes ILC3 Th1
do
echo -n "$i $j "
cat $i.celltype.txt|awk -F '\t' -v c=$j -v a=0  -v b=0 '{if($2==c){a++};b++}END{print a,b}'
done
done > ./data/cell.counts.txt
```

## Pre-step 7: calculate cell-type exrpression specificty from pCD data [R]
```
library(EWCE)

df<-read.table('./pCD.counts.txt')
cdf<-read.table('./pCD.celltype.txt',sep='\t')
annotLevels <- list(level1class = cdf[,2], level2class =cdf[,1])
df<-t(df[,-1])
#names(df)<-paste0("I",seq(ncol(df)))
gdf<-read.table('discovery.gene.txt')
row.names(df)<-gdf[,1]

ctd_file <- generate_celltype_data(
    exp=as.matrix(df),
    annotLevels=annotLevels,
    groupName='pCD',savePath='./data/'
)
```


## Pre-step 8: get counts from replication h5ad file format [R]
```
library(rhdf5)
library(Matrix)
for(h5file in c('replication.immune.h5ad','replication.epithelial.h5ad','replication.stromal.h5ad')){
  h5<-H5Fopen(h5file)
  gdf<-data.frame(symbol=as.character(factor(h5$'var/feature_name/codes',levels=0:(length(h5$'var/feature_name/codes')-1),labels=h5$'var/feature_name/categories'))
  cdf<-data.frame(
      Diagnosis=as.character(factor(h5$'obs/disease/codes',levels=0:(length(h5$'obs/disease/categories')-1),labels=h5$'obs/disease/categories')),
      Category=as.character(factor(h5$'obs/cell_type/codes',levels=0:(length(h5$'obs/cell_type/categories')-1),labels=h5$'obs/cell_type/categories')),
      Subtype=as.character(factor(h5$'obs/Celltype/codes',levels=0:(length(h5$'obs/Celltype/categories')-1),labels=h5$'obs/Celltype/categories')),
      Tissue=as.character(factor(h5$'obs/tissue/codes',levels=0:(length(h5$'obs/tissue/categories')-1),labels=h5$'obs/tissue/categories')),
      Status=as.character(factor(h5$'obs/Type/codes',levels=0:(length(h5$'obs/Type/categories')-1),labels=h5$'obs/Type/categories'))
  )
  indptr<-h5$'X/indptr'
  h5closeAll()
  write.table(cdf,file=paste0(h5file,'.celltype.txt'),quote=F,row.names=F,col.names=F,sep='\t')
  write.table(gdf,file=paste0(h5file,'.gene.txt'),quote=F,row.names=F,col.names=F,sep='\t')

  unlink(paste0(h5file,".counts.txt"))
  cells<-diff(indptr)
  marker<-seq(0,length(cells),by=1000)
  index<-1
  for(k in marker){
    start<-k+1
    end<-k+1000
    if(end>length(cells)){end<-length(cells)}
    rowN<-rep(1:(end-start+1),cells[start:end])
    num<-length(rowN)
    colN<-h5read(args[6], "/X/indices", index=list(index:(index+num-1)))
    counts<-h5read(args[6], "/X/data", index=list(index:(index+num-1)))
    index<-index+num
    df<-as.matrix(sparseMatrix(i = rowN, j = colN + 1,x = as.numeric(counts),dims=c(end-start+1,27345)))
    write.table(df,file=paste0(h5file,".counts.txt"),append=T,quote=F,row.names=F,col.names=F,sep=' ')
  }
}
```

## Pre-step 9: split replication data [shell]
```
cat replication.immune.h5ad.celltype.txt|tr ' ' '_'|paste - replication.immune.h5ad.counts.txt|awk -F '\t' '$2=="monocyte"{print $1,$4,$5,$6}'|gzip -f >  replication.monocyte.gz
cat replication.immune.h5ad.celltype.txt eplication.epithelial.h5ad.celltype.txt replication.stromal.h5ad.celltype.txt|tr ' ' '_' >  replication.celltype.txt
cat replication.immune.h5ad.counts.txt eplication.epithelial.h5ad.counts.txt replication.stromal.h5ad.counts.txt|paste replication.celltype.txt -|awk -F '\t' '$1=="Crohn_disease"{print $2,$3,$6}'|gzip -f > aCD.counts.txt.gz
```

## Pre-step 10: calculate cell-type exrpression specificty from aCD data [R]
```
library(EWCE)

df<-read.table(gzfile('./aCD.counts.txt.gz'))
annotLevels <- list(level1class = df[,2], level2class =df[,1])
df<-t(df[,-c(1,2)])
#names(df)<-paste0("I",seq(ncol(df)))
gdf<-read.table('replication.immune.h5ad.gene.txt')
row.names(df)<-gdf[,1]

ctd_file <- generate_celltype_data(
    exp=as.matrix(df),
    annotLevels=annotLevels,
    groupName='aCD',savePath='./data/'
)
```

## 2. IBD code script, you can start from here, step-by-step [R]
```
git clone --recursive https://github.com/liuhankui/IBD.git
cd IBD
R
```

## Require R packages
```
library(RColorBrewer) #Fig1-3
library(ggplot2)      #Fig1-3
library(ggpubr)       #Fig1C
library(geomtextpath) #Fig1DE
library(intergraph)   #Fig2
library(Matrix)       #Fig2
library(huge)         #Fig2
library(igraph)       #Fig2
library(network)      #Fig2
library(ggnetwork)    #Fig2
library(reshape2)     #Fig3A
library(ggrepel)      #Fig3C
library(dendextend)   #Fig4
library(polynom)      #Fig4
```

## Fig. 1A

```
source('./bin/bootstrap_enrichment_test.r')
source('./bin/cell_list_dist.r')
source('./bin/generate_controlled_bootstrap_geneset.r')
source('./bin/get_summed_proportions.r')
load('./data/ctd_pCD.rda')

gdf<-read.table('./data/ibd.gene')
x<-unique(gdf$V1)
bg<-rownames(ctd[[1]]$specificity)
hits<-x[x %in% bg]
set.seed(2023)
rdf<-bootstrap_enrichment_test(sct_data=ctd,
                               hits=hits,
                               bg=bg,
                               reps=10000,
                               annotLevel=1)
rdf$results$celltype<-row.names(rdf$results)
rdf$results$FDR<-p.adjust(rdf$results$p,method='fdr')
write.table(rdf$results,file='enrichment.txt',sep='\t',quote=F,row.names=F,col.names=T)

df<-read.table('enrichment.txt',sep='\t',head=T)
df$sign<-ifelse(df$FDR<0.05,'*',ifelse(df$FDR<0.01,'**',ifelse(df$FDR<0.005,'***',NA)))

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
df<-read.table('./data/cell.counts.txt')
ggplot(df,aes(V1,V3/V4*100)+
  geom_histogram(stat='identity',aes(fill=V2,group=V2),colour='black',width=0.8,position='dodge')+
  geom_text(aes(V1,label=V3,group=V2),position=position_dodge(width = 0.9),hjust=-0.25)+
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
df<-read.table(gzfile('./data/pCD.Monocytes.txt.gz'))
ss<-colSums(df)>0
df<-df[,ss]
df<-log(df+1)
sdf<-read.table('./data/discovery.gene.txt')
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
                                  c('HI','LOFT'),
                                  c('IBD','LOFT')),
              parse = T,y_position = c(0.8,1,1.2),tip_length=0.01,test = "t.test",
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

ibd_s<-mean(sdf$sd[sdf$V1 %in% IBD])
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
zdf<-read.table('./data/gene.zscore',head=T)
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
gdf<-read.table('./data/ibd.gene')
gene<-read.table('./data/discovery.gene.txt')$V1

#-----------network pCD----------------------
df<-read.table(gzfile('./data/pCD.Monocytes.txt.gz'))
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
adj<-as.matrix(mbOptRICGraph)
net<-igraph::graph_from_adjacency_matrix(adj,mode="undirected")
tmp<-as.matrix(sparseMatrix(i=mbOptRIC$refit@i+1,p=mbOptRIC$refit@p,x=mbOptRIC$refit@x,dims=c(ncol(df),ncol(df))),diag=T)
igraph::V(net)$name<-names(df)
igraph::V(net)$size<-as.numeric(exp)
igraph::E(net)$weight<-abs(cor.matrix[which(tmp==1 & lower.tri(tmp),arr.ind=TRUE)])
igraph::V(net)$colour<-as.character(factor(V(net)$name,levels=gdf$V1[gdf$V1 %in% V(net)$name],labels=gdf$V2[gdf$V1 %in% V(net)$name]))
gnet1<-net
ndf1<-ggnetwork(net)
ndf1$source<-'Paediatric CD'

#-----------network HP----------------------
df<-read.table(gzfile('./data/HP.Monocytes.txt.gz'))
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
adj<-as.matrix(mbOptRICGraph)
net<-igraph::graph_from_adjacency_matrix(adj,mode="undirected")
tmp<-as.matrix(sparseMatrix(i=mbOptRIC$refit@i+1,p=mbOptRIC$refit@p,x=mbOptRIC$refit@x,dims=c(ncol(df),ncol(df))),diag=T)
igraph::V(net)$name<-names(df)
igraph::V(net)$size<-as.numeric(exp)
igraph::E(net)$weight<-abs(cor.matrix[which(tmp==1 & lower.tri(tmp),arr.ind=TRUE)])
igraph::V(net)$colour<-as.character(factor(V(net)$name,levels=gdf$V1[gdf$V1 %in% V(net)$name],labels=gdf$V2[gdf$V1 %in% V(net)$name]))
gnet2<-net
ndf2<-ggnetwork(net)
ndf2$source<-'Paediatric Healthy'

#-----------network HF----------------------
df<-read.table(gzfile('./data/HF.Monocytes.txt.gz'))
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
adj<-as.matrix(mbOptRICGraph)
net<-igraph::graph_from_adjacency_matrix(adj,mode="undirected")
tmp<-as.matrix(sparseMatrix(i=mbOptRIC$refit@i+1,p=mbOptRIC$refit@p,x=mbOptRIC$refit@x,dims=c(ncol(df),ncol(df))),diag=T)
igraph::V(net)$name<-names(df)
igraph::V(net)$size<-as.numeric(exp)
igraph::E(net)$weight<-abs(cor.matrix[which(tmp==1 & lower.tri(tmp),arr.ind=TRUE)])
igraph::V(net)$colour<-as.character(factor(V(net)$name,levels=gdf$V1[gdf$V1 %in% V(net)$name],labels=gdf$V2[gdf$V1 %in% V(net)$name]))
gnet3<-net
ndf3<-ggnetwork(net)
ndf3$source<-'Fetal Healthy'

#-----------network HA----------------------
df<-read.table(gzfile('./data/HA.Monocytes.txt.gz'))
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
adj<-as.matrix(mbOptRICGraph)
net<-igraph::graph_from_adjacency_matrix(adj,mode="undirected")
tmp<-as.matrix(sparseMatrix(i=mbOptRIC$refit@i+1,p=mbOptRIC$refit@p,x=mbOptRIC$refit@x,dims=c(ncol(df),ncol(df))),diag=T)
igraph::V(net)$name<-names(df)
igraph::V(net)$size<-as.numeric(exp)
igraph::E(net)$weight<-abs(cor.matrix[which(tmp==1 & lower.tri(tmp),arr.ind=TRUE)])
igraph::V(net)$colour<-as.character(factor(V(net)$name,levels=gdf$V1[gdf$V1 %in% V(net)$name],labels=gdf$V2[gdf$V1 %in% V(net)$name]))
gnet4<-net
ndf4<-ggnetwork(net)
ndf4$source<-'Adult Healthy'

#----network merge----------------
ndf<-rbind(ndf1,ndf2,ndf3,ndf4)
ndf$source<-factor(ndf$source,levels=unique(ndf$source),order=T)
netA<-ggplot(ndf,aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(aes(linewidth=weight),colour='gray') +
  geom_nodes(aes(colour=colour,size=size))+
  geom_nodetext_repel(aes(label=name),size=2,box.padding = 0.1)+
  scale_color_brewer('Genes',palette = "Set2") +
  scale_linewidth_continuous('Correlations',range = c(0.1,1))+
  scale_size('Expression')+
  facet_wrap(~source,scale='free',nrow=2)+
  ggtitle('A')+
  theme_blank()+
  theme(panel.border = element_rect(colour='black',fill=NA))

#-------rep--------------
gdf<-read.table('ibd.gene')
gene<-read.table('replication.immune.h5ad.gene.txt')$V1
dfr<-read.table(gzfile('replication.monocyte.gz'))
names(dfr)<-c('V1','V2','V3',gene)
dfr<-dfr[,c(rep(T,3),gene %in% gdf$V1)]

# aCD Infl
df<-dfr[dfr$V1=="Crohn_disease" & dfr$V2=="sigmoid_colon" & dfr$V3=="Infl",]
df<-df[,colSums(df)>0]
df<-df[,colSums(df==0)/(nrow(df)-1)<0.9]
df<-log(df+1)
exp<-colSums(df)/nrow(df)
cor.matrix<- cor(df,use = "pairwise.complete.obs")
mbModel <- huge(as.matrix(df), method = "mb")
set.seed(2023)
mbOptRIC = huge.select(mbModel,criterion="stars")
mbOptRICGraph = mbOptRIC$refit
adj<-as.matrix(mbOptRICGraph)
net<-igraph::graph_from_adjacency_matrix(adj,mode="undirected")
tmp<-as.matrix(sparseMatrix(i=mbOptRIC$refit@i+1,p=mbOptRIC$refit@p,x=mbOptRIC$refit@x,dims=c(ncol(df),ncol(df))),diag=T)
igraph::V(net)$name<-names(df)
igraph::V(net)$size<-as.numeric(exp)
igraph::E(net)$weight<-abs(cor.matrix[which(tmp==1 & lower.tri(tmp),arr.ind=TRUE)])
igraph::V(net)$colour<-as.character(factor(V(net)$name,levels=gdf$V1[gdf$V1 %in% V(net)$name],labels=gdf$V2[gdf$V1 %in% V(net)$name]))
net1<-net
ndf1<-ggnetwork(net)
ndf1$source<-'Adult CD inflammation'

#aCD NonI---------------
df<-dfr[dfr$V1=="Crohn_disease" & dfr$V2=="right_colon" & dfr$V3=="NonI",]
df<-df[,colSums(df)>0]
df<-df[,colSums(df==0)/(nrow(df)-1)<0.9]
df<-log(df+1)
exp<-colSums(df)/nrow(df)
cor.matrix<- cor(df,use = "pairwise.complete.obs")
mbModel <- huge(as.matrix(df), method = "mb")
set.seed(2023)
mbOptRIC = huge.select(mbModel,criterion="stars")
mbOptRICGraph = mbOptRIC$refit
adj<-as.matrix(mbOptRICGraph)
net<-igraph::graph_from_adjacency_matrix(adj,mode="undirected")
tmp<-as.matrix(sparseMatrix(i=mbOptRIC$refit@i+1,p=mbOptRIC$refit@p,x=mbOptRIC$refit@x,dims=c(ncol(df),ncol(df))),diag=T)
igraph::V(net)$name<-names(df)
igraph::V(net)$size<-as.numeric(exp)
igraph::E(net)$weight<-abs(cor.matrix[which(tmp==1 & lower.tri(tmp),arr.ind=TRUE)])
igraph::V(net)$colour<-as.character(factor(V(net)$name,levels=gdf$V1[gdf$V1 %in% V(net)$name],labels=gdf$V2[gdf$V1 %in% V(net)$name]))
net2<-net
ndf2<-ggnetwork(net)
ndf2$source<-'Adult CD non-inflammation'

#-----------HC----------------------
df<-dfr[dfr$V1=="normal" & dfr$V2=='right_colon',]
df<-df[,colSums(df)>0]
df<-df[,colSums(df==0)/(nrow(df)-1)<0.9]
df<-log(df+1)
exp<-colSums(df)/nrow(df)
cor.matrix<- cor(df,use = "pairwise.complete.obs")
mbModel <- huge(as.matrix(df), method = "mb")
set.seed(2023)
mbOptRIC = huge.select(mbModel,criterion="stars")
mbOptRICGraph = mbOptRIC$refit
adj<-as.matrix(mbOptRICGraph)
net<-igraph::graph_from_adjacency_matrix(adj,mode="undirected")
tmp<-as.matrix(sparseMatrix(i=mbOptRIC$refit@i+1,p=mbOptRIC$refit@p,x=mbOptRIC$refit@x,dims=c(ncol(df),ncol(df))),diag=T)
igraph::V(net)$name<-names(df)
igraph::V(net)$size<-as.numeric(exp)
igraph::E(net)$weight<-abs(cor.matrix[which(tmp==1 & lower.tri(tmp),arr.ind=TRUE)])
igraph::V(net)$colour<-as.character(factor(V(net)$name,levels=gdf$V1[gdf$V1 %in% V(net)$name],labels=gdf$V2[gdf$V1 %in% V(net)$name]))
net3<-net
ndf3<-ggnetwork(net)
ndf3$source<-'Healthy adult'

ndf<-rbind(ndf1,ndf2,ndf3)
ndf$source<-factor(ndf$source,levels=unique(ndf$source),order=T)
netB<-ggplot(ndf,aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(aes(linewidth=weight),colour='gray') +
  geom_nodes(aes(colour=colour,size=size))+
  geom_nodetext_repel(aes(label=name),size=2,box.padding = 0.1)+
  scale_color_brewer('Genes',palette = "Set2") +
  scale_linewidth_continuous('Correlations',range = c(0.1,1))+
  scale_size('Expression')+
  facet_wrap(~source,scale='free',nrow=2)+
  ggtitle('B')+
  theme_blank()+
  theme(panel.border = element_rect(colour='black',fill=NA))

ggarrange(netA,netB,nrow=2,heights = c(1,1))

```


## Fig3. A
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

tmp1<-data.frame(V=length(V(gnet1)),E=length(E(gnet1)),D=graph.density(gnet1),T=transitivity(gnet1),C=cl(gnet1),S='pCD')
tmp2<-data.frame(V=length(V(gnet2)),E=length(E(gnet2)),D=graph.density(gnet2),T=transitivity(gnet2),C=cl(gnet2),S='PH')
tmp3<-data.frame(V=length(V(gnet3)),E=length(E(gnet3)),D=graph.density(gnet3),T=transitivity(gnet3),C=cl(gnet3),S='FH')
tmp4<-data.frame(V=length(V(gnet4)),E=length(E(gnet4)),D=graph.density(gnet4),T=transitivity(gnet4),C=cl(gnet4),S='AH')
tdf<-rbind(tmp1,tmp2,tmp3,tmp4)
tdf<-melt(tdf[,-4],id='S')

tdf$variable<-factor(tdf$variable,
                     levels=c('V','E','D','C'),
                     labels=c('Nodes','Edges','Density','Cluster'))

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
```

## Fig. 3B
```
tmp1<-data.frame(degree=igraph::degree(gnet1),center=betweenness(gnet1))
tmp1$gene<-row.names(tmp1)
tmp1$source<-'pCD'

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
```

## Fig. 3C
```
ldf<-ddf[ddf$center>180 & ddf$degree>5,]
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
```

## Fig. 4
```
gdf<-read.table('./data/ibd.gene')
gene<-read.table('./data/discovery.gene.txt')$V1

df<-read.table(gzfile('./data/HP.Monocytes.txt.gz'))
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
adj<-as.matrix(mbOptRICGraph)
colnames(adj)<-names(df)
rownames(adj)<-names(df)
ss<-rowSums(adj)>1
gg<-names(df)[ss]
adj<-adj[ss,ss]
net<-igraph::graph_from_adjacency_matrix(adj,mode="undirected")
Net_pHC<-net

df<-read.table(gzfile('./data/pCD.Monocytes.txt.gz'))
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
adj<-as.matrix(mbOptRICGraph)
colnames(adj)<-names(df)
rownames(adj)<-names(df)
ss<-names(df) %in% gg
adj<-adj[ss,ss]
net<-igraph::graph_from_adjacency_matrix(adj,mode="undirected")
Net_pCD<-net

kc1<-fastgreedy.community(Net_pCD)
kc2<-fastgreedy.community(Net_pHC)
dnd1 <- as.dendrogram(kc1)
dnd2 <- as.dendrogram(kc2)
dnd1 <- phylogram::ladder(dnd1)
dnd2 <- phylogram::ladder(dnd2)

treeplot<-tanglegram(rank_branches(dnd2), rank_branches(dnd1), edge.lwd = 2,
                  margin_inner = 5, type = "t", center = TRUE,
                  axes=F,
                  k_branches = 7)

plot(treea,margin_inner = 5)
```

## Fig. S1
```
source('./bin/bootstrap_enrichment_test.r')
source('./bin/cell_list_dist.r')
source('./bin/generate_controlled_bootstrap_geneset.r')
source('./bin/get_summed_proportions.r')

load('D:/IBD/ctd_aCD.rda')
gdf<-read.table('ibd.gene')
x<-unique(gdf$V1)
bg<-row.names(ctd[[1]]$specificity)
hits<-x[x %in% bg]
set.seed(1)
rdf<-bootstrap_enrichment_test(sct_data=ctd,
                               hits=hits,
                               bg=bg,
                               reps=10000,
                               annotLevel=1)
rdf$results$celltype<-row.names(rdf$results)
rdf$results$FDR<-p.adjust(rdf$results$p,method='fdr')
write.table(rdf$results,file='enrichment.txt',sep='\t',quote=F,row.names=F,col.names=T)

df<-read.table('enrichment.txt',sep='\t',head=T)
df<-df[order(df$p,decreasing = T),]
df$CellType<-factor(df$CellType,levels=df$CellType,order=T)
df$sign<-ifelse(df$p<0.005,'*',NA)

ggplot(df,aes(CellType,-log10(p)))+
  geom_histogram(stat='identity')+
  geom_text(aes(label=sign),hjust=0,vjust=0.75,size=5)+
  coord_flip()+
  xlab('')+
  scale_y_continuous(expression(paste(-log[10],"P-value")),
                     limits = c(0,5),expand=c(0,0))+
  theme_classic()+
  theme(legend.position = c(0.85,0.9),
        plot.title = element_text(size = 15,colour="black"),
        axis.text = element_text(size=12,colour="black"),
        axis.title = element_text(size=15,colour="black"),
        strip.text = element_text(size=12,colour="black"))
```



## Fig S2
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

tmp1<-data.frame(V=length(V(net1)),E=length(E(net1)),D=graph.density(net1),T=transitivity(net1),C=cl(net1),S='CDI')
tmp2<-data.frame(V=length(V(net2)),E=length(E(net2)),D=graph.density(net2),T=transitivity(net2),C=cl(net2),S='CDN')
tmp3<-data.frame(V=length(V(net3)),E=length(E(net3)),D=graph.density(net3),T=transitivity(net3),C=cl(net3),S='HA')
tdf<-rbind(tmp1,tmp2,tmp3)
tdf<-melt(tdf[,-4],id='S')
tdf$S<-factor(tdf$S,levels=c('CDI','CDN','HA'),
              labels=c("CD\ninflammation","CD non\ninflammation","Healthy\nadult"),
              order=T)
tdf$variable<-factor(tdf$variable,
                     levels=c('V','E','D','C'),
                     labels=c('Nodes','Edges','Density','Cluster'))

figs2A<-ggplot()+
  geom_histogram(data=tdf,aes(S,value,fill=variable),stat='identity',colour='black')+
  scale_fill_brewer(palette = "Set2")+
  xlab('')+ylab('# summary')+ggtitle('A')+
  facet_wrap(~variable,scale='free',nrow=2)+
  #coord_flip()+
  theme_classic()+
  theme(legend.position = 'none',
        plot.title = element_text(size = 15),
        axis.text = element_text(size=12,colour="black"),
        axis.text.x = element_text(angle=-45,hjust=0,size=8),
        axis.title = element_text(size=15,colour="black"),
        strip.text = element_text(size=12,colour="black"))



tmp1<-data.frame(degree=igraph::degree(net1),center=betweenness(net1))
tmp1$gene<-row.names(tmp1)
tmp1$source<-'CDI'

tmp2<-data.frame(degree=igraph::degree(net2),center=betweenness(net2))
tmp2$gene<-row.names(tmp2)
tmp2$source<-'CDN'

tmp3<-data.frame(degree=igraph::degree(net3),center=betweenness(net3))
tmp3$gene<-row.names(tmp3)
tmp3$source<-'HA'

ddf<-rbind(tmp1,tmp2,tmp3)
ddf$source<-factor(ddf$source,levels=c('CDI','CDN','HA'),
                   labels=c('CD inflammation','CD non-inflammation','Healthy adult'),
                   order=T)
figs2B<-
  ggplot()+
  geom_histogram(data=ddf,aes(degree,fill=source),binwidth = 1,colour='black')+
  scale_fill_brewer(palette = "Set2")+
  xlab('Node degree')+ylab('Counts')+ggtitle('B')+
  facet_wrap(~source,scale='free_y',nrow=3)+
  #coord_flip()+
  theme_classic()+
  theme(legend.position = 'none',
        plot.title = element_text(size = 15),
        axis.text = element_text(size=12,colour="black"),
        axis.title = element_text(size=15,colour="black"),
        strip.text = element_text(size=12,colour="black"))


ldf<-ddf[ddf$center>50 & ddf$degree>4,]
figs2C<-ggplot(ddf,aes(degree,center))+
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

figs2AB<-ggarrange(figsA,figsB,nrow=1,widths = c(1,1))
ggarrange(figs2AB,figs2C,ncol=1,heights  = c(1,1))
```

