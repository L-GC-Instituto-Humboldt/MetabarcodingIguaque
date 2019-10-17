---
title: "Routine_metabar_Serengeti2017_withLocalDB"
author: "Johan Pansu"
date: "12/12/2018"
output: html_document
---

Routine series for METABAR data quality checks and filtering
========================================================

*Framework authors: L. Zinger, J. Pansu, J. Chave, and others*

*`r date()`*

```{r, echo=F}
#chunk1
#for formatting purpose only: add captions in html format
#knit_hooks$set(htmlcap = function(before, options, envir) {
#  if(!before) {
#    paste('<p class="caption">',options$htmlcap,"</p>",sep="")
#    }
#    })
```


This is an R code framework for **METABAR data pre-analyses**. It consists in:
   - Monitoring initial characteristics of the data after obitools procedures (and eventually during)
   - Excluding putative contaminants
   - Detecting samples for which PCR did not work
  

Source characteristics and tips:
- Written in Markdown format (more here: http://www.rstudio.com/ide/docs/r_markdown, and here:http://support.iawriter.com/help/kb/general-questions/markdown-syntax-reference-guide)
- To be viewed/edited in Rstudio (available here: http://www.rstudio.com). 
- Requires knitr library to process R codes (to be installed in R, more here:http://yihui.name/knitr/ and here: https://github.com/yihui/knitr). KnitR should be set as the weave default program in R studio (see more here: http://yihui.name/knitr/demo/rstudio/)
- To be run on R 64 bit version (several intensive computation steps).  
- Produces html outputs including plain text, codes and figures
- Can be translated into pdf using the following command: pandoc file.md -t latex -o file.pdf --latex-engine=[set your pdflatex path] (requires pandoc, available here: http://johnmacfarlane.net/pandoc/ and lateX)



## 1/Libraries needed: ##


```{r, message=F, warning=F, comment=F}
#chunk2
library(vegan)
library(lattice)
library(ROBITools)
library(ROBITaxonomy)
library(plotrix)
library(ENmisc)
library(igraph)
library(pander)
library(stringr)
```



## 2/Inputs ##

### Set paths and file names: ###
```{r}
#chunk3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!!! MODIFY HERE !!!!
PATH="/home/henry/Documents/Research/MetabarcodingIguaque/results/testR" #working directory
OBJ="/home/henry/Documents/Research/MetabarcodingIguaque/results/testR/plants.tab" #obitool output
DB_N="/home/henry/Documents/Research/MetabarcodingIguaque/Reference_DB/embl_r134" # ecopcr database, used only for taxid manipulation purposes
replication=T
```

Inputs description:
  - obitool output: a text file obtained from the obitool pipeline:
      1. paired-end reads assembly (illuminapairedend)
      2. read assignment to a sample (ngsfilter)
      3. Removal of short reads (min length depends on the marker) and reads with "N"s (obigrep)
      4. Removal of reads  not assembled in step a/ (obigrep)
      5. Reads dereplication (obiuniq)
      6. Removal of singletons (sequences of count = 1 in the overall dataset)
      7. Small variants identification and affiliation (obiclean)
      8. Removal of small variants (obigrep (head+singleton)>internal for now) ***!!!!! TO BE MODIFIED AFTER ERIC's CHANGES ON OBICLEAN***
      9. Chimera status annotation (MOTHUR's chimera.perseus)
      10. Taxonomic assignment (ecoTag)
      11. Sequence clustering (optionnal, depends on the marker, sumatra+mcl FOR NOW)
      12. Fasta to table format (obitab), final obitool output
  - ecopcr database: raw embl database in ecopcr format  or enriched with local taxids (extentions rdx, ldx, ndx, tdx).
  - sample coordinates: a text file containing a table with geographic coordinates, i.e. at least 2 columns, named "x" and "y"
  - sample coordinates in PCR 96-wells (checks for problems in tags)
  - obitool data size check file: a text file containing the number of reads/otus remaining after each obitool filtering steps
  - colors chosen for plots.


### Initial formatting ###

#### Data import ####
```{r, cache=T, autodep=TRUE}
#chunk4
setwd(PATH) 
DB=read.taxonomy(DB_N)
OBI=import.metabarcoding.data(OBJ)


# Aggregate the sequence belonging to the same OTU: Insert Henry script
OBI_reads<-aggregate(t(OBI@reads), by=list(OBI@motus$infomapclust_d3), sum)
rownames(OBI_reads)<-OBI_reads[,1]
OBI_reads<-OBI_reads[,-1]

# Explore funtion split by lapply (OBI@reads) to sort OTUs base on their abundance and pick up the most abundant sequence as representative


```

Data format tip:  
Data (obitools final output) are imported in a metabarcoding format (R S4 object).  
Here, OBI is a S4 object containing a set of tables with the following structure (+ other stuffs not important and not used below):

read table  
R R R R R R R R R R R R-------S S S S  
R R R R R R R R R R R R-------S S S S  
R R R R R R R R R R R R-------S S S S  sample table  
R R R R R R R R R R R R-------S S S S  
R R R R R R R R R R R R-------S S S S  
R R R R R R R R R R R R-------S S S S  
| | | | | | | | | | | |  
M M M M M M M M M M M M  
M M M M M M M M M M M M   motu table  
M M M M M M M M M M M M   

 - read table (call: OBI@reads): corresponds to the community table, the main object used here
 - sample table (call: OBI@samples): stores sample characteristics (e.g. store xy coordinates, etc...)
 - motus table (call: OBI@motus): stores motus characteristics (e.g. taxonomic assignements, etc...)

S4 object class is very useful for data manipulation, e.g. when one wants to get rid of certain motus or samples: it allows keeping characteristics of the remaining samples/motus and therefore avoid cascading mistakes that may occur when one juggles with several tables.  

Examples: 
OBI[1,] selects the first line of OBI@reads and corresponding OBI@sample line.   
OBI[1:10,1:100] selects the 10 first rows in OBI@reads, OBI@sample, the first 100 columns in OBI@reads and the first 100 rows in OBI@motus.

#### Taxonomy formatting ####
At this point, one should standardize the taxonomic information amongst datasets (e.g. get always the same taxonomic levels, same columns names and so on).  
When several databases were used for taxonomic identification (e.g. local vs genebank), this step gives priority to taxonomic assignements based on a given database (usually the local, custom one) when the assignement score is satisfying (used here 0.97, lower may assign contaminant Salix sequences to Xylosma). Otherwise, the best taxonomic assignement is kept. This is done with the function decide_taxo. This function needs improvements and is not part of ROBITOOLS yet (too long to be shown here, will be attached to the script)

```{r, cache=T, autodep=TRUE}
#chunk5
if(length(grep("taxid_by_db", colnames(motus(OBI))))>1) {
  source('decide_taxo.R', chdir = TRUE) #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!!! MODIFY HERE !!!!
  OBI=decide_taxo(OBI,DB,"PNG",0.99999) #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!!! CUTOFF HERE !!!! PNG will be the default, the threshold reflects the level of similarity above which the Global DB will be used instead of the local DB. 
  #example 1: if threshold = 0.999, only sequences that matches with the local DB with a score higher than 0.999 will be assigned to the local DB (i.e. in practice, only perfect matches)
  #example 2: if threshold = 0.98, only sequences that matches with the local DB with a score higher than 0.98 will be assigned to the local DB, the others will be assigned to the global DB  
  OBI@motus$bid_ok=round(OBI@motus$bid_ok, 3) 
} else {
  OBI@motus$rank_ok=taxonomicrank(DB,OBI@motus$taxid)
  OBI@motus$species_name_ok=scientificname(DB, taxonatrank(DB,OBI@motus$taxid,"species"))
  OBI@motus$genus_name_ok=scientificname(DB, taxonatrank(DB,OBI@motus$taxid,"genus"))
  OBI@motus$family_name_ok=scientificname(DB, taxonatrank(DB,OBI@motus$taxid,"family"))
  OBI@motus$order_name_ok=scientificname(DB, taxonatrank(DB,OBI@motus$taxid,"order"))
  OBI@motus$class_name_ok=scientificname(DB, taxonatrank(DB,OBI@motus$taxid,"class"))
  OBI@motus$phylum_name_ok=scientificname(DB, taxonatrank(DB,OBI@motus$taxid,"phylum"))
  OBI@motus$kingdom_name_ok=scientificname(DB, taxonatrank(DB,OBI@motus$taxid,"kingdom"))
  OBI@motus$bid_ok=round(OBI@motus$best_identity, 3)
  OBI@motus$scientific_name_ok=OBI@motus$scientific_name
}


```




#### Sample formatting ####
Add any sample characteristics. (e.g. here geographic and plate sample coordinates, data remaining per sample at each obitools step)  
***!!!to be modified accordingly!!!***
```{r, cache=T, autodep=TRUE}
#chunk6
#plate coordinates #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!!! MODIFY HERE !!!!
#Not mandatory but help visualizing if errors occur in some plates specifically. Needs to be modified for each library

PLATE=c(rownames(OBI@reads)[-c(grep("BLK", rownames(OBI@reads)), grep("NEG", rownames(OBI@reads)), grep("MI", rownames(OBI@reads)), grep("MC", rownames(OBI@reads)), grep("MP", rownames(OBI@reads)))])

rownames(OBI@reads)[-c(grep("BLK", rownames(OBI@reads)), grep("NEG", rownames(OBI@reads)), grep("MI", rownames(OBI@reads)), grep("MC", rownames(OBI@reads)), grep("MP", rownames(OBI@reads)), grep("c", rownames(OBI@reads)), grep("b", rownames(OBI@reads)))]<-paste(c(rownames(OBI@reads)[-c(grep("BLK", rownames(OBI@reads)), grep("NEG", rownames(OBI@reads)), grep("MI", rownames(OBI@reads)), grep("MC", rownames(OBI@reads)), grep("MP", rownames(OBI@reads)), grep("c", rownames(OBI@reads)), grep("b", rownames(OBI@reads)))]), "a", sep="")

ngsfilt=read.table('/home/henry/Documents/Research/MetabarcodingIguaque/RawData/ngsfilterIguaquePlant_a.txt') #Read this ngsfilter

ngsfilt_rep1<-as.data.frame(as.matrix(ngsfilt)[grep('a', as.matrix(ngsfilt)[,2]),]) #Keep only plates with replicate 'a'
NB_PLATE_rep1=2 #Nb of plates with replicate 'a'
ind_rep1={1:nrow(ngsfilt_rep1)}
tmp_rep1=seq(0, 8*NB_PLATE_rep1, by=8)
x_rep1=NULL; for (i in 1:(length(tmp_rep1)-1)) {x_rep1=append(x_rep1,rep((tmp_rep1[i]+1):tmp_rep1[i+1],12))}
x_rep1<-x_rep1[1:nrow(ngsfilt_rep1)]
xLoc_rep1<-data.frame(as.vector(ngsfilt_rep1$V2), rep(rep(1:12, each=8), NB_PLATE_rep1)[c(1:length(as.vector(ngsfilt_rep1$V2)))])
colnames(xLoc_rep1)<-c('sample', 'column_in_plate')

ngsfilt_rep2<-as.data.frame(as.matrix(ngsfilt)[grep('_b', as.matrix(ngsfilt)[,2]),]) #Keep only plates with replicate 'b'
ngsfilt_rep2<-rbind(ngsfilt_rep2[-c(1:25),], ngsfilt_rep2[c(1:25),])
NB_PLATE_rep2=5
ind_rep2={1:nrow(ngsfilt_rep2)}
tmp_rep2=seq(NB_PLATE_rep2*8, NB_PLATE_rep2*8+8*NB_PLATE_rep2, by=8)
x_rep2=NULL; for (i in 1:(length(tmp_rep2)-1)) {x_rep2=append(x_rep2,rep((tmp_rep2[i]+1):tmp_rep2[i+1],12))}
x_rep2<-x_rep2[1:nrow(ngsfilt_rep2)]
xLoc_rep2<-data.frame(as.vector(ngsfilt_rep2$V2), rep(rep(1:12, each=8), NB_PLATE_rep2)[c(1:length(as.vector(ngsfilt_rep2$V2)))])
colnames(xLoc_rep2)<-c('sample', 'column_in_plate')

ngsfilt_rep3<-as.data.frame(as.matrix(ngsfilt)[grep('_c', as.matrix(ngsfilt)[,2]),])
ngsfilt_rep3<-rbind(ngsfilt_rep3[-c(1:25),], ngsfilt_rep3[c(1:25),])
NB_PLATE_rep3=5
ind_rep3={1:nrow(ngsfilt_rep3)}
tmp_rep3=seq(((NB_PLATE_rep2*8)+(NB_PLATE_rep3*8)), ((NB_PLATE_rep2*8)+(NB_PLATE_rep3*8))+8*NB_PLATE_rep3, by=8)
x_rep3=NULL; for (i in 1:(length(tmp_rep3)-1)) {x_rep3=append(x_rep3,rep((tmp_rep3[i]+1):tmp_rep3[i+1],12))}
x_rep3<-x_rep3[1:nrow(ngsfilt_rep3)]
xLoc_rep3<-data.frame(as.vector(ngsfilt_rep3$V2), rep(rep(1:12, each=8), NB_PLATE_rep3)[c(1:length(as.vector(ngsfilt_rep3$V2)))])
colnames(xLoc_rep3)<-c('sample', 'column_in_plate')

xLoc_all<-rbind(xLoc_rep1, xLoc_rep2, xLoc_rep3)
yLoc_all<-c(x_rep1, x_rep2, x_rep3)
xyLoc_all<-cbind(xLoc_all, yLoc_all)

xyLoc_all2<-xyLoc_all[which(xyLoc_all$sample %in% rownames(OBI@samples)),]
rownames(xyLoc_all2)<-xyLoc_all2[,1]
xyLoc_all3<-xyLoc_all2[order(rownames(xyLoc_all2)),]

#For each samples: add locations into plates 
OBI@samples$xPlate<-xyLoc_all3$column_in_plate
OBI@samples$yPlate<-xyLoc_all3$yLoc_all



save.image(paste(OBJ, ".Rdata", sep=""))
```




#### Initial data Export ####
An intermediate, standardized table (obitools raw outputs nb of columns and taxo settings may depends on the user).
```{r, cache=T, autodep=TRUE}
#chunk7
tmp=t(OBI@reads)
colnames(tmp)=paste("sample:", colnames(tmp), sep="")
write.table(data.frame(OBI@motus,tmp), paste(OBJ, "R1.txt", sep=""), row.names=F, col.names=T, sep="\t")
```

#### Controls and experimental replicates identification ####
***!!! modify here according to your replicate/negative control denomination!!!***  
*e.g.: _a,_b,_c... or _1,_2,_3, or r when duplicated, T for negative controls...*

```{r, cache=T, autodep=TRUE}
#chunk8
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!!! MODIFY HERE !!!!
#Identify all controls in the dataset
EMPTY=grep("BLK", rownames(OBI@reads)) #black
PCR_CTL=grep("NEG", rownames(OBI@reads)) #grey
POS_CTL=grep("M", rownames(OBI@reads)) #red
CONTROLS=c(PCR_CTL, EMPTY, POS_CTL)
CONTROLS=c(PCR_CTL, EMPTY)

#Identify the replicates
REPLICATES=as.factor(gsub("a", "", rownames(OBI@samples))) # factor, each level is a sample
REPLICATES=gsub("b", "", REPLICATES)
REPLICATES=gsub("c", "", REPLICATES)

#color settings
COL="lightgrey" 
BOR="grey"
COL.s="cyan3"
BOR.s="cyan4"
COL.neg='red'
BOR.neg="darkred"
COL.all=rep("cyan3",nrow(OBI)); COL.all[CONTROLS]="red"
BOR.all=rep("cyan4",nrow(OBI)); BOR.all[CONTROLS]="darkred"
```

## 3/ Overall data quality ##

### Sequencing depth and Richness per sample ###
```{r, cache=T, autodep=TRUE, fig.width=10, fig.height=7, htmlcap='Fig. 2: Reads/OTUs per sample and localization in PCR plates (dot size log10 scaled). Controls in red'}

#chunk10
ngsfilt<-ngsfilt[which(rownames(ngsfilt) %in% rownames(OBI@reads)),]
par(mfrow=c(2,1))
barplot(rowSums(OBI@reads), col=COL.all, border=BOR.all, xlab="Samples", ylab="Nb reads")
#rownames(OBI@reads[order(ngsfilt$V2),])

#log10(rowSums(OBI@reads[POS_CTL,]))
hist(log10(rowSums(OBI@reads[-CONTROLS,])), breaks=40,  col=rgb(0,0,1,0.5), main="", xlab="log10 nb Reads", ylab="Nb samples", xlim=c(0,6))
hist(log10(rowSums(OBI@reads[CONTROLS,])), breaks=40,  col=rgb(1,0,0,0.5), add=T)
thresh.seqdepth=2 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!!! CUTOFF HERE !!!!
abline(v=thresh.seqdepth, col="black", lty=2, lwd=2); mtext(side=3, paste("cutoff < ", thresh.seqdepth, sep=""), cex=0.8, font=3)
par(mfrow=c(1,1))
plot(OBI@samples$xPlate, -OBI@samples$yPlate, pch=21, col=BOR.all, bg=COL.all, xlab="x plate", ylab="y plate", xlim=c(-5,17))
points(OBI@samples$xPlate[which(log10(rowSums(OBI@reads))<thresh.seqdepth)]+0.2,
       -OBI@samples$yPlate[which(log10(rowSums(OBI@reads))<thresh.seqdepth)], pch=8, cex=0.5, lwd=1.5, col='yellow4')
abline(h=-seq(8.5,128.5,8), lty=2, col="grey")

par(mfrow=c(2,1))
barplot(specnumber(OBI@reads, MARGIN=1), col=COL.all, border=BOR.all, xlab="Samples", ylab="Nb OTUs", las=2)
hist(log10(specnumber(OBI@reads)), breaks=40, col=COL, main="", xlab="log10 nb Reads", ylab="Nb samples")
thresh.rich=1.8 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!!! CUTOFF HERE !!!!
abline(v=thresh.rich, col="black", lty=2, lwd=2); mtext(side=3, paste("cutoff < ", thresh.rich, sep=""), cex=0.8, font=3)
par(mfrow=c(1,1))
plot(OBI@samples$xPlate, -OBI@samples$yPlate, pch=21, col=BOR.all, bg=COL.all, xlab="x plate", ylab="y plate", xlim=c(-5,17))
points(OBI@samples$xPlate[which(log10(specnumber(OBI@reads))<thresh.rich)]+0.2,
       -OBI@samples$yPlate[which(log10(specnumber(OBI@reads))<thresh.rich)],
       pch=8, cex=0.5, lwd=1.5, col='yellow4')
abline(h=-seq(8.5,128.5,8), lty=2, col="grey")

SMALL=rownames(OBI)[-CONTROLS][which(log10(rowSums(OBI@reads)[-CONTROLS])<thresh.seqdepth)]
```


### Taxonomic Resolution ###
```{r, cache=T, fig.width=10, fig.height=10}
#chunk11
###============= function TaxoRes
TaxoRes=function(x,y,z, thresh){
  #x=metabarcoding object
  #y=column name for taxonomic rank
  #z=column name for identification score
  #thresh=threshold below which sequences are considered as not really identified
  
  #inital settings
  require(ROBITools)
  #a vector encompassing all possible taxonomic levels
  taxorank=c("superkingdom", "kingdom", "subkingdom", "superphylum", "phylum", "subphylum", "superclass", "class", "subclass",
             "superorder", "order", "suborder", "infraorder", "superfamily", "family", "subfamily", "supertribe", "tribe", "subtribe",
             "supergenus", "genus", "subgenus", "superspecies", "species", "subspecies", "varietas", "no rank")
  
  #nb of otus
  tmp=table(x@motus[y])
  taxores.otu=tmp[match(taxorank, names(tmp))]
  names(taxores.otu)=taxorank
  taxores.otu[which(is.na(taxores.otu)==T)]=0
  
  #nb of reads
  tmp=aggregate(x@motus$count, by=list(x@motus[,y]), sum)
  taxores.reads=tmp[match(taxorank,tmp[,1]),2]
  names(taxores.reads)=taxorank
  taxores.reads[which(is.na(taxores.reads))]=0
  
  #set below thresh to not assigned
  tmp=x@motus[,y]
  tmp[which(OBI@motus[,z]<thresh)]="not assigned"
  
  #nb of reads above thresh
  tmp2=aggregate(OBI@motus$count, by=list(tmp), sum)
  taxores.reads.t=tmp2[match(c(taxorank, "not assigned"),tmp2[,1]),2]
  names(taxores.reads.t)=c(taxorank, "not assigned")
  taxores.reads.t[which(is.na(taxores.reads.t))]=0
  
  #nb of otus above thresh
  tmp2=table(tmp)
  taxores.otu.t=tmp2[match(c(taxorank, "not assigned"), names(tmp2))]
  names(taxores.otu.t)=c(taxorank, "not assigned")
  taxores.otu.t[which(is.na(taxores.otu.t))]=0
  
  layout(matrix(c(1,2,3,1,4,5),3,2),heights=c(0.4,1,1))
  col.tmp=c(rainbow(length(taxorank)-1,start=0, end=0.5, alpha=0.6), "lightgrey", "darkgrey")
  par(mar=c(1,1,1,1), oma=c(0,0,2,0))
  frame()
  legend("bottom", names(taxores.otu.t), ncol=5, cex=0.7, fill=col.tmp)
  pie(taxores.otu, col=col.tmp, border="lightgrey", labels="", clockwise=T)
  mtext("All data", side=2, cex=0.8)
  mtext(expression(OTUs[100]), side=3, cex=0.8)
  pie(taxores.otu.t, col=col.tmp, border="lightgrey", labels="", clockwise=T)
  mtext(paste("Best identities >", thresh) , side=2, cex=0.8)
  pie(taxores.reads, col=col.tmp, border="lightgrey", labels="", clockwise=T)
  mtext("Reads", side=3, cex=0.8)
  pie(taxores.reads.t, col=col.tmp, border="lightgrey", labels="", clockwise=T)
  
  out=data.frame(otu=c(taxores.otu,0), reads=c(taxores.reads,0), otu.thresh=taxores.otu.t, reads.thresh=taxores.reads.t)
  rownames(out)[length(taxorank)+1]="not assigned"
  out
}
###============= end

raw.taxores.all=TaxoRes(OBI,"rank_ok", "bid_ok", 0.95); title("All data", outer=T)
```

### Contamination level ###
#### True Contaminant ####
Identification: OTU detected with max frequency in negative controls + relatives (clustering if done or best_matchID)
```{r, cache=T, fig.width=10, fig.height=10}

#chunk12
###============= function ContaSlayer
ContaSlayer=function(x,y,clust=NULL){
  #x a metabarcoding object
  #y a vector of samples names where conta are suspected to be (typically negative controls)
  #clust a column name in x@motus which indicates the clustering level used for aggregating taxa (facultative)
  require(vegan)
  x.fcol=decostand(x@reads, "total", 2)
  x.max=rownames(x.fcol[apply(x.fcol, 2, which.max),])
  conta1=colnames(x)[which(is.na(match(x.max,y))==F)]
  if (length(clust)!=0) {
    conta2=c(conta1[which(is.na(x@motus[conta1, clust])==T)],
             rownames(x@motus)[which(is.na(match(x@motus[,clust],unique(x@motus[conta1, clust][-which(is.na(x@motus[conta1,clust])==T)])))==F)])
  } else {
    conta2=conta1
  }
  conta2
}  
###============= end

CONTA=ContaSlayer(OBI,rownames(OBI)[CONTROLS], NULL)


CONTA=as.vector(unlist(MyList))
length(CONTA)
dim(OBI@reads[-CONTROLS,][,which(colnames(OBI@reads[-CONTROLS,]) %in% CONTA)])
colSums(OBI@reads[-CONTROLS,][,which(colnames(OBI@reads[-CONTROLS,]) %in% CONTA)])
max(colSums(OBI@reads[-CONTROLS,][,which(colnames(OBI@reads[-CONTROLS,]) %in% CONTA)]))

pandoc.table(data.frame(OBI@motus[CONTA,c(grep("best_identity", colnames(OBI@motus)),grep("best_match", colnames(OBI@motus)))],
                        OBI@motus[CONTA,c("bid_ok", "scientific_name_ok","count")],
                        do.call("rbind", lapply(CONTA, function(x) {
                          ind=which.max(OBI@reads[,x])
                          sample.max=names(ind)
                          count.max=OBI@reads[sample.max,x]
                          data.frame(sample.max, count.max)
                          })),
                        nb.sample.occ=unlist(lapply(CONTA, function(x) length(which(OBI@reads[,x]!=0))))
                        ),split.tables=Inf)

save.image('DataSavedAfterContaminantFiltering.Rdata')
```

Freq of contaminant in samples
```{r, cache=T, autodep=TRUE,  fig.width=15, fig.height=10}
#chunk13
OBI.freq=decostand(OBI@reads,"total",1)
OBI.freq.parse=c(seq(0,round(nrow(OBI@samples), digit=0),round(nrow(OBI@samples)/4, digit=0)), nrow(OBI@samples))

#Data visualization: not mandatory
par(mfrow=c(2,2), oma=c(0,0,2,0), mar=c(3,5,1,1))
for(i in 1:(length(OBI.freq.parse-1)-1)){
  image(-log10(OBI.freq[(OBI.freq.parse[i]+1):(OBI.freq.parse[i+1]),CONTA]), xaxt="n", yaxt="n")
  abline(v=seq(0,1,l=nrow(OBI@samples)/4)
         [grep("C",rownames(OBI.freq)[(OBI.freq.parse[i]+1):(OBI.freq.parse[i]+round(nrow(OBI.freq)/4, digit=0))])],
         col=COL.neg, lty=3)
  axis(side=1,at=seq(0,1,l=round(nrow(OBI@samples)/4, digit=0)),
       labels=rownames(OBI@samples)[(OBI.freq.parse[i]+1):(OBI.freq.parse[i]+round(nrow(OBI@samples)/4, digit=0))],
       las=2, cex.axis=0.3)
  axis(side=2, at=seq(0,1,l=length(CONTA)), labels=paste(CONTA, OBI@motus[CONTA,"sci_name_ok"]),
       cex.axis=0.3, las=2)
}
mtext(side=3, paste("Nb. conta (max count in negative controls + related sequences (cluster and or best match)) = ", length(CONTA)), outer=T, cex=0.7, font=3)
```


Heavily contaminated samples
```{r,  cache=T, autodep=TRUE,  fig.width=10, fig.height=5}
#chunk14
par(mar=c(4,4,1,1), mfrow=c(1,2))
hist(log10(rowSums(OBI.freq[-CONTROLS,CONTA])+1e-5), breaks=40, col=COL, xlim=log10(c(0,1)+1e-5), main="", xlab="log10 Prop conta reads (non filtered data)")

#most likely to be bad samples
thresh.contasamp=-2
CONTASAMP=rownames(OBI)[-CONTROLS][which(log10(rowSums(OBI.freq[-CONTROLS,CONTA])+1e-5)>=thresh.contasamp)]
abline(v=thresh.contasamp, lty=2, lwd=2); mtext(paste("cutoff >= ", thresh.contasamp, sep=""),side=3, line=0, cex=0.7, font=3)

plot(OBI@samples$xPlate, -OBI@samples$yPlate, pch=21, col=BOR.all, bg=COL.all, cex=1, xlab="x plate", ylab="y plate", xlim=c(-5,17))
abline(h=-seq(8.5,128.5,8), lty=2, col="grey")
points(OBI@samples$xPlate[match(CONTASAMP,rownames(OBI))]+0.2, -OBI@samples$yPlate[match(CONTASAMP,rownames(OBI))], pch=8, cex=0.7, lwd=2)
```

#### Odd OTUs with low identification scores ####
Identification
```{r,  cache=T, autodep=TRUE,  fig.width=10, fig.height=7}
#chunk15
par(mfrow=c(1,2), oma=c(0,0,2,0))
weighted.hist(OBI@motus[, "bid_ok"],OBI@motus[, "count"],breaks=20,col=COL, ylab="Nb reads", xlab="Ecotag scores")
thresh.oddotu=0.95 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!!! CUTOFF HERE !!!!
abline(v=thresh.oddotu, lty=2, lwd=2)
hist(OBI@motus[, "bid_ok"], breaks=40, col=COL, xlab="Ecotag scores", main="", ylab="Nb. OTUs")
abline(v=thresh.oddotu, lty=2, lwd=2)
mtext(side=3, paste("cutoff < ", thresh.oddotu, sep=""), cex=0.7, font=3, outer=T, line=1)

tmp=rownames(OBI@motus)[which(OBI@motus[, "bid_ok"]<thresh.oddotu)]
ODDOTU=if(length(grep("cluster", rownames(OBI@motus)))!=0) {
  c(tmp[which(is.na(OBI2@motus[tmp, "cluster95"])==T)],
  rownames(OBI2@motus)[which(is.na(match(OBI2@motus[,"cluster95"],unique(OBI2@motus[tmp, "cluster95"][-which(is.na(OBI2@motus[tmp,"cluster95"])==T)])))==F)])
  } else {
  tmp
}

mtext(side=3, paste("nb seq. with score < cutoff = ", length(ODDOTU)), cex=0.7, font=3, outer=T, line=0)

pandoc.table(data.frame(OBI@motus[ODDOTU,c(grep("best_identity", colnames(OBI@motus)),grep("best_match", colnames(OBI@motus)))],
##Viasualize data: not mandatory
OBI@motus[ODDOTU,c("bid_ok","scientific_name_ok","count")],
                        do.call("rbind", lapply(ODDOTU, function(x) {
                          ind=which.max(OBI@reads[,x])
                          sample.max=names(ind)
                          count.max=OBI@reads[sample.max,x]
                          data.frame(sample.max, count.max)
                        })),
                        nb.sample.occ=unlist(lapply(ODDOTU, function(x) length(which(OBI@reads[,x]!=0)))))
             [order(ODDOTU),],split.tables=Inf, split.cells=20)


dim(OBI@reads[-CONTROLS,][,which(colnames(OBI@reads[-CONTROLS,]) %in% ODDOTU)])
colSums(OBI@reads[-CONTROLS,][,which(colnames(OBI@reads[-CONTROLS,]) %in% ODDOTU)])
max(colSums(OBI@reads[-CONTROLS,][,which(colnames(OBI@reads[-CONTROLS,]) %in% ODDOTU)]))
```

Distribution of odd OTU in samples
```{r, cache=T, autodep=TRUE,  fig.width=10, fig.height=4}
#chunk16
par(mfrow=c(2,2), oma=c(0,0,2,0), mar=c(3,5,1,1))
for(i in 1:(length(OBI.freq.parse)-1)){
  image(log10(OBI.freq[(OBI.freq.parse[i]+1):(OBI.freq.parse[i+1]),ODDOTU]), xaxt="n", yaxt="n")
  abline(v=seq(0,1,l=nrow(OBI@samples)/4)
         [grep("C",rownames(OBI.freq)[(OBI.freq.parse[i]+1):(OBI.freq.parse[i]+round(nrow(OBI.freq)/4, digit=0))])],
         col=COL.neg, lty=3)
  axis(side=1,at=seq(0,1,l=round(nrow(OBI@samples)/4), digits=0),
       labels=rownames(OBI.freq)[(OBI.freq.parse[i]+1):(OBI.freq.parse[i]+round(nrow(OBI.freq)/4, digit=0))],
       las=2, cex.axis=0.3)
  axis(side=2, at=seq(0,1,l=length(ODDOTU)), labels= OBI@motus[ODDOTU,"sci_name_ok"],cex.axis=0.3, las=2)
}
mtext(side=3, paste("Nb. oddotu (max count in seq with low exotag score = ", length(ODDOTU)), outer=T, cex=0.7, font=3)
```


Sample with many odd OTUs
```{r,  cache=T, fig.width=10, fig.height=5}
#chunk17
par(mar=c(4,4,1,1), mfrow=c(1,2))

hist(log10(rowSums(OBI.freq[-CONTROLS,ODDOTU])+1e-5), breaks=20, col=COL, xlim=log10(c(0,1)+1e-5), main="", xlab="Prop conta2 reads")
#most likely to be bad samples
thresh.oddotusamp=-2  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!!! CUTOFF HERE !!!!
ODDOTUSAMP=rownames(OBI)[which(log10(rowSums(OBI.freq[,ODDOTU])+1e-5)>=thresh.oddotusamp)]
abline(v=thresh.oddotusamp, lty=2, lwd=2); mtext(paste("cutoff >= ", thresh.oddotusamp, sep=""),side=3, line=0, cex=0.7, font=3)

plot(OBI@samples$xPlate, -OBI@samples$yPlate, pch=21, col=BOR.all, bg=COL.all, cex=1, xlab="x plate", ylab="y plate", xlim=c(-5,17))
abline(h=-seq(8.5,95.5,8), lty=2, col="grey")
points(OBI@samples$xPlate[match(ODDOTUSAMP,rownames(OBI))]+0.2, -OBI@samples$yPlate[match(ODDOTUSAMP,rownames(OBI))], pch=8, cex=0.7, lwd=2)
```








### Data filtering 1 ###
Get rid off identified contaminant and odd OTUs and effects on dataset sequencing depth, taxo resolution
```{r, cache=T, fig.width=10, fig.height=5}
#chunk18
OBI2<-OBI
OBI2@reads=OBI@reads[,-match(unique(c(CONTA, ODDOTU)), colnames(OBI@reads))]
OBI2@motus=OBI2@motus[which(rownames(OBI2@motus) %in% colnames(OBI2@reads)),]
OBI2@samples=OBI2@samples[which(rownames(OBI2@samples) %in% rownames(OBI2@reads)),]
OBI2@samples$reads.contaodd=rowSums(OBI2@reads)
OBI2@samples$otu100.contaodd=specnumber(OBI2@reads)
OBI2@motus$count=colSums(OBI2@reads)


#Checks
#summary.metabarcoding(OBI2)
#filcontaodd.taxores.all=TaxoRes(OBI2,"rank_ok", "bid_ok", 0.95); title("Data -Conta/odd OTUs", outer=T)
```



STOP HERE!!!!! THEN WE WILL remove replicates with weird behaviour






















Similarity behaviour
```{r, cache=T,fig.width=10, fig.height=7}
#chunk19
###============= function vegsim
vegsim=function(x) {
  x.dist=vegdist(x, "bray")
  x.sim=1-as.matrix(x.dist)
  out=as.dist(x.sim, diag=F, upper=F)
  out
}
###============= end
OBI.sim=vegsim(OBI@reads); OBI.sim2=vegsim(OBI2@reads)

###============= function simlist: convert to similarity list 3 columns format
simlist=function(x) {
  out=data.frame(t(combn(labels(x),2)), as.numeric(x))
  names(out) <- c("obj1", "obj2", "sim")
  out
}
###============= end
OBI.simlist=simlist(OBI.sim); OBI.simlist2=simlist(OBI.sim2)


###============= functiont comp.type: annotates pairwise sim per comparison type (e.g. sample vs control, replicates vs replicates, etc..)
comp.type=function(x,y,z1,z2){
  #x a distance list table
  #y a metabarcoding object
  #z1 a vector indicating negative indices in the metabarcoding object
  #z2 a factor of sample names, repeats are indicated with the same sample name
  #controls intravar
  x$type=rep(NA, nrow(x))
  ind.negneg=intersect(which(is.na(match(as.vector(x[,1]), rownames(y)[z1]))==F), which(is.na(match(as.vector(x[,2]), rownames(y)[z1]))==F))
  x[ind.negneg,"type"]="neg.neg"
  #controls vs sample var
  ind.negsamp=c(intersect(which(is.na(match(as.vector(x[,1]), rownames(y)[z1]))==F), which(is.na(match(as.vector(x[,2]), rownames(y)[z1]))==T)),
                intersect(which(is.na(match(as.vector(x[,2]), rownames(y)[z1]))==F), which(is.na(match(as.vector(x[,1]), rownames(y)[z1]))==T)))
  x[ind.negsamp,"type"]="neg.samp"
  #replicates intravar
  repltech=names(table(z2))[which(table(z2)>1)]
  ind.repin=unlist(lapply(repltech, function(a) intersect(grep(a,x[,1]), grep(a, x[,2]))))
  x[ind.repin, "type"]="rep.in"
  #sample vs sample var
  x[-c(ind.negneg, ind.negsamp, ind.repin),"type"]="samp.samp"
  x
}

OBI.compsim=comp.type(OBI.simlist, OBI, CONTROLS, REPLICATES)
OBI.compsim2=comp.type(OBI.simlist2, OBI, CONTROLS, REPLICATES)

par(mfrow=c(1,2), mar=c(4,4,2,1))
boxplot(OBI.compsim$sim~OBI.compsim$type, col=COL, border=BOR, main="raw data", las=2, ylab="Pairwise sim")
boxplot(OBI.compsim2$sim~OBI.compsim2$type, col=COL, border=BOR, main="post- conta/odd OTU removal", las=2,  ylab="Pairwise sim")

save.image('PostCompSim.Rdata')
```


Additional filtering steps
```{r, cache=T, fig.width=10, fig.height=10, results="hide", warning=F}
#chunk20
#Examination of similarities behaviour with various filtering thresholds
###============= function DataSmooth a function that smooth reads count for both basal noise (typically cross sample chimeras) and detection level (min abundance taken into account)
###Quite long to run, not mandatory but can help defining detection threshold

DataSmooth=function(X,noise.cut,detect.cut){
  #OBI=a metabarcoding object
  #noise.cut= a vector containing thresholds for removing noise (0 to 1)
  #detect.cut= a vector containing minimum abundance values to be considered (otherwise set to 0)
  #returns a list of reads table
  X.cut=lapply(noise.cut, function(x){
    if(x==0) {X.x=X} else {X.x=threshold.set(X,2,x)}
    lapply(detect.cut, function(y) {
      if(y==0) {X.x@reads} else {apply(X.x@reads, 2, function(vec) {
        vec[vec<y]=0; vec
      })}
    })
  })
  names(X.cut)=paste("noise_",noise.cut, sep="")
  X.cut=unlist(lapply(X.cut, function(x) {names(x)=paste("detect_", detect.cut, sep=""); x}), recursive=F)
  X.cut=lapply(X.cut, function(x) x[,which(colSums(x)!=0)])
  X.cut
}
###=============end

noise.cut=c(0,0.99,0.97,0.95)
detect.cut=c(0,2,5,10,50)
OBI.cut2=DataSmooth(OBI2, noise.cut, detect.cut)
OBI.cut2.sim=lapply(OBI.cut2, function(x) simlist(vegsim(x)))
```


```{r, cache=T, results="hide", fig.width=15, fig.height=15}
#chunk21
###Visualization of the precedent step
par(mfrow=c(length(noise.cut), length(detect.cut)), mar=c(3,3,3,1))
lapply(1:length(OBI.cut2.sim), function(x)
  hist(OBI.cut2.sim[[x]][,"sim"], col=COL, border=BOR, main=names(OBI.cut2.sim)[x], xlab=" Sample Sim", cex.main=0.7)
)
par(mfrow=c(length(noise.cut), length(detect.cut)), mar=c(4,4,2,1))
lapply(1:length(OBI.cut2.sim), function(x)
  hist(log10(OBI.cut2.sim[[x]][,"sim"]), col=COL, border=BOR, main=names(OBI.cut2.sim)[x], xlab="log10 Sample Sim", cex.main=0.7)
)
par(mfrow=c(length(noise.cut), length(detect.cut)), mar=c(4,4,2,1))
lapply(1:length(OBI.cut2), function(x) boxplot(OBI.cut2.sim[[x]]$sim~OBI.compsim2$type, col=COL, border=BOR, main=names(OBI.cut2.sim)[x], las=2))
```


### Putative spurious PCRs ###
Detect outliers: Step 1
```{r, cache=T}
#chunk22
###============= simfun a function that aggregates similarities per sample with a given function
simfun=function(x,method){
  #x a simlist object
  #method a function to be applied to each sample (e.g. sum, mean)
  x.new=data.frame(sample=c(as.vector(x[,"obj1"]),as.vector(x[,"obj2"])), sim=rep(x[,"sim"],2))
  x.new$method=ave(x.new[,"sim"], x.new[, "sample"], FUN=method)
  x.new[match(unique(x.new[,"sample"]),x.new[,"sample"]),]
}
###============= end
OBI.simmean2=simfun(OBI.simlist2, mean); 
OBI.cut2.simmean=lapply(OBI.cut2.sim, simfun, mean)
```


Detect outliers: Step 2 (Option 1) => based on averaged similarities
```{r, cache=T, fig.width=15, fig.height=15, results="hide"}
#chunk23
par(mfrow=c(length(noise.cut), length(detect.cut)), mar=c(4,4,2,1))
lapply(1:length(OBI.cut2.simmean), function(x)
  hist(log10(OBI.cut2.simmean[[x]][,"sim"]), col=COL, border=BOR, main=names(OBI.cut2.simmean)[x], xlab="log10 Sample Sim", cex.main=0.7)
)
```

```{r, cache=T, fig.width=10, fig.height=7, results="hide"}
#chunk24
###Visualization of the precedent step
par(mfrow=c(1,2), oma=c(0,0,2,0))
hist(log10(OBI.simmean2[,"method"]), col=COL, border=BOR, main="", xlab="log10 Sample Avg Sim")
hist(OBI.simmean2[,"method"], col=COL, border=BOR, main="", xlab="log10 Sample Avg Sim")
thresh.outsim=-1.5  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!!! CUTOFF HERE !!!!
abline(v=thresh.outsim, lty=2, lwd=2); mtext(side=3, paste("cutoff > ", thresh.outsim), cex=0.7, font=3, outer=F)
OUTSIM=rownames(OBI)[which(log10(OBI.simmean2$method)<thresh.outsim)]
plot(OBI@samples$xPlate, -OBI@samples$yPlate, pch=21, col=BOR.all, bg=COL.all, cex=1, xlab="x plate", ylab="y plate", xlim=c(-5,17))
abline(h=-seq(8.5,95.5,8), lty=2, col="grey")
points(OBI@samples$xPlate[match(OUTSIM,rownames(OBI))]+0.2, -OBI@samples$yPlate[match(OUTSIM,rownames(OBI))], pch=8, cex=0.7, lwd=2)
```



Detect outliers: Step 2 (Option 2) => Using a clustering approach (recommended!)
```{r, cache=T, fig.width=10, fig.height=10}
#chunk25
###============= function quickg: graph construction
quickg=function(x,col, nod.fac, thresh){
  x.m=as.matrix(x)
  g=graph.adjacency(x.m, mode="lower", weighted=T, diag=F)
  V(g)$type=nod.fac
  V(g)$color=as.vector(col)
  sub.g=subgraph.edges(g, eids=E(g)[which(E(g)$weight>thresh)], delete.vertices=F)
  sub.g
}
###============= end

layout(matrix(c(1,2,2,3,2,2), 3,2))
hist(log10(OBI.simlist2[,"sim"]), col=COL, border=BOR, main="", xlab="Sample Sim")

#Get 2 lists: one for similarity among replicates, and one among samples
MyList_IntraSamp<-list()
kk<-1
MyList_InterSamp<-list()
ll<-1
for (i in 1:nrow(OBI.simlist2)) {
  if( substr(OBI.simlist2[i,1],1,10) == substr(OBI.simlist2[i,2],1,10)) {
    MyList_IntraSamp[[i]]<-OBI.simlist2[i,3]
    kk<-kk+1 
  } else {
    MyList_InterSamp[[i]]<-OBI.simlist2[i,3]
    ll<-ll+1
  }
}
MyList_IntraSamp_OK<-unlist(MyList_IntraSamp)
MyList_InterSamp_OK<-unlist(MyList_InterSamp)

#Visualize similarity scores intra and inter-sample (log 10 Sample sim)
hist(log10(MyList_InterSamp_OK),  col=rgb(0,0,1,0.5), main="", xlab="log10 Sample sim", ylab="Nb comparison", xlim=c(-5,0))
hist(log10(MyList_IntraSamp_OK),  col=rgb(1,0,0,0.5), add=T)

#Visualize similarity scores intra and inter-sample
hist(MyList_IntraSamp_OK, col=COL, border=BOR, main="", xlab="Sample Sim", breaks=40, xlim=c(0,1), freq=T)
thresh.graph=0.1  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!!! CUTOFF HERE !!!!
abline(v=thresh.graph, lty=2, lwd=2); mtext(side=3, paste("cutoff > ", thresh.graph), cex=0.7, font=3, outer=F)
hist(MyList_InterSamp_OK, col=COL, border=BOR, main="", xlab="Sample Sim", xlim=c(0,1), breaks=20, freq=T)
thresh.graph=0.1  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!!! CUTOFF HERE !!!!
abline(v=thresh.graph, lty=2, lwd=2); mtext(side=3, paste("cutoff > ", thresh.graph), cex=0.7, font=3, outer=F)

##Generate cluster
OBI.g=quickg(OBI.sim, COL.all, COL.all, thresh.graph)
OBI.g2=quickg(OBI.sim2, COL.all, COL.all, thresh.graph)
OBI.clust2=infomap.community(OBI.g2, nb.trials=100)

##Visualize clusters: not mandatory
col.g=as.factor(membership(OBI.clust2))
levels(col.g)=rainbow(max(membership(OBI.clust2))+10)[1:max(membership(OBI.clust2))][sample(1:max(membership(OBI.clust2)))]
coords=layout.fruchterman.reingold(OBI.g2)
plot(OBI.g2, layout=coords, vertex.size=2, vertex.label=rownames(OBI2), vertex.color=as.vector(col.g), vertex.label.cex=0.7, vertex.label.dist=-0.2, vertex.frame.color=as.vector(col.g))
plot(OBI@samples$xPlate, -OBI@samples$yPlate, pch=1, col=BOR.all, bg=COL.all, cex=1, xlab="x plate", ylab="y plate", xlim=c(-5,17))
abline(h=-seq(8.5,95.5,8), lty=2, col="grey")
points(OBI@samples$xPlate, -OBI@samples$yPlate, pch=8, cex=0.5, col=as.vector(col.g), lwd=2)
legend("left",paste("cluster", levels(as.factor(membership(OBI.clust2)))), pch=8, col=levels(col.g), cex=0.5, ncol=3)


######## Detect spurious samples + outliers detection ######## 

# 1) Detect samples that are clustered with TRUE control samples (Extration + PCR)
TRUE_CONTROLS<-rownames(OBI2)[CONTROLS][-grep('BLK', rownames(OBI2)[CONTROLS])]
tmp1=names(unlist(lapply(unique(membership(OBI.clust2)[TRUE_CONTROLS]),function(x)which(membership(OBI.clust2)[-CONTROLS]==x))))
unique(str_replace_all(tmp1, "[abc]", ""))
length(table(str_replace_all(tmp1, "[abc]", ""))[table(str_replace_all(tmp1, "[abc]", ""))>1])


# 2) Detect samples that are clustered with nothing else
thresh.outgraph=3  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!!! CUTOFF HERE !!!! Threshold to detect samples belonging to cluster containing less than 1 x samples
tmp2=unlist(lapply(names(which(table(membership(OBI.clust2))<=thresh.outgraph)), 
                    function(x) names(which(membership(OBI.clust2)==x))))
unique(str_replace_all(tmp2, "[abc]", ""))
table(str_replace_all(tmp2, "[abc]", ""))


# 3) Detect replicates that are not clustered with their relatives
tmp3<-list()
for (i in 1:length(unique(REPLICATES))) {
   xx<-membership(OBI.clust2)[which(str_replace_all(rownames(OBI@samples), "[abc]", "")==unique(REPLICATES)[i])]
  ##print(xx)
  if (length(table(xx))==1) {
    tmp3[[i]]<-NA 
  }
  if (length(table(xx))==2) {
    titi<-min(as.vector(table(xx)))
    tutu<-names(table(xx)[table(xx)==titi])
    tmp3[[i]]<-names(xx[xx==tutu])
  }
  if (length(table(xx))==3) {
    tmp3[[i]]<-names(xx)
  }
}

tmp3<-unlist(tmp3)
tmp3<-tmp3[!is.na(tmp3)]

tmp3
unique(str_replace_all(tmp3, "[abc]", ""))
table(str_replace_all(tmp3, "[abc]", ""))


###Choose the filtration level to use: I often either use option a) or d) 
#a) OUTGRAPH_level1_tmp1 => samples that are not clustered with their relatives
OUTGRAPH_level1_tmp1=unique(c(tmp3[-which(tmp3 %in% rownames(OBI2)[CONTROLS]==T)]))

#b) OUTGRAPH_level2_tmp1_tmp2 => samples that are not clustered with their relatives + samples that are clustered with nothing else
OUTGRAPH_level2_tmp1_tmp2=unique(c(tmp3[-which(tmp3 %in% rownames(OBI2)[CONTROLS]==T)],tmp2[-which(tmp2 %in% rownames(OBI2)[CONTROLS]==T)]))

#c) OUTGRAPH_level3_tmp1_tmp2_tmp3 => samples that are not clustered with their relatives + samples that are clustered with nothing else + samples that are clustered with TRUE control samples
OUTGRAPH_level3_tmp1_tmp2_tmp3=unique(c(tmp3[-which(tmp3 %in% rownames(OBI2)[CONTROLS]==T)],tmp2[-which(tmp2 %in% rownames(OBI2)[CONTROLS]==T)], tmp1))

#d) OUTGRAPH_level2_Alt_tmp1_tmp3 => samples that are not clustered with their relatives samples that are clustered with TRUE control samples 
OUTGRAPH_level2_Alt_tmp1_tmp3=unique(c(tmp3[-which(tmp3 %in% rownames(OBI2)[CONTROLS]==T)], tmp1))
```


### Final checks and formatting ###
Look at most abundant OTUS in samples identified as "odd"
```{r, cache=T}
#Visualize where are the samples filtered out at the precedent step
OUTGRAPH<-OUTGRAPH_level2_Alt_tmp1_tmp3 #Choose your filtering level here

#chunk26
pandoc.table(data.frame(
  OBI@motus[unlist(lapply(OUTGRAPH, function(x) names(sort(OBI@reads[x,], decreasing=T)[1:10]))),
            c(grep("best_identity", colnames(OBI@motus)),grep("best_matchID", colnames(OBI@motus)))],
  OBI@motus[unlist(lapply(OUTGRAPH, function(x) names(sort(OBI@reads[x,], decreasing=T)[1:10]))),
            c("bid_ok","scientific_name","count")],
  sample=rep(OUTGRAPH, each=10),
  is.contaodd=unlist(lapply(OUTGRAPH, function(x) names(sort(OBI@reads[x,], decreasing=T)[1:10]))) %in% c(CONTA, ODDOTU)),split.tables=Inf, split.cells=20)
```

Data formatting (removal of outliers + controls)
```{r, fig.width=10, fig.height=10}
OBI3=OBI2
#Removal of all "bad" samples: OUTGRAPH (again choose the level of filtering) + those with low number of reads
OBI3@reads=OBI3@reads[-c(match(OUTGRAPH_level2_Alt_tmp1_tmp3, rownames(OBI2@reads)),CONTROLS, match(SMALL, rownames(OBI2@reads))),]
OBI3@motus=OBI3@motus[which(rownames(OBI3@motus) %in% colnames(OBI3@reads)),]
OBI3@samples=OBI3@samples[which(rownames(OBI3@samples) %in% rownames(OBI3@reads)),]
OBI3@motus$count=colSums(OBI2@reads)
```

Data export
```{r}
tmp=t(OBI3@reads)
colnames(tmp)=paste("sample:", colnames(tmp), sep="")
if (any(colnames(OBI3@motus)=='species_list')==T) {
  OBI3@motus<-OBI3@motus[,-which(colnames(OBI3@motus) %in% c('species_list', 'count'))]
}
write.table(data.frame(OBI3@motus,tmp), "SER2017_Diet_CleanDataset_WithLocalDB_TrimLevel2_Alt_noAggregationPerSample.txt", row.names=F, col.names=T, sep="\t")

save.image('Filtering_Serengeti2017_withLocalDB_data.Rdata')

###Stop here if you do not want to aggregate replicates into a single sample
```









**Data importation for aggregating replicate**
```{r}
OBI_Trim<-import.metabarcoding.data('/Users/johan/Desktop/Princeton_projects/Serengeti/SER2017_Diet_CleanDataset_WithLocalDB_TrimLevel2_Alt_noAggregationPerSample.txt')  ##Be sure that sample names start by 'sample:' (and not 'sample.')


#Check the number of samples remaining after filtering
unique(str_replace_all(rownames(OBI_Trim@reads), "[abc]", ""))
length(unique(str_replace_all(rownames(OBI_Trim@reads), "[abc]", "")))
table(substr(unique(str_replace_all(rownames(OBI_Trim@reads), "[abc]", "")), 1, 12))
min(rowSums(OBI_Trim@reads))

##Remove Empty rows and columns in reads and motus tables
if(any(rowSums(OBI_Trim@reads)==0)) {
    OBI_Trim@reads<-OBI_Trim@reads[-which(rowSums(OBI_Trim@reads)==0),]
}
if(any(colSums(OBI_Trim@reads)==0)) {
    OBI_Trim@reads<-OBI_Trim@reads[,-which(colSums(OBI_Trim@reads)==0)]
}
OBI_Trim@motus<-OBI_Trim@motus[which(rownames(OBI_Trim@motus) %in% colnames(OBI_Trim@reads)),]

##Recalculate count and occurence of each sequences in the new dataset
OBI_Trim@motus$count<-colSums(OBI_Trim@reads)
Occurence<-specnumber(t(OBI_Trim@reads))
OBI_Trim@motus<-data.frame(OBI_Trim@motus, Occurence)
```


**Select filtering option and parameters**
```{r}
#Chunk 3

# FILTERING 1
SmoothingData<-F #smooth reads count for both basal noise (typically cross sample chimeras) and detection level (min abundance taken into account)
##  iF SmoothingData==T, specify also:
noise.cut=1 #a vector containing thresholds for removing noise (0 to 1)
detect.cut=0 #a vector containing minimum abundance values to be considered (otherwise set to 0)

# FILTERING 2
Maj_seq<-F  #for each sample, keep only sequences retrieved in at least 2 replicates: TRUE or FALSE

# FILTERING 3
Multi_Rep<-T #Have replicates of each samples been included in the analysis?
#If multiple replicates per sample have been included in the analysis, this filtering allows to aggregate replicates from a given sample. Two option for the aggregation: mean of relative read abundance (mean_RRA, data will be analyzed based on rra) or mean of reads (mean_reads, data can be analyzed based on other metrics as rarefied data).
Agg_Rep<-'mean_reads' #Can be either: mean_RRA or mean_reads
```


**Filtering option**
```{r}
#Chunk 4

####OPTIONAL FILTERING 1####
###============= function DataSmooth a function that smooth reads count for both basal noise (typically cross sample chimeras) and detection level (min abundance taken into account)
DataSmooth=function(X,noise.cut,detect.cut){
  #OBI=a metabarcoding object
  #noise.cut= a vector containing thresholds for removing noise (0 to 1)
  #detect.cut= a vector containing minimum abundance values to be considered (otherwise set to 0)
  #returns a list of reads table
  X.cut=lapply(noise.cut, function(x){
    if(x==0) {X.x=X} else {X.x=threshold.set(X,2,x)}
    lapply(detect.cut, function(y) {
      if(y==0) {X.x@reads} else {apply(X.x@reads, 2, function(vec) {
        vec[vec<y]=0; vec
      })}
    })
  })
  names(X.cut)=paste("noise_",noise.cut, sep="")
  X.cut=unlist(lapply(X.cut, function(x) {names(x)=paste("detect_", detect.cut, sep=""); x}), recursive=F)
  X.cut=lapply(X.cut, function(x) x[,which(colSums(x)!=0)])
  X.cut
}
###=============end

if (SmoothingData==T) {
  OBI.cut=DataSmooth(OBI_Trim, noise.cut, detect.cut)
  OBI_Trim@reads<-OBI.cut[[1]]
} 
####END OF OPTIONAL FILTERING 1####


####OPTIONAL FILTERING 2####
#### For each sample, keep only sequences retrieved in the majority of replicates (i.e. at least 2 replicates)
if (Maj_seq==T) {
  DataFormat<-OBI_Trim@reads
  SAMPLES=unique(substr(rownames(DataFormat), 1,10))
  for (i in 1:ncol(OBI_Trim@reads)) {
    for (j in 1:length(SAMPLES)) {
      xx<-OBI_Trim@reads[which(substr(rownames(OBI_Trim@reads), 1,10) ==SAMPLES[j]), i]
      xx2<-decostand(xx, 'pa')
      if (sum(xx2)<2) {
        DataFormat[which(substr(rownames(DataFormat), 1,10) ==SAMPLES[j]), i]<-0
      } else {
        DataFormat[which(substr(rownames(DataFormat), 1,10) == SAMPLES[j]), i]<-xx
      }
    }
  }
  if(any(rowSums(DataFormat)==0)) {
    DataFormat<-DataFormat[-which(rowSums(DataFormat)==0),]
  }
  if(any(colSums(DataFormat)==0)) {
    DataFormat<-DataFormat[,-which(colSums(DataFormat)==0)]
  }
  OBI_Trim@reads<-DataFormat
}
#### END OF THE OPTIONAL FILTERING 2


####OPTIONAL FILTERING 3####
if (Multi_Rep==T) {
  
  if (Agg_Rep=='mean_RRA') {
  #### Aggregate replicates fronm the same sample: mean of rra (to analyze data   from on relative read abundance)
    MyData<-decostand(OBI_Trim@reads, 'total', MARGIN=1)
    MyData2<-aggregate(MyData, by=list(substr(rownames(MyData),1,10)), FUN=mean)
    rownames(MyData2)<-MyData2[,1]
    MyData2<-MyData2[,-1]
    OBI_Trim@reads=as.matrix(MyData2)
  }
  ####
  
  if (Agg_Rep=='mean_reads') {
  #### Aggregate replicates from the same sample: mean of sample reads (to analyzed data from rarefied data)
    MyData<-OBI_Trim@reads
    MyData2<-aggregate(MyData, by=list(substr(rownames(MyData),1,10)), FUN=mean)
    rownames(MyData2)<-MyData2[,1]
    MyData2<-MyData2[,-1]
    OBI_Trim@reads=as.matrix(MyData2)
  }
}
#### END OF THE OPTIONAL FILTERING 3

```


#Remove sequences that represents less than 1% of the individual sample
OBI_Trim@reads[decostand(OBI_Trim@reads, 'total')<0.01]<-0
if (any(colSums(OBI_Trim@reads)==0)) {
  OBI_Trim@reads<-OBI_Trim@reads[,-which(colSums(OBI_Trim@reads)==0)]
}
if (any(rowSums(OBI_Trim@reads)==0)) {
  OBI_Trim@reads<-OBI_Trim@reads[-which(rowSums(OBI_Trim@reads)==0),]
}
dim(OBI_Trim@reads)
OBI_Trim@motus<-OBI_Trim@motus[which(rownames(OBI_Trim@motus) %in% colnames(OBI_Trim@reads)),]


**Generate the final table**
```{r}
##Note: Data can be rarefied at this stage if needed

TMP<-t(OBI_Trim@reads)
write.table(data.frame(OBI_Trim@motus,TMP), "SER2017_withLocalDB_Diet_Filtered.tab", row.names=F, col.names=T, sep="\t")
```
