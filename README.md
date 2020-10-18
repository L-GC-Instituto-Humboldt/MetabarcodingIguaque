# Metabarcoding assessment of high tropical Andes diversity along an elevational gradient

## Laboratorio de Genética de la Conservación, Instituto Humboldt

Project description and aims...


## Data generation

### Soil sampling

We sampled 29 plots along an altitudinal gradient (3,200-3,700 m.a.s.l.) at the SFF Iguaque (Colombia) in the Cordillera Oriental of the Northern Andes. Each sample consisted of 16 pooled soil-scoops taken 1m apart from each other in a 4m x 4m grid. Shovels were sterilized with alcohol and fire between sampling localities and soil samples were transported under refrigerated conditions. We thoroughly mixxed each sample and took three 15g-subsamples (Pansu 2017). 

### eDNA extraction

We extracted environmental DNA using a saturated phosphate buffer and the NucleoSpin Soil Kit (Macherey-Nagel, Germany) in a mobile laboratory in less than 16h from collection (Taberlet et al. 2018). DNA extraction was stopped before elution and stored in dry conditions to complete the extraction in a better equipped laboratory. 


### DNA amplification

Primers (*16S* bacteria, *18S* eukarya, *ITS* fungi, *trnL* plants, *16S* insects) and controls (negatives and mock community)...

PCR conditions, replicates, and spatial arrangement of samples in the plate...


### Libraries preparation and sequencing

Multiplexing and Illumina sequencing...


## Data filtering

The root names of the files called in this script are ```iguaque``` and ```insects```, however they should be replaced by the file name of the corresponding bacteria, eukarya, fungi, plants, and insects datasets. 

### Required softwares and packages

+ OBITools: 
+ ecoPCR: 
+ Sumaclust: 
+ R packages: vegan, lattice, ROBITools, ROBITaxonomy, plotrix, ENmisc, igraph, pander, stringr.


### Taxonomic database generation

Download EMBL reference database (i.e. embl_r134).

```
mkdir DB
./go_embl.bash
```

Create a reference database per barcode that preserve only the taxonomic information of the corresponding group and convert it to fasta format. Primers varies accoding the barcode as follow:

+ Bacteria: forward_primer=GGATTAGATACCCTGGTAGT; reverse_primer=CACGACACGAGCTGACG.
+ Eukarya: forward_primer=TCACAGACCTGTTATTGC; reverse_primer=TTTGTCTGCTTAATTSCG.
+ Fungi: forward_primer=CAAGAGATCCGTTGTTGAAAGTK; reverse_primer=GGAAGTAAAAGTCGTAACAAGG.
+ Plants: forward_primer=GGGCAATCCTGAGCCAA; reverse_primer=CCATTGAGTCTCTGCACCTATC.
+ Insects: forward_primer=TRRGACGAGAAGACCCTATA; reverse_primer=TCTTAATCCAACATCGAGGTC.

```
ecoPCR -d DB/embl_r140 -e 3 -l 5 -L 200 TRRGACGAGAAGACCCTATA TCTTAATCCAACATCGAGGTC > insects_db
obiconvert --fasta-output insects_db > insects_ref.fasta
```

Keep only those sequences that contain information at the species level.

```
obigrep -d DB/embl_r140 --require-rank=species --require-rank=genus --require-rank=family insects_ref.fasta > insects_ref_clean.fasta
```

Group and dereplicate sequences, keeping only sequence records that are unique.

```
obiuniq -d DB/embl_r140 insects_ref_clean.fasta > insects_ref_clean_uniq.fasta
awk '/^>/{gsub(";$", "", $1);print;}1' insects_ref_clean_uniq.fasta | obiannotate --uniq-id > insects_database_r134.fasta
```


### Sequence quality control and demultiplexing

Align reads from paired-end sequencing. ```--index-file``` points to the file containing the illumina index reads.

```
illuminapairedend --sanger --score-min=40 --index-file=index.fastq -r iguaque_reverse.fastq iguaque_forward.fastq > iguaque_align.fasta
```

Assign sequences to their original samples based on DNA tags and primers. ```-t``` specifies the file containing the samples description. ```-u``` specifies a filename to store the sequences unassigned to any sample. ```-e``` specifies the number of errors allowed for matching primers.

```
ngsfilter -t ngsfilter.txt -e 2 -u iguaque_unassigned.fastq iguaque_align.fastq > iguaque_align_filterE2.fastq
```

Dereplicate unique sequences, grouping together sequence records. ```-m sample```dereplicates per sample.

```
obiuniq -m sample iguaque_align_filterE2.fasta > iguaque_align_filterE2_uniq.fasta
```

Filter out sequences that were not aligned, have ambiguous bases, or are too short. The number between curly brackets ```{}```specifies the minimum sequence length, which varies according to the barcode as follow:

+ Bacteria 16S: 200 bp.
+ Eukarya 18S: 80 bp.
+ Fungi ITS: 100 bp.
+ Plants trnL: 7 bp.
+ Insects 16S: 75 bp.

```
obigrep -s '^[acgt]{75,}$' -a mode:alignment --fasta-output iguaque_align_filterE2_uniq > iguaque_align_filterE2_uniq_nl.fasta
```

Rename sequences (no mandatory).

```
obiannotate --set-identifier='"seq_%03d" % counter' iguaque_align_filterE2_uniq_nl.fasta > iguaque_align_filterE2_uniq_nl_setid.fasta
```

Filter out singletons and keep only interesting attributes. ```'count>10'``` or higher is advisable if taxonomic assignment takes an exaggerate amoung of time given the available computational capabilities.

```
obigrep -p 'count>10' iguaque_align_filterE2_uniq_nl_setid.fasta | obiannotate -k merged_sample -k count > iguaque_align_filterE2_uniq_nl_setid_c10.fasta
```

### Taxonomic assignment 

Make taxonomic assignment.

```
ecotag -d embl_r140 -R iguaque_database_r140.fasta iguaque_align_filterE2_uniq_nl_setid_c10.fasta > iguaque_align_filterE2_uniq_nl_setid_c10_assign_r140.fasta
```

Keep only the sequences assigned to the taxa of interest. ```-r``` specifies the taxid, which varies according to the barcode as follow:

+ Bacteria: 2.
+ Eukarya: 2759.
+ Fungi: 4751.
+ Plants (Spermatophyta): 58024.
+ Insects: 50557.

```
obigrep -d embl_r140 -r 50557 iguaque_align_filterE2_uniq_nl_setid_c1_assign_r140.fasta > GWM-iguaque_align_filterE2_uniq_nl_setid_c10_assign_r140_insects.fasta
```

### Sequence clustering

Cluster sequences by distance or similarity. ```-d``` indicates that the score threshold specified by ```-t``` is expressed as distance (bp), otherwise it is assumed it is expressed as similarity (%). The cluster center is the most abundant sequence.

+ 3 bp distance:  

```
sumaclust -d -r -t 3 iguaque_align_filterE2_uniq_nl_setid_c10_assign_insects.fasta > iguaque_align_filterE2_uniq_nl_setid_c10_assign_r140_insects_t3.fasta
```

+ 97% of similarity: 

```
sumaclust -t 0.97 iguaque_align_filterE2_uniq_nl_setid_c10_assign_r140_insects.fasta > iguaque_align_filterE2_uniq_nl_setid_c10_assign_insects_t97.fasta
```

+ 99% of similarity: 

```
sumaclust -t 0.99 iguaque_align_filterE2_uniq_nl_setid_c10_assign_r140_insects.fasta > iguaque_align_filterE2_uniq_nl_setid_c10_assign_insects_t99.fasta
```

Subsequent analyses were independently done on the three datasets above generated. For the sake of simplicity, only ```-t 3``` is called in the following lines.

Sort the sequences according to their count and export all the information in a tab-delimited file

```
obisort -r -k count iguaque_align_filterE2_uniq_nl_setid_c10_assign_r140_insects_t3.fasta | obitab -o > iguaque_align_filterE2_uniq_nl_setid_c10_assign_r140_insects_t3.tab
```

### MOTUs aggregation

The below post-OBITools filtering and visualisation steps are done in ```R``` based on an OBITools output object.

```
# OBITools output
OBI_OBJ <- "/media/henry/UNTITLED/Henry_06June2019/Iguaque/2020/GWM-841_align_filterE2_uniq_nl_setid_c10_assign_r140_Eukarya_t3.tab"

# ngsfilter file
ngsfilt=read.table('/media/henry/UNTITLED/Henry_06June2019/Iguaque/ngsfiltering/ngsfilterIguaqueEuca.txt')

# Working directory
PATH="[PATH]/Iguaque/"
setwd(PATH)

```

They require loading the following packages:

```
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

Aggregate sequences from the same cluster (MOTU). Keep only the informative columns: id (col. 1), cluster count (col. 8), taxonomic information (col. 3, 4, 10:15, -2:-12), samples (col. 16:-13), and sequence (col. -1). Negative numbers indicate column position from the last to the first.

```
tab <- read.csv(OBI_OBJ, sep="\t", header=T, check.names=F)
# Extract sequence ID, taxonomic DB matching score and sequence, and cluster size
match <- tab[,c(1,3,4,8)]
# Extract taxonomic information and sequence
taxo <- tab[,c(1,(ncol(tab)-11):ncol(tab))]
# Extract abundances per sample and pool them (by summing up) by cluster ID
samples <- aggregate(x=tab[,16:(ncol(tab)-12)], by=list(tab[,5]), FUN=sum)
colnames(samples)[1] <- colnames(tab)[1]
# Add taxonomic DB matching information to pooled abundances
tab <- merge(match, samples, by="id")
# Add taxonomic information and sequence of the "cluster center" ID to its corresponding pooled cluster abundancies
tab <- merge(tab, taxo, by="id")
# Add a missing "a" at the end of the name of the first replicate of localities and negative controls
a1 <- which(colnames(tab) %in% grep("sample:[[:digit:]].*[[:digit:]]$", colnames(tab), value=T)) # Samples missing an "a" in their names
a2 <- which(colnames(tab) %in% grep("sample:TNEGPART.*[[:digit:]]$", colnames(tab), value=T)) # Controls missing an "a" in their names
for (i in c(a1, a2)) {
  colnames(tab)[i] <- paste(colnames(tab)[i], "a", sep="")
}
# Rename count column
colnames(tab)[4] <- "count"
# Export the matrix as a tab-delimited table
write.table(tab, paste(stri_sub(OBI_OBJ, 1, -5), "_ag.tab", sep=""), quote=F, sep="\t", row.names=F)
rm(tab, match, taxo, samples, a1, a2)
```

### Taxonomic formatting

Standarise taxonomic information across all levels and export it as a new tab-delimited table. Note that ":" is replaced by "." in column names.

```
# ecopcr database, used only for taxid manipulation purposes
DB_N="/media/henry/UNTITLED/Henry_06June2019/Iguaque/Pre2020/Reference_DB/embl_r134" 
replication=T

# Import data
OBJ=paste(stri_sub(OBI_OBJ, 1, -5), "_ag.tab", sep="")
DB=read.taxonomy(DB_N)
OBI=import.metabarcoding.data(OBJ)
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

# Export matrix
tmp=t(OBI@reads)
colnames(tmp)=paste("sample:", colnames(tmp), sep="")
write.table(data.frame(OBI@motus,tmp), paste(stri_sub(OBJ, 1, -5), "_taxo.tab", sep=""), row.names=FALSE, col.names=T, quote=F, sep="\t")
```


## Data quality inspection

The following sections will output a serie of plots to visually inspect data quality based on presence and abundance of contaminant sequences and taxonomic assignment scores.

Identify and set different colors to samples and controls.

```
EMPTY=grep("BLK", rownames(OBI@reads)) #black
PCR_CTL=grep("NEG", rownames(OBI@reads)) #grey
POS_CTL=grep("M", rownames(OBI@reads)) #red
CONTROLS=c(PCR_CTL, EMPTY, POS_CTL)

COL="lightgrey" 
BOR="grey"
COL.s="cyan3"
BOR.s="cyan4"
COL.neg='red'
BOR.neg="darkred"
COL.all=rep("cyan3",nrow(OBI)); COL.all[CONTROLS]="red"
BOR.all=rep("cyan4",nrow(OBI)); BOR.all[CONTROLS]="darkred"
```

Index samples and controls, get theis position in the library preparation plates, and assign a color code. A ngsfilter file with coordinates of all the samples in the library preparation plates is needed per taxonomic group.

```
#For each samples: add locations into plates 
OBI@samples$xPlate<-rep("NA", length(OBI@samples$sample))
OBI@samples$yPlate<-rep("NA", length(OBI@samples$sample))

A<-1; B<-2; C<-3; D<-4; E<-5; F<-6; G<-7; H<-8
for (i in 1:length(OBI@samples$sample)) {
  PLATE <- as.numeric(substr(ngsfilt[which(ngsfilt[,2] == as.character(OBI@samples$sample[i])), 6], 14, 15))
  x_col <- as.numeric(substr(ngsfilt[which(ngsfilt[,2] == as.character(OBI@samples$sample[i])), 6], 17, 18))
  y_row <- get(substr(ngsfilt[which(ngsfilt[,2] == as.character(OBI@samples$sample[i])), 6], 19, 19)) + ((PLATE-1)*8)
  OBI@samples$xPlate[i] <- x_col
  OBI@samples$yPlate[i] <- y_row
}
rm(A,B,C,D,E,F,G,H)
```

### Reads count

Plot reads count per sample and across controls to set up a cut-off and indicate with an asterisk the position of the samples below such cut-off at the library plate.

```
ngsfilt<-ngsfilt[which(ngsfilt$V2 %in% rownames(OBI@reads)),]

par(mar=c(5,5,4,4)); layout(matrix(c(1,1,2,1,1,2,1,1,3,5,5,3,5,5,4,5,5,4), 6, 3, byrow = TRUE))
barplot(rowSums(OBI@reads), col=COL.all, border=BOR.all, xlab="Samples", ylab="Number of reads", main="Reads count")

hist(log10(rowSums(OBI@reads[-CONTROLS,])), breaks=40,  col=rgb(0,0,1,0.5), main="PCR CTL", xlab="log10 number of reads", ylab="Nb samples", xlim=c(0,6), lty="blank")
hist(log10(rowSums(OBI@reads[PCR_CTL,])), breaks=40,  col=rgb(1,0,0,0.5), add=T, lty="blank")
hist(log10(rowSums(OBI@reads[-CONTROLS,])), breaks=40,  col=rgb(0,0,1,0.5), main="EMPTY", xlab="log10 number of reads", ylab="Nb samples", xlim=c(0,6), lty="blank")
hist(log10(rowSums(OBI@reads[EMPTY,])), breaks=40,  col=rgb(1,0,0,0.5), add=T, lty="blank")
hist(log10(rowSums(OBI@reads[-CONTROLS,])), breaks=40,  col=rgb(0,0,1,0.5), main="POS CTL", xlab="log10 number of reads", ylab="Nb samples", xlim=c(0,6), lty="blank")
hist(log10(rowSums(OBI@reads[POS_CTL,])), breaks=40,  col=rgb(1,0,0,0.5), add=T, lty="blank")

# Set up cut-off value after visual inspection of the previous histograms
thresh.seqdepth=2 # Change accordingly after visual examination
abline(v=thresh.seqdepth, col="black", lty=2, lwd=2); mtext(side=4, paste("Cut-off < ", thresh.seqdepth, sep=""), cex=0.8, font=3)

plot(OBI@samples$yPlate, OBI@samples$xPlate, pch=21, col=BOR.all, bg=COL.all, xlab="", ylab="", main="Sample position at library plates")
points(as.numeric(OBI@samples$yPlate[which(log10(rowSums(OBI@reads))<thresh.seqdepth)]),
       as.numeric(OBI@samples$xPlate[which(log10(rowSums(OBI@reads))<thresh.seqdepth)])+0.2, 
       pch=8, cex=0.5, lwd=1.5, col='yellow4')
abline(v=seq(8.5,24.5,8), lty=2, col="grey")
```

![Alt text](GWM-841_align_filterE2_uniq_nl_setid_c10_assign_r140_Eukarya_t3_ag_taxo_readscount.jpeg?raw=true)

### MOTUs count

Plot MOTUs count per sample and across controls to set up a cut-off and indicate with an asterisk the position of the samples below such cut-off at the library plate.


```
layout(matrix(c(1,1,1,2,5,5,3,5,5,4,5,5), 4, 3, byrow = TRUE))
barplot(specnumber(OBI@reads, MARGIN=1), col=COL.all, border=BOR.all, xlab="Samples", ylab="Nb OTUs", cex.names=0.5, main="MOTUs count")

hist(log10(specnumber(OBI@reads)[-CONTROLS]), breaks=40, col=rgb(0,0,1,0.5), main="PCR CTL", xlab="log10 number of OTUs", ylab="Nb samples", xlim=c(0,max(log10(specnumber(OBI@reads)))+0.1), lty="blank")
hist(log10(specnumber(OBI@reads)[PCR_CTL]), breaks=40, col=rgb(1,0,0,0.5), add=T, lty="blank")
hist(log10(specnumber(OBI@reads)[-CONTROLS]), breaks=40, col=rgb(0,0,1,0.5), main="EMPTY", xlab="log10 number of OTUs", ylab="Nb samples", xlim=c(0,max(log10(specnumber(OBI@reads)))+0.1), lty="blank")
hist(log10(specnumber(OBI@reads)[EMPTY]), breaks=40, col=rgb(1,0,0,0.5), add=T, lty="blank")
hist(log10(specnumber(OBI@reads)[-CONTROLS]), breaks=40, col=rgb(0,0,1,0.5), main="POS CTL", xlab="log10 number of OTUs", ylab="Nb samples", xlim=c(0,max(log10(specnumber(OBI@reads)))+0.1), lty="blank")
hist(log10(specnumber(OBI@reads)[POS_CTL]), breaks=40, col=rgb(1,0,0,0.5), add=T, lty="blank")

# Set up cut-off value after visual inspection of the previous histograms
thresh.rich=1.2 # Change accordingly after visual examination
abline(v=thresh.rich, col="black", lty=2, lwd=2); mtext(side=4, paste("Cut-off < ", thresh.rich, sep=""), cex=0.8, font=3)

plot(OBI@samples$yPlate, OBI@samples$xPlate, pch=21, col=BOR.all, bg=COL.all, xlab="", ylab="", main="Samples position at library plates")
points(as.numeric(OBI@samples$yPlate[which(log10(specnumber(OBI@reads))<thresh.rich)]),
       as.numeric(OBI@samples$xPlate[which(log10(specnumber(OBI@reads))<thresh.rich)])+0.2,
       pch=8, cex=0.5, lwd=1.5, col='yellow4')
abline(v=seq(8.5,24.5,8), lty=2, col="grey")
```

![Alt text](GWM-841_align_filterE2_uniq_nl_setid_c10_assign_r140_Eukarya_t3_ag_taxo_motuscount.jpeg?raw=true)

### Taxonomic resolution

Plot taxonomic resolution of reads and MOTUs assignments both before and after filtering by a given best identity score threshold and export a table summarising the taxonomic assignment results.

```
# Input data
x = OBI; y = "rank_ok"; z = "bid_ok"
thresh = 0.95
  
# All possible taxonomic levels
taxorank=c("superkingdom", "kingdom", "subkingdom", "superphylum", "phylum", "subphylum", "superclass", "class", "subclass",
             "superorder", "order", "suborder", "infraorder", "superfamily", "family", "subfamily", "supertribe", "tribe", 
             "subtribe", "supergenus", "genus", "subgenus", "superspecies", "species", "subspecies", "varietas", "no rank")
  
# Number of MOTUs
tmp=table(x@motus[y])
taxores.otu=tmp[match(taxorank, names(tmp))]
names(taxores.otu)=taxorank
taxores.otu[which(is.na(taxores.otu)==T)]=0
  
# Number of reads
tmp=aggregate(x@motus$count, by=list(x@motus[,y]), sum)
taxores.reads=tmp[match(taxorank,tmp[,1]),2]
names(taxores.reads)=taxorank
taxores.reads[which(is.na(taxores.reads))]=0
  
# Set below threshold to not assigned
tmp=x@motus[,y]
tmp[which(OBI@motus[,z]<thresh)]="not assigned"
  
# Number of reads above threshold
tmp2=aggregate(OBI@motus$count, by=list(tmp), sum)
taxores.reads.t=tmp2[match(c(taxorank, "not assigned"),tmp2[,1]),2]
names(taxores.reads.t)=c(taxorank, "not assigned")
taxores.reads.t[which(is.na(taxores.reads.t))]=0
  
# Number of MOTUs above threshold
tmp2=table(tmp)
taxores.otu.t=tmp2[match(c(taxorank, "not assigned"), names(tmp2))]
names(taxores.otu.t)=c(taxorank, "not assigned")
taxores.otu.t[which(is.na(taxores.otu.t))]=0

# Plots taxonomic resolution
layout(matrix(c(2,4,1,3,5,1,6,7,7),3,3),heights=c(1,1,0.4))
col.tmp=c(rainbow(length(taxorank)-1,start=0, end=0.5, alpha=0.6), "lightgrey", "darkgrey")
par(mar=c(2,2,1,1), oma=c(0,0,2,0))
frame()
legend("bottom", names(taxores.otu.t), ncol=5, cex=0.7, fill=col.tmp, bty="n")
pie(taxores.reads, col=col.tmp, border="lightgrey", labels="", clockwise=T)
mtext("All data", side=3, cex=0.8)
mtext("Reads", side=2, cex=0.8)
pie(taxores.reads.t, col=col.tmp, border="lightgrey", labels="", clockwise=T)
mtext(paste("Best identities >", thresh) , side=3, cex=0.8)
pie(taxores.otu, col=col.tmp, border="lightgrey", labels="", clockwise=T)
mtext("MOTUs", side=2, cex=0.8)
pie(taxores.otu.t, col=col.tmp, border="lightgrey", labels="", clockwise=T)
  
## Plot reads and MOTUs with low identification scores
par(mar=c(1,4,1,2))
weighted.hist(OBI@motus[, "bid_ok"],OBI@motus[, "count"],breaks=20,col=COL, ylab="Nb reads")
par(mar=c(8,4,3,2))
hist(OBI@motus[, "bid_ok"], breaks=40, col=COL, xlab="Ecotag scores", main="", ylab="Nb. MOTUs")
abline(v=thresh, lty=2, lwd=2)
```

![Alt text](GWM-841_align_filterE2_uniq_nl_setid_c10_assign_r140_Eukarya_t3_ag_taxo_idscores.jpeg?raw=true)

### Contaminant MOTUs

Identify contaminant MOTUs as those that have a maximum frequency in controls. 

```
# Function ContaSlayer
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

# Identify MOTUs with a max. abundance in controls
CONTA=ContaSlayer(OBI,rownames(OBI)[CONTROLS], NULL)
```

Plot frequency of contaminant MOTUs in samples.

```
# Frequency of contaminant in samples
OBI.freq=decostand(OBI@reads,"total",1)

# Set up a four-panels plot
m4 <- round(nrow(OBI@samples)/4, digit=0)
if (m4*4!=nrow(OBI@samples)) {m4 <- m4+1}
OBI.freq.parse=c(seq(0,round(nrow(OBI@samples), digit=0),m4), nrow(OBI@samples))

#Data visualization
par(mfrow=c(2,2), oma=c(0,0,2,0), mar=c(3,5,1,1))
for(i in 1:(length(OBI.freq.parse-1)-1)){
  image(-log10(OBI.freq[(OBI.freq.parse[i]+1):(OBI.freq.parse[i+1]),CONTA]), xaxt="n", yaxt="n", col=rainbow(12))
  dm=dim(OBI.freq[(OBI.freq.parse[i]+1):(OBI.freq.parse[i+1]),CONTA])
  abline(v=seq(0,1,l=dm)
         [grep("BLK|NEG|M",rownames(OBI.freq)[(OBI.freq.parse[i]+1):(OBI.freq.parse[i]+m4)])],
         col=COL.neg, lty=3, lwd=0.5)
  axis(side=1,at=seq(0,1,l=dm),
       labels=rownames(OBI@samples)[(OBI.freq.parse[i]+1):(OBI.freq.parse[i]+dm)],
       las=2, cex.axis=0.3)
  axis(side=2, at=seq(0,1,l=length(CONTA)), labels=paste(CONTA, OBI@motus[CONTA,"sci_name_ok"]),
       cex.axis=0.3, las=2)
}
mtext(side=3, paste("Frequency of contaminant MOTUs in samples based on its abundance in controls. N =", length(CONTA)), outer=T, cex=1.5)
```

![Alt text](GWM-841_align_filterE2_uniq_nl_setid_c10_assign_r140_Eukarya_t3_ag_contseq.jpeg?raw=true)

### Odd MOTUs

Identify odd MOTUs with low identification scores and plot its frequency in samples.

```
# Identify odd MOTUs
tmp=rownames(OBI@motus)[which(OBI@motus[, "bid_ok"]<thresh)]
ODDOTU=if(length(grep("cluster", rownames(OBI@motus)))!=0) {
  c(tmp[which(is.na(OBI2@motus[tmp, "cluster95"])==T)],
    rownames(OBI2@motus)[which(is.na(match(OBI2@motus[,"cluster95"],unique(OBI2@motus[tmp, "cluster95"][-which(is.na(OBI2@motus[tmp,"cluster95"])==T)])))==F)])
} else {
  tmp
}

# Plot frequency of odd MOTUs in samples
par(mfrow=c(2,2), oma=c(0,0,2,0), mar=c(3,5,1,1))
for(i in 1:(length(OBI.freq.parse)-1)){
  image(log10(OBI.freq[(OBI.freq.parse[i]+1):(OBI.freq.parse[i+1]),ODDOTU]), xaxt="n", yaxt="n", col=rainbow(12))
  dm=dim(OBI.freq[(OBI.freq.parse[i]+1):(OBI.freq.parse[i+1]),ODDOTU])
  abline(v=seq(0,1,l=dm)
         [grep("BLK|NEG|M",rownames(OBI.freq)[(OBI.freq.parse[i]+1):(OBI.freq.parse[i]+round(nrow(OBI.freq)/4, digit=0))])],
         col=COL.neg, lty=3)
  axis(side=1,at=seq(0,1,l=dm, digits=0),
       labels=rownames(OBI.freq)[(OBI.freq.parse[i]+1):(OBI.freq.parse[i]+dm)],
       las=2, cex.axis=0.3)
  axis(side=2, at=seq(0,1,l=length(ODDOTU)), labels= OBI@motus[ODDOTU,"sci_name_ok"],cex.axis=0.3, las=2)
}
mtext(side=3, paste("Frequency of odd MOTUs in samples based on low ecotag scores. N =", length(ODDOTU)), outer=T, cex=1)
```

![Alt text](GWM-841_align_filterE2_uniq_nl_setid_c10_assign_r140_Eukarya_t3_ag_oddseq.jpeg?raw=true)

### Samples with contaminant and odd MOTUs

Plot samples with more than 1% of MOTUs either more abundant in the controls than in the samples or with a low ecotag score.

```
# Heavily contaminated samples based on MOTUs more abundant in controls
par(mar=c(5,5,2,1), mfrow=c(2,2), oma=c(0,0,2,0))
hist(log10(rowSums(OBI.freq[-CONTROLS,CONTA])+1e-5), breaks=40, col=COL, xlim=log10(c(0,1)+1e-5), main="", xlab="log10 proportion contaminant reads")

# Most likely to be bad samples
thresh.contasamp=-2
CONTASAMP=rownames(OBI)[-CONTROLS][which(log10(rowSums(OBI.freq[-CONTROLS,CONTA])+1e-5)>=thresh.contasamp)]
abline(v=thresh.contasamp, lty=2, lwd=2); mtext(paste("Cut-off >= ", thresh.contasamp, sep=""),side=3, line=0, cex=0.7, font=3)

plot(OBI@samples$yPlate, OBI@samples$xPlate, pch=21, col=BOR.all, bg=COL.all, cex=1, xlab="y plate", ylab="x plate")
abline(v=seq(8.5,24.5,8), lty=2, col="grey")
points(as.numeric(OBI@samples$yPlate[match(CONTASAMP,rownames(OBI))]), 
       as.numeric(OBI@samples$xPlate[match(CONTASAMP,rownames(OBI))])+0.2,
       pch=8, cex=0.7, lwd=2)


# Sample with many odd OTUs
hist(log10(rowSums(OBI.freq[-CONTROLS,ODDOTU])+1e-5), breaks=20, col=COL, xlim=log10(c(0,1)+1e-5), main="", xlab="log10 proportion odd MOTUs (low ecotag scores)")

# Most likely to be bad samples
thresh.oddotusamp=-2
ODDOTUSAMP=rownames(OBI)[which(log10(rowSums(OBI.freq[,ODDOTU])+1e-5)>=thresh.oddotusamp)]
abline(v=thresh.oddotusamp, lty=2, lwd=2); mtext(paste("Cut-off >= ", thresh.oddotusamp, sep=""),side=3, line=0, cex=0.7, font=3)

plot(OBI@samples$yPlate, OBI@samples$xPlate, pch=21, col=BOR.all, bg=COL.all, cex=1, xlab="y plate", ylab="x plate")
abline(v=seq(8.5,95.5,8), lty=2, col="grey")
points(as.numeric(OBI@samples$yPlate[match(ODDOTUSAMP,rownames(OBI))]), 
       as.numeric(OBI@samples$xPlate[match(ODDOTUSAMP,rownames(OBI))])+0.2,
       pch=8, cex=0.7, lwd=2)

mtext(side=3, "Samples with >1% of contaminant and odd MOTUs", outer=T)
```

![Alt text](GWM-841_align_filterE2_uniq_nl_setid_c10_assign_r140_Eukarya_t3_ag_conoddsamples.jpeg?raw=true)


### Effect of MOTUs filtering

Get rid of previously identified contaminant and odd MOTUs and check how this affects the similarity patterns between and within controls, samples, and replicates.

```
# Duplicate matrix and remove contaminant and odd MOTUs
OBI2<-OBI
OBI2@reads=OBI@reads[,-match(unique(c(CONTA, ODDOTU)), colnames(OBI@reads))]
COL.all2 <- COL.all[-which(rowSums(OBI2@reads)==0)]
OBI2@reads <- OBI2@reads[-which(rowSums(OBI2@reads)==0),]
OBI2@motus=OBI2@motus[which(rownames(OBI2@motus) %in% colnames(OBI2@reads)),]
OBI2@samples=OBI2@samples[which(rownames(OBI2@samples) %in% rownames(OBI2@reads)),]
OBI2@samples$reads.contaodd=rowSums(OBI2@reads)
OBI2@samples$otu100.contaodd=specnumber(OBI2@reads)
OBI2@motus$count=colSums(OBI2@reads)

# Calculate pairwise similarity distances
# Function vegsim
vegsim=function(x) {
  x.dist=vegdist(x, "bray")
  x.sim=1-as.matrix(x.dist)
  out=as.dist(x.sim, diag=F, upper=F)
  out
}
OBI.sim=vegsim(OBI@reads)
OBI.sim2=vegsim(OBI2@reads)

# Convert similarity matrix into a three-columns format list
# Function simlist
simlist=function(x) {
  out=data.frame(t(combn(labels(x),2)), as.numeric(x))
  names(out) <- c("obj1", "obj2", "sim")
  out
}
OBI.simlist=simlist(OBI.sim); OBI.simlist2=simlist(OBI.sim2)

# Annotate pairwise comparisons by type: sample vs. control, replicate vs. replicate, etc.
# Functiont comp.type
comp.type=function(x,y,z1,z2){
  #x a distance list table
  #y a metabarcoding object
  #z1 a vector indicating negative indices in the metabarcoding object
  #z2 a factor of sample names, repeats are indicated with the same sample name
  #controls intravar
  x$type=rep(NA, nrow(x))
  ind.negneg=intersect(which(is.na(match(as.vector(x[,1]), rownames(y)[z1]))==F), which(is.na(match(as.vector(x[,2]), rownames(y)[z1]))==F))
  x[ind.negneg,"type"]="CONTROL-CONTROL"
  #controls vs sample var
  ind.negsamp=c(intersect(which(is.na(match(as.vector(x[,1]), rownames(y)[z1]))==F), which(is.na(match(as.vector(x[,2]), rownames(y)[z1]))==T)),
                intersect(which(is.na(match(as.vector(x[,2]), rownames(y)[z1]))==F), which(is.na(match(as.vector(x[,1]), rownames(y)[z1]))==T)))
  x[ind.negsamp,"type"]="CONTROL-SAMPLE"
  #replicates intravar
  repltech=names(table(z2))[which(table(z2)>1)]
  ind.repin=unlist(lapply(repltech, function(a) intersect(grep(a,x[,1]), grep(a, x[,2]))))
  x[ind.repin, "type"]="REPLICATES"
  #sample vs sample var
  x[-c(ind.negneg, ind.negsamp, ind.repin),"type"]="SAMPLE-SAMPLE"
  x
}

# Identify the replicates
REPLICATES=as.factor(gsub("a", "", rownames(OBI@samples))) # factor, each level is a sample
REPLICATES=gsub("b", "", REPLICATES)
REPLICATES=gsub("c", "", REPLICATES)

# Plot
OBI.compsim=comp.type(OBI.simlist, OBI, CONTROLS, REPLICATES)
OBI.compsim2=comp.type(OBI.simlist2, OBI, CONTROLS, REPLICATES)
par(mfrow=c(1,2), mar=c(10,5,2,1))
boxplot(OBI.compsim$sim~OBI.compsim$type, col=COL, border=BOR, main="Raw data", las=2, ylab="Pairwise similarity")
boxplot(OBI.compsim2$sim~OBI.compsim2$type, col=COL, border=BOR, main="After conta/odd MOTUs filtering", las=2,  ylab="Pairwise similarity")
```

![Alt text](GWM-841_align_filterE2_uniq_nl_setid_c10_assign_r140_Eukarya_t3_ag_motusfiltering.jpeg?raw=true)

### Similarity across replicates

Plot similarity scores among replicates of the same sample and among replicates of different samples.

```
# Detect outliers using a clustering approach
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


#Get two lists: one for similarity among replicates, and one among samples
MyList_IntraSamp<-list()
kk<-1
MyList_InterSamp<-list()
ll<-1
for (i in 1:nrow(OBI.simlist2)) {
  if( substr(OBI.simlist2[i,1],1,str_length(OBI.simlist2[i,1])-1) == substr(OBI.simlist2[i,2],1,str_length(OBI.simlist2[i,2])-1)) {
    MyList_IntraSamp[[i]]<-OBI.simlist2[i,3]
    kk<-kk+1 
  } else {
    MyList_InterSamp[[i]]<-OBI.simlist2[i,3]
    ll<-ll+1
  }
}
MyList_IntraSamp_OK<-unlist(MyList_IntraSamp)
MyList_InterSamp_OK<-unlist(MyList_InterSamp)

# Plot
layout(matrix(c(2,1,1,3,1,1), 3,2))

#Visualize similarity scores intra and inter-sample (log10 of sample similarity)
hist(log10(MyList_InterSamp_OK),  col=rgb(0,0,1,0.5), main="", xlab="log10 sample similarity", ylab="Nb comparison", xlim=c(-5,0))
hist(log10(MyList_IntraSamp_OK),  col=rgb(1,0,0,0.5), add=T)

#Visualize similarity scores intra and inter-sample
hist(MyList_IntraSamp_OK, col=COL, border=BOR, main="", xlab="Intra-sample similarity", breaks=20, xlim=c(0,1), freq=T)
thresh.graph=0.1  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!!! CUTOFF HERE !!!!
abline(v=thresh.graph, lty=2, lwd=2); mtext(side=3, paste("cutoff = ", thresh.graph), cex=0.7, font=3, outer=F)
hist(MyList_InterSamp_OK, col=COL, border=BOR, main="", xlab="Inter-sample similarity", xlim=c(0,1), breaks=20, freq=T)
thresh.graph=0.1  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!!! CUTOFF HERE !!!!
abline(v=thresh.graph, lty=2, lwd=2); mtext(side=3, paste("cutoff = ", thresh.graph), cex=0.7, font=3, outer=F)

```

![Alt text](GWM-841_align_filterE2_uniq_nl_setid_c10_assign_r140_Eukarya_t3_ag_similarity.png?raw=true)

### Similarity clusters

Plot clusters of replicates based on their similarity scores.

```
##Generate cluster
OBI.g2=quickg(OBI.sim2, COL.all2, COL.all2, thresh.graph)
OBI.clust2=infomap.community(OBI.g2, nb.trials=100)

##Visualize clusters: not mandatory
par(mar=c(3,3,1,3))
layout(matrix(c(1,1,2,3), 2,2))
col.g=as.factor(membership(OBI.clust2))
levels(col.g)=rainbow(max(membership(OBI.clust2)))
coords=layout.fruchterman.reingold(OBI.g2)
plot(OBI.g2, layout=coords, vertex.size=2, vertex.label=rownames(OBI2), vertex.color=as.vector(col.g), vertex.label.cex=0.7, vertex.label.dist=-0.2, vertex.frame.color=as.vector(col.g))
legend("top", c("Sample", "Control"), pch=8, col=c("cyan4", "darkred"), cex=1, ncol=3)
legend("bottom",paste("cluster", levels(as.factor(membership(OBI.clust2)))), pch=8, col=levels(col.g), cex=1, ncol=3)
plot(OBI@samples$yPlate, OBI@samples$xPlate, pch=21, col=BOR.all, bg=COL.all, cex=1, xlab="y plate", ylab="x plate")
abline(v=seq(8.5,24.5,8), lty=2, col="grey")
plot(OBI@samples$yPlate, OBI@samples$xPlate, pch=1, col=BOR.all, bg=COL.all, cex=1, xlab="y plate", ylab="x plate")
abline(v=seq(8.5,24.5,8), lty=2, col="grey")
points(OBI@samples$yPlate, OBI@samples$xPlate, pch=8, cex=0.5, col=as.vector(col.g), lwd=2)

```

![Alt text](GWM-841_align_filterE2_uniq_nl_setid_c10_assign_r140_Eukarya_t3_ag_clusters.png?raw=true)



## Samples filtering

THERE ARE TWO APPROACHES...

### Approach A

Filter out samples that are not clustered with their relatives and samples that are clustered with control samples (BLK).

```
# Detect samples that are clustered with TRUE control samples (Extration + PCR)
TRUE_CONTROLS<-rownames(OBI)[CONTROLS][-grep('BLK', rownames(OBI)[CONTROLS])]
tmp1=names(unlist(lapply(unique(membership(OBI.clust2)[TRUE_CONTROLS]),function(x)which(membership(OBI.clust2)[-CONTROLS]==x))))

# Detect replicates that are not clustered with their relatives
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

#Samples that are not clustered with their relatives + samples that are clustered with TRUE control samples 
OUTGRAPH=unique(c(tmp3[-which(tmp3 %in% rownames(OBI)[CONTROLS]==T)], tmp1))

# Samples with few reads, similar to some controls. Threshold value was visually detected above.
SMALL=rownames(OBI)[-CONTROLS][which(log10(rowSums(OBI@reads)[-CONTROLS])<thresh.seqdepth)]

#Data formatting (removal of outliers + controls)
OBI3=OBI
#Removal of all "bad" samples: OUTGRAPH (again choose the level of filtering) + those with low number of reads
OBI3@reads=OBI3@reads[-na.omit(c(match(OUTGRAPH, rownames(OBI2@reads)),CONTROLS, match(SMALL, rownames(OBI2@reads)))),]
OBI3@motus=OBI3@motus[which(rownames(OBI3@motus) %in% colnames(OBI3@reads)),]
OBI3@samples=OBI3@samples[which(rownames(OBI3@samples) %in% rownames(OBI3@reads)),]
OBI3@motus$count=colSums(OBI@reads)

# Remove MOTUs with zero count
if (any(OBI3@motus$count==0)) {
  OBI3@motus = OBI3@motus[-which(OBI3@motus$count==0),]
  OBI3@reads = OBI3@reads[,which(colnames(OBI3@reads) %in% rownames(OBI3@motus))]
}

#Data export
tmp=t(OBI3@reads)
colnames(tmp)=paste("sample:", colnames(tmp), sep="")
write.table(data.frame(OBI3@motus,tmp), paste(stri_sub(OBI_OBJ, 1, -5),"_ag_taxo_cleanA.tab", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
```

### Approach B

Filter out non-replicating samples.

```
OBI2<-OBI
# Be sure that sample names start by 'sample:' (and not 'sample.')
# @reads contains the reads table
# @motus contains information about each motus (taxo, nb of reads etc.)
# @sample lists all samokes

#Before starting you need to create a metabarlist object
metabarlist<-list(reads=OBI2@reads, pcrs=data.frame(row.names=rownames(OBI2@reads), rownames(OBI2@reads)))

#You also need to group your replicates per sample
REPLICATES_CLEAN=as.factor(gsub("a", "", rownames(OBI2@samples))) # factor, each level is a sample
REPLICATES_CLEAN=gsub("b", "", REPLICATES_CLEAN)
REPLICATES_CLEAN=gsub("c", "", REPLICATES_CLEAN)


#' Identifying outlier replicates
#'
#' Identifying  non-replicating samples or controls in the sample x OTU table from a \code{\link{metabarlist}} object.
#' Process many iterations to compare the density of pairwise distances within replicates vs between samples.
#'
#' @param metabarlist a \code{metabarlist} object
#' @param FUN a function returning a distance matrix. The distance matrix should be an object of class `dist` that has the same length as table `reads`.
#' @param groups a vector containing the identifier of replicates. The vector must has the same length of the motu x sample table from a \code{\link{metabarlist}} object. Default = metabarlist$pcrs$sample_id
#' @param graphics a boolean value to plot  distances densities for each iteration. Default = FALSE
#'
#' @name pcr_outlayer
#'
#' @return a dataframe with the replicates groups and a column `replicating` or throws a stop
#'
#' @details
#'
#' This function identifies non-replicating samples or controls.
#'
#' The parameter `groups` defines groups of replicates (i.e. replicates from a same sample). The vector should be sorted like the table `PCRs` from the \code{\link{metabarlist}}.
#' Note: if the distance between replicates is too high compared to the distance between samples, the function cannot return result because all replicates are removed.
#' The parameter `FUN` defines the function used to compute the distance matrix. The function must return a object of class `dist` with the same length that its input table.
#' The default function uses the function `decostand` and `vegdist` of package `vegan` to perform a correspondance analysis of \code{\link{metabarlist}} table `reads`, and return a matrix distance.
#' Default function detail:
#' bray_function <- function(reads) {
#'   distance_matrix <- vegdist(decostand(reads, method = 'total'), method='bray')
#'   return(distance_matrix)
#' }
#'
#' When the parameter `graphics` is True, a graphic is plotted with the density of distance between replicates and between samples. The threshold is plotted as a vertical line at the intersection of two density curves.
#'
#' Note: when many projects are pulled in the same plate, you must process this function for each project. If you execute this function on many projects the variability between projects can disturb the computation of densities and then removed all samples for one project.

# recursive function to find non replicating PCRs or controls
filter_replicat <- function(sub_matrix, threshold) {
  replicat_to_remove <- c()
  if (any(sub_matrix > threshold)) {
    if (nrow(sub_matrix) == 2) {
      replicat_to_remove <- c(replicat_to_remove, rownames(sub_matrix))
    } else {
      replicat_to_remove <- c(
        colnames(sub_matrix)[which.max(colSums(sub_matrix))],
        filter_replicat(
          sub_matrix[
            -which.max(colSums(sub_matrix)),
            -which.max(colSums(sub_matrix))
          ],
          threshold
        )
      )
    }
  }
  return(replicat_to_remove)
}

# distance function with ade4 package and coa analysis
coa_function <- function(reads) {
  correspondence_analysis <- dudi.coa(sqrt(reads), scannf = FALSE, nf = 2)
  distance_matrix <- dist(correspondence_analysis$li)
  return(distance_matrix)
}

# distance function with vegan package and Bray-Curtis distance
bray_function <- function(reads) {
  distance_matrix <- vegdist(decostand(reads, method = "total"), method = "bray")
  return(distance_matrix)
}

# main function
pcr_outlier <- function(metabarlist,
                               FUN = bray_function,
                               groups = REPLICATES_CLEAN,
                               graphics = FALSE) {
    if (length(groups) != nrow(metabarlist$pcrs)) {
      stop("provided groups should match the number of pcrs")
    }

    subset_data <- data.frame(
      groups = groups, replicating = TRUE,
      row.names = rownames(metabarlist$pcrs)
    )

    print(subset_data)
    subset_data[which(rowSums(metabarlist$reads) == 0), "replicating"] <- FALSE

    iteration <- 0
    repeat {
      iteration <- iteration + 1
      print(paste("Iteration", iteration))

      # get only reads for replicating samples or controls
      matrix_with_replicate <- metabarlist$reads[
        rownames(subset_data),
      ][subset_data$replicating, ]

      # compute matrix distance
      function_result <- FUN(matrix_with_replicate)

      if (class(function_result) != "dist") {
        stop("The result of the provided function is incorrect! The function must return object 'dist'!")
      }

      if (length(labels(function_result)) != length(rownames(matrix_with_replicate))) {
        stop("The result of the provided function is incorrect! The dimension of the function output is incorrect!")
      }

      if (!all(labels(function_result) %in% rownames(matrix_with_replicate))) {
        stop("The result of the provided function is not correct! The labels of function output do not match with the data!")
      }

      distance_matrix <- as.matrix(function_result)

      replicates <- subset_data[rownames(distance_matrix), "groups"]
      within_replicates <- outer(replicates,
        replicates,
        FUN = "=="
      ) & upper.tri(distance_matrix)
      between_replicates <- outer(replicates,
        replicates,
        FUN = "!="
      ) & upper.tri(distance_matrix)

      if (length(distance_matrix[within_replicates]) < 2) {
        stop("Too many replicates have been removed!")
      }
      within_replicate_density <- density(distance_matrix[within_replicates],
        from = 0, to = max(distance_matrix),
        n = 1000
      )

      if (length(distance_matrix[between_replicates]) < 2) {
        stop("Too many replicates have been removed!")
      }
      between_replicate_density <- density(distance_matrix[between_replicates],
        from = 0, to = max(distance_matrix),
        n = 1000
      )

      threshold_distance <- between_replicate_density$x[
        min(which(within_replicate_density$y < between_replicate_density$y))
      ]

      if (graphics) {
        plot(within_replicate_density$x, within_replicate_density$y,
          type = "l", xlab = "Distances", ylab = "Density",
          main = paste("Distances densities iteration", iteration)
        )
        lines(between_replicate_density, col = "blue")
        abline(v = threshold_distance, col = "red")
      }

      need_to_be_checked <- unique(subset_data[
        rownames(which(
          (distance_matrix > threshold_distance) & within_replicates,
          arr.ind = T
        )),
        "groups"
      ])
      if (length(need_to_be_checked) > 0) {
        for (group in need_to_be_checked) {
          sub_matrix <-
            distance_matrix[
              subset_data[rownames(distance_matrix), "groups"] == group &
                subset_data[rownames(distance_matrix), "replicating"],
              subset_data[rownames(distance_matrix), "groups"] == group &
                subset_data[rownames(distance_matrix), "replicating"]
            ]

          non_replicating <- filter_replicat(
            sub_matrix,
            threshold_distance
          )
          subset_data$replicating[
            rownames(subset_data) %in% non_replicating
          ] <- FALSE
        }
      }
      else {
        break
      }
    }

    #### warning if more than 20% of replicates are removed
    if (dim(subset_data[subset_data$replicating == F, ])[1] / dim(subset_data)[1] > 0.2) {
      warning("More than 20% of replicates are removed !")
    }
    return(subset_data)
}

# Samples replicability
MyResult<-pcr_outlier(metabarlist, FUN=bray_function, groups = REPLICATES_CLEAN, graphics = T)
```

![Alt text](GWM-841_align_filterE2_uniq_nl_setid_c10_assign_r140_Eukarya_t3_ag_taxo_cleanB.png?raw=true)

```
# Remove non-replicating samples
OBI2@reads = OBI2@reads[MyResult[,2],]

# Redefine controls position
EMPTY=grep("BLK", rownames(OBI2@reads))
PCR_CTL=grep("NEG", rownames(OBI2@reads))
POS_CTL=grep("M", rownames(OBI2@reads))
CONTROLS=c(PCR_CTL, EMPTY, POS_CTL)

# Remove control samples
OBI2@reads = OBI2@reads[-CONTROLS,]
OBI2@motus = OBI2@motus[which(rownames(OBI2@motus) %in% colnames(OBI2@reads)),]
OBI2@samples = as.data.frame(OBI2@samples[which(rownames(OBI2@samples) %in% rownames(OBI2@reads)),])

# Remove MOTUs with zero count
OBI2@motus$count=colSums(OBI2@reads)
if (any(OBI2@motus$count==0)) {
  OBI2@motus = OBI2@motus[-which(OBI2@motus$count==0),]
  OBI2@reads = OBI2@reads[,which(colnames(OBI2@reads) %in% rownames(OBI2@motus))]
}

#Data export
tmp=t(OBI2@reads)
colnames(tmp)=paste("sample:", colnames(tmp), sep="")
write.table(data.frame(OBI2@motus,tmp), paste(str_sub(OBI_OBJ, 1, -5),"_ag_taxo_cleanB.tab", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
```

## Data analysis

### Required softwares and packages



### Alpha-diversity





### Beta-diversity






### Environmental associations





## References

+ Pansu J. 2017. Muestras ambientales. In: Gonzalez MA, Arenas-Castro H (Eds). Recolección de tejidos biológicas para análisis genéticos. Instituto Humboldt.
+ Taberlet P, Bonin A, Zinger L, Coissac E. 2018. Environmental DNA for biodiversity research and monitoring. Oxford University Press.
