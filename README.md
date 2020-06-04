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
illuminapairedend --sanger --score-min=40 --index-file=index.fastq -r iguaque_reverse.fastq iguaque_forward.fastq > iguaque_align.fastq
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


## Community matrix curation

Post-OBITools filtering is done in ```R```. It requires loading the following packages:

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

### MOTU's abundance pooling

Aggregate sequences from the same cluster (MOTU). Keep only the informative columns: id (col. 1), cluster count (col. 8), taxonomic information (col. 3, 4, 10:15, -2:-12), samples (col. 16:-13), and sequence (col. -1). Negative numbers indicate column position from the last to the first.

```
tab <- read.csv("iguaque_align_filterE2_uniq_nl_setid_c10_assign_r140_insects_t3.tab", sep="\t", header=T, check.names=F)
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
# Export the matrix as a tab-delimited table
write.table(tab, "iguaque_align_filterE2_uniq_nl_setid_c10_assign_r140_insects_t3_ag.tab", quote=F, sep="\t", row.names=F)
rm(tab, match, taxo, samples, a1, a2)
```

### Envirnoment setting

Set paths and file names.

```
# Working directory
PATH="[PATH]/Iguaque/"
# Tab-delimited community matrix
OBJ="/media/henry/UNTITLED/Henry_06June2019/Iguaque/2020/GWM-841_align_filterE2_uniq_nl_setid_c10_assign_r140_Eukarya_t3_ag.tab"
# ecopcr database, used only for taxid manipulation purposes
DB_N="/media/henry/UNTITLED/Henry_06June2019/Iguaque/Pre2020/Reference_DB/embl_r134" 
replication=T
```

Import data.

```
setwd(PATH) 
DB=read.taxonomy(DB_N)
OBI=import.metabarcoding.data(OBJ)
```

### Taxonomic formatting

Standarise taxonomic information across all levels and export it as a new tab-delimited table. Note that ":" is replaced by "." in column names.

```
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

### Data quality inspection

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
# Modify path accordingly
ngsfilt=read.table('/media/henry/UNTITLED/Henry_06June2019/Iguaque/ngsfiltering/ngsfilterIguaqueEuca.txt')

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

# Set cut-off value after visual inspection of the previous histograms
thresh.seqdepth=2.4 
abline(v=thresh.seqdepth, col="black", lty=2, lwd=2); mtext(side=4, paste("Cut-off < ", thresh.seqdepth, sep=""), cex=0.8, font=3)
plot(OBI@samples$yPlate, OBI@samples$xPlate, pch=21, col=BOR.all, bg=COL.all, xlab="", ylab="", main="Sample position at library plates")
points(as.numeric(OBI@samples$yPlate[which(log10(rowSums(OBI@reads))<thresh.seqdepth)]),
       as.numeric(OBI@samples$xPlate[which(log10(rowSums(OBI@reads))<thresh.seqdepth)])+0.2, 
       pch=8, cex=0.5, lwd=1.5, col='yellow4')
abline(v=seq(8.5,24.5,8), lty=2, col="grey")
```

![Alt text](GWM-841_align_filterE2_uniq_nl_setid_c10_assign_r140_Eukarya_t3_ag_taxo_readscount.jpeg?raw=true)




## Data analysis

### Required softwares and packages



### Alpha-diversity





### Beta-diversity






### Environmental associations





## References

+ Pansu J. 2017. Muestras ambientales. In: Gonzalez MA, Arenas-Castro H (Eds). Recolección de tejidos biológicas para análisis genéticos. Instituto Humboldt.
+ Taberlet P, Bonin A, Zinger L, Coissac E. 2018. Environmental DNA for biodiversity research and monitoring. Oxford University Press.
