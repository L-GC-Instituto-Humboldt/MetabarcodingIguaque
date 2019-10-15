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
+ Sumatra: 
+ Infomap: 
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
ecoPCR -d DB/embl_r134 -e 3 -l 5 -L 200 TRRGACGAGAAGACCCTATA TCTTAATCCAACATCGAGGTC > insects_db
obiconvert --fasta-output insects_db > insects_ref.fasta
```

Keep only those sequences that contain information at the species level.

```
obigrep -d DB/embl_r134 --require-rank=species --require-rank=genus --require-rank=family insects_ref.fasta > insects_ref_clean.fasta
```

Group and dereplicate sequences, keeping only sequence records that are unique.

```
obiuniq -d DB/embl_r134 insects_ref_clean.fasta > insects_ref_clean_uniq.fasta
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

Filter out sequences that were not aligned, have ambiguous bases, or are too short. ```-l```specifies the minimum sequence length, which varies according to the barcode as follow:

+ Bacteria 16S: 200 bp.
+ Eukarya 18S: 80 bp.
+ Fungi ITS: 100 bp.
+ Plants trnL: 7 bp.
+ Insects 16S: 75 bp.

```
obigrep -s '^[acgt]$' -l 75 -a mode:alignment --fasta-output iguaque_align_filterE2_uniq > iguaque_align_filterE2_uniq_nl.fasta
```

Rename sequences (no mandatory).

```
obiannotate --set-identifier='"seq_%03d" % counter' iguaque_align_filterE2_uniq_nl.fasta > iguaque_align_filterE2_uniq_nl_setid.fasta
end
```

Filter out singletons and keep only interesting attributes. ```'count>10'``` or higher is advisable if taxonomic assignment takes an exaggerate amoung of time given the available computational capabilities.

```
obigrep -p 'count>1' iguaque_align_filterE2_uniq_nl_setid.fasta | obiannotate -k merged_sample -k count > iguaque_align_filterE2_uniq_nl_setid_c1.fasta
```

### Taxonomic assignment 

Make taxonomic assignment.

```
ecotag -d embl_r134 -R iguaque_database_r134.fasta iguaque_align_filterE2_uniq_nl_setid_c1.fasta > iguaque_align_filterE2_uniq_nl_setid_c1_assign.fasta
```

Keep only the sequences assigned to the taxa of interest. ```-r``` specifies the taxid, which varies according to the barcode as follow:

+ Bacteria: 2.
+ Eukarya: 2759.
+ Fungi: 4751.
+ Plants (Spermatophyta): 58024.
+ Insects: 50557.

```
obigrep -d embl_r134 -r 50557 iguaque_align_filterE2_uniq_nl_setid_c1_assign.fasta > GWM-iguaque_align_filterE2_uniq_nl_setid_c1_assign_insects.fasta
```

Cluster sequences by: 
+ 3 bp...
+ 97%...
+ 99%...






### Community matrix curation





## Data analysis

### Required softwares and packages



### Alpha-diversity





### Beta-diversity






### Environmental associations





## References

+ Pansu J. 2017. Muestras ambientales. In: Gonzalez MA, Arenas-Castro H (Eds). Recolección de tejidos biológicas para análisis genéticos. Instituto Humboldt.
+ Taberlet P, Bonin A, Zinger L, Coissac E. 2018. Environmental DNA for biodiversity research and monitoring. Oxford University Press.
