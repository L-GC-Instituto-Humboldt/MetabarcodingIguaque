#!/bin/bash 


NBJOBS=3
DSTDIR="./DB"


wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

wget -nH --cut-dirs=4 -m -np ftp://ftp.ebi.ac.uk/pub/databases/embl/release/std/ --accept 'rel_std_*.dat.gz'


make -j $NBJOBS -f Makefile DSTDIR=$DSTDIR

cd $DSTDIR 

find . -name \* -type l -exec rm {} \;

ls  *.sdx | LC_ALL=en_US sort | gawk '{FILENAME=$1;split(FILENAME, a, "_");section=a[3]; part=a[4];release=a[5];print "mv", FILENAME, "embl_"release"_" a[3] "_" sprintf("%03d",a[4])".sdx"; print "ln -s","embl_"release"_" a[3] "_" sprintf("%03d",a[4])".sdx", "embl_"release"_"sprintf("%03d", ++nb)".sdx"; }' - | sh


release=`gawk '{split(FILENAME,a,"_");gsub(/\.adx$/,"",a[2]);print a[2];nextfile;}' embl_*.adx`


for f in `ls embl_${release}_???_???.sdx | gawk '{id=$1;gsub(/_[0-9]{3}\.sdx$/,"",id);print id;}' | uniq`
do

    ln -s embl_$release.rdx $f.rdx
    ln -s embl_$release.tdx $f.tdx
    ln -s embl_$release.ndx $f.ndx
    ln -s embl_$release.adx $f.adx

done
