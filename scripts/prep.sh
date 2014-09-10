#!bin/bash

# This script updates the chromosome differences database from wormbase.

wget ftp://ftp.sanger.ac.uk/pub2/wormbase/software/Remap-between-versions/remap.tar.bz2
tar jxf remap.tar.bz2

mv Remap-for-other-groups/CHROMOSOME_DIFFERENCES/ chrom_diff/

python load_chrom_diff.py

rm -f -r chrom_diff/
rm -f -r Remap-for-other-groups
rm remap*

mv chrom_diff.db ../chrom_diff.db

echo "Chromosome differences updated"