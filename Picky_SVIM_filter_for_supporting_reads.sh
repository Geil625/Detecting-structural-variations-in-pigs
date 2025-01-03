#!/bin/bash

autosomes_X=... # there are two columns in autosomes_X.txt, 1st column being chromosomes, 2nd column being lengths

### SVIM with full ONT data
awk '{if ($2 > 0 || $1 ~ /#/) print $0}' .../variants.vcf | \
 awk '{if ($7 == "PASS" || $1 ~ /#/)  print $0}' | \
 egrep $autosomes_X | grep -v 'SUPPORT=1;|SUPPORT=2;|SUPPORT=3;|SUPPORT=4;' \
 > .../svim_filtered.vcf

### Picky with full ONT data
awk '{if ($2 > 0 || $1 ~ /#/) print $0}' .../picky_sv.allsv.vcf | \
 awk '{if ($7 == "PASS" || $1 ~ /#/)  print $0}' | \
 egrep $autosomes_X | \
 egrep -v 'RE=1;|RE=2;|RE=3;|RE=4;' > .../Picky_filtered.vcf
