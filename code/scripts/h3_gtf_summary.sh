#! usr/bin/env bash

wget ftp://ftp.flybase.net/genomes/dmel/current/gtf/dmel-all-r6.36.gtf.gz
wget ftp://ftp.flybase.net/genomes/dmel/current/gtf/md5sum.txt

md5sum dmel-all-r6.36.gtf.gz > check3.txt

diff -a check3.txt md5sum.txt

gunzip dmel-all-r6.36.gtf.gz

gawk -F '\t' '{print $3}' dmel-all-r6.36.gtf | sort | uniq -c | sort -k1,1nr > unique_col3_sorted.txt

gawk -F '\t' '$3=="gene"{print $1}'  dmel-all-r6.36.gtf | grep -P '^[4XY23][LR]?$' | sort | uniq -c > chromosome_gene_count.txt
