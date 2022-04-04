#! usr/bin/env bash

wget ftp://ftp.flybase.net/genomes/dmel/current/fasta/dmel-all-chromosome-r6.36.fasta.gz
wget ftp://ftp.flybase.net/genomes/dmel/current/fasta/md5sum.txt

grep -E chromosome md5sum.txt > check1.txt
md5sum dmel-all-chromosome-r6.36.fasta.gz > check2.txt

diff -a check1.txt check2.txt

gunzip dmel-all-chromosome-r6.36.fasta.gz
faSize dmel-all-chromosome-r6.36.fasta
