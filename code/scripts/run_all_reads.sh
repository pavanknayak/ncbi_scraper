#!/bin/bash

for filename in *.fastq
do
echo  ${filename} >>number_of_Ns.txt
echo "Total number of reads:">>number_of_Ns.txt
T=`echo $(cat ${filename} |wc -l)/4|bc`
echo $T>>number_of_Ns.txt
echo "Number of reads with 10 consecutive N's">>number_of_Ns.txt
N=`echo $(grep NNNNNNNNNN ${filename} | wc -l)`
echo $N>>number_of_Ns.txt
echo "percentage of bad reads in ${filename}:">>number_of_Ns.txt
P=`echo $(echo "scale=7;($N/$T)*100"|bc)`
echo $P>>number_of_Ns.txt
done


