#!/bin/bash

echo ">betp" > out1.fasta
grep -v "*" out | grep -v "Align" |  awk '{FS=""; if($1!="") {if(a) printf "%s", a}; FS=" "} {a=$2} END {print}' >> out1.fasta
echo ">leut" > out2.fasta
grep -v "*" out | grep -v "Align" |  awk '{FS=""; if($1=="") {if(a) printf "%s", a}; FS=" "} {a=$2} END {print}' >> out2.fasta
cat out1.fasta out2.fasta > out.fasta
rm out
