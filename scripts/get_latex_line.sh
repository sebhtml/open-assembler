#!/bin/bash

averageLength=$(perl -e "printf '%u',$(getN50 4JoinedContigs.fasta|head -n2|tail -n1|awk '{print $2}')")
minimumCoverage=$( grep Covera Parameters.txt |grep Min|awk '{print $2}')
minutes=$(ruby -e "require 'time';a=Time.parse((File.open 'START').read);b=Time.parse((File.open 'END').read);puts ((b-a)/60).to_i.to_s")
#echo "Contigs Average N50 Largest Total"
echo "$(getN50 4JoinedContigs.fasta|head -n1|awk '{print $2}') & $averageLength & $(getN50 4JoinedContigs.fasta|tail -n1|awk '{print $2}') & $(getlengths 4JoinedContigs.fasta |awk '{print $2}'|sort -n|tail -n1) & $(($(cat 4JoinedContigs.fasta|grep -v '>'|wc |awk '{print $3}')-$(cat 4JoinedContigs.fasta|grep -v '>'|wc |awk '{print $1}'))) & $minimumCoverage &  $minutes \\\\"
