rm -rf D*

for i in $(ls *.fa.fa)
do
	(
	dna_DeBruijnAssembler  -assemblyDirectory D$i $i  > D$i.log
	cd D$i
	bash CreateBank.sh
	bash Merge.sh
	cd ..
	blat $(echo $i|sed 's/.fa$//') D$i/Merger.fasta D$i/blat.psl -fastMap 
	blat $(echo $i|sed 's/.fa$//') D$i/Merger.fasta D$i/blat.blast -fastMap -out=blast
	)&
done
