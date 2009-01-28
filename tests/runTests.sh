rm -rf D*

for i in $(ls *.fa.fa)
do
	(
	dna_DeBruijnAssembler  -assemblyDirectory D$i $i 
	dna_Merger D$i/Contigs.fasta D$i/Assembly.fa > D$i/Merger.log
	blat $(echo $i|sed 's/.fa$//') D$i/Assembly.fa D$i/blat.psl -fastMap 
	blat $(echo $i|sed 's/.fa$//') D$i/Assembly.fa D$i/blat.blast -fastMap -out=blast
	)&
done
