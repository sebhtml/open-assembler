rm -rf D*

for i in $(ls *.fa.fa)
do
	(
	dna_DeBruijnAssembler -wordSize 15 -assemblyDirectory D$i $i 
	dna_Merger D$i/Contigs.fa D$i/Assembly.fa > D$i/Merger.log
	blat $(echo $i|sed 's/.fa$//') D$i/Assembly.fa D$i/blat.psl -fastMap 
	)&
done
