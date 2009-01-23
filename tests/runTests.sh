rm -rf D*

for i in $(ls *.fa.fa)
do
	(
	dna_DeBruijnAssembler -assemblyDirectory D$i $i 
	blat $(echo $i|sed 's/.fa$//') D$i/Contigs.fa D$i/blat.psl -fastMap 
	)&
done
