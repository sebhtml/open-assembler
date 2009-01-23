rm -rf D*

for i in $(ls *.fa.fa)
do
	dna_DeBruijnAssembler -assemblyDirectory D$i $i &
done
