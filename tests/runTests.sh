rm -rf D*

for i in $(ls *.fa.fa)
do
	(
	dna_DeBruijnAssembler  -assemblyDirectory D$i $i 
	cp $i D$i/fasta.reads
	tarchive2amos -s D$i/fasta.reads -o D$i/reads.afg
	bank-transact  -m D$i/reads.afg -b D$i/bank -c
	bank-transact -m D$i/Assembly.afg -b D$i/bank
	dna_Merger D$i/fasta.contigs D$i/Assembly.fa > D$i/Merger.log
	#blat $(echo $i|sed 's/.fa$//') D$i/Assembly.fa D$i/blat.psl -fastMap 
	#blat $(echo $i|sed 's/.fa$//') D$i/Assembly.fa D$i/blat.blast -fastMap -out=blast
	)&
done
