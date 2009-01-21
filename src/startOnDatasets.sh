killall 454dna
for i in $(find Datasets/ -name fastq)
do
	echo $i
	cd $i
	rm nohup.out
	nohup 454dna *fastq -assemblyDirectory 20081212_Assembly
	cd /home/boiseb01/bin/denovoassembler/trunk/src
done
