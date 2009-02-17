cpu=30

while true
do
	for i in $(ls .|grep -v C)
	do
		if test -f $i/DONE
		then
			echo 1 > /dev/null
		else
			cd $i
			if $(ps aux|grep dna|wc -l) -lt $cpu
			then
				dna_DeBruijnAssembler *.fasta > 1 &
				touch DONE
			fi
		fi
	done
	sleep 1
done
