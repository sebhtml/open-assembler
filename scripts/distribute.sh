#!/bin/bash

# distribute computation on a multi-core machine.

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
			processes=$(ps aux|grep dna|wc -l)
			if test $processes -lt $cpu
			then	
				echo "Starting $i ($processes / $cpu)"
				dna_DeBruijnAssembler *.fasta > 1 &
				touch DONE
				sleep 0.01
			fi
			cd ..
		fi
	done
	sleep 1
done
