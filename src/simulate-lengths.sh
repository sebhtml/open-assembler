for i in 30 36 50 100 200 230 300 400 500 600 800 1000 1200 1500 1700  2000
do 
	#ruby ../src/simulate-reads.rb Streptococcus-pneumoniae-R6.fasta $i > LeNGTH_$i.fasta  
	#DNA.rb LeNGTH_$i.fasta -directory l$i &> /dev/null
	echo "$i $(getlengths l$i/3MergedContigs.fasta |wc -l)"
done
