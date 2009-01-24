for i in $(ls|grep D)
do
	echo $i
	cat $i/Assembly.fa|grep '>'|wc -l
	echo ""
done
