for i in $(ls|grep D)
do
	echo $i
	cat $i/Log.txt|ruby parse.rb
	echo ""
done
