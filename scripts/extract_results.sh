
echo "Contigs Average N50 Largest Total MinimumC Minutes "

for i in $(find . -name jeudi)
do
	path=$(pwd)
	cd $i
	name=$(echo $i|sed 's/jeudi//'|sed 's/\///g'|sed 's/\.//g')
	echo  "$name & $(bash ~/denovoassembler/trunk/scripts/get_latex_line.sh) "
	cd $path
done
