# 230000000
# 230000000/1024/1024 ->  219 GB
# 606 fastq files

ulimit -v  230000000

nohup dna -buckets 1000000000 -assemblyDirectory HumanGenome $(ls ~/Datasets/SRA000271/*fastq|head -n100) > /dev/null &
