# 454 S. pneumoniae
nohup dna -assemblyDirectory 1020 ~/Datasets/SRA001020/sff/ETJITFZ02.sff > /dev/null &
# 454 S. cerevisiae
nohup dna -assemblyDirectory 257 ~/Datasets/SRA000257/sff/*.sff > /dev/null &
# Solexa S. cerevisiae
nohup dna -assemblyDirectory 1177 ~/Datasets/SRA001177/SRR003681.*fastq* > /dev/null &
# 454 L. tarentolae
nohup dna -assemblyDirectory tar ~/Datasets/tar-454/sff/*.sff > /dev/null &
# 454 r6 S pneumoniae
nohup dna -assemblyDirectory r6  ~/Datasets/Marc/r6/sff/E*.sff > /dev/null &
