
llvmc -O5 $(ls *.cpp|grep -v spli|grep -v sca |grep -v mer|grep -v sff|grep -v fasta)  -o dna_gdb
