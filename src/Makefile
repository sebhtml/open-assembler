all: 454dna  454dna.static

454dna.static: 454_De_Bruijn_De_Novo_Assembler_main.o Read.o DeBruijnAssembler.o Logger.o
	g++  454_De_Bruijn_De_Novo_Assembler_main.o Read.o DeBruijnAssembler.o Logger.o -O4 -Wall -static -o 454dna.static


454dna: 454_De_Bruijn_De_Novo_Assembler_main.o Read.o DeBruijnAssembler.o Logger.o
	g++  454_De_Bruijn_De_Novo_Assembler_main.o Read.o DeBruijnAssembler.o Logger.o -O4 -Wall -o 454dna


gdb: 454_De_Bruijn_De_Novo_Assembler_main.cpp Read.cpp Read.h DeBruijnAssembler.cpp DeBruijnAssembler.h  Logger.cpp Logger.h
	g++ *cpp   -Wall -g  -o gdb.out 
	gdb --exec=gdb.out --command=gdb.conf

valgrind:
	g++ 454_De_Bruijn_De_Novo_Assembler_main.cpp Read.cpp DeBruijnAssembler.cpp   -Wall -g  -o gdb.out 
	valgrind ./gdb.out -assemblyDirectory test454 fastq/Run73_FKKSCFA01.sff.fastq
test: 454dna
	./454dna fastq/ESWYO7U01.sff.fastq
	#./454dna -assemblyDirectory test454 fastq/Run73_FKKSCFA01.sff.fastq
	#./454dna -assemblyDirectory test454 fastq/EQYPFFW01.sff.fastq

simulate:
	ruby simulate-454.rb Datasets/L_major_simulation/LmjFwholegenome_20070731_V5.2.fasta > L.major.txt

dist: 
	bash dist.sh


