make
as BinarySearch.s -o BinarySearch.o
g++  -O3 -Wall   -O3 -Wall  -o dna_DeBruijnAssembler De_Bruijn_De_Novo_Assembler_main.o Read.o DeBruijnAssembler.o Loader.o SffLoader.o VertexData.o SequenceDataFull.o CoverageDistribution.o SortedList.o GraphData.o BinarySearch.o Merger.o
