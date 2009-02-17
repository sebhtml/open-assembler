bin_PROGRAMS = dna_DeBruijnAssembler dna_sffPairedToFasta dna_Merger dna_Scaffolder dna_fastaToAMOS 
#dna_DeBruijnSplitter 

dna_DeBruijnAssembler_SOURCES = De_Bruijn_De_Novo_Assembler_main.cpp  DeBruijnAssembler.h  Read.cpp \
DeBruijnAssembler.cpp                               Read.h Loader.cpp Loader.h  SffLoader.h SffLoader.cpp \
VertexData.h VertexData.cpp  SequenceDataFull.h SequenceDataFull.cpp CoverageDistribution.h CoverageDistribution.cpp \
SortedList.h SortedList.cpp GraphData.h GraphData.cpp BinarySearch.h BinarySearch.cpp

#dna_DeBruijnSplitter_SOURCES = splitter_main.cpp Loader.cpp Read.cpp DeBruijnAssembler.cpp VertexData.cpp SequenceDataFull.cpp SffLoader.cpp \
#LightVertex.h LightVertex.cpp CoverageDistribution.h CoverageDistribution.cpp  \
#SortedList.h SortedList.cpp GraphData.h GraphData.cpp BinarySearch.h BinarySearch.cpp

dna_Merger_SOURCES = merger_main.cpp  Read.cpp DeBruijnAssembler.cpp Loader.cpp SffLoader.cpp VertexData.cpp SequenceDataFull.cpp  CoverageDistribution.h CoverageDistribution.cpp \
SortedList.h SortedList.cpp GraphData.h GraphData.cpp BinarySearch.h BinarySearch.cpp

dna_Scaffolder_SOURCES =   CoverageDistribution.h CoverageDistribution.cpp scaffolder_main.cpp  Read.cpp DeBruijnAssembler.cpp Loader.cpp SffLoader.cpp  SequenceDataFull.cpp \
 SortedList.h SortedList.cpp VertexData.h VertexData.cpp   GraphData.h GraphData.cpp BinarySearch.h BinarySearch.cpp

dna_sffPairedToFasta_SOURCES = SffLoader.cpp SffLoader.h Read.h Read.cpp DeBruijnAssembler.cpp DeBruijnAssembler.h sffPairedToFasta_main.cpp  \
CustomMap.hpp  Loader.cpp Loader.h SequenceDataFull.h SequenceDataFull.cpp  CoverageDistribution.h CoverageDistribution.cpp \
SortedList.h SortedList.cpp VertexData.cpp VertexData.h  GraphData.h GraphData.cpp BinarySearch.h BinarySearch.cpp

dna_fastaToAMOS_SOURCES = fasta_to_afg_main.cpp  Read.cpp DeBruijnAssembler.cpp Loader.cpp SffLoader.cpp VertexData.cpp SequenceDataFull.cpp  CoverageDistribution.h CoverageDistribution.cpp \
SortedList.h SortedList.cpp  GraphData.h GraphData.cpp BinarySearch.h BinarySearch.cpp
