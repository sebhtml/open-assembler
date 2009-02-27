bin_PROGRAMS =   dna_ConvertSffToFasta dna_MergeContigs dna_BuildGraph dna_ScaffoldContigs  dna_ExtractContigs  dna_KeepLargeContigs \
	dna_ConvertFastaToAfg dna_JoinContigs # dna_SplitReads dna_ConvertFastaToSql

# dna_ConvertFastaToSql_SOURCES =  dump_sql_main.cpp  common_functions.h common_functions.cpp

dna_BuildGraph_SOURCES = De_Bruijn_De_Novo_Assembler_main.cpp  DeBruijnAssembler.h  Read.cpp \
DeBruijnAssembler.cpp                               Read.h Loader.cpp Loader.h  SffLoader.h SffLoader.cpp \
VertexData.h VertexData.cpp  SequenceDataFull.h SequenceDataFull.cpp CoverageDistribution.h CoverageDistribution.cpp \
SortedList.h SortedList.cpp GraphData.h GraphData.cpp BinarySearch.h BinarySearch.cpp Merger.h Merger.cpp    common_functions.h common_functions.cpp

dna_ExtractContigs_SOURCES = extract_contigs_main.cpp  DeBruijnAssembler.h  Read.cpp \
DeBruijnAssembler.cpp                               Read.h Loader.cpp Loader.h  SffLoader.h SffLoader.cpp \
VertexData.h VertexData.cpp  SequenceDataFull.h SequenceDataFull.cpp CoverageDistribution.h CoverageDistribution.cpp \
SortedList.h SortedList.cpp GraphData.h GraphData.cpp BinarySearch.h BinarySearch.cpp Merger.h Merger.cpp    common_functions.h common_functions.cpp



dna_KeepLargeContigs_SOURCES = keepLargeContigs_main.cpp Read.h Read.cpp Loader.cpp Loader.h common_functions.h common_functions.cpp SffLoader.h  SffLoader.cpp

dna_SplitReads_SOURCES =  LightVertex.h LightVertex.cpp  \
SortedList.h SortedList.cpp GraphData.h GraphData.cpp BinarySearch.h BinarySearch.cpp GraphDataLight.h GraphDataLight.cpp \
splitter_main.cpp Loader.cpp Read.cpp DeBruijnAssembler.cpp VertexData.cpp SequenceDataFull.cpp SffLoader.cpp \
CoverageDistribution.h CoverageDistribution.cpp  

dna_MergeContigs_SOURCES = merger_main.cpp  Read.cpp DeBruijnAssembler.cpp Loader.cpp SffLoader.cpp VertexData.cpp SequenceDataFull.cpp  CoverageDistribution.h CoverageDistribution.cpp \
SortedList.h SortedList.cpp GraphData.h GraphData.cpp BinarySearch.h BinarySearch.cpp Merger.h Merger.cpp   common_functions.h common_functions.cpp


dna_ScaffoldContigs_SOURCES =   CoverageDistribution.h CoverageDistribution.cpp scaffolder_main.cpp  Read.cpp DeBruijnAssembler.cpp Loader.cpp SffLoader.cpp  SequenceDataFull.cpp \
 SortedList.h SortedList.cpp VertexData.h VertexData.cpp   GraphData.h GraphData.cpp BinarySearch.h BinarySearch.cpp   common_functions.h common_functions.cpp


dna_ConvertSffToFasta_SOURCES = SffLoader.cpp SffLoader.h Read.h Read.cpp DeBruijnAssembler.cpp DeBruijnAssembler.h sffPairedToFasta_main.cpp  \
CustomMap.hpp  Loader.cpp Loader.h SequenceDataFull.h SequenceDataFull.cpp  CoverageDistribution.h CoverageDistribution.cpp \
SortedList.h SortedList.cpp VertexData.cpp VertexData.h  GraphData.h GraphData.cpp BinarySearch.h BinarySearch.cpp   common_functions.h common_functions.cpp


dna_ConvertFastaToAfg_SOURCES = fasta_to_afg_main.cpp  Read.cpp DeBruijnAssembler.cpp Loader.cpp SffLoader.cpp VertexData.cpp SequenceDataFull.cpp  CoverageDistribution.h CoverageDistribution.cpp \
SortedList.h SortedList.cpp  GraphData.h GraphData.cpp BinarySearch.h BinarySearch.cpp   common_functions.h common_functions.cpp


dna_JoinContigs_SOURCES =  DeBruijnAssembler.h  Read.cpp module_join_main.cpp \
DeBruijnAssembler.cpp                               Read.h Loader.cpp Loader.h  SffLoader.h SffLoader.cpp \
VertexData.h VertexData.cpp  SequenceDataFull.h SequenceDataFull.cpp CoverageDistribution.h CoverageDistribution.cpp \
SortedList.h SortedList.cpp GraphData.h GraphData.cpp BinarySearch.h BinarySearch.cpp Merger.h Merger.cpp   common_functions.h common_functions.cpp


