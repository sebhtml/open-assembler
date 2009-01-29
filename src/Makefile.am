bin_PROGRAMS = dna_DeBruijnAssembler dna_sffPairedToFasta dna_DeBruijnSplitter dna_Merger dna_Scaffolder dna_fastaToAMOS

dna_DeBruijnAssembler_SOURCES = De_Bruijn_De_Novo_Assembler_main.cpp  DeBruijnAssembler.h  Read.cpp \
DeBruijnAssembler.cpp                               Read.h Loader.cpp Loader.h  SffLoader.h SffLoader.cpp \
VertexData.h VertexData.cpp  SequenceDataFull.h SequenceDataFull.cpp

dna_DeBruijnSplitter_SOURCES = splitter_main.cpp Loader.cpp Read.cpp DeBruijnAssembler.cpp VertexData.cpp SequenceDataFull.cpp SffLoader.cpp \
LightVertex.h LightVertex.cpp

dna_Merger_SOURCES = merger_main.cpp  Read.cpp DeBruijnAssembler.cpp Loader.cpp SffLoader.cpp VertexData.cpp SequenceDataFull.cpp

dna_Scaffolder_SOURCES = scaffolder_main.cpp

dna_sffPairedToFasta_SOURCES = SffLoader.cpp SffLoader.h Read.h Read.cpp DeBruijnAssembler.cpp DeBruijnAssembler.h sffPairedToFasta_main.cpp VertexData.h VertexData.cpp CustomMap.hpp  Loader.cpp Loader.h SequenceDataFull.h SequenceDataFull.cpp

dna_fastaToAMOS_SOURCES = fasta_to_afg_main.cpp  Read.cpp DeBruijnAssembler.cpp Loader.cpp SffLoader.cpp VertexData.cpp SequenceDataFull.cpp
