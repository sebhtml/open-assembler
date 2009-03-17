/*
	dna: De Novo Assembler
    Copyright (C) 2008, 2009 Sébastien Boisvert
	$Id$

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _DeBruijnAssembler
#define _DeBruijnAssembler

#define MAP_TYPE map

#define SOFTWARE_VERSION "dna-0.0.0-pre-alpha-dev"

#include"common_functions.h"
#include"SequenceDataFull.h"
#include<hash_set>
#include"VertexData.h"
#include"Read.h"
#include<map>
#include<set>
#include"GraphData.h"
#include<cstring>
#include<hash_map>
#include<vector>
#include<PairedRead.h>
#include<sstream>
using namespace __gnu_cxx;
using namespace std;

#define AMOS_FILE_NAME "contigs-amos.afg"
#define FASTA_FILE_NAME "contigs.fasta"
#define COVERAGE_FILE_NAME "contigs-coverage.txt"



class DeBruijnAssembler{
	int m_coverage_mean;
	bool m_Solexa_detected;
	string m_graphFile;
	string m_minimumCoverageParameter;
	string m_onlyFirstMer;
	string m_onlyOneStrand;
	int m_minimumCoverage;
	hash_map<uint32_t,PairedRead,hash<uint32_t> > m_paired_reads;
	int m_REPEAT_DETECTION;
	int m_wordSize;
	bool m_DEBUG;
	//CustomMap<VertexData>*m_data;
	GraphData m_data;

	SequenceDataFull*m_sequenceData;
	// map edges to reads


	// minimu contig size
	string m_assemblyDirectory;

	void build_From_Scratch(SequenceDataFull*sequenceData);
	void writeGraph();

	string pathToDNA(vector<VERTEX_TYPE>*path);

	void version2_Walker(uint64_t  a,vector<uint64_t>*b,
			vector<map<int,map<char,int> > >* currentReadPositions,
			vector<int>*repeatAnnotations);

	void Walk_In_GRAPH();

	void indexReadStrand(int readId,char strand,SequenceDataFull*sequenceData);

	void writeContig_fasta(vector<VERTEX_TYPE>*path,ofstream*file,int i);
	void writeContig_Amos(vector<map<int,map<char,int> > >*currentReadPositions,vector<VERTEX_TYPE>*path,ofstream*file,int i);



	void ThreadReads(VertexData*aData,hash_set<int>*usedReads,hash_set<int>*readsInRange,
		hash_map<int,int>*readsReadPosition,hash_map<int,int>*readsContigPositions,hash_map<int,char>*readsReadStrands,
		vector<uint64_t>*contig,
		map<uint64_t,int>*sumScores,map<uint64_t,vector<int> >*annotationsForEach,
		vector<uint64_t>*children,bool skipThoroughtCheck,
		vector<map<int,map<char,int> > >* currentReadPositions);

	void addAnnotations(VertexData*aData,hash_set<int>*usedReads,hash_set<int>*readsInRange,
		hash_map<int,int>*readsReadPosition,hash_map<int,int>*readsContigPositions,hash_map<int,char>*readsReadStrands,
vector<uint64_t>*contig,
		vector<map<int,map<char,int> > >* currentReadPositions);

	void writeContig_Coverage(vector<map<int,map<char,int> > >*currentReadPositions,vector<VERTEX_TYPE>*path,ofstream*file,int i);

	void writeContig_RepeatAnnotation(vector<int>*repeatAnnotations,int i,ofstream*file,vector<VERTEX_TYPE>*path);


	int visitVertices(VERTEX_TYPE a,set<VERTEX_TYPE>*nodes,int maxDepth,bool parents);


	bool getDebug();
	void removeDebug();
	void addDebug();
	void updateDebug();

public:
	DeBruijnAssembler();
	void setWordSize(int k);
	void setSequenceData(SequenceDataFull*sequenceData);
	void buildGraph();
	void setAssemblyDirectory(string assemblyDirectory);
	void setStrandUsage(string onlyOneStrand);
	void setMerUsage(string onlyFirstMer);
	void outputContigs();
	void debug();
	void loadParameters();
	void load_graphFrom_file();
	void loadPairedInformation();
	~DeBruijnAssembler();
	void Algorithm_Assembler_20090121();
	void setMinimumCoverage(string coverage);
};

#endif
