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

#include"SequenceDataFull.h"
#include"VertexData.h"
#include"Read.h"
#include<map>
#include"CustomMap.hpp"
#include<set>
#include<cstring>
#include<vector>
#include<sstream>
using namespace std;



class DeBruijnAssembler{
	string m_graphFile;
	string m_minimumCoverageParameter;
	int m_minimumCoverage;
	int m_wordSize;
	double m_threshold;
	int m_Coverage_From_DepletionCurve;
	bool m_pairedAvailable;
	uint64_t m_buckets;
	string m_pairedInfoFile;
	bool m_DEBUG;
	uint64_t m_solidMers;
	CustomMap<VertexData>*m_data;

	SequenceDataFull*m_sequenceData;
	// map edges to reads


	// minimu contig size
	ostream*m_cout;
	string m_assemblyDirectory;

	void load_graphFrom_file();
	void build_From_Scratch(SequenceDataFull*sequenceData);
	void writeGraph();

	void contig_From_SINGLE(vector<map<int,map<char,int> > >*currentReadPositions,vector<VERTEX_TYPE>*path,vector<VERTEX_TYPE>*newSources);
	vector<VERTEX_TYPE> getWalk(VERTEX_TYPE prefix,vector<VERTEX_TYPE>*path,int length,vector<map<int,map<char,int > > >*currentReadPositions);

	char getLastSymbol(VERTEX_TYPE i);
	string pathToDNA(vector<VERTEX_TYPE>*path);

	vector<VERTEX_TYPE> nextVertices(vector<VERTEX_TYPE>*path,vector<map<int,map<char,int> > >*currentReadPositions,vector<VERTEX_TYPE>*newSources);
	bool DETECT_BUBBLE(vector<VERTEX_TYPE>*path,VERTEX_TYPE a,VERTEX_TYPE b);

	VERTEX_TYPE reverseComplement_VERTEX(VERTEX_TYPE a);

	//bool addNewContig(vector<vector<VERTEX_TYPE> >*newContigs,vector<VERTEX_TYPE>*newContig,int currentContigId,int otherContigId,set<int>*contigsProcessed);

	void Walk_In_GRAPH();


	void indexReadStrand(int readId,char strand,SequenceDataFull*sequenceData,CustomMap<int>*solidMers);

	void writeContig_fasta(vector<VERTEX_TYPE>*path,ofstream*file,int i);
	void writeContig_Amos(vector<map<int,map<char,int> > >*currentReadPositions,vector<VERTEX_TYPE>*path,ofstream*file,int i);
public:
	DeBruijnAssembler(ostream*m_cout);
	void setPairedInfo(string a);
	void setWordSize(int k);
	void buildGraph(SequenceDataFull*sequenceData);
	void setAssemblyDirectory(string assemblyDirectory);
	void setMinimumContigSize(int minimumContigSize);
	void outputContigs();
	void debug();

	~DeBruijnAssembler();
	void setBuckets(uint64_t buckets);
	void Algorithm_Assembler_20090121();
	void setMinimumCoverage(string coverage);

	static string reverseComplement(string a);
	static VERTEX_TYPE wordId(const char*a);
	static string idToWord(VERTEX_TYPE i,int wordSize);
	static char complement(char a);
	static int m_WordSize;
	static int m_NUCLEOTIDE_A;
	static int m_NUCLEOTIDE_T;
	static int m_NUCLEOTIDE_C;
	static int m_NUCLEOTIDE_G;
};

#endif
