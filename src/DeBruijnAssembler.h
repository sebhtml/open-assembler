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
	int m_minimumCoverage_for_walk;
	int m_wordSize;
	double m_threshold;
	bool m_pairedAvailable;
	uint64_t m_buckets;
	int m_longReadMode_threshold;
	bool m_longReadAvailable;
	string m_pairedInfoFile;
	int m_default_window;
	uint64_t m_solidMers;
	CustomMap<VertexData>*m_data;

	SequenceDataFull*m_sequenceData;
	// map edges to reads

	// stock directly a vector of vertices
	vector<vector<VERTEX_TYPE> > m_contig_paths;

	// minimu contig size
	int m_minimumContigSize;
	ostream*m_cout;
	string m_assemblyDirectory;

	void load_graphFrom_file();
	void build_From_Scratch(SequenceDataFull*sequenceData);
	void writeGraph();

	vector<VERTEX_TYPE> optimizedNextVertices(vector<VERTEX_TYPE>*path,map<VERTEX_TYPE,int>*visits,int C,int l);
	void contig_From_SINGLE(map<int,map<int,map<char,int> > >*currentReadPositions,vector<VERTEX_TYPE>*path,vector<VERTEX_TYPE>*newSources,bool repeat_aware);
	vector<VERTEX_TYPE> getWalk(VERTEX_TYPE prefix,vector<VERTEX_TYPE>*path,int length,map<int,map<int,map<char,int > > >*currentReadPositions);
	vector<VERTEX_TYPE> removeBubblesAndTips(vector<VERTEX_TYPE> vertices,vector<VERTEX_TYPE>*path,map<int,map<int,map<char,int> > >*currentReadPositions);
	char getLastSymbol(VERTEX_TYPE i);
	string pathToDNA(vector<VERTEX_TYPE>*path);

	vector<VERTEX_TYPE> nextVertices(vector<VERTEX_TYPE>*path,map<int,map<int,map<char,int> > >*currentReadPositions,bool repeat_aware);
	vector<VERTEX_TYPE> nextVertices_SAFE(vector<VERTEX_TYPE>*path,map<int,map<int,map<char,int> > >*currentReadPositions,bool repeat_aware);
	int recThreading(VERTEX_TYPE prefix,VERTEX_TYPE suffix,map<int,int>*allowedReads);


	VERTEX_TYPE reverseComplement_VERTEX(VERTEX_TYPE a);

	vector<int> not_Processed(vector<int>contigs,set<int>processed,int self);

	bool addNewContig(vector<vector<VERTEX_TYPE> >*newContigs,vector<VERTEX_TYPE>*newContig,int currentContigId,int otherContigId,set<int>*contigsProcessed);

	vector<vector<VERTEX_TYPE> >Remove_Small_Contigs(vector<vector<VERTEX_TYPE> > contigs);
	vector<vector<VERTEX_TYPE> > ExtendReverseComplements(vector<vector<VERTEX_TYPE> > contigs);
	vector<vector<VERTEX_TYPE> >Filter_Remove_Smaller_Duplicates(vector<vector<VERTEX_TYPE> > contigs);


	vector<vector<VERTEX_TYPE> > Filter_Remove_Smaller_Duplicates_Cached(vector<vector<VERTEX_TYPE > > largeContigs);
	void Walk_In_GRAPH();

	vector<AnnotationElement>annotationsWithCurrent(vector<AnnotationElement>*elements,map<int,map<int,map<char,int> > >*currentReadPositions);

	void indexReadStrand(int readId,char strand,SequenceDataFull*sequenceData,CustomMap<int>*solidMers);


	vector<VERTEX_TYPE> nextVertices_OLD(vector<VERTEX_TYPE>*path,vector<map<char,int> >*currentReadPositions);


	void writeContig_fasta(vector<VERTEX_TYPE>*path,ofstream*file,int i);
	void writeContig_Amos(map<int,map<int,map<char,int> > >*currentReadPositions,vector<VERTEX_TYPE>*path,ofstream*file,int i);
public:
	DeBruijnAssembler(ostream*m_cout);
	void setPairedInfo(string a);
	void setWordSize(int k);
	void buildGraph(SequenceDataFull*sequenceData);
	void setAssemblyDirectory(string assemblyDirectory);
	void setMinimumContigSize(int minimumContigSize);
	void outputContigs();

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
