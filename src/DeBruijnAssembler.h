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
#include"VertexData.h"
#include"Read.h"
#include<map>
#include<set>
#include"GraphData.h"
#include<cstring>
#include<vector>
#include<sstream>
using namespace std;

#define AMOS_FILE_NAME "contigs-amos.afg"
#define FASTA_FILE_NAME "contigs.fasta"
#define COVERAGE_FILE_NAME "contigs-coverage.txt"




class DeBruijnAssembler{
	int m_coverage_mean;
	int m_carry_forward_offset;
	string m_graphFile;
	string m_minimumCoverageParameter;
	int m_minimumCoverage;
	int m_REPEAT_DETECTION;
	int m_wordSize;
	bool m_DEBUG;
	//CustomMap<VertexData>*m_data;
	GraphData m_data;

	SequenceDataFull*m_sequenceData;
	// map edges to reads


	// minimu contig size
	ostream*m_cout;
	string m_assemblyDirectory;



	bool is_d_Threading(AnnotationElement*annotation,vector<map<int,map<char,int> > >*currentReadPositions,vector<VERTEX_TYPE>*path,map<int,int>*usedReads,bool beforeAdding);

	void build_From_Scratch(SequenceDataFull*sequenceData);
	void writeGraph();

	void contig_From_SINGLE(vector<map<int,map<char,int> > >*currentReadPositions,vector<VERTEX_TYPE>*path,vector<VERTEX_TYPE>*newSources,vector<int>*repeatAnnotations,VERTEX_TYPE source);
	vector<VERTEX_TYPE> getWalk(VERTEX_TYPE prefix,vector<VERTEX_TYPE>*path,int length,vector<map<int,map<char,int > > >*currentReadPositions);

	string pathToDNA(vector<VERTEX_TYPE>*path);

	vector<VERTEX_TYPE> nextVertices(vector<VERTEX_TYPE>*path,vector<map<int,map<char,int> > >*currentReadPositions,vector<VERTEX_TYPE>*newSources,map<int,int>*usedReads,map<int,char>*readStrands,
		set<int>*theReadsThatAreRelevant);
	bool DETECT_BUBBLE(vector<VERTEX_TYPE>*path,VERTEX_TYPE a,VERTEX_TYPE b);

	VERTEX_TYPE reverseComplement_VERTEX(VERTEX_TYPE a);

	//bool addNewContig(vector<vector<VERTEX_TYPE> >*newContigs,vector<VERTEX_TYPE>*newContig,int currentContigId,int otherContigId,set<int>*contigsProcessed);

	void Walk_In_GRAPH();


	void indexReadStrand(int readId,char strand,SequenceDataFull*sequenceData);

	void writeContig_fasta(vector<VERTEX_TYPE>*path,ofstream*file,int i);
	void writeContig_Amos(vector<map<int,map<char,int> > >*currentReadPositions,vector<VERTEX_TYPE>*path,ofstream*file,int i);
	int DFS_watch(VERTEX_TYPE a,int color);

	void writeContig_Coverage(vector<map<int,map<char,int> > >*currentReadPositions,vector<VERTEX_TYPE>*path,ofstream*file,int i);

	void writeContig_RepeatAnnotation(vector<int>*repeatAnnotations,int i,ofstream*file,vector<VERTEX_TYPE>*path);


	int visitVertices(VERTEX_TYPE a,set<VERTEX_TYPE>*nodes,int maxDepth,bool parents);
public:
	DeBruijnAssembler(ostream*m_cout);
	void setWordSize(int k);
	void setSequenceData(SequenceDataFull*sequenceData);
	void buildGraph();
	void setAssemblyDirectory(string assemblyDirectory);
	void outputContigs();
	void debug();
	void loadParameters();
	void load_graphFrom_file();
	~DeBruijnAssembler();
	void Algorithm_Assembler_20090121();
	void setMinimumCoverage(string coverage);

	static int m_WordSize;
};

#endif
