/*
	454 De Novo Assembler
    Copyright (C) 2008 Sébastien Boisvert

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

#include"Read.h"
#include<map>
#include<set>
#include<list>
#include<cstring>
#include<vector>
#include"Logger.hpp"
#include<sstream>
using namespace std;



class DeBruijnAssembler{
	string m_graphFile;
	int m_windowSize;

	static int m_NUCLEOTIDE_A;
	static int m_NUCLEOTIDE_T;
	static int m_NUCLEOTIDE_C;
	static int m_NUCLEOTIDE_G;
	int m_minimumCoverage;
	int m_wordSize;

	vector<Read*>*m_reads;

	// map edges to reads
	map<VERTEX_TYPE,map<VERTEX_TYPE, vector<int > > > m_graph;

	string m_useCache;
	// parents..
	map<VERTEX_TYPE,set<VERTEX_TYPE> > m_vertex_parents;

	set<VERTEX_TYPE> m_solidMers;

	// stock directly a vector of vertices
	vector<vector<VERTEX_TYPE> > m_contig_paths;

	// minimu contig size
	int m_minimumContigSize;
	Logger*m_cout;
	string m_assemblyDirectory;

	void load_graphFrom_file();
	void build_From_Scratch();
	void writeGraph();

	vector<VERTEX_TYPE> contig_From_SINGLE(VERTEX_TYPE prefix,map<VERTEX_TYPE,int> visits,vector<VERTEX_TYPE> path);
	vector<VERTEX_TYPE> getWalk(VERTEX_TYPE prefix,vector<VERTEX_TYPE>path,int length);
	vector<VERTEX_TYPE> removeBubblesAndTips(vector<VERTEX_TYPE> vertices,vector<VERTEX_TYPE>path);
	char getLastSymbol(VERTEX_TYPE i);
	string pathToDNA(vector<VERTEX_TYPE> path);

	string idToWord(VERTEX_TYPE i,int wordSize);

	vector<VERTEX_TYPE> nextVertices(VERTEX_TYPE prefix,vector<VERTEX_TYPE> path,int l,int C);
	bool passFilterCoverage(vector<VERTEX_TYPE> path,int l,int C);
	vector<vector<VERTEX_TYPE> >Filter_Remove_Smaller_Duplicates(vector<vector<VERTEX_TYPE> > contigs);
public:
	DeBruijnAssembler(Logger*m_cout);
	void setMinimumCoverage(int m);
	void setWordSize(int k);
	void buildGraph(vector<Read*>*reads);
	void setAssemblyDirectory(string assemblyDirectory);
	void setMinimumContigSize(int minimumContigSize);
	void outputContigs();
	void setWindowSize(int windowSize);
	void setUseCache(string useCache);
	void run_New_Algorithm_Assembler_20090102();


	static string reverseComplement(string a);
	static unsigned long int wordId(const char*a);
	static char complement(char a);
};

#endif
