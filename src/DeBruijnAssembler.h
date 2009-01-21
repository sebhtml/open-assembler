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

#include"Read.h"
#include<map>
#include"Logger.hpp"
#include<set>
#include<list>
#include <utility>
#include<cstring>
#include<vector>
#include<sstream>
using namespace std;

#ifndef _DeBruijnAssembler
#define _DeBruijnAssembler

// TODO: clean obsolete functions

class DeBruijnAssembler{
	string m_graphFile;
	// for removed k-mers
	// initial state
	// TODO: remove
 
	int m_COLOR_IN_PROGRESS_ASSEMBLER;
	int m_COLOR_NOT_ASSEMBLED;
	int m_COLOR_IN_PROGRESS_SOLVER;
	int m_COLOR_DISCARDED;
	int m_COLOR_REPEAT;
	int m_windowSize;

	int m_NUCLEOTIDE_A;
	int m_NUCLEOTIDE_T;
	int m_NUCLEOTIDE_C;
	int m_NUCLEOTIDE_G;
	int m_minimumCoverage;
	int m_minimumQuality;
	int m_wordSize;

	// TODO: remove, deprecated...
	// stat of all (k+1) mers (edges)
	// used to eliminate edges
	map<unsigned long int,map<unsigned long int,int > > m_edge_states;

	vector<Read*>*m_reads;
	// contig id for the next one

	// map (k+1) mers to reads
	map<unsigned long int,vector<int> > m_mer_to_read_table;

	// map edges to reads
	map<unsigned long int,map<unsigned long int , vector<int > > > m_graph;

	string m_useCache;
	// parents..
	map<unsigned long int,vector<unsigned long int> > m_vertex_parents;

	// k-mer in graph
	set<unsigned long int> m_graph_mers;

	// stock directly a vector of vertices
	vector<vector<unsigned long int> > m_contig_paths;

	// minimu contig size
	int m_minimumContigSize;
	Logger*m_cout;
	ofstream*m_solver_log;
	string m_assemblyDirectory;
	set<unsigned long int> m_solidMers;

	void load_graphFrom_file();
	void build_From_Scratch();
	void writeGraph();

	// TODO: cleanup private methods when the code will work..
	vector<unsigned long int> contig_From_SINGLE(unsigned long int prefix,map<unsigned long int,int> visits,vector<unsigned long int> path);
	vector<unsigned long int> getWalk(unsigned long int prefix,vector<unsigned long int>path,int length);
	vector<unsigned long int> removeBubbles(vector<unsigned long int> vertices,vector<unsigned long int>path);
	bool passFilterZ(vector<unsigned long int>contig,int l,int C,int maxZ);
	char getLastSymbol(unsigned long int i);
	string pathToDNA(vector<unsigned long int> path);
	void showPath(vector<unsigned long int> path);
	void makeContigPath(unsigned long int prefix,int contig);
	string getTrivialPath(unsigned long int pathVertex,unsigned long int vertex);
	vector<unsigned long int> getHighQualityMers(int readId);
	void solveMultiPath(unsigned long int vertex,int depth);
	string reverseComplement(string a);
	char complement(char a);
	unsigned long int wordId(const char*a);

	string idToWord(unsigned long int i,int wordSize);
	vector<unsigned long int> getNeighbours(unsigned long int vertex);
	unsigned long int toKPlusOneMer(unsigned long int prefix,unsigned long int suffix);

	vector<vector<unsigned long int> > contigs_From(unsigned long int prefix,map<unsigned long int,int> visits,int depth,vector<unsigned long int>previousEdges,map<int,map<unsigned long int, int> > readUsages);
	bool passFilter(vector<unsigned long int>,int l,int C,int Z);
	vector<unsigned long int> nextVertices(unsigned long int prefix,vector<unsigned long int> path,int l,int C);
	bool passFilterCoverage(vector<unsigned long int> path,int l,int C);
	vector<vector<unsigned long int> >removeSmallContigs(vector<vector<unsigned long int> > contigs);
	vector<vector<unsigned long int> >Filter_Remove_Smaller_Duplicates(vector<vector<unsigned long int> > contigs);
	unsigned long int toPrefix(unsigned long int edge);
	unsigned long int toSuffix(unsigned long int edge);
	void run_New_Algorithm_Assembler();
	void run_Assembler();
public:
	DeBruijnAssembler(Logger*m_cout);
	void setMinimumCoverage(int m);
	void setMinimumQuality(int q);
	void setWordSize(int k);
	void buildGraph(vector<Read*>*reads);
	void setAssemblyDirectory(string assemblyDirectory);
	void setMinimumContigSize(int minimumContigSize);
	void outputContigs();
	void setWindowSize(int windowSize);
	void setUseCache(string useCache);
	void run_New_Algorithm_Assembler_20090102();
};

#endif
