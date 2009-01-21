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
#include"Logger.h"
#include<set>
#include<list>
#include <utility>
#include<cstring>
#include<vector>
#include<sstream>
using namespace std;

#ifndef _DeBruijnAssembler
#define _DeBruijnAssembler


class DeBruijnAssembler{
	// for removed k-mers
	int m_COLOR_DISCARDED; 
	// initial state
	int m_COLOR_NOT_ASSEMBLED;
	int m_COLOR_IN_PROGRESS_ASSEMBLER;
	int m_COLOR_IN_PROGRESS_SOLVER;
	int m_NUCLEOTIDE_A;
	int m_NUCLEOTIDE_T;
	int m_NUCLEOTIDE_C;
	int m_NUCLEOTIDE_G;
	int m_minimumCoverage;
	int m_minimumQuality;
	int m_wordSize;
	int m_last_vertices_size;
	set<unsigned long int> m_repeat_vertices;
	unsigned long int m_repeat_vertex;
	bool m_DETECTED_REPEAT;
	vector<Read*>*m_reads;
	// contig id for the next one
	int m_contig_id;
	// contig for each k-mer
	map<unsigned long int,int> m_colors;
	// the graph
	map<unsigned long int,set<unsigned long int> > m_graph;
	// index (k+1)-mer to reads
	map<unsigned long int,vector<int> > m_mer_to_read_table;
	// k-mer in graph
	set<unsigned long int> m_graph_mers;
	// where a contig begins TODO: stock directly a vector of vertices...
	map<int,list<unsigned long int> > m_contig_paths;
	// minimu contig size
	int m_minimumContigSize;
	Logger*m_cout;

	char getLastSymbol(unsigned long int i);

	string getTrivialPath(unsigned long int pathVertex,unsigned long int vertex);
	vector<unsigned long int> getHighQualityMers(int readId);
	void perform_Assembly(unsigned long int vertex);
	void solveMultiPath(unsigned long int vertex,int depth);
	string reverseComplement(string a);
	char complement(char a);
	unsigned long int wordId(const char*a);

	string idToWord(unsigned long int i,int wordSize);
	void printBinary(unsigned long long int i,int wordSize);
	vector<unsigned long int> getNeighbours(unsigned long int vertex);
public:
	DeBruijnAssembler(Logger*m_cout);
	void setMinimumCoverage(int m);
	void setMinimumQuality(int q);
	void setWordSize(int k);
	void buildGraph(vector<Read*>*reads);
	void run_Assembler();
	void setMinimumContigSize(int minimumContigSize);
	void outputContigs(ofstream * f);

};

#endif
