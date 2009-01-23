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

#ifndef _SequenceData
#define _SequenceData

#include<string>
#include<map>
#include<vector>
#include"Read.h"
#include<stdint.h>
using namespace std;

// this class gives access to sequence data without loading all the files
// only 2 files are loaded at any point
// therefore, fore paired end reads it is necessary to have the exact file order..
class SequenceData{
	map<string,int > m_file_first;
	map<string,int > m_file_last;
	vector<string> m_inputFiles;
	ostream*m_cout;
	uint64_t m_bases;
	int m_reads;

	int m_file_1;
	vector<Read*> m_reads_1;
	int m_start_1;
	int m_end_1;

	int m_file_2;
	vector<Read*> m_reads_2;
	int m_start_2;
	int m_end_2;


	void loadFiles();
public:
	// Public interface:
	//
	SequenceData(vector<string>*files,ostream*logger);
	int size();
	uint64_t bases();
	Read*at(int i);
	int getFirst(string a);
	int getLast(string a);
	bool hasFile(string a);
};

#endif
