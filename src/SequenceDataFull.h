/*
	dna: De Novo Assembler
    Copyright (C) 2008, 2009 Sébastien Boisvert
	$Id: SequenceData.h 11 2009-01-23 00:57:32Z boiseb01 $

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

#ifndef _SequenceDataFull
#define _SequenceDataFull

#include<string>
#include<map>
#include<vector>
#include"Read.h"
#include<stdint.h>
using namespace std;

class SequenceDataFull{
	map<string,int > m_file_first;
	map<string,int > m_file_last;
	vector<Read*> m_reads_vector;
	vector<string> m_inputFiles;
	ostream*m_cout;
	uint64_t m_bases;
	int m_reads;


	void loadFiles();
public:
	// Public interface:
	//
	SequenceDataFull(vector<string>*files,ostream*logger);
	int size();
	uint64_t bases();
	Read*at(int i);
	int getFirst(string a);
	int getLast(string a);
	bool hasFile(string a);
};

#endif
