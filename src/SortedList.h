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



#ifndef _SortedList
#define _SortedList

#include<list>
#include<vector>
#include"Read.h"
#include<map>
using namespace std;

class SortableElement{
	uint64_t m_kmer;
	uint16_t m_count;
public:
	uint64_t getKMer();
	uint16_t getCount();
	SortableElement(uint64_t m_kmer,uint16_t m_count);
	bool operator()(const SortableElement&a,const SortableElement&b);
};




class SortedList{
	vector<SortableElement> *m_list;
public:
	void add(VERTEX_TYPE a);
	vector<VERTEX_TYPE> elementsWithALeastCCoverage(int c);
	map<int,int>getDistributionOfCoverage();
	void sort();
	SortedList();
	void clear();
};

#endif
