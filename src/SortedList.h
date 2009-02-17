/*
	dna: De Novo Assembler
    Copyright (C) 2008, 2009 Sébastien Boisvert
	$Id: De_Bruijn_De_Novo_Assembler_main.cpp 116 2009-02-16 21:19:41Z boiseb01 $

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

class SortedList{
	vector<VERTEX_TYPE> m_list;
public:
	void add(VERTEX_TYPE a);
	vector<VERTEX_TYPE> elementsWithALeastCCoverage(int c);
	map<int,int>getDistributionOfCoverage();
	void sort();
	void clear();
};

#endif
