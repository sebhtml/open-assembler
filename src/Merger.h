/*
	dna: De Novo Assembler
    Copyright (C) 2008, 2009 Sébastien Boisvert
	$Id: De_Bruijn_De_Novo_Assembler_main.cpp 122 2009-02-17 13:24:50Z boiseb01 $

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



#ifndef _Merger
#define _Merger

#include<vector>
#include"Read.h"
using namespace std;

class Merger{
	bool m_DEBUG;
	int m_wordSize;
	vector<Read*> reverseMerge(vector<Read*> contigSequences);
	vector<Read*> forwardOverlapTailToHead(vector<Read*> contigSequences);
	vector<Read*> reverseOverlapTailToTail(vector<Read*> contigSequences);
	vector<Read*> reverseOverlapHeadToHead(vector<Read*> contigSequences);
	vector<Read*> mergeForward(vector<Read*> contigSequences);
public:
	Merger();
	vector<Read*>mergeContigs(vector<Read*>contigs);
};

#endif