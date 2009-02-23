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

#ifndef _Loader
#define _Loader

#include<vector>
#include"Read.h"
#include<fstream>
#include<string>
using namespace std;

class Loader{
	int m_total;
	int m_bases;
	void add(vector<Read*>*reads,string*id,ostringstream*sequence,ostringstream*quality);
public:
	Loader();
	void load(string file,vector<Read*>*reads);
	int getBases();
};

#endif
