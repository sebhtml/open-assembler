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

#ifndef _Read
#define _Read

#include"common_functions.h"
#include<string>
#include<vector>
using namespace std;

class Read{
	char*m_id;
	char*m_sequence;
	int m_startForward;
	int m_startReverse;
public:
	Read(const char*id,const char*sequence);
	char*getId();
	~Read();
	char*getSeq();
	int getStartForward();
	int getStartReverse();
	void setStartForward(int i);
	int length();
	void setStartReverse(int i);
};

#endif
