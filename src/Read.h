/*
	dna: De Novo Assembler
    Copyright (C) 2008, 2009 Sébastien Boisvert
	$Id: Read.h 274 2009-01-13 23:18:48Z sebhtml $

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

/*
 *  unsigned long int
 *			for 64-bit machines, allows wordSize<=31
 *  unsigned int
 *  			for 32-bit machines, allows wordSize<=15
 */

#define VERTEX_TYPE unsigned long int

#include<string>
#include<vector>
using namespace std;

class Read{
	char*m_sequence;
public:
	Read(const char*sequence);
	~Read();
	char*getSeq();
	vector<VERTEX_TYPE>getHighQualityMers(int wordSize);
	bool isValidDNA(string*x);
};

#endif
