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




#ifndef _common_functions
#define _common_functions

/*
 *  unsigned long int
 *			for 64-bit machines, allows wordSize<=31
 *  unsigned int
 *  			for 32-bit machines, allows wordSize<=15
 */

//#define VERTEX_TYPE unsigned long int
//#define VERTEX_TYPE uint32_t
#define VERTEX_TYPE uint64_t

#include<string>
using namespace std;

string reverseComplement(string a);
VERTEX_TYPE wordId(const char*a);
string idToWord(VERTEX_TYPE i,int wordSize);
char complement(char a);
void CommonHeader(ostream*out);

bool isValidDNA(const char*x);

char getLastSymbol(VERTEX_TYPE i,int w);

void coutBIN(uint64_t a);

#endif
