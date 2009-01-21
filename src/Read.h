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

#ifndef _Read
#define _Read

/*
 *  unsigned long int
 *			for 64-bit machines, allows wordSize<=31
 *  unsigned int
 *  			for 32-bit machines, allows wordSize<=15
 */

#define VERTEX_TYPE unsigned long int

#include<vector>
using namespace std;

class Read{
	char*m_id;
	char*m_sequence;
	char*m_quality;
	vector<VERTEX_TYPE> m_highQualityMers;
public:
	Read(const char*id,const char*sequence,const char*quality);
	~Read();
	char*getQual();
	char*getId();
	char*getSeq();
	vector<VERTEX_TYPE>*getHighQualityMers(int wordSize);
};

#endif
