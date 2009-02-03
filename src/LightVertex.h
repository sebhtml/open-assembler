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

#ifndef _LightVertex
#define _LightVertex

#include"Read.h"
#include<vector>
using namespace std;

class LightVertex{
	uint8_t m_children;
	uint8_t m_parents;
	uint32_t m_color;
public:
	LightVertex();
	vector<VERTEX_TYPE> getChildren(VERTEX_TYPE prefix,int wordSize);
	vector<VERTEX_TYPE> getParents(VERTEX_TYPE prefix,int wordSize);
	void addChild(VERTEX_TYPE suffix,int wordSize);
	void addParent(VERTEX_TYPE v,int wordSize);
	void setColor(uint32_t c);
	uint32_t getColor();
};

#endif

