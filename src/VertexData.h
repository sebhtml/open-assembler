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

#ifndef _VertexData
#define _VertexData

#include<map>
#include<vector>
#include"Read.h"
#include"AnnotationElement.h"
#include<stdint.h>
using namespace std;

class GraphData;

class VertexData{
	static vector<AnnotationElement> m_empty_vector;
	// bit    nucleotide
	// 0      A
	// 1      T
	// 2      C
	// 3      G
	// 4	Unused
	// 5	Unused
	// 6	Unused
	// 7	Unused
	
	bool m_isEliminated;
	bool m_assembled;
	uint8_t m_parents;
	uint32_t m_color;

	uint8_t m_children;
	map<char,vector<AnnotationElement> >*m_annotations;
	
public:
	vector<VERTEX_TYPE> getChildren(VERTEX_TYPE prefix);
	vector<VERTEX_TYPE> getParents(VERTEX_TYPE prefix,GraphData*m_data);
	vector<AnnotationElement>*getAnnotations(VERTEX_TYPE suffix);
	void addAnnotation(VERTEX_TYPE suffix,uint32_t read,POSITION_TYPE position,uint8_t strand);
	void addParent(VERTEX_TYPE a);
	void addChild(VERTEX_TYPE a);
	bool hasChild(VERTEX_TYPE suffix);
	VertexData();
	~VertexData();
	bool IsEliminated();
	void eliminateNow();
	uint32_t getColor();
	void setColor(uint32_t c);
	void assemble();
	bool IsAssembled();
	bool NotTrivial(VERTEX_TYPE a);
};

#endif
