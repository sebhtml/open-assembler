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
	
	//bool m_isEliminated;
	//map<VERTEX_TYPE,vector<int> > m_positionInContig;
	//uint32_t m_color;
	vector<AnnotationElement>*m_annotations;
	bool m_assembled;
	bool m_trivial;
	uint8_t m_children;
	uint8_t m_parents;
	
public:
	void set_topology_1_1();
	bool Is_1_1();
	vector<VERTEX_TYPE> getChildren(VERTEX_TYPE prefix,int w);
	vector<VERTEX_TYPE> getParents(VERTEX_TYPE prefix,int w);
	vector<AnnotationElement>*getAnnotations();
	void addAnnotation(uint32_t read,POSITION_TYPE position,uint8_t strand);
	void addParent(VERTEX_TYPE a,int w);
	void addChild(VERTEX_TYPE a,int w);
	bool hasChild(VERTEX_TYPE suffix);
	VertexData();
	~VertexData();
	bool IsAssembled();
	void assemble();
	bool hasManyChildren(VERTEX_TYPE a,int w);
};

#endif
