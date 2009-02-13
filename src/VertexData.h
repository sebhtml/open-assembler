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
#include "CustomMap.hpp"
#include"Read.h"
#include"AnnotationElement.h"
#include<stdint.h>
using namespace std;


class VertexData{
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
	uint8_t m_parents;
	vector<int> m_distances;
	vector<VERTEX_TYPE> m_kmer_apart;

	vector<AnnotationElement> A_elements;
	vector<AnnotationElement> T_elements;
	vector<AnnotationElement> C_elements;
	vector<AnnotationElement> G_elements;
	
public:
	vector<VERTEX_TYPE> getChildren(VERTEX_TYPE prefix);
	vector<VERTEX_TYPE> getParents(VERTEX_TYPE prefix,CustomMap<VertexData>*m_data);
	vector<AnnotationElement>*getAnnotations(VERTEX_TYPE suffix);
	void addParent(VERTEX_TYPE parent);
	void addAnnotation(VERTEX_TYPE suffix,uint32_t read,POSITION_TYPE position,uint8_t strand);
	bool hasChild(VERTEX_TYPE suffix);
	void addPaired(VERTEX_TYPE kmer_apart,int distance);
	VertexData();
	~VertexData();
	bool IsEliminated();
	void eliminateNow();
};

#endif
