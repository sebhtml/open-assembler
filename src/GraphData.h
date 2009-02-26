/*
	dna: De Novo Assembler
    Copyright (C) 2008, 2009 Sébastien Boisvert
	$Id: DeBruijnAssembler.cpp 116 2009-02-16 21:19:41Z boiseb01 $

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



#ifndef _GraphData
#define _GraphData

#include"VertexData.h"
#include<vector>
using namespace std;


class GraphData{
	int m_size;
	vector<VERTEX_TYPE> m_nodes;
	VertexData*m_node_data;
public:
	GraphData();
	~GraphData();
	VertexData*get(VERTEX_TYPE a);
	vector<VERTEX_TYPE>*getNodes();
	VertexData*getNodeData();
	int size();
	void add(VERTEX_TYPE a);
	bool hasNode(VERTEX_TYPE a);
	void makeMemory();
};

#endif
