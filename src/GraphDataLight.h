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



#ifndef _GraphDataLight
#define _GraphDataLight

#include"LightVertex.h"
#include<vector>
using namespace std;


class GraphDataLight{
	vector<VERTEX_TYPE> m_nodes;
	vector<LightVertex> m_node_data;
public:
	LightVertex*get(VERTEX_TYPE a);
	vector<VERTEX_TYPE>*getNodes();
	vector<LightVertex>*getNodeData();
	int size();
	void add(VERTEX_TYPE a);
	bool find(VERTEX_TYPE a);
};

#endif
