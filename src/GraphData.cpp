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

#include"GraphData.h"
#include<iostream>
#include<stdlib.h>
#include"BinarySearch.h"
using namespace std;

VertexData*GraphData::get(VERTEX_TYPE a){
	int index=BinarySearch(&m_nodes,a);
	if(index==-1){
		cout<<"Error, not found (should not happen...)"<<endl;
		exit(0);
	}
	return &(m_node_data.at(index));
}

vector<VERTEX_TYPE>*GraphData::getNodes(){
	return &m_nodes;
}

vector<VertexData>*GraphData::getNodeData(){
	return &m_node_data;
}

int GraphData::size(){
	return m_nodes.size();
}

void GraphData::add(VERTEX_TYPE a){
	VertexData b;
	m_nodes.push_back(a);
	m_node_data.push_back(b);
}
