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

#include"GraphData.h"
#include<iostream>
#include<stdlib.h>
#include"BinarySearch.h"
using namespace std;

VertexData*GraphData::get(VERTEX_TYPE a){
	int index=BinarySearch(&m_nodes,a,m_size);
	if(index==-1){
		cout<<"Error, not found (should not happen...)"<<endl;
		cout<<idToWord(a,21)<<endl;
		cout<<a<<endl;
		exit(0);
	}
	return &(m_node_data[index]);
}

bool GraphData::hasNode(VERTEX_TYPE a){
	return BinarySearch(&m_nodes,a,m_size)!=-1;
}

vector<VERTEX_TYPE>*GraphData::getNodes(){
	return &m_nodes;
}

VertexData*GraphData::getNodeData(){
	return m_node_data;
}

int GraphData::size(){
	return m_size;
}

void GraphData::add(VERTEX_TYPE a){
	m_nodes.push_back(a);
}

void GraphData::makeMemory(){
	m_size=m_nodes.size();
	m_node_data=new VertexData[m_size];
}

GraphData::GraphData(){
	m_node_data=NULL;
}
GraphData::~GraphData(){
	if(m_node_data!=NULL){
		delete [] m_node_data;
		m_node_data=NULL;
	}
}
