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

#include "VertexData.h"
#include<iostream>
#include<stdint.h>
#include<string>
#include<stdlib.h>

using namespace std;

vector<AnnotationElement> VertexData::m_empty_vector;

VertexData::VertexData(){
	m_parents=0;
	m_children=0;
	m_annotations=NULL;
	m_assembled=false;
	m_trivial=false;
	m_deleted=false;
	m_count=0;
}

void VertexData::addToCount(uint16_t a){
	uint16_t newValue=m_count+a;
	if(m_count<newValue)
		m_count+=a;
}

void VertexData::assemble(){
	m_assembled=true;
}

void VertexData::addAnnotation(uint32_t read,POSITION_TYPE position,uint8_t strand){
	if(m_annotations==NULL)
		m_annotations=new vector<AnnotationElement>;
	AnnotationElement element(read,position,strand);
	m_annotations->push_back(element);
}

/*
 * 63 62 61 60 ... 1 0
 *
 */
void VertexData::addChild(VERTEX_TYPE child,int m_wordSize){
	m_children=m_children|(1<<(child<<(64-2*m_wordSize)>>62));
}



void VertexData::addParent(VERTEX_TYPE parent,int m_WordSize){
	//string a=DeBruijnAssembler::idToWord(parent,m_WordSize);
	//char symbol=a[0];
	//cout<<"Parent "<<parent<<endl;
	//cout<<hex<<parent<<endl;
	parent=parent<<62;
	//cout<<hex<<parent<<endl;
	parent=parent>>62;
	//cout<<hex<<parent<<endl;
	int position=parent;
	//cout<<symbol<<" "<<position<<endl;
	uint8_t toAdd=1;
	toAdd=toAdd<<position;
	m_parents=m_parents|toAdd;
	//cout<<"Parents "<<(int)m_parents<<endl;
}

/*
 * TODO: here, pass the list of reads to build the annotations here...
 *
 */
vector<AnnotationElement>*VertexData::getAnnotations(){
	if(m_annotations==NULL)
		return &m_empty_vector;
	return m_annotations;
}





vector<VERTEX_TYPE> VertexData::getChildren(VERTEX_TYPE prefix,int m_wordSize){
	vector<VERTEX_TYPE> output;
	for(uint64_t i=0;i<4;i++){
		uint8_t toCheck=m_children;
		toCheck=(toCheck<<(7-i));
		toCheck=toCheck>>7;
		if(toCheck==1){
			VERTEX_TYPE dataNode=(prefix>>2)|(i<<(2*(m_wordSize-1)));
			output.push_back(dataNode);
		}
	}
	return output;
}

vector<VERTEX_TYPE> VertexData::getParents(VERTEX_TYPE prefix,int m_WordSize){
	vector<VERTEX_TYPE> output;
	for(uint64_t i=0;i<4;i++){
		uint8_t toCheck=m_parents;
		toCheck=(toCheck<<(7-i));
		toCheck=toCheck>>7;
		if(toCheck==1){
			/*
 				63 62 61 60 ... 5 4 3 2 1 0

			*/
			VERTEX_TYPE dataNode=((prefix<<(64-2*m_WordSize+2))>>(64-2*m_WordSize))|i;
			output.push_back(dataNode);
		}
	}
	return output;
}




VertexData::~VertexData(){
	if(m_annotations!=NULL){
		delete m_annotations;
		m_annotations=NULL;
	}
}

bool VertexData::IsAssembled(){
	//return m_positionInContig.size()>0;
	return m_assembled;
}

void VertexData::set_topology_1_1(){
	m_trivial=true;
}

bool VertexData::Is_1_1(){
	return m_trivial;
}

void VertexData::Delete(){
	m_deleted=true;
}

bool VertexData::Deleted(){
	return m_deleted;
}

int VertexData::getCount(){
	return m_count;
}


