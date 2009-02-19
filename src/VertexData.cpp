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
#include"GraphData.h"
#include<iostream>
#include"DeBruijnAssembler.h"
#include<stdint.h>
#include<string>
#include<stdlib.h>

using namespace std;

vector<AnnotationElement> VertexData::m_empty_vector;

VertexData::VertexData(){
	m_isEliminated=false;
	m_assembled=false;
	m_parents=0;
	m_children=0;
	m_color=-1;
	m_annotations=NULL;
}


void VertexData::addAnnotation(VERTEX_TYPE suffix,uint32_t read,POSITION_TYPE position,uint8_t strand){
	if(m_annotations==NULL){
		m_annotations=new map<char,vector<AnnotationElement> >;
	}
	string a=DeBruijnAssembler::idToWord(suffix,DeBruijnAssembler::m_WordSize);
	char symbol=a[a.length()-1];
	AnnotationElement element;
	element.readId=read;
	element.readPosition=position;
	element.readStrand=strand;
	(*m_annotations)[symbol].push_back(element);
	//cout<<(*m_annotations)[symbol].size()<<" annotations"<<endl;
}

void VertexData::addChild(VERTEX_TYPE child){
	//string a=DeBruijnAssembler::idToWord(child,DeBruijnAssembler::m_WordSize);
	//char symbol=a[DeBruijnAssembler::m_WordSize-1];
	//cout<<symbol<<endl;
	//cout<<hex<<child<<endl;
	child=child<<(62-DeBruijnAssembler::m_WordSize-DeBruijnAssembler::m_WordSize);
	child=child>>62;
	//cout<<hex<<child<<endl;
	//cout<<dec;
	int position=child;
	uint8_t toAdd=1;
	toAdd=toAdd<<position;
	m_children=m_children|toAdd;
	//cout<<"Children "<<(int)m_children<<endl;
}



void VertexData::addParent(VERTEX_TYPE parent){
	//string a=DeBruijnAssembler::idToWord(parent,DeBruijnAssembler::m_WordSize);
	//char symbol=a[0];
	//cout<<"Parent "<<parent<<endl;
	//cout<<hex<<parent<<endl;
	parent=parent<<60;
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
vector<AnnotationElement>*VertexData::getAnnotations(VERTEX_TYPE suffix){
	string a=DeBruijnAssembler::idToWord(suffix,DeBruijnAssembler::m_WordSize);
	char symbol=a[DeBruijnAssembler::m_WordSize-1];
	if(m_annotations==NULL){
		return &m_empty_vector;
	}
	return &((*m_annotations)[symbol]);
}





vector<VERTEX_TYPE> VertexData::getChildren(VERTEX_TYPE prefix){
	string a=DeBruijnAssembler::idToWord(prefix,DeBruijnAssembler::m_WordSize);
	//cout<<"Prefix "<<a<<endl;
	vector<VERTEX_TYPE> output;
	//cout<<a<<endl;
	//cout<<(int)m_parents<<endl;
	//cout<<"parents " <<m_parents<<endl;
	for(int i=0;i<4;i++){
		uint8_t toCheck=m_children;
		//00001000
		//10000000
		toCheck=(toCheck<<(7-i));
		toCheck=toCheck>>7;
		if(toCheck==1){
			char symbol='A';
			if(i==0){
				symbol='A';
			}else if(i==1){
				symbol='T';
			}else if(i==2){
				symbol='C';
			}else if(i==3){
				symbol='G';
			}
			string sequence=a.substr(1,DeBruijnAssembler::m_WordSize-1)+symbol;
			//cout<<sequence<<endl;
			VERTEX_TYPE dataNode=DeBruijnAssembler::wordId(sequence.c_str());
			output.push_back(dataNode);
			//cout<<sequence<<endl;
		}
	}
	return output;
}

vector<VERTEX_TYPE> VertexData::getParents(VERTEX_TYPE prefix,GraphData*m_data){
	string a=DeBruijnAssembler::idToWord(prefix,DeBruijnAssembler::m_WordSize);
	vector<VERTEX_TYPE> output;
	//cout<<a<<endl;
	//cout<<(int)m_parents<<endl;
	//cout<<"parents " <<m_parents<<endl;
	for(int i=0;i<4;i++){
		uint8_t toCheck=m_parents;
		//00001000
		//10000000
		toCheck=(toCheck<<(7-i));
		toCheck=toCheck>>7;
		if(toCheck==1){
			char symbol='A';
			if(i==0){
				symbol='A';
			}else if(i==1){
				symbol='T';
			}else if(i==2){
				symbol='C';
			}else if(i==3){
				symbol='G';
			}
			string sequence=symbol+a.substr(0,DeBruijnAssembler::m_WordSize-1);
			//cout<<sequence<<endl;
			VERTEX_TYPE dataNode=DeBruijnAssembler::wordId(sequence.c_str());
			output.push_back(dataNode);
		}
	}
	vector<VERTEX_TYPE> eliminated;
	for(vector<VERTEX_TYPE>::iterator i=output.begin();i!=output.end();i++){
		if(m_data!=NULL&&m_data->get(*i)->IsEliminated())
			continue;
		eliminated.push_back(*i);
	}
	return eliminated;
}




VertexData::~VertexData(){
	if(m_annotations!=NULL){
		delete m_annotations;
		m_annotations=NULL;
	}
}

void VertexData::eliminateNow(){
	m_isEliminated=true;
}

bool VertexData::IsEliminated(){
	return m_isEliminated;
}

void VertexData::assemble(){
	m_assembled=true;
}

bool VertexData::IsAssembled(){
	return m_assembled;
}


uint32_t VertexData::getColor(){
	return m_color;
}

void VertexData::setColor(uint32_t c){
	m_color=c;
}

bool VertexData::NotTrivial(VERTEX_TYPE a){
	return getParents(a,NULL).size()>1||getChildren(a).size()>1;
}
