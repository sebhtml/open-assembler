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
	m_parents=0;
	m_children=0;
	m_color=-1;
}


void VertexData::addAnnotation(uint32_t read,POSITION_TYPE position,uint8_t strand){
	AnnotationElement element;
	element.readId=read;
	element.readPosition=position;
	element.readStrand=strand;
	m_annotations.push_back(element);
}

void VertexData::addChild(VERTEX_TYPE child,int m_wordSize){
	//string a=DeBruijnAssembler::idToWord(child,m_wordSize);
	//char symbol=a[m_WordSize-1];
	//cout<<symbol<<endl;
	//cout<<hex<<child<<endl;
	child=child<<(62-m_wordSize-m_wordSize);
	child=child>>62;
	//cout<<hex<<child<<endl;
	//cout<<dec;
	int position=child;
	uint8_t toAdd=1;
	toAdd=toAdd<<position;
	m_children=m_children|toAdd;
	//cout<<"Children "<<(int)m_children<<endl;
}



void VertexData::addParent(VERTEX_TYPE parent,int m_WordSize){
	//string a=DeBruijnAssembler::idToWord(parent,m_WordSize);
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
vector<AnnotationElement>*VertexData::getAnnotations(){
	return &m_annotations;
}





vector<VERTEX_TYPE> VertexData::getChildren(VERTEX_TYPE prefix,int m_WordSize){
	string a=idToWord(prefix,m_WordSize);
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
			string sequence=a.substr(1,m_WordSize-1)+symbol;
			//cout<<sequence<<endl;
			VERTEX_TYPE dataNode=wordId(sequence.c_str());
			output.push_back(dataNode);
			//cout<<sequence<<endl;
		}
	}
	return output;
}

vector<VERTEX_TYPE> VertexData::getParents(VERTEX_TYPE prefix,GraphData*m_data,int m_WordSize){
	string a=idToWord(prefix,m_WordSize);
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
			string sequence=symbol+a.substr(0,m_WordSize-1);
			//cout<<sequence<<endl;
			VERTEX_TYPE dataNode=wordId(sequence.c_str());
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
}

void VertexData::eliminateNow(){
	m_isEliminated=true;
}

bool VertexData::IsEliminated(){
	return m_isEliminated;
}


bool VertexData::IsAssembled(){
	return m_positionInContig.size()>0;
}


uint32_t VertexData::getColor(){
	return m_color;
}

void VertexData::setColor(uint32_t c){
	m_color=c;
}

bool VertexData::NotTrivial(VERTEX_TYPE a,int m_wordSize){
	return getParents(a,NULL,m_wordSize).size()>1||getChildren(a,m_wordSize).size()>1;
}


void VertexData::addPositionInContig(VERTEX_TYPE a,int b){
	m_positionInContig[a].push_back(b);
}

map<VERTEX_TYPE,vector<int> >*VertexData::getPositions(){
	return &m_positionInContig;
}

void VertexData::printPositions(){
	for(map<VERTEX_TYPE,vector<int> >::iterator i=m_positionInContig.begin();i!=m_positionInContig.end();i++){
		for(vector<int>::iterator j=i->second.begin();j!=i->second.end();j++){
			cout<<i->first<<" "<<*j<<endl;
		}
	}
}
