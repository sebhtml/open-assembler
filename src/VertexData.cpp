/*
	dna: De Novo Assembler
    Copyright (C) 2008, 2009 Sébastien Boisvert
	$Id: VertexData.cpp 290 2009-01-15 14:23:28Z sebhtml $

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
#include"DeBruijnAssembler.h"
#include<stdint.h>
#include<string>
#include<stdlib.h>

using namespace std;

VertexData::VertexData(){
	m_parents=0;
}


void VertexData::addAnnotation(VERTEX_TYPE suffix,uint32_t read,uint16_t position,uint8_t strand){
	string a=DeBruijnAssembler::idToWord(suffix,DeBruijnAssembler::m_WordSize);
	char symbol=a[a.length()-1];
	AnnotationElement element;
	element.readId=read;
	element.readPosition=position;
	element.readStrand=strand;

	if(symbol=='A'){
		A_elements.push_back(element);
	}else if(symbol=='T'){
		T_elements.push_back(element);
	}else if(symbol=='C'){
		C_elements.push_back(element);
	}else if(symbol=='G'){
		G_elements.push_back(element);
	}
}

void VertexData::addParent(VERTEX_TYPE parent){
	string a=DeBruijnAssembler::idToWord(parent,DeBruijnAssembler::m_WordSize);
	char symbol=a[0];
	int position=0;
	if(symbol=='A'){
		position=DeBruijnAssembler::m_NUCLEOTIDE_A;
	}else if(symbol=='T'){
		position=DeBruijnAssembler::m_NUCLEOTIDE_T;
	}else if(symbol=='C'){
		position=DeBruijnAssembler::m_NUCLEOTIDE_C;
	}else if(symbol=='G'){
		position=DeBruijnAssembler::m_NUCLEOTIDE_G;
	}
	uint8_t toAdd=1;
	toAdd=toAdd<<position;
	m_parents=m_parents|toAdd;
	//cout<<(int)m_parents<<endl;
}

vector<AnnotationElement>*VertexData::getAnnotations(VERTEX_TYPE suffix){
	string a=DeBruijnAssembler::idToWord(suffix,DeBruijnAssembler::m_WordSize);
	char symbol=a[DeBruijnAssembler::m_WordSize-1];
	if(symbol=='A'){
		return &A_elements;
	}else if(symbol=='T'){
		return &T_elements;
	}else if(symbol=='C'){
		return &C_elements;
	}else if(symbol=='G'){
		return &G_elements;
	}
	cout<<"Error "<<a<<" is not a child"<<endl;
	exit(0);
	return NULL;
}

bool VertexData::hasChild(VERTEX_TYPE suffix){
	string a=DeBruijnAssembler::idToWord(suffix,DeBruijnAssembler::m_WordSize);
	char symbol=a[DeBruijnAssembler::m_WordSize-1];
	if(symbol=='A'){
		return A_elements.size()>0;
	}else if(symbol=='T'){
		return T_elements.size()>0;
	}else if(symbol=='C'){
		return C_elements.size()>0;
	}else if(symbol=='G'){
		return G_elements.size()>0;
	}
	return false;
}

vector<VERTEX_TYPE> VertexData::getChildren(VERTEX_TYPE prefix){
	string a=DeBruijnAssembler::idToWord(prefix,DeBruijnAssembler::m_WordSize);
	vector<VERTEX_TYPE> output;

	if(A_elements.size()>0){
		string sequence=a.substr(1)+"A";
		output.push_back(DeBruijnAssembler::wordId(sequence.c_str()));
	}
	if(T_elements.size()>0){
		string sequence=a.substr(1)+"T";
		output.push_back(DeBruijnAssembler::wordId(sequence.c_str()));
	}
	if(C_elements.size()>0){
		string sequence=a.substr(1)+"C";
		output.push_back(DeBruijnAssembler::wordId(sequence.c_str()));
	}
	if(G_elements.size()>0){
		string sequence=a.substr(1)+"G";
		output.push_back(DeBruijnAssembler::wordId(sequence.c_str()));
	}


	return output;
}

vector<VERTEX_TYPE> VertexData::getParents(VERTEX_TYPE prefix){
	string a=DeBruijnAssembler::idToWord(prefix,DeBruijnAssembler::m_WordSize);
	vector<VERTEX_TYPE> output;
	//cout<<a<<endl;
	//cout<<(int)m_parents<<endl;
	for(int i=0;i<4;i++){
		uint8_t toCheck=m_parents;
		toCheck=toCheck>>i;
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
			output.push_back(DeBruijnAssembler::wordId(sequence.c_str()));
		}
	}
	return output;
}


void VertexData::addPaired(VERTEX_TYPE kmer_apart,int distance){
	m_distances.push_back(distance);
	m_kmer_apart.push_back(kmer_apart);
}

VertexData::~VertexData(){
}
