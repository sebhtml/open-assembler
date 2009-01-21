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


void VertexData::addRead(VERTEX_TYPE suffix,int read){
	string a=DeBruijnAssembler::idToWord(suffix,DeBruijnAssembler::m_WordSize);
	char symbol=a[DeBruijnAssembler::m_WordSize-1];
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
	//cout<<"Adding read "<<symbol<<endl;
	m_reads[position].push_back(read);
}

void VertexData::addParent(VERTEX_TYPE parent){
	string a=DeBruijnAssembler::idToWord(parent,DeBruijnAssembler::m_WordSize);
	//cout<<a<<endl;
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

vector<int> VertexData::getReads(VERTEX_TYPE suffix){
	string a=DeBruijnAssembler::idToWord(suffix,DeBruijnAssembler::m_WordSize);
	char symbol=a[DeBruijnAssembler::m_WordSize-1];
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
	if(m_reads.count(position)==0){
		cout<<"Error "<<a<<" is not a child"<<endl;
		exit(0);
	}
	return m_reads[position];
}

bool VertexData::hasChild(VERTEX_TYPE suffix){
	string a=DeBruijnAssembler::idToWord(suffix,DeBruijnAssembler::m_WordSize);
	char symbol=a[DeBruijnAssembler::m_WordSize-1];
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
	return m_reads.count(position)>0;
}

vector<VERTEX_TYPE> VertexData::getChildren(VERTEX_TYPE prefix){
	string a=DeBruijnAssembler::idToWord(prefix,DeBruijnAssembler::m_WordSize);
	if(a=="ATTGTAACAAATTCTCCTGCCTCTG"){
		cout<<"[getChildren] "<<m_reads.size()<<endl;
	}
	vector<VERTEX_TYPE> output;
	//cout<<a<<endl;
	//cout<<m_reads.size()<<" children"<<endl;
	for(int i=0;i<4;i++){
		if(m_reads.count(i)>0){
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
			string sequence=a.substr(1)+symbol;
			output.push_back(DeBruijnAssembler::wordId(sequence.c_str()));
		}
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
