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

#include"Read.h"
#include<cstdlib>
#include<iostream>
#include<cstring>
#include"DeBruijnAssembler.h"

using namespace  std;

Read::Read(const char*id,const char*sequence){
	m_sequence=(char*)malloc(strlen(sequence)+1);
	m_id=(char*)malloc(strlen(id)+1);
	strcpy(m_sequence,sequence);
	strcpy(m_id,id);
	//cout<<strlen(sequence)<<endl;
	m_startForward=0;
	m_startReverse=0;
}

Read::~Read(){
	free(m_sequence);
	m_sequence=NULL;
}

char*Read::getId(){
	return m_id;
}

char*Read::getSeq(){
	return m_sequence;
}

// return (k+1)-mers from read, without Ns
vector<VERTEX_TYPE>Read::getHighQualityMers(int wordSize){
	vector<VERTEX_TYPE> highQualityMers;

	string sequence=m_sequence;

	for(int j=0;j<(int)sequence.length()-(wordSize);j++){
		string wordFoward=sequence.substr(j,wordSize+1);
		string wordReverse;
		if(!isValidDNA(wordFoward.c_str()))
			continue;
		wordReverse=DeBruijnAssembler::reverseComplement(wordFoward);
		highQualityMers.push_back(DeBruijnAssembler::wordId(wordFoward.c_str()));
		highQualityMers.push_back(DeBruijnAssembler::wordId(wordReverse.c_str()));
	}

	return highQualityMers;
}

bool Read::isValidDNA(const char*x){
	int len=strlen(x);
	for(int i=0;i<len;i++){
		char a=x[i];
		if(!(a=='A'||a=='T'||a=='C'||a=='G'))
			return false;
	}
	return true;
}

int Read::getStartForward(){
	return m_startForward;
}

int Read::getStartReverse(){
	return m_startReverse;
}

void Read::setStartForward(int i){
	m_startForward=i;
}

void Read::setStartReverse(int i){
	m_startReverse=i;
}
