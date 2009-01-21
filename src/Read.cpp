/*
	454 De Novo Assembler
    Copyright (C) 2008 Sébastien Boisvert

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


Read::Read(const char*id,const char*sequence,const char*quality){
	m_id=(char*)malloc(strlen(id)+1);
	m_sequence=(char*)malloc(strlen(sequence)+1);
	m_quality=(char*)malloc(strlen(quality)+1);
	strcpy(m_id,id);
	strcpy(m_sequence,sequence);
	strcpy(m_quality,quality);
	//cout<<m_id<<endl;
}

Read::~Read(){
	free(m_id);
	free(m_sequence);
	free(m_quality);
	m_id=NULL;
	m_sequence=NULL;
	m_quality=NULL;
}

char*Read::getQual(){
	return m_quality;
}

char*Read::getId(){
	return m_id;
}

char*Read::getSeq(){
	return m_sequence;
}

// return (k+1)-mers from read
// check that the k+1 nucleotides at position i are good quality and that they are not the N
vector<VERTEX_TYPE>*Read::getHighQualityMers(int wordSize){
	if(m_highQualityMers.size()>0)
		return &m_highQualityMers;

	string sequence=m_sequence;

	for(int j=0;j<(int)sequence.length()-(wordSize+1);j++){
		string wordFoward=sequence.substr(j,wordSize+1);
		string wordReverse;
		bool ok=true;
		for(int p=0;p<wordSize+1;p++){
			if(wordFoward[p]=='N'){
				ok=false;
				break;
			}
		}
		if(!ok)
			continue;
		wordReverse=DeBruijnAssembler::reverseComplement(wordFoward);
		m_highQualityMers.push_back(DeBruijnAssembler::wordId(wordFoward.c_str()));
		m_highQualityMers.push_back(DeBruijnAssembler::wordId(wordReverse.c_str()));
	}

	return &m_highQualityMers;
}

