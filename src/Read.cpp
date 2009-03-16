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

#include"common_functions.h"
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
	free(m_id);
	m_id=NULL;
	free(m_sequence);
	m_sequence=NULL;
}

char*Read::getId(){
	return m_id;
}

char*Read::getSeq(){
	return m_sequence;
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

int Read::length(){
	return strlen(m_sequence);
}

/*                      
 *           -----------------------------------
 *           -----------------------------------
 *                     p p-1 p-2               0
 */
uint64_t Read::Vertex(int pos,int w,char strand){
	uint64_t key=0;
	if(strand=='F'){
		for(int i=0;i<w;i++){
			char a=m_sequence[pos+i];
			uint64_t mask=0;
			if(a=='A'){
				mask=0;
			}else if(a=='T'){
				mask=1;
			}else if(a=='C'){
				mask=2;
			}else if(a=='G'){
				mask=3;
			}
			key=key|(mask<<(2*i));
		}
	}else{
		for(int i=0;i<w;i++){
			char a=m_sequence[strlen(m_sequence)-1-i-pos];
			uint64_t mask=0;
			if(a=='A'){
				mask=1;
			}else if(a=='T'){
				mask=0;
			}else if(a=='C'){
				mask=3;
			}else if(a=='G'){
				mask=2;
			}
			key=key|(mask<<(2*i));
		}
	}

	return key;
}
