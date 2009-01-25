/*
	dna: De Novo Assembler
    Copyright (C) 2008, 2009 Sébastien Boisvert
	$Id: splitter_main.cpp 11 2009-01-23 00:57:32Z boiseb01 $

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


#include"DeBruijnAssembler.h"
#include"LightVertex.h"

LightVertex::LightVertex(){
	m_children=0;
	m_color=-1;
}

vector<VERTEX_TYPE> LightVertex::getChildren(VERTEX_TYPE prefix,int wordSize){
	vector<VERTEX_TYPE> output;
	string self=DeBruijnAssembler::idToWord(prefix,wordSize);
	string thePart=self.substr(1);
	for(int i=0;i<4;i++){
		//     7 6 5 4 3 2 1 0
		
		int v=((m_children<<(7-i))>>7);
		if(v==1){
			char a='A';
			if(i==0){
				a='A';
			}else if(i==1){
				a='T';
			}else if(i==2){
				a='C';
			}else if(i==3){
				a='G';
			}
			string b=thePart+a;
			output.push_back(DeBruijnAssembler::wordId(b.c_str()));
		}
	}
	return output;
}

void LightVertex::addChild(VERTEX_TYPE suffix,int wordSize){
	string v=DeBruijnAssembler::idToWord(suffix,wordSize);
	char a=v[wordSize-1];
	uint8_t toAdd=1;
	if(a=='A'){
		toAdd<<0;
	}else if(a=='T'){
		toAdd<<1;
	}else if(a=='C'){
		toAdd<<2;
	}else if(a=='G'){
		toAdd<<3;
	}
	m_children|=toAdd;
}


void LightVertex::setColor(uint32_t c){
	m_color=c;
}

uint32_t LightVertex::getColor(){
	return m_color;
}
