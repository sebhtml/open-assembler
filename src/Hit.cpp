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



#include<Hit.h>
#include<iostream>
using namespace std;

char Hit::getReadStrand(){
	return m_readStrand;
}

int Hit::getReadNumber(){
	return m_readNumber;
}

int Hit::getReadPosition(){
	return m_readPosition;
}

Hit::Hit(){
}

int Hit::getContigPosition(){
	return m_contigPosition;
}

int Hit::getContigNumber(){
	return m_contigNumber;
}

char Hit::getContigStrand(){
	return m_contigStrand;
}

Hit::Hit(int contigNumber,int contigPosition,char contigStrand,
		int readNumber,int readPosition,char readStrand){
	//cout<<"Adding a Hit"<<endl;
	m_contigNumber=contigNumber;
	m_contigPosition=contigPosition;
	m_contigStrand=contigStrand;
	m_readNumber=readNumber;
	m_readPosition=readPosition;
	m_readStrand=readStrand;
}

void Hit::show(){
	cout<<m_contigNumber<<" "<<m_contigStrand<<" "<<m_contigPosition<<" ~ "<<m_readNumber<<" "<<m_readStrand<<" "<<m_readPosition<<endl;
}

void Hit::setL_or_R(char a){
	m_L_or_R=a;
}

char Hit::side(){
	return m_L_or_R;
}

