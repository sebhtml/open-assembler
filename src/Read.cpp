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


