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



#include"HitPair.h"
#include<iostream>
using namespace std;

Hit*HitPair::getLeft(){
	return &m_left;
}

Hit*HitPair::getRight(){
	return &m_right;
}

HitPair::HitPair(Hit*left,Hit*right){
	m_left=*left;
	m_right=*right;
}

HitPair::HitPair(){
}

bool HitPair::valid(){
	if(m_left.getContigPosition()==0&&m_right.getContigPosition()==0)
		return false;
	if(//m_left.getContigStrand()==m_right.getContigStrand()&&
		((m_left.getContigPosition()==0&&m_right.getContigPosition()==0)||
		(m_left.getContigPosition()!=0&&m_right.getContigPosition()!=0)))
		return false;

	return true;
}

void HitPair::show(){
	cout<<"Hitpair"<<endl;
	m_left.show();
	m_right.show();
}