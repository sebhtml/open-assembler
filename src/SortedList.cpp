/*
	dna: De Novo Assembler
    Copyright (C) 2008, 2009 Sébastien Boisvert
	$Id: De_Bruijn_De_Novo_Assembler_main.cpp 116 2009-02-16 21:19:41Z boiseb01 $

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




#include"SortedList.h"
#include<iostream>
#include<algorithm>
using namespace std;


void SortedList::add(VERTEX_TYPE a){
	m_list.push_back(a);
}


vector<VERTEX_TYPE> SortedList::elementsWithALeastCCoverage(int c){
	vector<VERTEX_TYPE> output;
	vector<VERTEX_TYPE>::iterator i=m_list.begin();
	while(i!=m_list.end()){
		int currentCount=1;
		VERTEX_TYPE currentValue=*i;
		i++;
		while(i!=m_list.end()&&*i==currentValue){
			currentCount++;
			i++;
		}
		if(currentCount>=c)
			output.push_back(currentValue);
	}
	cout<<"Elements: "<<output.size()<<", c="<<c<<endl;
	return output;
}


map<int,int> SortedList::getDistributionOfCoverage(){
	map<int,int> m_coverageDistribution;
	vector<VERTEX_TYPE>::iterator i=m_list.begin();
	while(i!=m_list.end()){
		VERTEX_TYPE currentValue=*i;
		int currentCount=1;
		i++;
		while(i!=m_list.end()&&*i==currentValue){
			currentCount++;
			i++;
		}
		m_coverageDistribution[currentCount]++;
	}

	return m_coverageDistribution;
}

void SortedList::sort(){
	cout<<"Sorting "<<m_list.size()<<" elements"<<endl;
	std::sort(m_list.begin(),m_list.end());
}

void SortedList::clear(){
	m_list.clear();
}

