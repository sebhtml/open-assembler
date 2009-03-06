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




#include"SortedList.h"
#include<iostream>
#include<algorithm>
using namespace std;


SortedList::SortedList(){
	m_list=new vector<SortableElement>;
}

void SortedList::add(VERTEX_TYPE a){
	SortableElement element(a,1);
	m_list->push_back(element);
}


vector<VERTEX_TYPE> SortedList::elementsWithALeastCCoverage(int c){
	vector<VERTEX_TYPE> output;
	vector<SortableElement>::iterator i=m_list->begin();
	while(i!=m_list->end()){
		int currentCount=(*i).getCount();
		VERTEX_TYPE currentValue=(*i).getKMer();
		i++;
		while(i!=m_list->end()&&(*i).getKMer()==currentValue){
			currentCount+=(*i).getCount();
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
	vector<SortableElement>::iterator i=m_list->begin();
	while(i!=m_list->end()){
		VERTEX_TYPE currentValue=(*i).getKMer();
		int currentCount=(*i).getCount();
		i++;
		while(i!=m_list->end()&&(*i).getKMer()==currentValue){
			currentCount+=(*i).getCount();
			i++;
		}
		m_coverageDistribution[currentCount]++;
	}

	return m_coverageDistribution;
}

void SortedList::sort(){
	cout<<"Sorting "<<m_list->size()<<" elements"<<endl;
	SortableElement e(0,0);
	std::sort(m_list->begin(),m_list->end(),e);
	vector<SortableElement>*newList=new vector<SortableElement>;
	
	vector<SortableElement>::iterator i=m_list->begin();
	while(i!=m_list->end()){
		int currentCount=1;
		VERTEX_TYPE currentValue=(*i).getKMer();
		i++;
		while(i!=m_list->end()&&(*i).getKMer()==currentValue){
			currentCount+=(*i).getCount();
			i++;
		}
		SortableElement newElement(currentValue,currentCount);
		newList->push_back(newElement);
	}
	cout<<"Transformation: "<<m_list->size()<<" -> "<<newList->size()<<endl;
	delete m_list;
	m_list=newList;
}

void SortedList::clear(){
	m_list->clear();
	delete m_list;
	m_list=NULL;
}


bool SortableElement::operator()(const SortableElement&a,const SortableElement&b){
	return a.m_kmer<b.m_kmer;
}

uint64_t SortableElement::getKMer(){
	return m_kmer;
}

uint16_t SortableElement::getCount(){
	return m_count;
}

SortableElement::SortableElement(uint64_t kmer,uint16_t count){
	m_kmer=kmer;
	m_count=count;
}
