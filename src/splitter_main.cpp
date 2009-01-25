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

#include"CustomMap.hpp"
#include"LightVertex.h"
#include<vector>
#include"DeBruijnAssembler.h"
#include"Read.h"
#include"Loader.h"
#include<map>
#include<string>
#include<stdlib.h>
#include<iostream>
using namespace std;

void applyColor(VERTEX_TYPE v,CustomMap<LightVertex>*graph,int color,int wordSize){
	cout<<"Coloring "<<v <<" "<<color<<endl;
	graph->get(v).setColor(color);
	vector<VERTEX_TYPE>parents=graph->get(v).getParents(v,wordSize);
	vector<VERTEX_TYPE>children=graph->get(v).getChildren(v,wordSize);
	//cout<<parents.size()<<" parents"<<endl;
	//cout<<children.size()<<" children"<<endl;
	for(int k=0;k<parents.size();k++){
		if(graph->get(parents[k]).getColor()==color)
			continue;
		applyColor(parents[k],graph,color,wordSize);
	}
	
	for(int k=0;k<children.size();k++){
		if(graph->get(children[k]).getColor()==color)
			continue;
		applyColor(children[k],graph,color,wordSize);
	}
}

int main(int argc,char*argv[]){
	cout<<"usage"<<endl;
	cout<<"dna_DeBruijnSplitter [-outputDirectory parts] [-wordSize 21] [-buckets 100000000] [-minimumCoverage auto] <sequence files>"<<endl;
	string outputDirectory="parts";
	int wordSize=21;
	int buckets=100000000;
	string minimumCoverage="auto";
	vector<string> inputFiles;

	// collect arguments
	for(int i=1;i<argc;i++){
		string option=argv[i];
		if(option=="-outputDirectory"){
			i++;
			outputDirectory=argv[i];
		}else if(option=="-minimumCoverage"){
			i++;
			minimumCoverage=argv[i];
		}else if(option=="-buckets"){
			i++;
			buckets=atoi(argv[i]);
		}else if(option=="-wordSize"){
			i++;
			if(wordSize>31){
				cout<<"Wordsize cannot be greater than 31"<<endl;
			}else{
				wordSize=atoi(argv[i]);
			}
		}else{
			inputFiles.push_back(argv[i]);
		}
	}

	CustomMap<int> wordCount(buckets);
	
	for(int i=0;i<inputFiles.size();i++){
		vector<Read*> reads;
		Loader loader(&cout);
		loader.load(inputFiles[i],&reads);
		for(int j=0;j<reads.size();j++){
			vector<VERTEX_TYPE> mers=reads.at(j)->getHighQualityMers(wordSize);
			for(int k=0;k<mers.size();k++){
				if(!wordCount.find(mers[k]))
					wordCount.add(mers[k],0);
				wordCount.set(mers[k],wordCount.get(mers[k])+1);
			}
		}
	}


	// depletion curve analysis
	//
	//
	

	(cout)<<"Quality analysis (Depletion curve), k+1 = "<<wordSize+1<<endl;
	(cout)<<"MinimumCoverage NumberOfSolidMers NormalizedNumberOfSolidMers Change"<<endl;
	map<int,int> solidCounts;
	for(int minimumCoverage=1;minimumCoverage<=32;minimumCoverage++){
		solidCounts[minimumCoverage]=0;
	}
	for(CustomMap<int>::iterator i=wordCount.begin();i!=wordCount.end();i++){
		for(int minimumCoverage=1;minimumCoverage<=32;minimumCoverage++){
			if(i.second()>=minimumCoverage){
				solidCounts[minimumCoverage]++;
			}
		}
	}

	map<int,double> solidRatio;
	
	for(int minimumCoverage=1;minimumCoverage<=32;minimumCoverage++){
		solidRatio[minimumCoverage]=solidCounts[minimumCoverage]/(0.0+solidCounts[1]);
	}
	map<int,double> solidDiff;
	for(int minimumCoverage=2;minimumCoverage<=32;minimumCoverage++){
		solidDiff[minimumCoverage]=solidRatio[minimumCoverage-1]-solidRatio[minimumCoverage];
	}

	int m_minimumCoverage=5;
	double minimum=1;
	for(int minimumCoverage=1;minimumCoverage<=32;minimumCoverage++){
		cout<<minimumCoverage<<" "<<solidCounts[minimumCoverage]<<" "<<solidRatio[minimumCoverage];
		if(minimumCoverage!=1)
			cout<<" "<<solidDiff[minimumCoverage];
		cout<<endl;
	}
	cout<<endl;
	double cutOff=0.015;
	cout<<"Cutoff: "<<cutOff<<endl;
	for(int minimumCoverage=2;minimumCoverage<=32;minimumCoverage++){
		if(solidDiff[minimumCoverage]<minimum){
			m_minimumCoverage=minimumCoverage;
			minimum=solidDiff[minimumCoverage];
		}
		if(solidDiff[minimumCoverage]>minimum)
			break;
		if(solidDiff[minimumCoverage]<cutOff){
			m_minimumCoverage=minimumCoverage;
			cout<<"Best Coverage <- "<<m_minimumCoverage<<endl;
			break;
		}
	}

	if(minimumCoverage=="auto"){
		cout<<"Using depletion curve (-minimumCoverage auto)"<<endl;
		cout<<"m_minimumCoverage <- "<<m_minimumCoverage<<endl;
	}else{
		cout<<"Not using the depletion curve"<<endl;
		cout<<"m_minimumCoverage <- "<<atoi(minimumCoverage.c_str())<<endl;
		m_minimumCoverage=atoi(minimumCoverage.c_str());
	}

	CustomMap<LightVertex> graphWithoutData(buckets);
	for(CustomMap<int>::iterator i=wordCount.begin();i!=wordCount.end();i++){
		if(i.second()>=m_minimumCoverage){
			VERTEX_TYPE mer=i.first();
			string merString=DeBruijnAssembler::idToWord(mer,wordSize+1);
			string prefix=merString.substr(0,wordSize);
			string suffix=merString.substr(1,wordSize);
			VERTEX_TYPE prefixInteger=DeBruijnAssembler::wordId(prefix.c_str());
			VERTEX_TYPE suffixInteger=DeBruijnAssembler::wordId(suffix.c_str());
			if(!graphWithoutData.find(prefixInteger)){
				LightVertex vertex;
				graphWithoutData.add(prefixInteger,vertex);
			}
			if(!graphWithoutData.find(suffixInteger)){
				LightVertex vertex;
				graphWithoutData.add(suffixInteger,vertex);
			}
			graphWithoutData.get(prefixInteger).addChild(suffixInteger,wordSize);
			graphWithoutData.get(suffixInteger).addParent(prefixInteger,wordSize);
		}
	}
	
	// find connected components, but how?
	int color=1;
	for(CustomMap<LightVertex>::iterator i=graphWithoutData.begin();i!=graphWithoutData.end();i++){
		if(graphWithoutData.get(i.first()).getColor()!=-1)
			continue;
		applyColor(i.first(),&graphWithoutData,color,wordSize);
		color++;
	}
	cout<<color<<" colors"<<endl;
	return 0;
}

