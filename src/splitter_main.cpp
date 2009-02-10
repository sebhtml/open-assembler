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
#include<sstream>
#include"LightVertex.h"
#include"CoverageDistribution.h"
#include<vector>
#include"DeBruijnAssembler.h"
#include"Read.h"
#include"Loader.h"
#include<map>
#include<string>
#include<stdlib.h>
#include<queue>
#include<iostream>
using namespace std;

void applyColor(VERTEX_TYPE v2,CustomMap<LightVertex>*graph,int color,int wordSize){
	queue<VERTEX_TYPE> stackOfVertices;
	stackOfVertices.push(v2);
	while(stackOfVertices.size()>0){
		VERTEX_TYPE v=stackOfVertices.front();
		stackOfVertices.pop();
		graph->get(v).setColor(color);
		vector<VERTEX_TYPE>parents=graph->get(v).getParents(v,wordSize);
		for(vector<VERTEX_TYPE>::iterator i=parents.begin();i!=parents.end();i++){
			if(graph->get(*i).getColor()==-1){
				stackOfVertices.push(*i);
			}
		}
		vector<VERTEX_TYPE>children=graph->get(v).getChildren(v,wordSize);
		for(vector<VERTEX_TYPE>::iterator i=children.begin();i!=children.end();i++){
			if(graph->get(*i).getColor()==-1){
				stackOfVertices.push(*i);
			}
		}
	}
}

int main(int argc,char*argv[]){
	DeBruijnAssembler::CommonHeader(&cout);
	cout<<"usage"<<endl;
	cout<<"dna_DeBruijnSplitter [-outputDirectory parts] [-wordSize 21] [-buckets 100000000] [-minimumCoverage auto] <sequence files>"<<endl;
	string outputDirectory="parts";
	int wordSize=21;
	int buckets=100000000;
	string m_minimumCoverageParameter="2";
	vector<string> inputFiles;

	// collect arguments
	for(int i=1;i<argc;i++){
		string option=argv[i];
		if(option=="-outputDirectory"){
			i++;
			outputDirectory=argv[i];
		}else if(option=="-minimumCoverage"){
			i++;
			m_minimumCoverageParameter=argv[i];
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
	if(inputFiles.size()==0){
		cout<<"No files provided."<<endl;
		return 0;
	}


	cout<<"Creating "<<outputDirectory<<endl;
	string command = "mkdir -p "+outputDirectory;
	system(command.c_str());
	CustomMap<int> words(buckets);
	
	cout<<"********** Loading mers from files..."<<endl;
	for(int i=0;i<inputFiles.size();i++){
		vector<Read*> reads;
		Loader loader(&cout);
		loader.load(inputFiles[i],&reads);
		for(int j=0;j<reads.size();j++){
			vector<VERTEX_TYPE> mers=reads.at(j)->getHighQualityMers(wordSize);
			for(int k=0;k<mers.size();k++){
				if(!words.find(mers[k]))
					words.add(mers[k],0);
				words.set(mers[k],words.get(mers[k])+1);
			}
		}
	}




	// depletion curve analysis
	//
	//
	string m_assemblyDirectory=outputDirectory;
	int m_coverage_mean=0;
	int m_minimumCoverage=0;
	//
	//

	CoverageDistribution coverageDistributionObject(&words,outputDirectory);
	m_minimumCoverage=coverageDistributionObject.getMinimumCoverage();
	m_coverage_mean=coverageDistributionObject.getMeanCoverage();

	if(m_minimumCoverageParameter!="auto"){
		m_minimumCoverage=atoi(m_minimumCoverageParameter.c_str());
		cout<<"Setting minimumCoverage <- "<<m_minimumCoverage<<endl;
	}


	cout<<"********** Building graph..."<<endl;

	CustomMap<LightVertex> graphWithoutData(buckets);
	int edges=0;
	for(CustomMap<int>::iterator i=words.begin();i!=words.end();i++){
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
			edges++;
		}
	}
	cout<<"Total: "<<edges<<" edges."<<endl;
	// find connected components, but how?
	int color=1;
	cout<<"********** Finding connected parts..."<<endl;
	for(CustomMap<LightVertex>::iterator i=graphWithoutData.begin();i!=graphWithoutData.end();i++){
		if(graphWithoutData.get(i.first()).getColor()!=-1)
			continue;
		applyColor(i.first(),&graphWithoutData,color,wordSize);
		color++;
	}
	color--;
	cout<<color<<" colors"<<endl;
	map<int,int> colorSizes;
	for(CustomMap<LightVertex>::iterator i=graphWithoutData.begin();i!=graphWithoutData.end();i++){
		colorSizes[i.second().getColor()]++;
	}
	for(map<int,int>::iterator i=colorSizes.begin();i!=colorSizes.end();i++){
		//cout<<i->first<<" "<<i->second<<endl;
	}
	string createFiles=outputDirectory+"/CreateDirectories.sh";
	ofstream streamCreate(createFiles.c_str());
	int Threshold=100;
	for(int i=1;i<=color;i++){
		if(colorSizes[i]<Threshold)
			continue;
		streamCreate<<"mkdir -p "<<outputDirectory<<"/";
		streamCreate<<i<<endl;
	}
	streamCreate.close();
	command="bash "+createFiles;
	system(command.c_str());

	map<string,FILE*> fileStreams;

	cout<<"********** Spitting sequences..."<<endl;
	for(int i=0;i<inputFiles.size();i++){
		vector<Read*> reads;
		Loader loader(&cout);
		loader.load(inputFiles[i],&reads);
		for(int j=0;j<reads.size();j++){
			string sequence=reads.at(j)->getSeq();
			string word=sequence.substr(0,wordSize);
			if(!reads.at(j)->isValidDNA(&word))
				continue;
			VERTEX_TYPE mer=DeBruijnAssembler::wordId(sequence.substr(0,wordSize).c_str());
			VERTEX_TYPE revMer=DeBruijnAssembler::wordId(DeBruijnAssembler::reverseComplement(DeBruijnAssembler::idToWord(mer,wordSize)).c_str());
			string baseName=inputFiles[i].substr(inputFiles[i].find_last_of("/")+1);
			if(graphWithoutData.find(mer)&&colorSizes[graphWithoutData.get(mer).getColor()]>=Threshold){
				ostringstream file;
				file<<outputDirectory<<"/"<<graphWithoutData.get(mer).getColor()<<"/"<<baseName<<".fasta";				
				if(fileStreams.count(file.str())==0){
					//cout<<"adding "<<file.str()<<endl;
					FILE*fp=fopen(file.str().c_str(),"w+");
					fileStreams[file.str()]=fp;
				}
				//cout<<fileStreams[file.str()]<<endl;
				//cout<<"baseName "<<baseName<<endl;
				fprintf(fileStreams[file.str()],">%s\n%s\n",reads.at(j)->getId(),reads.at(j)->getSeq());
			}
			if(graphWithoutData.find(revMer)&&colorSizes[graphWithoutData.get(revMer).getColor()]>=Threshold){
				ostringstream file;
				file<<outputDirectory<<"/"<<graphWithoutData.get(revMer).getColor()<<"/"<<baseName<<".fasta";
				if(fileStreams.count(file.str())==0){
					FILE*fp=fopen(file.str().c_str(),"w+");
					fileStreams[file.str()]=fp;
				}
				if(fileStreams.count(file.str())==0)
					cout<<"Error: no file found."<<endl;

				FILE*filePTR=fileStreams[file.str()];
				if(filePTR==NULL){
					cout<<"Error: "<<file.str()<<" leads to NULL object."<<endl;
				}
				fprintf(filePTR,">%s\n%s\n",reads.at(j)->getId(),reads.at(j)->getSeq());
			}
		}
	}

	for(map<string,FILE*>::iterator i=fileStreams.begin();i!=fileStreams.end();i++){
		fclose(i->second);
	}

	return 0;
}

