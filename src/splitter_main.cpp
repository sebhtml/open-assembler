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
#include"GraphDataLight.h"
#include<sstream>
#include"LightVertex.h"
#include"SortedList.h"
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

void applyColor(VERTEX_TYPE v2,GraphDataLight*graph,int color,int wordSize){
	queue<VERTEX_TYPE> stackOfVertices;
	stackOfVertices.push(v2);
	while(stackOfVertices.size()>0){
		VERTEX_TYPE v=stackOfVertices.front();
		stackOfVertices.pop();
		graph->get(v)->setColor(color);
		vector<VERTEX_TYPE>parents=graph->get(v)->getParents(v,wordSize);
		for(vector<VERTEX_TYPE>::iterator i=parents.begin();i!=parents.end();i++){
			if(graph->get(*i)->getColor()!=color){
				stackOfVertices.push(*i);
			}
		}
		vector<VERTEX_TYPE>children=graph->get(v)->getChildren(v,wordSize);
		for(vector<VERTEX_TYPE>::iterator i=children.begin();i!=children.end();i++){
			if(graph->get(*i)->getColor()!=color){
				stackOfVertices.push(*i);
			}
		}
	}
}

int main(int argc,char*argv[]){
	DeBruijnAssembler::CommonHeader(&cout);
	cout<<"usage"<<endl;
	string m_minimumCoverageParameter="2";
	cout<<"dna_DeBruijnSplitter [-outputDirectory parts] [-wordSize 21]  [-minimumCoverage "<<m_minimumCoverageParameter<<"] <sequence files>"<<endl;
	string outputDirectory="parts";
	int wordSize=21;
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
	string command = "rm -rf "+outputDirectory+"; mkdir -p "+outputDirectory;
	system(command.c_str());
	SortedList myList;
	
	cout<<endl;
	cout<<"********** Loading mers from files..."<<endl;
	cout<<endl;
	for(int i=0;i<inputFiles.size();i++){
		cout<<endl;
		cout<<"File: "<<i+1<<" / "<<inputFiles.size()<<endl;
		vector<Read*> reads;
		Loader loader(&cout);
		loader.load(inputFiles[i],&reads);
		for(int j=0;j<reads.size();j++){
			if(j%10000==0){
				cout<<"Read: "<<j+1<<" / "<<reads.size()<<endl;
			}
			vector<VERTEX_TYPE> mers=reads.at(j)->getHighQualityMers(wordSize);
			for(int k=0;k<mers.size();k++){
		/*
				if(!words.find(mers[k]))
					words.add(mers[k],0);
				words.set(mers[k],words.get(mers[k])+1);
	*/
				myList.add(mers[k]);
			}
		}

		cout<<"Read: "<<reads.size()<<" / "<<reads.size()<<endl;
	}
	cout<<endl;
	myList.sort();


	// depletion curve analysis
	//
	//
	string m_assemblyDirectory=outputDirectory;
	int m_coverage_mean=0;
	int m_minimumCoverage=0;
	//
	//
	cout<<endl;

	CoverageDistribution coverageDistributionObject(myList.getDistributionOfCoverage(),outputDirectory);
	m_minimumCoverage=coverageDistributionObject.getMinimumCoverage();
	m_coverage_mean=coverageDistributionObject.getMeanCoverage();

	if(m_minimumCoverageParameter!="auto"){
		m_minimumCoverage=atoi(m_minimumCoverageParameter.c_str());
		cout<<"Setting minimumCoverage <- "<<m_minimumCoverage<<endl;
	}

	cout<<endl;
	cout<<"********** Building graph..."<<endl;
	cout<<endl;
	cout<<"Extracting edges"<<endl;
	vector<VERTEX_TYPE> solidMers=myList.elementsWithALeastCCoverage(m_minimumCoverage);
	GraphDataLight graphWithoutData;
	SortedList graphNodesList;
	for(vector<VERTEX_TYPE>::iterator i=solidMers.begin();i!=solidMers.end();i++){
		VERTEX_TYPE node=*i;
		string wordString=DeBruijnAssembler::idToWord(node,wordSize+1);
		VERTEX_TYPE prefix=DeBruijnAssembler::wordId(wordString.substr(0,wordSize).c_str());
		VERTEX_TYPE suffix=DeBruijnAssembler::wordId(wordString.substr(1,wordSize).c_str());
		graphNodesList.add(prefix);
		graphNodesList.add(suffix);
	}
	graphNodesList.sort();
	cout<<"Extracting vertices"<<endl;
	vector<VERTEX_TYPE> nodes=graphNodesList.elementsWithALeastCCoverage(1);
	graphNodesList.clear();
	for(vector<VERTEX_TYPE>::iterator i=nodes.begin();i!=nodes.end();i++){
		graphWithoutData.add(*i);
	}

	
	int edges=0;
	for(vector<VERTEX_TYPE>::iterator myIterator=solidMers.begin();myIterator!=solidMers.end();myIterator++){
		VERTEX_TYPE mer=*myIterator;
		string merString=DeBruijnAssembler::idToWord(mer,wordSize+1);
		string prefix=merString.substr(0,wordSize);
		string suffix=merString.substr(1,wordSize);
		VERTEX_TYPE prefixInteger=DeBruijnAssembler::wordId(prefix.c_str());
		VERTEX_TYPE suffixInteger=DeBruijnAssembler::wordId(suffix.c_str());
		graphWithoutData.get(prefixInteger)->addChild(suffixInteger,wordSize);
		graphWithoutData.get(suffixInteger)->addParent(prefixInteger,wordSize);
		edges++;
		//}
	}
	cout<<"Total: "<<edges<<" edges."<<endl;
	// find connected components, but how?
	int color=1;
	cout<<endl;
	cout<<"********** Finding connected parts..."<<endl;
	cout<<endl;
	vector<VERTEX_TYPE>*graphNodes=graphWithoutData.getNodes();
	vector<LightVertex>*graphNodeData=graphWithoutData.getNodeData();
	for(int i=0;i<graphNodes->size();i++){
		if(i%1000==0){
			cout<<"Vertex: "<<i+1<<" / "<<graphNodes->size()<<endl;
		}
		if(graphNodeData->at(i).getColor()!=-1)
			continue;
		applyColor(graphNodes->at(i),&graphWithoutData,color,wordSize);
		color++;
	}

	cout<<"Vertex: "<<graphNodes->size()<<" / "<<graphNodes->size()<<endl;
	color--;
	//cout<<color<<" colors"<<endl;
	map<int,int> colorSizes;
	for(int i=0;i<graphNodes->size();i++){
		colorSizes[graphNodeData->at(i).getColor()]++;
	}
	cout<<colorSizes.size()<<" colors"<<endl;
	string partSizes=outputDirectory+"/ConnectedParts.txt";
	ofstream aStreamForParts(partSizes.c_str());
	for(map<int,int>::iterator i=colorSizes.begin();i!=colorSizes.end();i++){
		aStreamForParts<<i->first<<" "<<i->second<<endl;
	}
	aStreamForParts.close();

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

	cout<<endl;
	cout<<"********** Spitting sequences..."<<endl;
	cout<<endl;
	for(int i=0;i<inputFiles.size();i++){
		vector<Read*> reads;
		Loader loader(&cout);
		loader.load(inputFiles[i],&reads);
		string baseName=inputFiles[i].substr(inputFiles[i].find_last_of("/")+1);
		for(int j=0;j<reads.size();j++){
			if(j%1000==0){
				cout<<"Read: "<<j+1<<" / "<<reads.size()<<endl;
			}
			string sequence=reads.at(j)->getSeq();
			set<int> colorsForRead;
			bool gotIt=false;
			for(int k=0;k<sequence.length();k++){
				string word=sequence.substr(k,wordSize);
				if(word.length()!=wordSize)
					continue;
				if(!reads.at(j)->isValidDNA(word.c_str())){
					//cout<<"Invalid start"<<endl;
					continue;
				}
				VERTEX_TYPE mer=DeBruijnAssembler::wordId(word.c_str());
				VERTEX_TYPE revMer=DeBruijnAssembler::wordId(DeBruijnAssembler::reverseComplement(word.c_str()).c_str());
				if(graphWithoutData.find(mer)&&colorSizes[graphWithoutData.get(mer)->getColor()]>=Threshold){
					colorsForRead.insert(graphWithoutData.get(mer)->getColor());
					gotIt=true;
				}
				if(graphWithoutData.find(revMer)&&colorSizes[graphWithoutData.get(revMer)->getColor()]>=Threshold){
					colorsForRead.insert(graphWithoutData.get(revMer)->getColor());
					gotIt=true;
				}
			}
			if(gotIt==false){
				cout<<"Read is nowhere, "<<reads.at(j)->getId()<<endl;
			}

			for(set<int>::iterator colorIterator=colorsForRead.begin();colorIterator!=colorsForRead.end();colorIterator++){
				ostringstream file;
				file<<outputDirectory<<"/"<<*colorIterator<<"/"<<baseName<<".fasta";				
				FILE*fp=fopen(file.str().c_str(),"a+");
				fprintf(fp,">%s\n%s\n",reads.at(j)->getId(),reads.at(j)->getSeq());
				fclose(fp);
			}
		}

		cout<<"Read: "<<reads.size()<<" / "<<reads.size()<<endl;
	}
	cout<<endl;
	cout<<"Done."<<endl;
	return 0;
}

