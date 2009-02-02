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

#include<iostream>
#include"Read.h"
#include<set>
#include<map>
#include<string>
#include"Loader.h"
#include"DeBruijnAssembler.h"

#define wordSize 21

using namespace std;

void applyColor(int i,int color,map<int,int>*contigToColor,map<int,set<int> >*contigs_graph){
	if((*contigToColor)[i]!=-1)
		return;
	(*contigToColor)[i]=color;
	set<int> children=(*contigs_graph)[i];
	for(set<int>::iterator j=children.begin();j!=children.end();j++){
		applyColor(*j,color,contigToColor,contigs_graph);
	}
}

vector<string> overlapper(vector<string> contigSequences){
	bool reducing=true;
	cout<<contigSequences.size()<<" contigs"<<endl;
	int round=1;
	while(reducing){
		cout<<endl;
		cout<<"Round: "<<round<<endl;
		round++;
		map<VERTEX_TYPE,vector<int> > indexOfWords;
		map<VERTEX_TYPE,vector<int> > indexOfRevWords;
		cout<<"Indexing"<<endl;
		for(int i=0;i<contigSequences.size();i++){
			if(i%1000==0)
				cout<<i<<" / "<<contigSequences.size()<<endl;
			for(int j=0;j<contigSequences[i].length();j+=wordSize){
				string word=contigSequences[i].substr(j,wordSize);
				if(word.length()!=wordSize)
					continue;
				string revWord=DeBruijnAssembler::reverseComplement(word);
				indexOfWords[DeBruijnAssembler::wordId(word.c_str())].push_back(i);
				indexOfRevWords[DeBruijnAssembler::wordId(revWord.c_str())].push_back(i);
			}
		}
		cout<<contigSequences.size()<<" / "<<contigSequences.size()<<endl;
		vector<string> nextGeneration;
		set<int> contigsDone;
		cout<<"Merging"<<endl;
		for(int i=0;i<contigSequences.size();i++){
			if(i%100==0)
				cout<<i<<" / "<<contigSequences.size()<<endl;

			if(contigsDone.count(i)!=0)
				continue;

			// case 1: 
			//
			//       -------------------------------------------->       current contig
			//       <--------------------------------------------
			//
			//       					-------------------------------------->
			//       					<-------------------------------------- other contig

			VERTEX_TYPE seed=DeBruijnAssembler::wordId(contigSequences[i].substr(contigSequences[i].length()-wordSize,wordSize).c_str());
			vector<int> otherContigs=indexOfRevWords[seed];
			if(otherContigs.size()==1&&contigsDone.count(otherContigs[0])==0){
				int numberOfMersFound=0;
				int numberOfMersNotFound=0;
				set<VERTEX_TYPE> currentIndex;
				for(int j=0;j<contigSequences[i].length();j++){
					string word=contigSequences[i].substr(j,wordSize);
					if(word.length()!=wordSize)
						continue;
					currentIndex.insert(DeBruijnAssembler::wordId(word.c_str()));
				}
				for(int j=0;j<contigSequences[otherContigs[0]].length();j++){
					string word=contigSequences[otherContigs[0]].substr(j,wordSize);
					if(word.length()!=wordSize)
						continue;
					string revWord=DeBruijnAssembler::reverseComplement(word);
					VERTEX_TYPE theSeed=DeBruijnAssembler::wordId(revWord.c_str());
					if(theSeed==seed){
						numberOfMersFound++;
						break;
					}else if(currentIndex.count(theSeed)>0){
						numberOfMersFound++;
					}else{
						numberOfMersNotFound++;
					}
				}
				if(numberOfMersFound>numberOfMersNotFound){
					ostringstream newSequence;
					string firstToWatch=DeBruijnAssembler::reverseComplement(contigSequences[otherContigs[0]].substr(0,wordSize));
					for(int j=0;j<contigSequences[i].length();j++){
						string word=contigSequences[i].substr(j,wordSize);
						if(word==firstToWatch)
							break;
						newSequence<< contigSequences[i][j];
					}
					newSequence<< DeBruijnAssembler::reverseComplement(contigSequences[otherContigs[0]]);
					nextGeneration.push_back(newSequence.str());
					contigsDone.insert(i);
					contigsDone.insert(otherContigs[0]);
					cout<<"Overlap"<<endl;
				}
			}

		}

		cout<<contigSequences.size()<<" / "<<contigSequences.size()<<endl;
		for(int i=0;i<contigSequences.size();i++){
			if(contigsDone.count(i)==0)
				nextGeneration.push_back(contigSequences[i]);
		}

		reducing=!(nextGeneration.size()==contigSequences.size());
		contigSequences=nextGeneration;
		cout<<contigSequences.size()<<" contigs"<<endl;
	}
	return contigSequences;

}

vector<int> merge(vector<string> contigSequences){
	map<int,set<int> > contigs_graph;
	cout<<contigSequences.size()<<" contigs"<<endl;
	map<VERTEX_TYPE,vector<int> > indexOfWords;
	map<VERTEX_TYPE,vector<int> > indexOfRevWords;
	cout<<"Indexing"<<endl;
	for(int i=0;i<contigSequences.size();i++){
		contigs_graph[i].size();// insert a node
		if(i%1000==0)
			cout<<i<<" / "<<contigSequences.size()<<endl;
		for(int j=0;j<contigSequences[i].length();j+=wordSize){
			string word=contigSequences[i].substr(j,wordSize);
			if(word.length()!=wordSize)
				continue;
			string revWord=DeBruijnAssembler::reverseComplement(word);
			indexOfWords[DeBruijnAssembler::wordId(word.c_str())].push_back(i);
			indexOfRevWords[DeBruijnAssembler::wordId(revWord.c_str())].push_back(i);
		}
	}
	int maxNotFound=2*wordSize;
	//maxNotFound=0;
	cout<<contigSequences.size()<<" / "<<contigSequences.size()<<endl;
	cout<<"Building a graph of contigs."<<endl;
	for(int i=0;i<contigSequences.size();i++){
		if(i%100==0)
			cout<<i<<" / "<<contigSequences.size()<<endl;

		set<int>otherContigs;
		set<int>otherRevContigs;
		for(int j=contigSequences[i].length()/2-3*wordSize;j<contigSequences[i].length()/2+3*wordSize;j+=1){
			string word=contigSequences[i].substr(j,wordSize);
			if(word.length()!=wordSize)
				continue;
			vector<int> otherContigs2=indexOfWords[DeBruijnAssembler::wordId(word.c_str())];
			for(int k=0;k<otherContigs2.size();k++)
				otherContigs.insert(otherContigs2[k]);

			vector<int> otherRevContigs2=indexOfRevWords[DeBruijnAssembler::wordId(word.c_str())];
			for(int k=0;k<otherRevContigs2.size();k++)
				otherRevContigs.insert(otherRevContigs2[k]);
		}
		for(set<int>::iterator matchContig=otherRevContigs.begin();matchContig!=otherRevContigs.end();matchContig++){
			if(i!=*matchContig&&
		contigSequences[i].length()<=contigSequences[*matchContig].length()
		){
				//cout<<"Possible match (Reverse Complement)"<<endl;
				set<VERTEX_TYPE> otherIndex;
				for(int k=0;k<contigSequences[*matchContig].length();k++){
					string word=contigSequences[*matchContig].substr(k,wordSize);
					if(word.length()!=wordSize)
						continue;
					otherIndex.insert(DeBruijnAssembler::wordId(word.c_str()));
				}
				int notFound=0;
				for(int k=0;k<contigSequences[i].length();k++){
					if(notFound>maxNotFound)
						break;
					string word=DeBruijnAssembler::reverseComplement(contigSequences[i].substr(k,wordSize));
					if(word.length()!=wordSize)
						continue;
					if(otherIndex.count(DeBruijnAssembler::wordId(word.c_str()))==0)
						notFound++;
				}
				if(notFound<=maxNotFound){
					//contigs_graph[i].insert(*matchContig);
					contigs_graph[*matchContig].insert(i);
				}
			}
		}

		for(set<int>::iterator matchContig=otherContigs.begin();matchContig!=otherContigs.end();matchContig++){
			if(i!=*matchContig&&
		contigSequences[i].length()<=contigSequences[*matchContig].length()
		){
				//cout<<"Possible match (Foward)"<<endl;
				set<VERTEX_TYPE> otherIndex;
				for(int k=0;k<contigSequences[*matchContig].length();k++){
					string word=contigSequences[*matchContig].substr(k,wordSize);
					if(word.length()!=wordSize)
						continue;
					otherIndex.insert(DeBruijnAssembler::wordId(word.c_str()));
				}
				int notFound=0;
				for(int k=0;k<contigSequences[i].length();k++){
					if(notFound>maxNotFound)
						break;
					string word=(contigSequences[i].substr(k,wordSize));
					if(word.length()!=wordSize)
						continue;
					if(otherIndex.count(DeBruijnAssembler::wordId(word.c_str()))==0)
						notFound++;
				}
				//cout<<"Not found"<<endl;
				if(notFound<=maxNotFound){
					//contigs_graph[i].insert(*matchContig);
					contigs_graph[*matchContig].insert(i);
				}
			}
		}

	}
	cout<<"The graph is ready."<<endl;

	map<int,int> contigToColor;
	
	for(int i=0;i<contigSequences.size();i++){
		contigToColor[i]=-1;
	}
	int color=1;
	
	for(int i=0;i<contigSequences.size();i++){
		applyColor(i,color,&contigToColor,&contigs_graph);
		color++;
	}
	map<int,vector<int> > colorToContigs;
	for(int i=0;i<contigSequences.size();i++){
		colorToContigs[contigToColor[i]].push_back(i);
	}
	
	vector<int> theBestContigs;
	for(map<int,vector<int> >::iterator i=colorToContigs.begin();i!=colorToContigs.end();i++){
		vector<int> contigs=i->second;
		cout<<"Color: "<<i->first<<endl;
		cout<<"   Contigs: ";
		
		int best=contigs[0];
		for(int j=0;j<contigs.size();j++){
			//cout<<contigsWithNames->at(contigs[j])->getId()<<" ";
			cout<<contigs[j]<<" ";
			if(contigSequences[contigs[j]].length()>contigSequences[best].length()){
				best=contigs[j];
			}
		}
		cout<<endl;
		theBestContigs.push_back(best);
	}
	return theBestContigs;
}

int main(int argc,char*argv[]){
	cout<<"This is dna_Merger."<<endl;
	if(argc!=3){
		cout<<"usage"<<endl;
		cout<<"dna_Merger contigs.fa mergedContig.fa"<<endl;
		return 0;
	}
	int minimumContigSize=500;
	string contigsFile=argv[1];
	string outputFile=argv[2];
	vector<Read*> contigs;
	Loader loader(&cout);
	loader.load(contigsFile,&contigs);
	vector<string> contigSequences;
	for(int i=0;i<contigs.size();i++){
		contigSequences.push_back(contigs[i]->getSeq());
	}

	vector<int> theBestContigs=merge(contigSequences);
	//contigSequences=overlapper(contigSequences);
	cout<<endl;

	ofstream f(outputFile.c_str());
	for(int j=0;j<theBestContigs.size();j++){
		int i=theBestContigs[j];
		f<<">"<<contigs.at(i)->getId()<<" "<<contigSequences[i].length()<<" nucleotides"<<endl;
		cout<<">"<<contigs.at(i)->getId()<<" "<<contigSequences[i].length()<<" nucleotides"<<endl;
		int columns=70;
		int j=0;
		while(j<contigSequences[i].length()){
			f<<contigSequences[i].substr(j,columns)<<endl;
			j+=columns;
		}
	}
	f.close();
	return 0;
}

