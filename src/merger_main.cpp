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

#define wordSize 15

using namespace std;

vector<string> merge(vector<string> contigSequences){
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
			for(int j=0;j<contigSequences[i].length();j++){
				string word=contigSequences[i].substr(j,wordSize);
				if(word.length()!=wordSize)
					continue;
				string revWord=DeBruijnAssembler::reverseComplement(word);
				indexOfWords[DeBruijnAssembler::wordId(word.c_str())].push_back(i);
				indexOfRevWords[DeBruijnAssembler::wordId(revWord.c_str())].push_back(i);
			}
		}
		vector<string> nextGeneration;
		set<int> contigsDone;
		cout<<"Merging"<<endl;
		for(int i=0;i<contigSequences.size();i++){
			set<int>otherContigs;
			set<int>otherRevContigs;
			for(int j=0;j<contigSequences[i].length();j++){
				if(contigsDone.count(i)!=0)
					break;
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
			&&contigsDone.count(*matchContig)==0 &&
			contigsDone.count(i)==0){
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
						string word=DeBruijnAssembler::reverseComplement(contigSequences[i].substr(k,wordSize));
						if(word.length()!=wordSize)
							continue;
						if(otherIndex.count(DeBruijnAssembler::wordId(word.c_str()))==0)
							notFound++;
					}
					if(notFound<2*wordSize){
						contigsDone.insert(i);
						contigsDone.insert(*matchContig);
						nextGeneration.push_back(contigSequences[*matchContig]);
						cout<<"Merging (notFound="<<notFound<<")"<<endl;
					}
				}
			}

			for(set<int>::iterator matchContig=otherContigs.begin();matchContig!=otherContigs.end();matchContig++){
				if(i!=*matchContig&&
			contigSequences[i].length()<=contigSequences[*matchContig].length()
			&&contigsDone.count(*matchContig)==0 &&
			contigsDone.count(i)==0){
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
						string word=(contigSequences[i].substr(k,wordSize));
						if(word.length()!=wordSize)
							continue;
						if(otherIndex.count(DeBruijnAssembler::wordId(word.c_str()))==0)
							notFound++;
					}
					if(notFound<2*wordSize){
						contigsDone.insert(i);
						contigsDone.insert(*matchContig);
						nextGeneration.push_back(contigSequences[*matchContig]);
						cout<<"Merging (notFound="<<notFound<<")"<<endl;
					}
				}
			}

		}

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

int main(int argc,char*argv[]){
	cout<<"This is dna_Merger."<<endl;
	if(argc!=3){
		cout<<"usage"<<endl;
		cout<<"dna_Merger contigs.fa mergedContig.fa"<<endl;
		return 0;
	}
	string contigsFile=argv[1];
	string outputFile=argv[2];
	vector<Read*> contigs;
	Loader loader(&cout);
	loader.load(contigsFile,&contigs);
	vector<string> contigSequences;
	for(int i=0;i<contigs.size();i++)
		contigSequences.push_back(contigs[i]->getSeq());

	contigSequences=merge(contigSequences);
	cout<<endl;

	ofstream f(outputFile.c_str());
	for(int i=0;i<contigSequences.size();i++){
		f<<">Contig"<<i+1<<" "<<contigSequences[i].length()<<" nucleotides"<<endl;
		cout<<">Contig"<<i+1<<" "<<contigSequences[i].length()<<" nucleotides"<<endl;
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

