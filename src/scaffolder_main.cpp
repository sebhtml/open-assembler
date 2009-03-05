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

#include"Loader.h"
#include"Hit.h"
#include<iostream>
#include"HitPair.h"
#include<stack>
#include"Read.h"
#include<stdlib.h>
#include"DeBruijnAssembler.h"
using namespace std;


int main(int argc,char*argv[]){
	CommonHeader(&cout);
	if(argc!=5){
		cout<<"usage"<<endl;
		cout<<argv[0]<<" left.fasta right.fasta contigs.fasta scaffolds.fasta"<<endl;
		return 0;
	}
	string leftReadsFile=argv[1];
	string rightReadsFile=argv[2];
	string contigsFile=argv[3];
	int wordSize=21;
	vector<Read*> contigs;
	Loader loader;
	loader.load(contigsFile,&contigs);
	vector<Read*> leftReads;
	vector<Read*> rightReads;
	Loader loaderReads;
	loaderReads.load(leftReadsFile,&leftReads);
	loaderReads.load(rightReadsFile,&rightReads);

	// word           read     position
	map<VERTEX_TYPE,map<int,int> > leftIndex;
	map<VERTEX_TYPE,map<int,int> > rightIndex;
	
	for(int i=0;i<leftReads.size();i++){
		if(i%10000==0){
			cout<<i<<" / "<<leftReads.size()<<endl;
		}
		string sequenceLeft=leftReads[i]->getSeq();
		string sequenceRight=rightReads[i]->getSeq();
		for(int j=0;j<sequenceLeft.length();j++){
			string word=sequenceLeft.substr(j,wordSize);
			if(word.length()!=wordSize||!leftReads[i]->isValidDNA(word.c_str()))
				continue;
			VERTEX_TYPE leftWord=wordId(word.c_str());
			leftIndex[leftWord][i]=j;
		}
		for(int j=0;j<sequenceRight.length();j++){
			string word=sequenceRight.substr(j,wordSize);
			if(word.length()!=wordSize||!rightReads[i]->isValidDNA(word.c_str()))
				continue;
			VERTEX_TYPE rightWord=wordId(word.c_str());
			rightIndex[rightWord][i]=j;
		}

	}
	cout<<leftReads.size()<<" / "<<leftReads.size()<<endl;

	cout<<"Done indexing paired mates., "<<leftIndex.size()<<" left seeds. "<<rightIndex.size()<<" right seeds."<<endl;

	cout<<"Processing contigs now.."<<endl;

	vector<Hit> allHits;

	for(int contigNumber=0;contigNumber<contigs.size();contigNumber++){
		if(contigNumber%100==0){
			cout<<contigNumber<<" / "<<contigs.size()<<endl;
		}
		string sequence=contigs[contigNumber]->getSeq();
		for(int sequencePosition=0;sequencePosition<sequence.length();sequencePosition++){
			string word=sequence.substr(sequencePosition,wordSize);
			if(word.length()!=wordSize||!leftReads[0]->isValidDNA(word.c_str()))
				continue;
			VERTEX_TYPE forwardWord=wordId(word.c_str());
			if(leftIndex.count(forwardWord)>0){
				for(map<int,int>::iterator i=leftIndex[forwardWord].begin();i!=leftIndex[forwardWord].end();i++){
					int readNumber=i->first;
					int readPosition=i->second;
					Hit hit(contigNumber,sequencePosition,'F',readNumber,readPosition,'F');
					hit.setL_or_R('L');
					allHits.push_back(hit);
				}
			}
			if(rightIndex.count(forwardWord)>0){
				for(map<int,int>::iterator i=rightIndex[forwardWord].begin();i!=rightIndex[forwardWord].end();i++){
					int readNumber=i->first;
					int readPosition=i->second;
					Hit hit(contigNumber,sequencePosition,'F',readNumber,readPosition,'F');
					hit.setL_or_R('R');
					allHits.push_back(hit);
				}
			}
		}
		sequence=reverseComplement(sequence);
		for(int sequencePosition=0;sequencePosition<sequence.length();sequencePosition++){
			string word=sequence.substr(sequencePosition,wordSize);
			if(word.length()!=wordSize||!leftReads[0]->isValidDNA(word.c_str()))
				continue;
			VERTEX_TYPE forwardWord=wordId(word.c_str());
			if(leftIndex.count(forwardWord)>0){
				for(map<int,int>::iterator i=leftIndex[forwardWord].begin();i!=leftIndex[forwardWord].end();i++){
					int readNumber=i->first;
					int readPosition=i->second;
					Hit hit(contigNumber,sequencePosition,'R',readNumber,readPosition,'F');
					hit.setL_or_R('L');
					allHits.push_back(hit);
				}
			}
			if(rightIndex.count(forwardWord)>0){
				for(map<int,int>::iterator i=rightIndex[forwardWord].begin();i!=rightIndex[forwardWord].end();i++){
					int readNumber=i->first;
					int readPosition=i->second;
					Hit hit(contigNumber,sequencePosition,'R',readNumber,readPosition,'F');
					hit.setL_or_R('R');
					allHits.push_back(hit);
				}
			}
		}
	}

	cout<<contigs.size()<<" / "<<contigs.size()<<endl;

	cout<<"Done with contigs"<<endl;
	cout<<"Analyzing hits ("<<allHits.size()<<")"<<endl;

	map<int,map<char,vector<int> > > readHits;
	for(int hitId=0;hitId<allHits.size();hitId++){
		readHits[allHits[hitId].getReadNumber()][allHits[hitId].side()].push_back(hitId);
	}

	cout<<readHits.size()<<" paired reads with hits"<<endl;

	vector<HitPair> hitPairs;
	int bothSidesNotSame=0;
	int bothMapped=0;
	for(map<int,map<char,vector<int> > >::iterator i=readHits.begin();i!=readHits.end();i++){
		int readNumber=i->first;
		vector<int> leftHits;
		if(i->second.count('L')>0)
			leftHits=i->second['L'];
		
		vector<int> rightHits;
		if(i->second.count('R')>0)
			rightHits=i->second['R'];
		if(leftHits.size()>0&&rightHits.size()>0){
			bothMapped++;
			int leftContig=allHits[leftHits[0]].getContigNumber();
			for(int hh=0;hh<leftHits.size();hh++){
				if(leftContig!=allHits[leftHits[hh]].getContigNumber()){
					leftContig=-1;
					break;
				}
			}
			// TODO accept mixed hit pools
			int rightContig=allHits[rightHits[0]].getContigNumber();
			for(int hh=0;hh<rightHits.size();hh++){
				if(rightContig!=allHits[rightHits[hh]].getContigNumber()){
					rightContig=-1;
					break;
				}
			}
			if(leftContig!=-1&&rightContig!=-1&&leftContig!=rightContig){
				bothSidesNotSame++;
				HitPair aPair(&(allHits[leftHits[0]]),&(allHits[rightHits[0]]));
				hitPairs.push_back(aPair);
			}
		}
	}

	cout<<"Total read pairs: "<<leftReads.size()<<endl;
	cout<<"Total read pairs with both mapped parts: "<<bothMapped<<endl;
	cout<<"Total read pairs with both mapped parts on different contigs: "<<bothSidesNotSame<<endl;

	ofstream graphFile("graph.txt");
	graphFile<<"digraph{"<<endl;
	for(vector<HitPair>::iterator i=hitPairs.begin();i!=hitPairs.end();i++){
		HitPair pair=*i;
		graphFile<<"c"<<pair.getLeft()->getContigNumber()<<" -> c"<<pair.getRight()->getContigNumber()<<" [ label=\"("<<pair.getLeft()->getContigStrand()<<","<<pair.getRight()->getContigStrand()<<")\" ]"<<endl;
	}
	graphFile<<"}"<<endl;
	graphFile.close();
	
	cout<<"Scaffolding now!"<<endl;

	return 0;
}

