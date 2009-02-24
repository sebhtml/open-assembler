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
#include<iostream>
#include<stack>
#include"Read.h"
#include<stdlib.h>
#include"DeBruijnAssembler.h"
using namespace std;

int main(int argc,char*argv[]){
	DeBruijnAssembler::CommonHeader(&cout);
	if(argc!=7){
		cout<<"usage"<<endl;
		cout<<"dna_Scaffolder left.fasta right.fasta insert stddev contigs.fasta scaffolds.fasta"<<endl;
		return 0;
	}
	int MINIMUM_TILING_WINDOWS=3;
	string leftReadsFile=argv[1];
	string rightReadsFile=argv[2];
	int insertSize=atoi(argv[3]);
	int stddev=atoi(argv[4]);
	string contigsFile=argv[5];
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
			VERTEX_TYPE leftWord=DeBruijnAssembler::wordId(word.c_str());
			leftIndex[leftWord][i]=j;
		}
		for(int j=0;j<sequenceRight.length();j++){
			string word=sequenceRight.substr(j,wordSize);
			if(word.length()!=wordSize||!rightReads[i]->isValidDNA(word.c_str()))
				continue;
			VERTEX_TYPE rightWord=DeBruijnAssembler::wordId(word.c_str());
			rightIndex[rightWord][i]=j;
		}

	}
	cout<<leftReads.size()<<" / "<<leftReads.size()<<endl;

	cout<<"Done indexing paired mates., "<<leftIndex.size()<<" left seeds. "<<rightIndex.size()<<" right seeds."<<endl;

	cout<<"Processing contigs now.."<<endl;

	//  read,  contig  position
	map<int,map<int,vector<int> > > leftForwardHits;
	map<int,map<int,vector<int> > >rightForwardHits;
	map<int,map<int,vector<int> > >leftReverseHits;
	map<int,map<int,vector<int> > >rightReverseHits;

	for(int contigNumber=0;contigNumber<contigs.size();contigNumber++){
		if(contigNumber%100==0){
			cout<<contigNumber<<" / "<<contigs.size()<<endl;
		}
		string sequence=contigs[contigNumber]->getSeq();
		for(int sequencePosition=0;sequencePosition<sequence.length();sequencePosition++){
			string word=sequence.substr(sequencePosition,wordSize);
			if(word.length()!=wordSize||!leftReads[0]->isValidDNA(word.c_str()))
				continue;
			string reverseWordSeq=DeBruijnAssembler::reverseComplement(word);
			VERTEX_TYPE forwardWord=DeBruijnAssembler::wordId(word.c_str());
			VERTEX_TYPE reverseWord=DeBruijnAssembler::wordId(reverseWordSeq.c_str());
			if(leftIndex.count(forwardWord)>0){
				for(map<int,int>::iterator i=leftIndex[forwardWord].begin();i!=leftIndex[forwardWord].end();i++){
					int readNumber=i->first;
					leftForwardHits[readNumber][contigNumber].push_back(sequencePosition);
				}
			}
			if(leftIndex.count(reverseWord)>0){
				for(map<int,int>::iterator i=leftIndex[reverseWord].begin();i!=leftIndex[reverseWord].end();i++){
					int readNumber=i->first;
					leftReverseHits[readNumber][contigNumber].push_back(sequencePosition);
				}
			}
			if(rightIndex.count(forwardWord)>0){
				for(map<int,int>::iterator i=rightIndex[forwardWord].begin();i!=rightIndex[forwardWord].end();i++){
					int readNumber=i->first;
					rightForwardHits[readNumber][contigNumber].push_back(sequencePosition);
				}
			}
			if(rightIndex.count(reverseWord)>0){
				for(map<int,int>::iterator i=rightIndex[reverseWord].begin();i!=rightIndex[reverseWord].end();i++){
					int readNumber=i->first;
					rightReverseHits[readNumber][contigNumber].push_back(sequencePosition);
				}
			}
		}
	}

	cout<<contigs.size()<<" / "<<contigs.size()<<endl;


	cout<<"Scaffolding now!"<<endl;
	map<int,set<int> > theScaffolderGraphChildren;
	map<int,set<int> > theScaffolderGraphParents;

	for(int contigNumber=0;contigNumber<contigs.size();contigNumber++){
		set<int> anEntry;
		theScaffolderGraphParents[contigNumber]=anEntry;
		theScaffolderGraphChildren[contigNumber]=anEntry;
	}

	int activePairedReads=0;
	int notOnTheSame=0;
	for(int readNumber=0;readNumber<leftReads.size();readNumber++){
		cout<<"Read: "<<readNumber<<endl;
		vector<int> validOnLeftForward;
		vector<int> validOnLeftReverse;
		vector<int> validOnRightForward;
		vector<int> validOnRightReverse;
		if(leftForwardHits.count(readNumber)>0){
			cout<<"leftForwardHits"<<endl;
			for(map<int,vector<int> >::iterator i=leftForwardHits[readNumber].begin();i!=leftForwardHits[readNumber].end();i++){
				if(i->second.size()>=MINIMUM_TILING_WINDOWS)
					validOnLeftForward.push_back(i->first);
				cout<<"Contig "<<i->first<<endl;
				cout<<"Positions: ";
				for(int j=0;j<i->second.size();j++){
					cout<<" "<<i->second[j];
				}
				cout<<endl;
			}
		}
		if(rightForwardHits.count(readNumber)>0){
			cout<<"rightForwardHits"<<endl;
			for(map<int,vector<int> >::iterator i=rightForwardHits[readNumber].begin();i!=rightForwardHits[readNumber].end();i++){
				if(i->second.size()>=MINIMUM_TILING_WINDOWS)
					validOnRightForward.push_back(i->first);
				cout<<"Contig "<<i->first<<endl;
				cout<<"Positions: ";
				for(int j=0;j<i->second.size();j++){
					cout<<" "<<i->second[j];
				}
				cout<<endl;
			}
		}
		if(leftReverseHits.count(readNumber)>0){
			cout<<"leftReverseHits"<<endl;
			for(map<int,vector<int> >::iterator i=leftReverseHits[readNumber].begin();i!=leftReverseHits[readNumber].end();i++){
				if(i->second.size()>=MINIMUM_TILING_WINDOWS)
					validOnLeftReverse.push_back(i->first);
				cout<<"Contig "<<i->first<<endl;
				cout<<"Positions: ";
				for(int j=0;j<i->second.size();j++){
					cout<<" "<<i->second[j];
				}
				cout<<endl;
			}
		}
		if(rightReverseHits.count(readNumber)>0){
			cout<<"rightReverseHits"<<endl;
			for(map<int,vector<int> >::iterator i=rightReverseHits[readNumber].begin();i!=rightReverseHits[readNumber].end();i++){
				if(i->second.size()>=MINIMUM_TILING_WINDOWS)
					validOnRightReverse.push_back(i->first);
				cout<<"Contig "<<i->first<<endl;
				cout<<"Positions: ";
				for(int j=0;j<i->second.size();j++){
					cout<<" "<<i->second[j];
				}
				cout<<endl;
			}
		}
		cout<<"validOnLeftForward "<<validOnLeftForward.size()<<endl;
		cout<<"validOnLeftReverse "<<validOnLeftReverse.size()<<endl;
		cout<<"validOnRightForward "<<validOnRightForward.size()<<endl;
		cout<<"validOnRightReverse "<<validOnRightReverse.size()<<endl;
		if(validOnLeftForward.size()+validOnLeftReverse.size()==1&&validOnRightForward.size()+validOnRightReverse.size()==1){
			activePairedReads++;
			int leftContig;
			int rightContig;
			if(validOnLeftForward.size()==1)
				leftContig=validOnLeftForward[0];
			else
				leftContig=validOnLeftReverse[0];
			if(validOnRightForward.size()==1)
				rightContig=validOnRightForward[0];
			else
				rightContig=validOnRightReverse[0];
			if(leftContig!=rightContig)
				notOnTheSame++;
			theScaffolderGraphChildren[leftContig].insert(rightContig);
			theScaffolderGraphParents[rightContig].insert(leftContig);
		}
	}
	cout<<"activePairedReads "<<activePairedReads<<endl;
	cout<<"not on the same "<<notOnTheSame<<endl;
	
	map<int,int> scaffoldColors;
	int currentColor=1;
	for(int contigNumber=0;contigNumber<contigs.size();contigNumber++){
		if(scaffoldColors.count(contigNumber)==0){
			stack<int> toDoList;
			toDoList.push(contigNumber);
			while(toDoList.size()>0){
				int aContigNumber=toDoList.top();
				toDoList.pop();
				scaffoldColors[aContigNumber]=currentColor;
				set<int> children=theScaffolderGraphChildren[aContigNumber];
				set<int> parents=theScaffolderGraphParents[aContigNumber];
				for(set<int>::iterator i=children.begin();i!=children.end();i++){
					toDoList.push(*i);
				}
				for(set<int>::iterator i=parents.begin();i!=parents.end();i++){
					toDoList.push(*i);
				}
			}
			currentColor++;
		}
	}

	cout<<"Scaffolds: "<<currentColor-1<<endl;

	return 0;
}

