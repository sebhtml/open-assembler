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
#include"DeBruijnAssembler.h"
#include<stdlib.h>
#include<string>
#include<vector>
#include"Read.h"
#include<iostream>
#include<sstream>
#include"HitPair.h"
#include"Hit.h"
using namespace std;





int main(int argc,char*argv[]){
	CommonHeader(&cout);

	if(argc<5){
		cout<<"Usage"<<endl;
		cout<<argv[0]<<" <wordSize> <contigs> <superContigs> <sequence files>"<<endl;
		return 0;
	}

/*
	f_start						f_end
	r_end						r_start
		--------------------------------------->
		<---------------------------------------

*/

	map<uint64_t,vector<int> > f_end_index;
	map<uint64_t,vector<int> > f_start_index;
	map<uint64_t,vector<int> > r_start_index;
	map<uint64_t,vector<int> > r_end_index;

	int wordSize=atoi(argv[1]);
	vector<Read*> contigs;
	string contigsFile=argv[2];

	Loader loader;
	loader.load(contigsFile,&contigs);

	map<int,set<int> > theGraph;
	map<int,map<int,HitPair> > edgeAnnotations;
	for(int contigNumber=0;contigNumber<contigs.size();contigNumber++){
		set<int> a;
		theGraph[contigNumber]=a;
		string contigSequence=contigs[contigNumber]->getSeq();
		string f_start=contigSequence.substr(0,wordSize);
		// 0 1 2 3 4
		//
		string f_end=contigSequence.substr(contigSequence.length()-wordSize,wordSize);
		string r_end=reverseComplement(f_start);
		string r_start=reverseComplement(f_end);
		f_end_index[wordId(f_end.c_str())].push_back(contigNumber);
		f_start_index[wordId(f_start.c_str())].push_back(contigNumber);
		r_end_index[wordId(r_end.c_str())].push_back(contigNumber);
		r_start_index[wordId(r_start.c_str())].push_back(contigNumber);
	}

	cout<<"Done indexing contigs"<<endl;


	Loader readsLoader;
	vector<Read*> theReads;
	for(int i=4;i<argc;i++){
		string readsFile=argv[i];
		readsLoader.load(readsFile,&theReads);
	}

	//system("sleep 100");
	//string outputRead=argv[4];
	//ofstream outputStream(outputRead.c_str());
	for(int readNumber=0;readNumber<theReads.size();readNumber++){
		if(readNumber%10000==0){
			cout<<readNumber<<" / "<<theReads.size()<<endl;
		}
		string readSequence=theReads[readNumber]->getSeq();
		vector<Hit> forwardHits;
		for(int readPosition=0;readPosition<readSequence.length();readPosition++){
			string wordAtPosition=readSequence.substr(readPosition,wordSize);
			if(wordAtPosition.length()!=wordSize){
				continue;
			}
			// f_end
			if(f_end_index.count(wordId(wordAtPosition.c_str()))>0){
				uint64_t keyValue=wordId(wordAtPosition.c_str());
				vector<int> f_end_vector=f_end_index[keyValue];
				for(vector<int>::iterator i=f_end_vector.begin();i!=f_end_vector.end();i++){
					int contigNumber=*i;
					string contigSequence=contigs[contigNumber]->getSeq();
					Hit myHit(contigNumber,contigSequence.length()-wordSize,'F',
						readNumber,readPosition,'F');
					forwardHits.push_back(myHit);
				}
			}

			// f_start
			if(f_start_index.count(wordId(wordAtPosition.c_str()))>0){
				for(vector<int>::iterator i=f_start_index[wordId(wordAtPosition.c_str())].begin();i!=f_start_index[wordId(wordAtPosition.c_str())].end();i++){
					int contigNumber=*i;
					string contigSequence=contigs[contigNumber]->getSeq();
					Hit myHit(contigNumber,0,'F',
						readNumber,readPosition,'F');
					forwardHits.push_back(myHit);
				}
			}
			// r_end
			if(r_end_index.count(wordId(wordAtPosition.c_str()))>0){
				for(vector<int>::iterator i=r_end_index[wordId(wordAtPosition.c_str())].begin();i!=r_end_index[wordId(wordAtPosition.c_str())].end();i++){
					int contigNumber=*i;
					string contigSequence=contigs[contigNumber]->getSeq();
					Hit myHit(contigNumber,contigSequence.length()-wordSize,'R',
						readNumber,readPosition,'F');
					forwardHits.push_back(myHit);
				}
			}
			// r_start
			if(r_start_index.count(wordId(wordAtPosition.c_str()))>0){
				for(vector<int>::iterator i=r_start_index[wordId(wordAtPosition.c_str())].begin();i!=r_start_index[wordId(wordAtPosition.c_str())].end();i++){
					int contigNumber=*i;
					string contigSequence=contigs[contigNumber]->getSeq();
					Hit myHit(contigNumber,0,'R',
						readNumber,readPosition,'F');
					forwardHits.push_back(myHit);
				}
			}
		}
		if(forwardHits.size()==2&&forwardHits[0].getContigNumber()!=forwardHits[1].getContigNumber()){
			HitPair hitPair(&(forwardHits[0]),&(forwardHits[1]));
			if(hitPair.valid()){
				cout<<endl;
				int leftContig=forwardHits[0].getContigNumber();
				int rightContig=forwardHits[1].getContigNumber();
				if(forwardHits[0].getContigPosition()==0&&forwardHits[1].getContigPosition()!=0){
					HitPair hitPair2(&forwardHits[1],&forwardHits[0]);
					hitPair=hitPair2;
					int t=leftContig;
					leftContig=rightContig;
					rightContig=t;
				}
				hitPair.show();
				theGraph[leftContig].insert(rightContig);
				edgeAnnotations[leftContig][rightContig]=hitPair;
				//outputStream<<">"<<theReads[readNumber]->getId()<<endl;
				//outputStream<<theReads[readNumber]->getSeq()<<endl;
			}
		}


		string reverseSequence=reverseComplement(readSequence);

		vector<Hit> reverseHits;
		for(int readPosition=0;readPosition<reverseSequence.length();readPosition++){
			string wordAtPosition=reverseSequence.substr(readPosition,wordSize);
			if(wordAtPosition.length()!=wordSize){
				continue;
			}
			// f_end
			if(f_end_index.count(wordId(wordAtPosition.c_str()))>0){
				for(vector<int>::iterator i=f_end_index[wordId(wordAtPosition.c_str())].begin();i!=f_end_index[wordId(wordAtPosition.c_str())].end();i++){
					int contigNumber=*i;
					string contigSequence=contigs[contigNumber]->getSeq();
					Hit myHit(contigNumber,contigSequence.length()-wordSize,'F',
						readNumber,readPosition,'R');
					reverseHits.push_back(myHit);
				}
			}
			// f_start
			if(f_start_index.count(wordId(wordAtPosition.c_str()))>0){
				for(vector<int>::iterator i=f_start_index[wordId(wordAtPosition.c_str())].begin();i!=f_start_index[wordId(wordAtPosition.c_str())].end();i++){
					int contigNumber=*i;
					string contigSequence=contigs[contigNumber]->getSeq();
					Hit myHit(contigNumber,0,'F',
						readNumber,readPosition,'R');
					reverseHits.push_back(myHit);
				}
			}
			// r_end
			if(r_end_index.count(wordId(wordAtPosition.c_str()))>0){
				for(vector<int>::iterator i=r_end_index[wordId(wordAtPosition.c_str())].begin();i!=r_end_index[wordId(wordAtPosition.c_str())].end();i++){
					int contigNumber=*i;
					string contigSequence=contigs[contigNumber]->getSeq();
					Hit myHit(contigNumber,contigSequence.length()-wordSize,'R',
						readNumber,readPosition,'R');
					reverseHits.push_back(myHit);
				}
			}
			// r_start
			if(r_start_index.count(wordId(wordAtPosition.c_str()))>0){
				for(vector<int>::iterator i=r_start_index[wordId(wordAtPosition.c_str())].begin();i!=r_start_index[wordId(wordAtPosition.c_str())].end();i++){
					int contigNumber=*i;
					string contigSequence=contigs[contigNumber]->getSeq();
					Hit myHit(contigNumber,0,'R',
						readNumber,readPosition,'R');
					reverseHits.push_back(myHit);
				}
			}
		}
		if(reverseHits.size()==2&&reverseHits[0].getContigNumber()!=reverseHits[1].getContigNumber()){
			HitPair hitPair(&(reverseHits[0]),&(reverseHits[1]));
			if(hitPair.valid()){
				cout<<endl;
				int leftContig=reverseHits[0].getContigNumber();
				int rightContig=reverseHits[1].getContigNumber();
				if(reverseHits[0].getContigPosition()==0&&reverseHits[1].getContigPosition()!=0){
					HitPair hitPair2(&reverseHits[1],&reverseHits[0]);
					hitPair=hitPair2;
					int t=leftContig;
					leftContig=rightContig;
					rightContig=t;
				}
				hitPair.show();
				theGraph[leftContig].insert(rightContig);
				edgeAnnotations[leftContig][rightContig]=hitPair;
				//outputStream<<">"<<theReads[readNumber]->getId()<<endl;
				//outputStream<<theReads[readNumber]->getSeq()<<endl;
			}
		}
	}
	//outputStream.close();

	ofstream fGraphViz("Graphviz.txt");

		

	fGraphViz<<"digraph{"<<endl;
	
	string outputFile=argv[3];
	ofstream fStream(outputFile.c_str());
	set<int> _heads; // where to start.
	set<int> _mixed; // mixed vertices
	set<int> _doneContig;
	for(int contigNumber=0;contigNumber<contigs.size();contigNumber++){
		if(theGraph[contigNumber].size()>2)
			_mixed.insert(contigNumber);
	}

	for(int contigNumber=0;contigNumber<contigs.size();contigNumber++){
		if(theGraph[contigNumber].size()==1||theGraph[contigNumber].size()==0)
			_heads.insert(contigNumber);
		if(theGraph[contigNumber].size()==2){
			set<int>::iterator anIterator=theGraph[contigNumber].begin();
			int theFirst=*anIterator;
			anIterator++;
			int theSecond=*anIterator;
			if(_mixed.count(theFirst)>0||_mixed.count(theSecond)>0)
				_heads.insert(contigNumber);
		}
	}

	for(int contigNumber=0;contigNumber<contigs.size();contigNumber++){
		if(_doneContig.count(contigNumber)>0)
			continue;
		set<int> links=theGraph[contigNumber];
		fGraphViz<<"c"<<contigNumber<<endl;
		for(set<int>::iterator i=links.begin();i!=links.end();i++){
			fGraphViz<<"c"<<contigNumber<<" -> c"<<*i<<endl;
		}
		if(_heads.count(contigNumber)>0||_mixed.count(contigNumber)>0){
			cout<<"Building new Contig"<<endl;
			ostringstream contigName;
			contigName<<">SuperContig_";
			int currentContig=contigNumber;
			int contigStartToPrintAtPosition=0;
			ostringstream contigSequence;
			char contigStrand='F';
			bool canChange=true;
			while(currentContig!=-1&&_doneContig.count(currentContig)==0){
				_doneContig.insert(currentContig);
				// show the contig
				int nextContigToGet=-1;
				char nextStrandToGet='F';
				int nextContigStartPosition=0;
				bool hasNext=false;
				string sequenceFromRead="";
				if(theGraph[currentContig].size()>0&&_mixed.count(currentContig)==0){
					for(set<int>::iterator k=theGraph[currentContig].begin();k!=theGraph[currentContig].end();k++){
						int nextContig=*k;
						HitPair aHitPair=edgeAnnotations[currentContig][nextContig];
						if(aHitPair.getLeft()->getContigStrand()!=contigStrand&&canChange){
							canChange=false;
							contigStrand=aHitPair.getLeft()->getContigStrand();
						}
						if(aHitPair.getLeft()->getContigStrand()==contigStrand){//found it yes
							nextContigToGet=aHitPair.getRight()->getContigNumber();
							nextStrandToGet=aHitPair.getRight()->getContigStrand();
							int theDistance=aHitPair.getRight()->getReadPosition()-(aHitPair.getLeft()->getReadPosition()+wordSize-1)-1;
							cout<<"Distance <- "<<theDistance<<endl;
							hasNext=true;
							string readSequenceForGap=theReads[aHitPair.getLeft()->getReadNumber()]->getSeq();
							if(aHitPair.getLeft()->getReadStrand()=='R'){
								readSequenceForGap=reverseComplement(readSequenceForGap);
							}
							sequenceFromRead=readSequenceForGap.substr(aHitPair.getLeft()->getReadPosition()+wordSize,theDistance);
							if(theDistance>=0){
								nextContigStartPosition=0;
							}else{
								nextContigStartPosition=-theDistance;
							}
							break;
						}
					}
				}
				cout<<"Contig: "<<currentContig<<" "<<contigStrand<<endl;
				contigName<<"_"<<contigs[currentContig]->getId();
				string sequence=contigs[currentContig]->getSeq();
				if(contigStrand=='R')
					sequence=reverseComplement(sequence);

				contigSequence<<sequence.substr(contigStartToPrintAtPosition);
				
				if(currentContig!=-1){
					//contigSequence<<"NNNNNNNNNNNNNN";
				}
				if(hasNext){
					if(nextContigStartPosition==0){
						//contigSequence<<"_INSERT_";
						contigSequence<<sequenceFromRead;
					}else{
						//contigSequence<<"_OVERLAP_";
					}
				}
				currentContig=nextContigToGet;
				contigStrand=nextStrandToGet;
				contigStartToPrintAtPosition=nextContigStartPosition;
				if(_mixed.count(currentContig)>0)
					currentContig=-1;
			}
			fStream<<contigName.str()<<endl;
			int pp=0;
			int columns=60;
			string superContigSequence=contigSequence.str();
			while(pp<superContigSequence.length()){
				fStream<<superContigSequence.substr(pp,columns)<<endl;
				pp+=columns;
			}
		}
	}
	fStream.close();
	fGraphViz<<"}"<<endl;
	fGraphViz.close();
	return 0;
}

