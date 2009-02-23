#include"Loader.h"
#include"DeBruijnAssembler.h"
#include<stdlib.h>
#include<string>
#include<vector>
#include"Read.h"
#include<iostream>
using namespace std;


int main(int argc,char*argv[]){
	DeBruijnAssembler::CommonHeader(&cout);

	if(argc!=4){
		cout<<"Usage"<<endl;
		cout<<"module_join <wordSize>A <contigs> <reads>"<<endl;
		return 0;
	}
	Loader readLoader;
	Loader contigsLoader;
	int wordSize=atoi(argv[1]);
	vector<Read*> contigs;
	vector<Read*> reads;
	string contigsFile=argv[2];
	string readsFile=argv[3];
	readLoader.load(readsFile,&reads);
	contigsLoader.load(contigsFile,&contigs);

	cout<<endl;
	cout<<"Reads: "<<reads.size()<<endl;
	cout<<"Contigs: "<<contigs.size()<<endl;
	cout<<"WordSize: "<<wordSize<<endl;

	map<VERTEX_TYPE,vector<int> >contigHeads;
	map<VERTEX_TYPE,vector<int> >contigTails;

	for(int contigNumber=0;contigNumber<contigs.size();contigNumber++){
		string sequenceOfContig=contigs.at(contigNumber)->getSeq();
		string contigHeadString=sequenceOfContig.substr(0,wordSize);
		string contigTailString=sequenceOfContig.substr(sequenceOfContig.length()-1-wordSize,wordSize);
		VERTEX_TYPE headInInt=DeBruijnAssembler::wordId(contigHeadString.c_str());
		VERTEX_TYPE tailInInt=DeBruijnAssembler::wordId(contigTailString.c_str());
		contigHeads[headInInt].push_back(contigNumber);
		contigTails[tailInInt].push_back(contigNumber);
	}
	
	cout<<endl;
	cout<<"Done indexing contigs"<<endl;

	map<VERTEX_TYPE,map<int,int> > readIndexForward;
	map<VERTEX_TYPE,map<int,int> > readIndexReverse;

	for(int readNumber=0;readNumber<reads.size();readNumber++){
		if(readNumber%100000==0){
			cout<<readNumber<<" / "<<reads.size()<<endl;
		}
		string readSequence=reads[readNumber]->getSeq();
		for(int readPosition=0;readPosition<readSequence.length();readPosition++){
			string wordSequence=readSequence.substr(readPosition,wordSize);
			if(wordSequence.length()!=wordSize)
				continue;
			VERTEX_TYPE wordInInt=DeBruijnAssembler::wordId(wordSequence.c_str());
			if(contigHeads.count(wordInInt)>0||contigTails.count(wordInInt)>0){
				readIndexForward[wordInInt][readNumber]=readPosition;
			}
			
			string wordSequenceReverse=DeBruijnAssembler::reverseComplement(wordSequence);
			VERTEX_TYPE wordInIntReverse=DeBruijnAssembler::wordId(wordSequenceReverse.c_str());
			if(contigHeads.count(wordInIntReverse)>0||contigTails.count(wordInIntReverse)>0){
				readIndexReverse[wordInIntReverse][readNumber]=readPosition;
			}
		}
	}

	cout<<reads.size()<<" / "<<reads.size()<<endl;

	cout<<"Done indexing reads"<<endl;
	return 0;
}

