#include"Loader.h"
#include"DeBruijnAssembler.h"
#include<stdlib.h>
#include<string>
#include<vector>
#include"Read.h"
#include<iostream>
using namespace std;


class Hit{
	int m_contigNumber;
	int m_contigPosition;
	char m_contigStrand;
	
	int m_readNumber;
	int m_readPosition;
	char m_readStrand;
public:
	Hit(int contigNumber,int contigPosition,char contigStrand,
		int readNumber,int readPosition,char readStrand);
	void show();
	int getContigNumber();
	int getContigPosition();
	char getContigStrand();
};


int Hit::getContigPosition(){
	return m_contigPosition;
}

int Hit::getContigNumber(){
	return m_contigNumber;
}

char Hit::getContigStrand(){
	return m_contigStrand;
}

Hit::Hit(int contigNumber,int contigPosition,char contigStrand,
		int readNumber,int readPosition,char readStrand){
	m_contigNumber=contigNumber;
	m_contigPosition=contigPosition;
	m_contigStrand=contigStrand;
	m_readNumber=readNumber;
	m_readPosition=readPosition;
	m_readStrand=readStrand;
}

void Hit::show(){
	cout<<m_contigNumber<<" "<<m_contigStrand<<" "<<m_contigPosition<<" ~ "<<m_readNumber<<" "<<m_readStrand<<" "<<m_readPosition<<endl;
}


class HitPair{
	Hit*m_left;
	Hit*m_right;
public:
	HitPair(Hit*left,Hit*right);
	HitPair();
	bool valid();
	void show();
};

HitPair::HitPair(Hit*left,Hit*right){
	m_left=left;
	m_right=right;
}

HitPair::HitPair(){
}

bool HitPair::valid(){
	if(m_left->getContigPosition()==0&&m_right->getContigPosition()==0)
		return false;
	if(m_left->getContigStrand()==m_right->getContigStrand()&&
		((m_left->getContigPosition()==0&&m_right->getContigPosition()==0)||
		(m_left->getContigPosition()!=0&&m_right->getContigPosition()!=0)))
		return false;

	return true;
}

void HitPair::show(){
	cout<<"Hitpair"<<endl;
	m_left->show();
	m_right->show();
}

int main(int argc,char*argv[]){
	CommonHeader(&cout);

	if(argc!=5){
		cout<<"Usage"<<endl;
		cout<<"module_join <wordSize> <contigs> <reads> <output.fasta>"<<endl;
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
		string f_end=contigSequence.substr(contigSequence.length()-1-wordSize,wordSize);
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
	string readsFile=argv[3];
	readsLoader.load(readsFile,&theReads);
	string outputRead=argv[4];
	ofstream outputStream(outputRead.c_str());
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
			for(vector<int>::iterator i=f_end_index[wordId(wordAtPosition.c_str())].begin();i!=f_end_index[wordId(wordAtPosition.c_str())].end();i++){
				int contigNumber=*i;
				string contigSequence=contigs[contigNumber]->getSeq();
				Hit myHit(contigNumber,contigSequence.length()-1-wordSize,'F',
					readNumber,readPosition,'F');
				forwardHits.push_back(myHit);
			}
			// f_start
			for(vector<int>::iterator i=f_start_index[wordId(wordAtPosition.c_str())].begin();i!=f_start_index[wordId(wordAtPosition.c_str())].end();i++){
				int contigNumber=*i;
				string contigSequence=contigs[contigNumber]->getSeq();
				Hit myHit(contigNumber,0,'F',
					readNumber,readPosition,'F');
				forwardHits.push_back(myHit);
			}
			// r_end
			for(vector<int>::iterator i=r_end_index[wordId(wordAtPosition.c_str())].begin();i!=r_end_index[wordId(wordAtPosition.c_str())].end();i++){
				int contigNumber=*i;
				string contigSequence=contigs[contigNumber]->getSeq();
				Hit myHit(contigNumber,contigSequence.length()-1-wordSize,'R',
					readNumber,readPosition,'F');
				forwardHits.push_back(myHit);
			}
			// r_start
			for(vector<int>::iterator i=r_start_index[wordId(wordAtPosition.c_str())].begin();i!=r_start_index[wordId(wordAtPosition.c_str())].end();i++){
				int contigNumber=*i;
				string contigSequence=contigs[contigNumber]->getSeq();
				Hit myHit(contigNumber,0,'R',
					readNumber,readPosition,'F');
				forwardHits.push_back(myHit);
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
				outputStream<<">"<<theReads[readNumber]->getId()<<endl;
				outputStream<<theReads[readNumber]->getSeq()<<endl;
			}
		}
		//string reverseSequence=reverseComplement(readSequence);
	}
	outputStream.close();
	return 0;
}

