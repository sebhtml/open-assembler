#include"Loader.h"
#include"DeBruijnAssembler.h"
#include<stdlib.h>
#include<string>
#include<vector>
#include"Read.h"
#include<iostream>
#include<sstream>
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
	Hit();
	void show();
	int getContigNumber();
	int getContigPosition();
	char getContigStrand();
	int getReadPosition();
	int getReadNumber();
	char getReadStrand();
};

char Hit::getReadStrand(){
	return m_readStrand;
}

int Hit::getReadNumber(){
	return m_readNumber;
}

int Hit::getReadPosition(){
	return m_readPosition;
}

Hit::Hit(){
}

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
	Hit m_left;
	Hit m_right;
public:
	HitPair(Hit*left,Hit*right);
	HitPair();
	bool valid();
	void show();
	Hit*getLeft();
	Hit*getRight();
};

Hit*HitPair::getLeft(){
	return &m_left;
}

Hit*HitPair::getRight(){
	return &m_right;
}

HitPair::HitPair(Hit*left,Hit*right){
	m_left=*left;
	m_right=*right;
}

HitPair::HitPair(){
}

bool HitPair::valid(){
	if(m_left.getContigPosition()==0&&m_right.getContigPosition()==0)
		return false;
	if(//m_left.getContigStrand()==m_right.getContigStrand()&&
		((m_left.getContigPosition()==0&&m_right.getContigPosition()==0)||
		(m_left.getContigPosition()!=0&&m_right.getContigPosition()!=0)))
		return false;

	return true;
}

void HitPair::show(){
	cout<<"Hitpair"<<endl;
	m_left.show();
	m_right.show();
}

int main(int argc,char*argv[]){
	CommonHeader(&cout);

	if(argc!=5){
		cout<<"Usage"<<endl;
		cout<<"module_join <wordSize> <contigs> <reads> <superContigs>"<<endl;
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
	map<int,set<int> > theParents;
	map<int,map<int,HitPair> > edgeAnnotations;
	for(int contigNumber=0;contigNumber<contigs.size();contigNumber++){
		set<int> a;
		theParents[contigNumber]=a;
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
	string readsFile=argv[3];
	readsLoader.load(readsFile,&theReads);
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
			for(vector<int>::iterator i=f_end_index[wordId(wordAtPosition.c_str())].begin();i!=f_end_index[wordId(wordAtPosition.c_str())].end();i++){
				int contigNumber=*i;
				string contigSequence=contigs[contigNumber]->getSeq();
				Hit myHit(contigNumber,contigSequence.length()-wordSize,'F',
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
				Hit myHit(contigNumber,contigSequence.length()-wordSize,'R',
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
				theParents[rightContig].insert(leftContig);
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
			for(vector<int>::iterator i=f_end_index[wordId(wordAtPosition.c_str())].begin();i!=f_end_index[wordId(wordAtPosition.c_str())].end();i++){
				int contigNumber=*i;
				string contigSequence=contigs[contigNumber]->getSeq();
				Hit myHit(contigNumber,contigSequence.length()-wordSize,'F',
					readNumber,readPosition,'R');
				reverseHits.push_back(myHit);
			}
			// f_start
			for(vector<int>::iterator i=f_start_index[wordId(wordAtPosition.c_str())].begin();i!=f_start_index[wordId(wordAtPosition.c_str())].end();i++){
				int contigNumber=*i;
				string contigSequence=contigs[contigNumber]->getSeq();
				Hit myHit(contigNumber,0,'F',
					readNumber,readPosition,'R');
				reverseHits.push_back(myHit);
			}
			// r_end
			for(vector<int>::iterator i=r_end_index[wordId(wordAtPosition.c_str())].begin();i!=r_end_index[wordId(wordAtPosition.c_str())].end();i++){
				int contigNumber=*i;
				string contigSequence=contigs[contigNumber]->getSeq();
				Hit myHit(contigNumber,contigSequence.length()-wordSize,'R',
					readNumber,readPosition,'R');
				reverseHits.push_back(myHit);
			}
			// r_start
			for(vector<int>::iterator i=r_start_index[wordId(wordAtPosition.c_str())].begin();i!=r_start_index[wordId(wordAtPosition.c_str())].end();i++){
				int contigNumber=*i;
				string contigSequence=contigs[contigNumber]->getSeq();
				Hit myHit(contigNumber,0,'R',
					readNumber,readPosition,'R');
				reverseHits.push_back(myHit);
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
				theParents[rightContig].insert(leftContig);
				edgeAnnotations[leftContig][rightContig]=hitPair;
				//outputStream<<">"<<theReads[readNumber]->getId()<<endl;
				//outputStream<<theReads[readNumber]->getSeq()<<endl;
			}
		}
	}
	//outputStream.close();

	ofstream fGraphViz("Graphviz.txt");

		

	fGraphViz<<"digraph{"<<endl;
	
	string outputFile=argv[4];
	ofstream fStream(outputFile.c_str());
	set<int> _heads; // where to start.
	set<int> _mixed; // mixed vertices
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
			while(currentContig!=-1){
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

