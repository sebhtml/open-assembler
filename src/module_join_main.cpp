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
		cout<<"module_join <wordSize> <contigs> <reads>"<<endl;
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
	map<VERTEX_TYPE,vector<int> >contigHeadsReverse;
	map<VERTEX_TYPE,vector<int> >contigTailsReverse;


	for(int contigNumber=0;contigNumber<contigs.size();contigNumber++){
		string sequenceOfContig=contigs.at(contigNumber)->getSeq();
		string contigHeadString=sequenceOfContig.substr(0,wordSize);
		string contigTailString=sequenceOfContig.substr(sequenceOfContig.length()-1-wordSize,wordSize);
		VERTEX_TYPE headInInt=DeBruijnAssembler::wordId(contigHeadString.c_str());
		VERTEX_TYPE tailInInt=DeBruijnAssembler::wordId(contigTailString.c_str());
		contigHeads[headInInt].push_back(contigNumber);
		contigTails[tailInInt].push_back(contigNumber);
		string contigHeadStringReverse=DeBruijnAssembler::reverseComplement(contigHeadString.c_str());
		VERTEX_TYPE headInIntReverse=DeBruijnAssembler::wordId(contigHeadStringReverse.c_str());
		string contigTailStringReverse=DeBruijnAssembler::reverseComplement(contigTailString.c_str());
		VERTEX_TYPE tailInIntReverse=DeBruijnAssembler::wordId(contigTailStringReverse.c_str());
		contigHeadsReverse[headInIntReverse].push_back(contigNumber);
		contigTailsReverse[tailInIntReverse].push_back(contigNumber);
	}
	
	cout<<endl;
	cout<<"Done indexing contigs"<<endl;


	map<int,map<int,vector<int> > > forwardForward;

	bool done=false;
	while(done==false){
		vector<Read*>stepForwardTailForwardHead_contigs;
		set<int>contigStat;
		for(int readNumber=0;readNumber<reads.size();readNumber++){
			if(readNumber%10000==0){
				cout<<readNumber<<" / "<<reads.size()<<endl;
			}
			//cout<<readNumber<<endl;
			vector<int> positionsHeads;
			vector<int> positionsTails;
			vector<int> positionsHeadsReverse;
			vector<int> positionsTailsReverse;
			string readSequence=reads[readNumber]->getSeq();
			for(int readPosition=0;readPosition<readSequence.length();readPosition++){
				string wordSequence=readSequence.substr(readPosition,wordSize);
				if(wordSequence.length()!=wordSize)
					continue;
				VERTEX_TYPE wordInInt=DeBruijnAssembler::wordId(wordSequence.c_str());
				
				if(contigHeads.count(wordInInt)>0){
					positionsHeads.push_back(readPosition);
				}
				if(contigTails.count(wordInInt)>0){
					positionsTails.push_back(readPosition);
				}
				if(contigHeadsReverse.count(wordInInt)>0){
					positionsHeadsReverse.push_back(readPosition);
				}
				if(contigTailsReverse.count(wordInInt)>0){
					positionsTailsReverse.push_back(readPosition);
				}
			}
			if(positionsHeads.size()>0&&positionsTails.size()>0){
				cout<<"Matching ForwardTail-ForwardHead "<<positionsTails.size()<<" "<<positionsHeads.size()<<endl;
				if(positionsHeads.size()==1&&positionsTails.size()==1){
					VERTEX_TYPE head=DeBruijnAssembler::wordId(readSequence.substr(positionsHeads[0],wordSize).c_str());
					VERTEX_TYPE tail=DeBruijnAssembler::wordId(readSequence.substr(positionsTails[0],wordSize).c_str());
					vector<int> headContigs=contigHeads[head];
					vector<int> tailContigs=contigTails[tail];
					/*
					cout<<headContigs.size()<<" Heads"<<endl;
				cout<<headContigs[0]<<endl;
				cout<<tailContigs.size()<<" Tails"<<endl;
				cout<<tailContigs[0]<<endl;
				*/
					if(headContigs.size()==1&&tailContigs.size()==1&&contigStat.count(tailContigs[0])==0&&
						contigStat.count(headContigs[0])==0&&
							tailContigs[0]!=headContigs[0]){
						Read*firstContig=contigs[tailContigs[0]];
						Read*lastContig=contigs[headContigs[0]];
						cout<<"Read "<<reads[readNumber]->getId()<<" suggests that "<<firstContig->getId()<<" and "<<lastContig->getId()<<" should be linked"<<endl;
						cout<<"read"<<endl;
						cout<<readSequence<<endl;
						cout<<"tail"<<endl;
						string firstSequence=firstContig->getSeq();
						cout<<firstSequence.substr(firstSequence.length()-400-1,400)<<endl;
						cout<<"head"<<endl;
						string lastSequence=lastContig->getSeq();
						cout<<lastSequence.substr(0,400)<<endl;
						forwardForward[tailContigs[0]][headContigs[0]].push_back(readNumber);
						
						// 2 cases:
					// 1) there is a gap (>=0)
					// 2) there is an overlap
						int tailPositionInRead=positionsTails[0];
						int headPositionInRead=positionsHeads[0];
						int difference=headPositionInRead-tailPositionInRead;
						cout<<"Difference: "<<difference<<endl;
						int leftHits=0;
						int rightHits=0;
						int i=0;
						while(i+headPositionInRead<readSequence.length()){
							if(lastSequence[i]==readSequence[headPositionInRead+i])
								rightHits++;
							i++;
						}
						double rightPercentage=rightHits/(i+0.0);
						cout<<"rightHits "<<rightHits<<"/"<<rightPercentage<<endl;
						i=0;
						while(tailPositionInRead-i+wordSize>=0){
							//cout<<firstSequence[firstSequence.length()-1-i]<<" "<<readSequence[tailPositionInRead-i+wordSize]<<endl;
							if(firstSequence[firstSequence.length()-1-i]==readSequence[tailPositionInRead-i+wordSize])
								leftHits++;
							i++;
						}
						double leftPercentage=leftHits/(i+0.0);
	
						cout<<"leftHits "<<leftHits<<"/"<<leftPercentage<<endl;
						if(rightHits>=30&&leftHits>=30&&leftPercentage>=0.9&&rightPercentage>=0.9){
							int overlap=wordSize-difference;
							if(overlap<0)
								overlap=0;
							ostringstream newContigSequence;
							newContigSequence<<firstSequence;
							for(int i=tailPositionInRead+wordSize+1;i<headPositionInRead;i++){
								newContigSequence<<readSequence[i];
							}
							newContigSequence<<lastSequence.substr(overlap)<<endl;
							cout<<"overlap "<<overlap<<endl;
							cout<<"merge"<<endl;
							//cout<<newContigSequence.str()<<endl;
							contigStat.insert(tailContigs[0]);
							contigStat.insert(headContigs[0]);
							ostringstream aName;
							aName<<firstContig->getId()<<"_"<<lastContig->getId();
							Read*read=new Read(aName.str().c_str(),newContigSequence.str().c_str());

							stepForwardTailForwardHead_contigs.push_back(read);
						}
				
					}else{
						cout<<"Not 1-1"<<endl;
					}
				}
			}
			if(positionsHeads.size()>0&&positionsTailsReverse.size()>0){
				cout<<"Matching ReverseTail-ForwardHead"<<endl;
			}
			if(positionsHeadsReverse.size()>0&&positionsTails.size()>0){
				cout<<"Matching ForwardTail-ReverseHead"<<endl;
			}
			if(positionsHeadsReverse.size()>0&&positionsTailsReverse.size()>0){
				cout<<"Matching ReverseTail-ReverseHead"<<endl;
			}
		}
		for(int contigNumber=0;contigNumber<contigs.size();contigNumber++){
			if(contigStat.count(contigNumber)==0){
				stepForwardTailForwardHead_contigs.push_back(contigs[contigNumber]);
			}
		}
		cout<<contigs.size()<<" -> "<<stepForwardTailForwardHead_contigs.size()<<endl;
		done=contigs.size()==stepForwardTailForwardHead_contigs.size();
		contigs=stepForwardTailForwardHead_contigs;
	}

/*
	cout<<"Statistics"<<endl;
	cout<<"ForwardTail-ForwardHead"<<endl;
	for(map<int,map<int,vector<int> > >::iterator i=forwardForward.begin();i!=forwardForward.end();i++){
		for(map<int,vector<int> >::iterator j=i->second.begin();j!=i->second.end();j++){
			cout<<contigs[i->first]->getId()<<" "<<contigs[j->first]->getId()<<" "<<j->second.size()<<endl;
		}
	}
*/


	cout<<reads.size()<<" / "<<reads.size()<<endl;

	cout<<"Done indexing reads"<<endl;
	return 0;
}

