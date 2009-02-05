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

#define wordSize 31
#define m_DEBUG false

using namespace std;

void applyColor(int i,int color,map<int,int>*contigToColor,map<int,set<int> >*contigs_graph){
	if((*contigToColor)[i]==color)
		return;
	(*contigToColor)[i]=color;
	set<int> children=(*contigs_graph)[i];
	for(set<int>::iterator j=children.begin();j!=children.end();j++){
		applyColor(*j,color,contigToColor,contigs_graph);
	}
}

/*
				||||
			-----------------------------------------------> Query
			<-----------------------------------------------
	
					     |||||
	----------------------------------------------->
	<----------------------------------------------- Subject



*/
vector<Read*> reverseOverlapHeadToHead(vector<Read*> contigSequences){
	cout<<"reverseOverlapHeadToHead"<<endl;
	int offset=50;
	bool reducing=true;
	while(reducing){
		map<string,vector<int> > index;
		for(int i=0;i<contigSequences.size();i++){
			if(i%1000==0){
				cout<<i<<" / "<<contigSequences.size()<<endl;
			}
			string sequence=contigSequences[i]->getSeq();
			if(sequence.length()<100)
				continue;
			string revWord=DeBruijnAssembler::reverseComplement(sequence.substr(offset,wordSize));
			index[revWord].push_back(i);
		}
		cout<<contigSequences.size()<<" / "<<contigSequences.size()<<endl;

		vector<Read*> newContigs;
		set<int> contigUsage;
		for(int i=0;i<contigSequences.size();i++){
			if(contigUsage.count(i)!=0)
				continue;
			string sequence=contigSequences[i]->getSeq();
			if(sequence.length()<100)
				continue;
			vector<int> hits;
			for(int j=0;j<sequence.length();j++){
				string word=sequence.substr(j,wordSize);
				if(word.length()<wordSize)
					continue;
				vector<int> localHits=index[word];
				for(int k=0;k<localHits.size();k++){
					if(localHits[k]==i)
						continue;
					if(contigUsage.count(localHits[k])!=0)
						continue;
					hits.push_back(localHits[k]);
				}
			}
			if(m_DEBUG)
				cout<<hits.size()<<" hits"<<endl;

			for(int j=0;j<hits.size();j++){
				string QuerySequence=sequence;
				string subjectSequenceRaw=contigSequences[hits[j]]->getSeq();
				string SubjectSequenceRev=DeBruijnAssembler::reverseComplement(subjectSequenceRaw);
				string startingSeed=QuerySequence.substr(offset,wordSize);
				string endingSeed=DeBruijnAssembler::reverseComplement(subjectSequenceRaw.substr(offset,wordSize));
				int AlignmentStartInQuery=0;
				int AlignmentStartInSubject=0;
				int AlignmentEndInQuery=0;
				int AlignmentEndInSubject=0;
				while(AlignmentStartInQuery<QuerySequence.length()&&QuerySequence.substr(AlignmentStartInQuery,wordSize)!=startingSeed)
					AlignmentStartInQuery++;
				while(AlignmentEndInQuery<QuerySequence.length()&&QuerySequence.substr(AlignmentEndInQuery,wordSize)!=endingSeed)
					AlignmentEndInQuery++;
				while(AlignmentStartInSubject<SubjectSequenceRev.length()&&SubjectSequenceRev.substr(AlignmentStartInSubject,wordSize)!=startingSeed)
					AlignmentStartInSubject++;
				while(AlignmentEndInSubject<SubjectSequenceRev.length()&&SubjectSequenceRev.substr(AlignmentEndInSubject,wordSize)!=endingSeed)
					AlignmentEndInSubject++;


/*
 			------------------------------------->
	------------------------------>

*/
				int found=0;
				set<string> currentIndex;
				if(AlignmentEndInQuery>AlignmentStartInQuery&&AlignmentEndInSubject>AlignmentEndInSubject){
					for(int iteratorI=AlignmentStartInQuery;iteratorI<=AlignmentEndInQuery;iteratorI++){
						currentIndex.insert(QuerySequence.substr(iteratorI,wordSize));
					}
					for(int iteratorI=AlignmentStartInSubject;iteratorI<=AlignmentEndInSubject;iteratorI++){
						if(currentIndex.count(SubjectSequenceRev.substr(iteratorI,wordSize))>0)
							found++;
					}
				}
				bool isTheSame=found/(AlignmentEndInSubject-AlignmentStartInSubject+0.0)>0.90;
				//cout<<AlignmentEndInQuery-AlignmentStartInQuery<<endl;
				//cout<<AlignmentEndInSubject-AlignmentStartInSubject<<endl;
				if(isTheSame&&contigUsage.count(i)==0&&contigUsage.count(hits[j])==0){
					//cout<<"Alignment! "<<AlignmentEndInQuery-AlignmentStartInQuery<<endl;
					ostringstream newName;
					newName<<"Contig_"<<contigSequences[hits[j]]->getId()<<"_R"<<contigSequences[i]->getId()<<"_F_";
					
					string sequence=SubjectSequenceRev.substr(0,AlignmentStartInSubject)+QuerySequence.substr(AlignmentStartInQuery);
					Read*read=new Read(newName.str().c_str(),sequence.c_str());
					newContigs.push_back(read);
					contigUsage.insert(i);
					contigUsage.insert(hits[j]);
				}
			}
		}

		for(int i=0;i<contigSequences.size();i++){
			if(contigUsage.count(i)==0)
				newContigs.push_back(contigSequences[i]);
		}
		reducing=contigSequences.size()>newContigs.size();
		contigSequences=newContigs;
	}
	return contigSequences;
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




vector<Read*> mergeForward(vector<Read*> contigSequences){
	map<int,set<int> > contigs_graph;
	cout<<contigSequences.size()<<" contigs"<<endl;
	map<string,vector<int> > indexOfWords;
	cout<<"[Forward to Forward]"<<endl;
	cout<<"Indexing"<<endl;
	for(int i=0;i<contigSequences.size();i++){
		contigs_graph[i].size();// insert a node
		if(i%10000==0)
			cout<<i+1<<" / "<<contigSequences.size()<<endl;
		string sequence=contigSequences[i]->getSeq();
		if(sequence.length()<wordSize)
			continue;
		//cout<<sequence.length()<<endl;
		string word=sequence.substr(sequence.length()-wordSize,wordSize);
		indexOfWords[word].push_back(i);
	}
	cout<<contigSequences.size()<<" / "<<contigSequences.size()<<endl;
	cout<<"Building a graph of contigs."<<endl;
	for(int i=0;i<contigSequences.size();i++){
		if(i%1000==0)
			cout<<i+1<<" / "<<contigSequences.size()<<endl;

		bool used=false;
		string sequence=contigSequences[i]->getSeq();
		//cout<<sequence.length()<<endl;

		if(sequence.length()<wordSize)
			continue;
		string word=sequence.substr(sequence.length()-wordSize,wordSize);
		vector<int>otherContigs=indexOfWords[word];
		if(m_DEBUG){
			cout<<"Forward hits: "<<otherContigs.size()<<endl;
			cout<<otherContigs.size()<< "hits"<<endl;
		}
		// forward hits
		for(vector<int>::iterator matchContig=otherContigs.begin();matchContig!=otherContigs.end();matchContig++){
			string shortContig=contigSequences[i]->getSeq();
			string longContig=contigSequences[*matchContig]->getSeq();
			if(i!=*matchContig&&
		shortContig.length()<=longContig.length()&&
		used==false
			){
			
				int ok=0;
				for(int k=0;k<shortContig.length();k++){
					if(shortContig[shortContig.length()-1-k]==longContig[longContig.length()-1-k])
						ok++;
					else
						break;
				}
				if(m_DEBUG)
					cout<<ok<<" "<<shortContig.length()<<endl;
				bool theSame=ok>(shortContig.length()/2);
				if(theSame){
					if(m_DEBUG)
						cout<<"Forward the same"<<endl;
					contigs_graph[i].insert(*matchContig);
					contigs_graph[*matchContig].insert(i);
					//used=true;
				}
			}
		}
	}

	cout<<contigSequences.size()<<" / "<<contigSequences.size()<<endl;
	cout<<"The graph is ready."<<endl;

	map<int,int> contigToColor;
	
	for(int i=0;i<contigSequences.size();i++){
		contigToColor[i]=-1;
	}
	int color=1;
	
	for(int i=0;i<contigSequences.size();i++){
		if(i%1000==0)
			cout<<i+1<<" / "<<contigSequences.size()<<endl;
		applyColor(i,color,&contigToColor,&contigs_graph);
		color++;
	}
	cout<<contigSequences.size()<<" / "<<contigSequences.size()<<endl;

/*
	ofstream graphFile("graph.graphviz");
	graphFile<<"digraph G{"<<endl;
	for(map<int,int>::iterator i=contigToColor.begin();i!=contigToColor.end();i++){
		//graphFi
	}
	for(map<int,set<int> >::iterator i=contigs_graph.begin();i!=contigs_graph.end();i++){
		for(set<int>::iterator j=i->second.begin();j!=i->second.end();j++){
			graphFile<<i->first<<" -> "<<*j<<endl;
		}
	}
	graphFile<<"}"<<endl;
	graphFile.close();
*/
	map<int,vector<int> > colorToContigs;
	for(int i=0;i<contigSequences.size();i++){
		colorToContigs[contigToColor[i]].push_back(i);
	}
	
	vector<Read*> theBestContigs;
	int colors=0;
	for(map<int,vector<int> >::iterator i=colorToContigs.begin();i!=colorToContigs.end();i++){
		vector<int> contigs=i->second;
		if(colors%10000==0)
			cout<<"Progress: "<<colors<<endl;
		colors++;
		//cout<<"Color: "<<i->first<<endl;
		//cout<<contigs.size()<<" contigs: ";
		
		int best=contigs[0];
		for(int j=0;j<contigs.size();j++){
			//cout<<contigsWithNames->at(contigs[j])->getId()<<" ";
			//cout<<contigs[j]<<" ";
			string a=contigSequences[contigs[j]]->getSeq();
			string b=contigSequences[best]->getSeq();
			if(a.length()>b.length()){
				best=contigs[j];
			}
		}
		//cout<<endl;
		theBestContigs.push_back(contigSequences[best]);
	}
	return theBestContigs;
}


vector<Read*> reverseMerge(vector<Read*> contigSequences){
	map<int,set<int> > contigs_graph;
	cout<<contigSequences.size()<<" contigs"<<endl;
	map<VERTEX_TYPE,vector<int> > indexOfRevWords;
	cout<<"Reverse to Reverse now."<<endl;
	cout<<"Indexing"<<endl;
	for(int i=0;i<contigSequences.size();i++){
		contigs_graph[i].size();// insert a node
		if(i%10000==0)
			cout<<i+1<<" / "<<contigSequences.size()<<endl;
		string sequence=contigSequences[i]->getSeq();
		for(int j=0;j<sequence.length();j+=wordSize){
			string word=sequence.substr(j,wordSize);
			if(word.length()!=wordSize)
				continue;
			string revWord=DeBruijnAssembler::reverseComplement(word);
			indexOfRevWords[DeBruijnAssembler::wordId(revWord.c_str())].push_back(i);
		}
	}
	int maxNotFound=2*wordSize;
	maxNotFound=0;
	cout<<contigSequences.size()<<" / "<<contigSequences.size()<<endl;
	cout<<"Building a graph of contigs."<<endl;
	for(int i=0;i<contigSequences.size();i++){
		if(i%100==0)
			cout<<i+1<<" / "<<contigSequences.size()<<endl;

		bool used=false;
		set<int>otherRevContigs;
		string sequence=contigSequences[i]->getSeq();
		for(int j=sequence.length()/2-3*wordSize;j<sequence.length()/2+3*wordSize;j+=1){
			string word=sequence.substr(j,wordSize);
			if(word.length()!=wordSize)
				continue;

			vector<int> otherRevContigs2=indexOfRevWords[DeBruijnAssembler::wordId(word.c_str())];
			for(int k=0;k<otherRevContigs2.size();k++)
				otherRevContigs.insert(otherRevContigs2[k]);
		}


		if(m_DEBUG){
			cout<<"Reverse hits: "<<otherRevContigs.size()<<endl;
		}
		// reverse complement hits
		for(set<int>::iterator matchContig=otherRevContigs.begin();matchContig!=otherRevContigs.end();matchContig++){
			string sequenceSmall=contigSequences[i]->getSeq();
			string sequenceLong=contigSequences[*matchContig]->getSeq();
			if(i!=*matchContig&&
		sequenceSmall.length()<=sequenceLong.length()&&
				used==false
			){
				string sequenceLong_Rev=DeBruijnAssembler::reverseComplement(sequenceLong);
				string shortContig=sequenceSmall;
				string longContig=sequenceLong_Rev;
				set<VERTEX_TYPE> longMap;
				for(int o=0;o<longContig.length();o++){
					string word=longContig.substr(o,wordSize);
					if(word.length()!=wordSize)
						continue;
					longMap.insert(DeBruijnAssembler::wordId(word.c_str()));
				}

				int notFound=0;
				for(int o=0;o<shortContig.length();o++){
					string word=shortContig.substr(o,wordSize);
					if(word.length()!=wordSize)
						continue;
					if(longMap.count(DeBruijnAssembler::wordId(word.c_str()))==0){
						//cout<<"P "<<o<<endl;
						notFound++;
					}
				}
				bool theSame=(notFound/(wordSize+0.0))/shortContig.length()<0.01;
				cout<<"Not found "<<notFound/wordSize<<", "<<shortContig.length()<<endl;
				if(theSame){
					contigs_graph[i].insert(*matchContig);
					contigs_graph[*matchContig].insert(i);
					//used=true;
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

	ofstream graphFile("graph.graphviz");
	graphFile<<"digraph G{"<<endl;
	for(map<int,int>::iterator i=contigToColor.begin();i!=contigToColor.end();i++){
		//graphFi
	}
	for(map<int,set<int> >::iterator i=contigs_graph.begin();i!=contigs_graph.end();i++){
		for(set<int>::iterator j=i->second.begin();j!=i->second.end();j++){
			graphFile<<i->first<<" -> "<<*j<<endl;
		}
	}
	graphFile<<"}"<<endl;
	graphFile.close();

	map<int,vector<int> > colorToContigs;
	for(int i=0;i<contigSequences.size();i++){
		colorToContigs[contigToColor[i]].push_back(i);
	}
	
	vector<Read*> theBestContigs;
	for(map<int,vector<int> >::iterator i=colorToContigs.begin();i!=colorToContigs.end();i++){
		vector<int> contigs=i->second;
		cout<<"Color: "<<i->first<<endl;
		cout<<contigs.size()<<" contigs: ";
		
		int best=contigs[0];
		for(int j=0;j<contigs.size();j++){
			//cout<<contigsWithNames->at(contigs[j])->getId()<<" ";
			cout<<contigs[j]<<" ";
			string sequenceA=contigSequences[contigs[j]]->getSeq();
			string sequenceB=contigSequences[best]->getSeq();
			if(sequenceA.length()>sequenceB.length()){
				best=contigs[j];
			}
		}
		cout<<endl;
		theBestContigs.push_back(contigSequences[best]);
	}
	return theBestContigs;
}

int main(int argc,char*argv[]){
	DeBruijnAssembler::CommonHeader(&cout);
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

	//for(int i=0;i<contigs.size();i++)
		//cout<<strlen(contigs[i]->getSeq())<<endl;
	vector<Read*> theBestContigs=(reverseOverlapHeadToHead(reverseMerge(mergeForward(contigs))));

	ofstream f(outputFile.c_str());
	for(int j=0;j<theBestContigs.size();j++){
		string sequence=theBestContigs[j]->getSeq();
		f<<">"<<theBestContigs.at(j)->getId()<<" "<<sequence.length()<<" nucleotides"<<endl;
		cout<<">"<<theBestContigs.at(j)->getId()<<" "<<sequence.length()<<" nucleotides"<<endl;
		int columns=70;
		int k=0;
		while(k<sequence.length()){
			f<<sequence.substr(k,columns)<<endl;
			k+=columns;
		}
	}
	f.close();
	return 0;
}

