/*
	454 De Novo Assembler
    Copyright (C) 2008 Sébastien Boisvert

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

//#include<multiset>
#include<cstdlib>
#include"Logger.h"
#include<iostream>
#include"DeBruijnAssembler.h"
#include<fstream>
#include<vector>

#define endl '\n'
using namespace std;

// TODO: add paired end reads ?
// with 454 it is 3k

// TODO expand repeats using reads

DeBruijnAssembler::DeBruijnAssembler(Logger*m_cout){
	this->m_cout=m_cout;
	m_COLOR_NOT_ASSEMBLED=-2;
	m_COLOR_IN_PROGRESS_ASSEMBLER=-3;
	m_COLOR_IN_PROGRESS_SOLVER=-5;
	m_COLOR_DISCARDED=-99;
	m_contig_id=1;
	m_NUCLEOTIDE_A=0;
	m_NUCLEOTIDE_T=1;
	m_NUCLEOTIDE_C=2;
	m_NUCLEOTIDE_G=3;
}



void DeBruijnAssembler::setMinimumQuality(int q){
	m_minimumQuality=q;
}

void DeBruijnAssembler::setWordSize(int k){
	m_wordSize=k;
}

void DeBruijnAssembler::setMinimumCoverage(int m){
	m_minimumCoverage=m;
}


char DeBruijnAssembler::complement(char a){
	Logger&m_cout=*(this->m_cout);
	if(a=='A')
		return 'T';
	if(a=='T')
		return 'A';
	if(a=='C')
		return 'G';
	if(a=='G')
		return 'C';
	m_cout<<"Error: symbol not allowed ("<<a<<")"<<endl;
	exit(0);
	return 'E';
}

string DeBruijnAssembler::reverseComplement(string a){
	ostringstream i;
	for(int p=a.length()-1;p>=0;p--)
		i<< complement(a[p]);
	return i.str();
}

// return (k+1)-mers from read
// check that the k+1 nucleotides at position i are good quality and that they are not the N
vector<unsigned long int> DeBruijnAssembler::getHighQualityMers(int readId){
	vector<unsigned long int> highQualityMers;
	Read*read=m_reads->at(readId);
	char*sequence=read->getSeq();
	char*quality=read->getQual();
	int l=strlen(sequence);

	char*wordFoward=(char*)malloc(m_wordSize+2);
	char*wordReverse=(char*)malloc(m_wordSize+2);
	// skip tag TCAG... 
	for(int j=0;j<=l-m_wordSize+1;j++){
		bool ok=true;
		for(int p=0;p<=m_wordSize+1;p++){
			if(sequence[j+p]=='N'||quality[j+p]<m_minimumQuality){
				ok=false;
				break;
			}
		}
		if(!ok)
			continue;
		for(int p=0;p<m_wordSize+1;p++){
			wordFoward[p]=sequence[j+p];
			wordReverse[p]=complement(sequence[j+m_wordSize-p]);
		}
		wordFoward[m_wordSize+1]='\0';
		wordReverse[m_wordSize+1]='\0';
		highQualityMers.push_back(wordId(wordFoward));
		highQualityMers.push_back(wordId(wordReverse));
	}

	free(wordFoward);
	free(wordReverse);
	return highQualityMers;
}

void DeBruijnAssembler::buildGraph(vector<Read*>*reads){
	Logger&m_cout=*(this->m_cout);
	m_reads=reads;


	m_cout<<endl;
	m_cout<<"Collecting high-quality mers from reads"<<endl;
	map<unsigned long int,int> words;
	for(int i=0;i<(int)reads->size();i++){
		if(i%10000==0){
			m_cout<<"Reads: "<<i<<" / "<<reads->size()<<endl;
		}
		vector<unsigned long int>highQualityMers=getHighQualityMers(i);
		for(vector<unsigned long int>::iterator iteratorMer=highQualityMers.begin();iteratorMer!=highQualityMers.end();iteratorMer++){
			//m_cout<<endl;
			//m_cout<<wordFoward<<endl;
			//m_cout<<wordReverse<<endl;
			//m_cout<<idToWord(wordId(wordFoward))<<endl;
			words[*iteratorMer]++;
			if((int)words.size()%1000000==0&&(int)words.size()!=m_last_vertices_size){
				m_cout<<"High-quality mers: "<<words.size()<<endl;
				m_last_vertices_size=words.size();
			}
		}
	}	


	m_cout<<"High-quality mers: "<<words.size()<<endl;
	m_cout<<"Reads: "<<reads->size()<<" / "<<reads->size()<<endl;
	int processed=0;
	int solid=0;




	m_cout<<endl;
	m_cout<<"Collecting solid mers"<<endl;
	set<unsigned long int> solidMers;
	for(map<unsigned long int,int>::iterator i=words.begin();i!=words.end();i++){
		processed++;
		if(i->second>=m_minimumCoverage){
			solid++;
			unsigned long int w=i->first;
			solidMers.insert(w);
		}
		if(processed%1000==0){
			//m_cout<<"Processed words: "<<processed<<" / "<<words.size()<<endl;
		}
		if(solid%1000==0){
			//m_cout<<"Solid words: "<<solid<<" / "<<words.size()<<endl;
		}
	}
	m_cout<<"Processed mers: "<<processed<<" / "<<words.size()<<endl;
	m_cout<<"Solid mers: "<<solid<<" / "<<words.size()<<" ---> "<<((solid+0.0)/words.size()+0.0)*100.0<<"%"<<endl;
	m_cout<<" (this should be roughly twice the genome size)"<<endl;
	m_cout<<"Not-solid mers: "<<processed-solid<<" / "<<words.size()<<" ---> "<<(processed-solid+0.0)/words.size()*100.0<<"%"<<endl;

	m_cout<<endl;
	words.clear();
	//m_cout<<"Clearing mers"<<endl;
	//m_cout<<endl;




	int solid_mers_processed=0;
	m_cout<<"Indexing solid mers in reads."<<endl; // <-------
	for(int readId=0;readId<(int)m_reads->size();readId++){
		if(readId%10000==0)
			m_cout<<readId<<" / "<<m_reads->size()<<endl;
		vector<unsigned long int> highQualityMers=getHighQualityMers(readId);
		for(vector<unsigned long int>::iterator merIterator=highQualityMers.begin();
			merIterator!=highQualityMers.end();merIterator++){
			if(solidMers.count(*merIterator)>0)
				m_mer_to_read_table[*merIterator].push_back(readId);
		}
	}

	m_cout<<m_reads->size()<<" / "<<m_reads->size()<<endl;




	m_cout<<"Building de Bruijn graph with solid mers"<<endl;
	for(set<unsigned long int>::iterator i=solidMers.begin();i!=solidMers.end();i++){
		solid_mers_processed++;
		//m_cout<<endl;
		if(solid_mers_processed%100000==0){
			m_cout<<"Building graph: "<<solid_mers_processed<<" / "<<solidMers.size()<<endl;
		}
		string wholeWord=idToWord(*i,m_wordSize+1);
		//m_cout<<"w   "<<wholeWord<<endl;
		//unsigned long long int prefix=(w*4);
		//m_cout<<sizeof(unsigned long long int)<<endl;
		//printBinary(w,m_wordSize+1);
		
		//m_cout<<"w<<2"<<endl;
		//printBinary(prefix,m_wordSize+1);
		//prefix=prefix/4;
		//m_cout<<"w<<2>>2"<<endl;
		//printBinary(prefix,m_wordSize+1);
		/*
-----------prefix-------
		-----------suffix-------
		*/
		//m_cout<<endl;
		//m_cout<<i->second<<" X "<<endl;
		unsigned long int prefix=wordId(wholeWord.substr(0,m_wordSize).c_str());
		unsigned long int suffix=wordId(wholeWord.substr(1,m_wordSize).c_str());
		//m_cout<<"pre "<<prefix<<endl;
		//printBinary(prefix,m_wordSize+1);
		//m_cout<<"suf  "<<suffix<<endl;
		//m_cout<<"  ";
		//printBinary(suffix,m_wordSize+1);
		/*
		m_cout<<idToWord(prefix,m_wordSize)<<endl;
		m_cout<<idToWord(suffix,m_wordSize)<<endl;
		m_cout<<idToWord(prefixRev,m_wordSize)<<endl;
		m_cout<<idToWord(suffixRev,m_wordSize)<<endl;
		*/
		/*
		unsigned long int prefixRev=wordId(reverseComplement(wholeWord.substr(0,m_wordSize)).c_str());
		unsigned long int suffixRev=wordId(reverseComplement(wholeWord.substr(1,m_wordSize)).c_str());
		TODO: how to not duplicate everything?
		if(m_graph_mers.count(prefixRev)>0||m_graph_mers.count(suffixRev)>0){
			m_graph_mers.insert(prefixRev);
			m_graph_mers.insert(suffixRev);
			m_graph[suffixRev].insert(prefixRev);
		}else{
			m_graph_mers.insert(prefix);
			m_graph_mers.insert(suffix);
			m_graph[prefix].insert(suffix);
		}
		*/
		m_graph_mers.insert(prefix);
		m_graph_mers.insert(suffix);
		m_graph[prefix].insert(suffix);
	}

	m_cout<<"Building graph: "<<solidMers.size()<<" / "<<solidMers.size()<<endl;
	m_cout<<"the graph is ready"<<endl;
	m_cout<<"Clearing solid mers"<<endl;
	solidMers.clear();
}


void DeBruijnAssembler::perform_Assembly(unsigned long int vertex){
	Logger&m_cout=*(this->m_cout);
	//m_cout<<"[perform_Assembly] "<<idToWord(vertex,m_wordSize)<<endl;
	if(m_colors[vertex]>=1){
		//m_cout<<"Already assembled. "<<m_colors[vertex]<<endl;
		return;
	}
	if(m_colors[vertex]==m_COLOR_IN_PROGRESS_ASSEMBLER){// already done
		//m_cout<<"Already in progress."<<endl;
		return;
	}
	m_colors[vertex]=m_COLOR_IN_PROGRESS_ASSEMBLER;
	vector<unsigned long int> nextVertices=getNeighbours(vertex);
	if(nextVertices.size()==1){
		//m_cout<<"One next, going in."<<endl;
		unsigned long int nextVertex=nextVertices[0];
		if(m_colors[nextVertex]==m_COLOR_IN_PROGRESS_ASSEMBLER){
			m_cout<<"REPEAT detected?"<<endl;
			m_repeat_vertices.insert(nextVertex);
			m_colors[vertex]=m_contig_id;
			m_contig_paths[m_colors[vertex]].push_front(vertex);
			m_contig_id++;
		}else{
			perform_Assembly(nextVertex);
			m_colors[vertex]=m_colors[nextVertex];
			m_contig_paths[m_colors[vertex]].push_front(vertex);
		}
	}else{
		//m_cout<<"Leaf, assigning contig "<<m_contig_id<<endl;
		m_colors[vertex]=m_contig_id;
		m_contig_paths[m_colors[vertex]].push_front(vertex);
		m_contig_id++;
	}
	if(m_colors[vertex]==m_COLOR_IN_PROGRESS_ASSEMBLER){
		m_cout<<"Error: color not set "<<idToWord(vertex,m_wordSize)<<endl;
	}
}

void DeBruijnAssembler::run_Assembler(){
	Logger&m_cout=*(this->m_cout);
	m_cout<<"setting colors"<<endl;
	for(set<unsigned long int>::iterator i=m_graph_mers.begin();i!=m_graph_mers.end();i++){
		m_colors[*i]=m_COLOR_NOT_ASSEMBLED;
	}

	m_cout<<"fixing graph"<<endl;
	int progressIndicator=0;
	int recursionDepth=1000;
	for(set<unsigned long int>::iterator i=m_graph_mers.begin();i!=m_graph_mers.end();i++){
		if(progressIndicator%100000==0)
			m_cout<<progressIndicator<<" / "<<m_graph_mers.size()<<endl;
		solveMultiPath(*i,recursionDepth);
		progressIndicator++;
	}
	m_cout<<m_graph_mers.size()<<" / "<<m_graph_mers.size()<<endl;


	progressIndicator=0;
	m_cout<<"Performing assembly"<<endl;
	for(set<unsigned long int>::iterator i=m_graph_mers.begin();i!=m_graph_mers.end();i++){
		if(progressIndicator%100000==0)
			m_cout<<progressIndicator<<" / "<<m_graph_mers.size()<<endl;
		perform_Assembly(*i);
		progressIndicator++;
	}
	m_cout<<m_graph_mers.size()<<" / "<<m_graph_mers.size()<<endl;
	
}

void DeBruijnAssembler::setMinimumContigSize(int minimumContigSize){
	m_minimumContigSize=minimumContigSize;
}

void DeBruijnAssembler::outputContigs(ofstream*output){
	Logger&m_cout=*(this->m_cout);
	// write contigs
	// TODO: don't duplicate each strand..

	int columns=60;
	vector<string> contigs;
	for(map<int,list<unsigned long int> >::iterator iteratorContigFirstVertex=m_contig_paths.begin();
		iteratorContigFirstVertex!=m_contig_paths.end();iteratorContigFirstVertex++){
		if((int)iteratorContigFirstVertex->second.size()+m_wordSize<m_minimumContigSize)
			continue;
		ostringstream contigSequence;
		for(list<unsigned long int>::iterator i=iteratorContigFirstVertex->second.begin();
			i!=iteratorContigFirstVertex->second.end();i++){
			if(i==iteratorContigFirstVertex->second.begin()){
				contigSequence<< idToWord(*i,m_wordSize);
			}else{
				contigSequence<<getLastSymbol(*i);
			}
		}
		
		contigs.push_back(contigSequence.str());
	}


	// remove duplicates
	set<int> duplicatedContig;

	for(int contigIndex=0;contigIndex<contigs.size();contigIndex++){
		//cout<<contigIndex<<endl;
		if(duplicatedContig.count(contigIndex)>0)// if contig is not already tagged as duplicated
			continue;
		for(int otherContigIndex=0;otherContigIndex<contigs.size();otherContigIndex++){
			//cout<<"  "<<otherContigIndex<<endl;
			if(otherContigIndex==contigIndex)// is the same
				continue;
			if(duplicatedContig.count(otherContigIndex)>0)// already tagged as duplicated
				continue;
			if(contigs[contigIndex].length()!=contigs[otherContigIndex].length())
				continue;
			string revComp=reverseComplement(contigs[contigIndex]);
			if(revComp==contigs[otherContigIndex]){
				//m_cout<<"duplicate"<<endl;
				duplicatedContig.insert(otherContigIndex);
			}else{
				//m_cout<<revComp.substr(0,70)<<endl;
				//m_cout<<contigs[otherContigIndex].substr(0,70)<<endl;
			}
		}
	}
	m_cout<<contigs.size()-duplicatedContig.size()<<" contigs ("<<contigs.size()<<" in the graph)"<<endl;	
	//m_cout<<contigs.size()-duplicatedContig.size()<<" after mergin of complement contigs (the graph contains both strands)"<<endl;
	int totalLength=0;
	int contigId=1;
	multiset<int> contigSizes;
	unsigned long int stat_minimumContigSize=9999999999;
	int stat_maximumContigSize=0;
	for(int i=0;i<contigs.size();i++){
		if(duplicatedContig.count(i)>0){
			continue;
		}

		string sequenceDNA=contigs[i];
		totalLength+=sequenceDNA.length();
		contigSizes.insert(sequenceDNA.length());
		if(sequenceDNA.length()<stat_minimumContigSize)
			stat_minimumContigSize=sequenceDNA.length();
		if(sequenceDNA.length()>stat_maximumContigSize)
			stat_maximumContigSize=sequenceDNA.length();

		m_cout<<sequenceDNA.length()<<" nucleotides"<<endl;
		//m_cout<<"contig "<<iteratorContigFirstVertex->first<<endl;
		int j=0;
		(*output)<<">Contig";
		(*output)<<contigId;
		(*output)<<endl;
		contigId++;
		while(j<(int)sequenceDNA.length()){
			(*output)<<sequenceDNA.substr(j,columns);
			(*output)<<endl;
			j+=columns;
		}

	}
	int averageContigLength=totalLength/(contigs.size()-duplicatedContig.size());
	m_cout<<"Asssembly statistics"<<endl;
	m_cout<<"Estimated genome size (total contig bases): "<<totalLength<<endl;
	m_cout<<"Average contig size: "<<averageContigLength<<endl;
	m_cout<<"Smallest contig size: "<<stat_minimumContigSize<<endl;
	m_cout<<"Largest contig size: "<<stat_maximumContigSize<<endl;
	m_cout<<"Contigs: "<<contigs.size()-duplicatedContig.size()<<endl;
	int cumulativeSize=0;
	set<int> stat_nx_done;
	for(multiset<int>::reverse_iterator i=contigSizes.rbegin();i!=contigSizes.rend();i++){
		cumulativeSize+=*i;
		for(int j=50;j<=90;j+=10){
			if(cumulativeSize>=totalLength*j/100.0&&stat_nx_done.count(j)==0){
				stat_nx_done.insert(j);
				cout<<"N"<<j<<" size: "<<*i<<endl;
				break;
			}
		}
	}
	m_cout<<"See http://www.cbcb.umd.edu/research/castats.shtml#N50 for the definition of N50."<<endl;
}



vector<unsigned long int> DeBruijnAssembler::getNeighbours(unsigned long int vertex){
	vector<unsigned long int> neighboursVector;


	for(set<unsigned long int>::iterator i=m_graph[vertex].begin();i!=m_graph[vertex].end();i++){
		unsigned long int currentNeighbour=*i;
		if(m_colors[currentNeighbour]!=m_COLOR_DISCARDED&&currentNeighbour!=vertex){
		//if(m_colors[currentNeighbour]==m_colors[vertex]&&vertex!=currentNeighbour){
			neighboursVector.push_back(currentNeighbour);
		}
	}


	return neighboursVector;
}



// returns a trivial for a vertex
string DeBruijnAssembler::getTrivialPath(unsigned long pathVertex,unsigned long int vertex){
	ostringstream path;
	m_DETECTED_REPEAT=false;
	path<< idToWord(vertex,m_wordSize);
	int maxIterations=10000;
	unsigned long int firstNeighbour=pathVertex;
	path<< getLastSymbol(firstNeighbour);
	set<unsigned long int> verticesInPath;
	verticesInPath.insert(vertex);
	verticesInPath.insert(pathVertex);
	int iterations=0;
	while(iterations<=maxIterations){
		vector<unsigned long int> localNeighbours1=getNeighbours(firstNeighbour);
		//m_cout<<iterations<<endl;
		if(localNeighbours1.size()==1){
			firstNeighbour=localNeighbours1[0];
			if(verticesInPath.count(firstNeighbour)>0){
				m_DETECTED_REPEAT=true;
				m_repeat_vertex=firstNeighbour;
				break;
			}
			path<< getLastSymbol(firstNeighbour);
			verticesInPath.insert(firstNeighbour);
		}else if(localNeighbours1.size()==2){
			// TODO code repeats here
			// if 2 neighbours, one is the repeat, the other is the path
			break;
		}else{
			break;
		}
		iterations++;
	}
	return path.str();
}

// 1. try majority vote if more than 1 neighbour
// 2. if that fails, try to find bubbles
// 3. if that fails, ... work in progress
// TODO: maybe split this function in parts, (not very important)
void DeBruijnAssembler::solveMultiPath(unsigned long int vertex,int depth){
	Logger&m_cout=*(this->m_cout);
	if(depth==0||m_colors[vertex]==m_COLOR_IN_PROGRESS_SOLVER)
		return ;
	if(m_colors[vertex]==m_COLOR_DISCARDED)
		return;

	m_colors[vertex]=m_COLOR_IN_PROGRESS_SOLVER;

	vector<unsigned long int> neighboursVector=getNeighbours(vertex);
	for(vector<unsigned long int>::iterator solverIterator=neighboursVector.begin();
		solverIterator!=neighboursVector.end();solverIterator++){
		if(*solverIterator==vertex){
			m_colors[*solverIterator]=m_COLOR_DISCARDED;
		}
		solveMultiPath(*solverIterator,depth-1);
	}
	// majority vote
	if(neighboursVector.size()>1){
		// TODO: put debug here only.
		// new algo here
		//m_cout<<"new algo"<<endl;
		map<int,map<int,int> > votes;
		int maximumThreshold=20;
		for(int i=0;i<(int)neighboursVector.size();i++){
			for(int threshold=0;threshold<=maximumThreshold;threshold++)
				votes[threshold][i]=0;
			// TODO: do majority vote in binary space here.
			// this will give the same results faster
			string path=getTrivialPath(neighboursVector[i],vertex);
			//m_cout<<"Path "<<path<<endl;
			// count the reads that agree..
			//map<unsigned long int,int> treeForPath;
			// build hash
			map<int,int> readVotes;
			for(int iteratorPath=0;iteratorPath<(int)path.length()-m_wordSize;iteratorPath++){
				unsigned long int kMerInPath=wordId(path.substr(iteratorPath,m_wordSize+1).c_str());
				//treeForPath[kMerInPath]=iteratorPath;
				//m_cout<<m_mer_to_read_table[kMerInPath].size()<<" reads for "<<idToWord(kMerInPath,m_wordSize+1)<<endl;
				for(vector<int>::iterator readItIterator=m_mer_to_read_table[kMerInPath].begin();
					readItIterator!=m_mer_to_read_table[kMerInPath].end();readItIterator++){
					readVotes[*readItIterator]++;
				}
			}
			for(int threshold=0;threshold<=maximumThreshold;threshold++){
				for(map<int,int>::iterator voteIterator=readVotes.begin();voteIterator!=readVotes.end();voteIterator++){
					if(voteIterator->second>=threshold)
						votes[threshold][i]++;
				}
			}
		}
		//m_cout<<"votes "<<endl;

		int winner=-1;
		int maxVotes=0;
		for(int threshold=maximumThreshold;threshold>=0;threshold--){
			for(map<int,int>::iterator voteIterator=votes[threshold].begin();voteIterator!=votes[threshold].end();voteIterator++){	
				//m_cout<<voteIterator->first<<" "<<voteIterator->second<<endl;
				if(voteIterator->second>maxVotes){
					winner=voteIterator->first;
					maxVotes=voteIterator->second;
				}else if(voteIterator->second==maxVotes){
					winner=-1;
				}
			}
			if(winner!=-1){
				// make sure that all other are 0...
				for(map<int,int>::iterator voteIterator=votes[threshold].begin();voteIterator!=votes[threshold].end();voteIterator++){	
					if(voteIterator->first!=winner){
						if(voteIterator->second!=0){
							//m_cout<<" not-0 vote found. fallback"<<endl;
							winner=-1;
							break;
						}
					}
				}
			}
			if(winner!=-1){
				//m_cout<<"Found a winner ;P, threshold="<<threshold<<endl;
				break; // found a winner
			}
		}
		if(winner!=-1){
			//m_cout<<winner<<" wins, has "<<maxVotes<<" other 0."<<endl;
			for(int iteratorOnNext=0;iteratorOnNext<(int)neighboursVector.size();iteratorOnNext++){
				if(iteratorOnNext!=winner){
					m_colors[neighboursVector[iteratorOnNext]]=m_COLOR_DISCARDED;
				}
			}
		}else{
			// TODO: add this as a parameter?
			// bubbles are caused by 454's homopolymers that are cursed by wizards
			int maximumDifferenceInBubbles=3;
			//m_cout<<"No winner..."<<endl;
			bool allTheSameLength=true;
			for(int u=0;u<(int)neighboursVector.size();u++){
				if(getTrivialPath(neighboursVector[u],vertex).length()!=getTrivialPath(neighboursVector[0],vertex).length()){
					allTheSameLength=false;
				}
			}
			if((int)(allTheSameLength&&getTrivialPath(neighboursVector[0],vertex).length())==m_wordSize+1){
				m_cout<<endl;
				m_cout<<"All too short, skipping."<<endl;
				for(int u=0;u<(int)neighboursVector.size();u++){
					m_cout<<getTrivialPath(neighboursVector[u],vertex)<<endl;
				}
			}else if(neighboursVector.size()==2){
				//m_cout<<"2 alternative paths"<<endl;
				string path1=getTrivialPath(neighboursVector[0],vertex);
				bool path1IsRepeat=m_DETECTED_REPEAT;
				unsigned long int repeatPath1=m_repeat_vertex;
				string path2=getTrivialPath(neighboursVector[1],vertex);
				unsigned long int repeatPath2=m_repeat_vertex;
				bool path2IsRepeat=m_DETECTED_REPEAT;
				string lastKMerPath1=path1.substr(path1.length()-m_wordSize,m_wordSize);
				string lastKMerPath2=path2.substr(path2.length()-m_wordSize,m_wordSize);
				if(lastKMerPath1==lastKMerPath2){// bubble
					int diff=path1.length()-path2.length();
					if(diff<0)
						diff=-diff;
					if(diff<=maximumDifferenceInBubbles){
						m_cout<<endl;
						// TODO: maybe chop it all
						m_cout<<"Bubble found!, removing it... "<<path1.length()<<" "<<path2.length()<<endl;
						// TODO SNP detection goes here ?? or it can be sequencing errors also..
						m_colors[neighboursVector[1]]=m_COLOR_DISCARDED;
						m_cout<<path1<<endl;
						m_cout<<path2<<endl;
					}else{
						m_cout<<endl;
						m_cout<<"Non-trivial ("<<diff<<") Bubble found!, skipping it... "<<path1.length()<<" "<<path2.length()<<endl;
						m_cout<<path1<<endl;
						m_cout<<path2<<endl;
					}
				}else{
					//TODO here
					int mismatches=0;
					int positionInPaths=0;
					int minimumPathSize=path1.length();
					if((int)minimumPathSize>(int)path2.length())
						minimumPathSize=(int)path2.length();
					while(positionInPaths<minimumPathSize){
						positionInPaths++;
						if(path1[positionInPaths]!=path2[positionInPaths])
							mismatches++;
					}
					if(mismatches<=3&&path1.length()!=path2.length()){
						m_cout<<endl;
						m_cout<<mismatches<<" mismatches in first "<<minimumPathSize<< " bases (2 things, not a bubble): taking the longuest"<<endl;
						m_cout<<path1<<endl;
						m_cout<<path2<<endl;
						if((int)minimumPathSize==(int)path1.length()){
							m_colors[neighboursVector[0]]=m_COLOR_DISCARDED;
						}else{
							m_colors[neighboursVector[1]]=m_COLOR_DISCARDED;
						}
					}else{
						m_cout<<endl;
						m_cout<<"Working on it (2 things, not a bubble), "<<mismatches<<" mismatches in the first "<<minimumPathSize<<" bases"<<endl;
						if(path1IsRepeat){
							m_cout<<" Path 1 is a REPEAT ("<<getTrivialPath(getNeighbours(repeatPath1)[0],repeatPath1)<<")"<<endl;
						}
						if(path2IsRepeat){
							m_cout<<" Path 2 is a REPEAT ("<<getTrivialPath(getNeighbours(repeatPath2)[0],repeatPath2)<<")"<<endl;
						}
						m_cout<<path1<<endl;
						m_cout<<path2<<endl;
					}
				}
			}else if(neighboursVector.size()==3){
				// TODO here
				m_cout<<endl;
				m_cout<<"3 alternative paths"<<endl;
				for(int u=0;u<(int)neighboursVector.size();u++){
					m_cout<<getTrivialPath(neighboursVector[u],vertex)<<endl;
				}
			}else if(neighboursVector.size()>3){
				// TODO here
				m_cout<<endl;
				m_cout<<"> 3 alternative paths"<<endl;
				for(int u=0;u<(int)neighboursVector.size();u++){
					m_cout<<getTrivialPath(neighboursVector[u],vertex)<<endl;
				}
			}

		}
	}

}





// concert k-mer to unsigned long int
// algo works with (k+1)-mer
// TODO: replace *= by <<
unsigned long int DeBruijnAssembler::wordId(const char*a){
	unsigned long int i=0;
	for(int j=0;j<(int)strlen(a);j++){
		unsigned long int k=0;
		if(a[j]=='A'){
			// binary=00
			// dec = 0
			k=m_NUCLEOTIDE_A;
		}else if(a[j]=='T'){
			// binary=01
			//  dec = 1
			k=m_NUCLEOTIDE_T;
		}else if(a[j]=='C'){
			// binary=10
			// dec = 2
			k=m_NUCLEOTIDE_C;
		}else if(a[j]=='G'){
			// binary=11
			// dec = 3
			k=m_NUCLEOTIDE_G;
		}
		for(int l=0;l<=j;l++){
			k*=4; // right shift two times two positions
		}
		i+=k;
	}
	//m_cout<<"[DeBruijnAssembler,wordId] "<<a<<" -> "<<i<<endl;
	//string word=idToWord(i);
	//m_cout<<"Conversion <"<<a<<"> <"<<word<<">"<<endl;
	//m_cout<<"getLastSymbol "<<getLastSymbol(i)<<endl;
	return i;
}






// TODO: replace *= by << and /= by >>
string DeBruijnAssembler::idToWord(unsigned long int i,int wordSize){
	Logger&m_cout=*(this->m_cout);
	string a="";
	int maxSize=sizeof(long int)*8/2; // 32
	for(int p=0;p<wordSize;p++){
		unsigned long int j=i;
		for(int k=0;k<(maxSize-p-2);k++){
			j*=4;
		}
		for(int k=0;k<(maxSize-1);k++){
			j/=4;
		}
		if(j==0){
			a+='A';
		}else if(j==1){
			a+='T';
		}else if(j==2){
			a+='C';
		}else if(j==3){
			a+='G';
		}else{
			m_cout<<"j "<<j<<endl;
		}
		

	}
	return a;
}

// TODO: test this, actually it is only for debuging purposes.
void DeBruijnAssembler::printBinary(unsigned long long int i,int wordSize){
	Logger&m_cout=*(this->m_cout);
	int maxSize=sizeof(unsigned long long int)*8/2; // 32
	for(int p=0;p<wordSize;p++){
		unsigned long long int j=i;
		for(int k=0;k<(maxSize-p-2);k++){
			j*=4;
		}
		for(int k=0;k<(maxSize-1);k++){
			j/=4;
		}
		if(j==0){
			m_cout<< "00";
		}else if(j==1){
			m_cout<<"01";
		}else if(j==2){
			m_cout<<"10";
		}else if(j==3){
			m_cout<<"11";
		}else{
		}
		

	}

	m_cout<<endl;
}


// get the last symbol
// TODO: replace /= by >>1
char DeBruijnAssembler::getLastSymbol(unsigned long int i){
        unsigned long int j=i;
        for(int i=0;i<m_wordSize;i++){
                j/=2;
                j/=2;
        }

        if((int)j==m_NUCLEOTIDE_A)
                return 'A';
        if((int)j==m_NUCLEOTIDE_T)
                return 'T';
        if((int)j==m_NUCLEOTIDE_C)
                return 'C';
        if((int)j==m_NUCLEOTIDE_G)
                return 'G';
        return 'E';
}



// THE END IS HERE  ------------> .
