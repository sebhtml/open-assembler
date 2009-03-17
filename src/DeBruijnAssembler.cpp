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

#include"CoverageDistribution.h"
#include<cmath>
#include<cstdlib>
#include<map>
#include<hash_set>
#include<stack>
#include<queue>
#include"BinarySearch.h"
#include<hash_map>
#include"SortedList.h"
#include<iostream>
#include"DeBruijnAssembler.h"
#include<fstream>
#include<vector>
#include"GraphData.h"
#include<string>
#include<vector>
using namespace std;


int min(int a,int b){
	if(a<b)
		return a;
	return b;
}

int abs_f(int a){
	if(a<0)
		return -a;
	return a;
}


DeBruijnAssembler::DeBruijnAssembler(){
	m_DEBUG=false;
	m_Solexa_detected=false;
}

void DeBruijnAssembler::setWordSize(int k){
	m_wordSize=k;
	m_minimumCoverage=5;
}






void DeBruijnAssembler::build_From_Scratch(SequenceDataFull*sequenceData){
	cout<<endl;
	cout<<"********** Collecting mers from reads..."<<endl;
	cout<<endl;
	cout<<"k+1 = "<<m_wordSize+1<<endl;
	SortedList myList;
	int last_vertices_size=-1;
	int seq36=0;
	for(int i=0;i<(int)sequenceData->size();i++){
		if(i%100000==0){
			cout<<"Reads: "<<i<<" / "<<sequenceData->size()<<endl;
		}
		if(i%1000000==0)
			myList.sort();
		string readSequence=sequenceData->at(i)->getSeq();
		if(readSequence.length()<50)
			seq36++;
		for(int p=0;p<readSequence.length();p++){
			string word=readSequence.substr(p,m_wordSize+1);
			if(word.length()==m_wordSize+1&&isValidDNA(word.c_str())){
				myList.add(wordId(word.c_str()));
				myList.add(wordId(reverseComplement(word).c_str()));
			}
		}
	}	

	if(seq36>=1000000)
		m_Solexa_detected=true;
	cout<<"Reads: "<<sequenceData->size()<<" / "<<sequenceData->size()<<endl;
	myList.sort();

	int processed=0;
	int solid=0;

	cout<<endl;

	//
	//
	CoverageDistribution coverageDistributionObject(myList.getDistributionOfCoverage(),m_assemblyDirectory);
	m_minimumCoverage=coverageDistributionObject.getMinimumCoverage();
	m_coverage_mean=coverageDistributionObject.getMeanCoverage();

	if(m_minimumCoverageParameter!="auto"){
		m_minimumCoverage=atoi(m_minimumCoverageParameter.c_str());
		cout<<"Setting minimumCoverage <- "<<m_minimumCoverage<<endl;
	}
	m_REPEAT_DETECTION=5*m_coverage_mean;
	if(m_minimumCoverage>m_coverage_mean)
		m_REPEAT_DETECTION=5*m_minimumCoverage;
	(cout)<<"REPEAT_DETECTION_COVERAGE =  "<<m_REPEAT_DETECTION<<endl;

	uint64_t total_bases=0;
	//uint64_t solid_bases=0;
	
	vector<VERTEX_TYPE>*solidMers=new vector<VERTEX_TYPE>;
	*solidMers=myList.elementsWithALeastCCoverage(m_minimumCoverage);
	myList.clear();
	cout<<"c-confident mers: "<<solidMers->size()<<", c="<<m_minimumCoverage<<endl;
	if(solidMers->size()==0){
		cout<<"Error: mers are depleted..."<<endl;
		exit(0);
	}
	cout<<endl;
	//words.clear();


	cout<<"********** Creating vertices..."<<endl;
	cout<<endl;
	SortedList graphNodesList;
	for(vector<VERTEX_TYPE>::iterator i=solidMers->begin();i!=solidMers->end();i++){
		VERTEX_TYPE node=*i;
		string wordString=idToWord(node,m_wordSize+1);
		VERTEX_TYPE prefix=wordId(wordString.substr(0,m_wordSize).c_str());
		VERTEX_TYPE suffix=wordId(wordString.substr(1,m_wordSize).c_str());
		graphNodesList.add(prefix);
		graphNodesList.add(suffix);
	}
	graphNodesList.sort();
	vector<VERTEX_TYPE>*nodes=new vector<VERTEX_TYPE>;
	*nodes=graphNodesList.elementsWithALeastCCoverage(1);
	graphNodesList.clear();
	for(vector<VERTEX_TYPE>::iterator i=nodes->begin();i!=nodes->end();i++){
		m_data.add(*i);
	}
	cout<<nodes->size()<<" vertices"<<endl;
	nodes->clear();
	delete nodes;
	nodes=NULL;
	cout<<endl;
	m_data.makeMemory();

	cout<<"********** Creating edges..."<<endl;
	cout<<endl;
	string edgesFile=m_assemblyDirectory+"/Edges.txt";
	ofstream edgesStream(edgesFile.c_str());
	for(vector<VERTEX_TYPE>::iterator i=solidMers->begin();i!=solidMers->end();i++){
		VERTEX_TYPE node=*i;
		string wordString=idToWord(node,m_wordSize+1);
		edgesStream<<wordString<<endl;
		VERTEX_TYPE prefix=wordId(wordString.substr(0,m_wordSize).c_str());
		VERTEX_TYPE suffix=wordId(wordString.substr(1,m_wordSize).c_str());
		m_data.get(prefix)->addChild(suffix,m_wordSize);
		m_data.get(suffix)->addParent(prefix,m_wordSize);
	}
	edgesStream.close();
	cout<<solidMers->size()<<" edges"<<endl;
	solidMers->clear();
	delete solidMers;
	solidMers=NULL;
	cout<<endl;


	cout<<"********** Indexing solid mers in reads..."<<endl; // <-------
	cout<<endl;

	for(int readId=0;readId<(int)sequenceData->size();readId++){
		if(readId%10000==0)
			cout<<"Reads: "<<readId<<" / "<<sequenceData->size()<<endl;
		indexReadStrand(readId,'F',sequenceData);
		indexReadStrand(readId,'R',sequenceData);
	
	}
	cout<<"Reads: "<<sequenceData->size()<<" / "<<sequenceData->size()<<endl;
	cout<<endl;

}



void DeBruijnAssembler::buildGraph(){
	bool debug=m_DEBUG;

	build_From_Scratch(m_sequenceData);
	writeGraph();
	string parametersFile=m_assemblyDirectory+"/Parameters.txt";
	ofstream f(parametersFile.c_str());
	f<<"WordSize "<<m_wordSize<<endl;
	f<<"MinimumCoverage "<<m_minimumCoverage<<"\n";
	f<<"PeakCoverage "<<m_coverage_mean<<"\n";
	f<<"RepeatDetectionCoverage "<<m_REPEAT_DETECTION<<"\n";
	f.close();
}

void DeBruijnAssembler::loadParameters(){
	string parametersFile=m_assemblyDirectory+"/Parameters.txt";
	ifstream f(parametersFile.c_str());
	string buffer;
	f>>buffer>>m_wordSize>>buffer>>m_minimumCoverage>>buffer>>m_coverage_mean>>buffer>>m_REPEAT_DETECTION;
	f.close();
}


void DeBruijnAssembler::load_graphFrom_file(){
	cout<<"********** Reading graph file."<<endl;
	cout<<endl;
	ifstream f(m_graphFile.c_str());
	string version;
	string buffer;
	f>>buffer;
	int n;
	
	f>>n;
	for(int i=0;i<n;i++){
		if(i%1000000==0){
			cout<<"Loading vertices: "<<i<<" / "<<n<<endl;
		}
		VERTEX_TYPE a;
		f>>a;
		m_data.add(a);
	}
	cout<<"Loading vertices: "<<n<<" / "<<n<<endl;
	int total_Edges=0;
	m_data.makeMemory();
	VertexData*dataPointer=m_data.getNodeData();
	f>>buffer>>buffer;
	cout<<endl;
	for(int i=0;i<n;i++){
		if(i%1000000==0){
			cout<<"Loading edges: "<<i<<" / "<<n<<endl;
		}
		VERTEX_TYPE a;
		f>>a;
		int nAnnotations;
		f>>nAnnotations;
		VertexData*dataPointerVertex=&(dataPointer[i]);
		for(int k=0;k<nAnnotations;k++){
			uint32_t readId;
			uint8_t readStrand;
			uint16_t readPosition;
			f>>readId>>readStrand>>readPosition;
			dataPointerVertex->addAnnotation(readId,readPosition,readStrand);
		}

		int childrenCount;
		f>>childrenCount;
		for(int j=0;j<childrenCount;j++){
			VERTEX_TYPE b;
			f>>b;
			VertexData*childContainer=m_data.get(b);
			childContainer->addParent(a,m_wordSize);
			dataPointerVertex->addChild(b,m_wordSize);
			total_Edges++;
		}
	}

	cout<<"Loading edges: "<<n<<" / "<<n<<endl;
	f.close();
}

void DeBruijnAssembler::writeGraph(){
	cout<<"********** Writing graph file."<<endl;

	ofstream f(m_graphFile.c_str());
	
	vector<VERTEX_TYPE>*nodes=m_data.getNodes();
	f<<"Vertices: "<<m_data.size()<<"\n";
	int k=0;
	for(int i=0;i<m_data.size();i++){
		if(k%100000==0){
			cout<<"Vertices: "<<k<<" / "<<m_data.size()<<"\n";
		}
		f<<(*nodes)[i]<<"\n";
		k++;
	}

	cout<<"Vertices: "<<m_data.size()<<" / "<<m_data.size()<<"\n";
	f<<"Data: "<<m_data.size()<<"\n";
	VertexData*theData=m_data.getNodeData();
	k=0;
	for(int i=0;i<m_data.size();i++){
		if(k%100000==0){
			cout<<"Edges: "<<k<<" / "<<m_data.size()<<endl;
		}
		k++;
		VERTEX_TYPE prefix=(*nodes)[(i)];
		f<<prefix<<"\n";
		VertexData*prefixData=&(theData[(i)]);
		vector<VERTEX_TYPE>children=prefixData->getChildren(prefix,m_wordSize);
		vector<AnnotationElement>*annotations=prefixData->getAnnotations();
		if(annotations!=NULL){
			f<<annotations->size()<<"\n";
			for(vector<AnnotationElement>::iterator u=annotations->begin();u!=annotations->end();u++){
				f<<(*u).readId<<" "<<(*u).readStrand<<" "<<(*u).readPosition<<"\n";
			}
		}

		f<<children.size()<<"\n";
		for(vector<VERTEX_TYPE>::iterator k=children.begin();k!=children.end();k++){
			f<<*k<<"\n";
		}
	}
	cout<<"Edges: "<<m_data.size()<<" / "<<m_data.size()<<endl;

	f.close();
}

string DeBruijnAssembler::pathToDNA(vector<VERTEX_TYPE>*path){
	ostringstream contigSequence;
	for(vector<VERTEX_TYPE>::iterator i=path->begin();i!=path->end();i++){
		if(i==path->begin()){
			contigSequence<< idToWord(*i,m_wordSize);
		}else{
			contigSequence<<getLastSymbol(*i,m_wordSize);
		}
	}
	return contigSequence.str();
}




 
void DeBruijnAssembler::setAssemblyDirectory(string assemblyDirectory){
	m_assemblyDirectory=assemblyDirectory;
	ostringstream name;
	name<<m_assemblyDirectory<<"/graph";
	m_graphFile=name.str();
}



void DeBruijnAssembler::Walk_In_GRAPH(){
	cout<<endl;
	vector<VERTEX_TYPE>*theNodes=m_data.getNodes();

	cout<<endl;
	vector<VERTEX_TYPE> withoutParents;
	(cout)<<"********* Inspecting the graph"<<endl;
	cout<<endl;
	map<int,map<int,int> > stats_parents_children;
	for(int i=0;i<m_data.size();i++){
		if(i%1000000==0){
			cout<<"Inspecting: "<<i<<" / "<<m_data.size()<<endl;
		}
		VERTEX_TYPE vertex=(*theNodes)[i];
		VertexData*dataNode=m_data.get(vertex);
		vector<uint64_t> theParents=dataNode->getParents(vertex,m_wordSize);

		int parents=theParents.size();
		int children=dataNode->getChildren(vertex,m_wordSize).size();
		stats_parents_children[parents][children]++;
		if(parents==1&&children==1){
			dataNode->set_topology_1_1();
		}
		if(parents==1&&children==1&&m_data.get(theParents[0])->getChildren(theParents[0],m_wordSize).size()<2)
			continue;
		if(children==0)
			continue;
		withoutParents.push_back(vertex);
	}

	cout<<"Inspecting: "<<m_data.size()<<" / "<<m_data.size()<<endl;

	cout<<endl;
	for(map<int,map<int,int> >::iterator i=stats_parents_children.begin();i!=stats_parents_children.end();i++){
		for(map<int,int>::iterator j=i->second.begin();j!=i->second.end();j++){
			cout<<i->first<<" parents, "<<j->first<<" children: "<<j->second<<" vertices"<<endl;
		}
	}
	cout<<endl;

	vector<VERTEX_TYPE>*nodes=m_data.getNodes();
	VertexData*nodeData=m_data.getNodeData();
	(cout)<<endl;
	cout<<endl;

	(cout)<<"Done..., "<<withoutParents.size()<<" sources."<<endl;
	
	(cout)<<endl;
	vector<VERTEX_TYPE> sources=withoutParents;
	string assemblyAmos=m_assemblyDirectory+"/"+AMOS_FILE_NAME;
	string contigsFile=m_assemblyDirectory+"/"+FASTA_FILE_NAME;
	string coverageFile=m_assemblyDirectory+"/"+COVERAGE_FILE_NAME;
	string repeatAnnotationFile=m_assemblyDirectory+"/contigs-repeats.txt";
	ofstream amosFile(assemblyAmos.c_str());
	ofstream contigsFileStream(contigsFile.c_str());
	ofstream coverageStream(coverageFile.c_str());
	ofstream repeatAnnotation(repeatAnnotationFile.c_str());
	int contigId=1;
	int round=1;

	cout<<"********** Assembling contigs..."<<endl;

	while(sources.size()>0){
		cout<<endl;
		vector<VERTEX_TYPE> newSources;
		for(int i=0;i<(int)sources.size();i++){
			if(i%1000==0){
				updateDebug();
			}
			VERTEX_TYPE prefix=sources[i];
			if(m_data.get(prefix)->IsAssembled())
				continue;
			cout<<endl;
			cout<<"Source "<<i+1<<" / "<<sources.size()<<endl;//" REPEAT MODE"<<endl;
			cout<<"From: "<<idToWord(prefix,m_wordSize)<<endl;
		//cout<<m_contig_paths.size()<<" contigs"<<endl;
			vector<VERTEX_TYPE>path;
			vector<map<int,map<char,int> > > currentReadPositions;
			vector<int> repeatAnnotations;
			vector<VERTEX_TYPE> localNewSources;
			version2_Walker(prefix,&path,&currentReadPositions,&repeatAnnotations);
			cout<<path.size()<<" vertices"<<endl;
			if(path.size()<=2)
				continue;
			//writeContig_Amos(&currentReadPositions,&path,&amosFile,contigId);
			writeContig_fasta(&path,&contigsFileStream,contigId);
			//writeContig_Coverage(&currentReadPositions,&path,&coverageStream,contigId);
			//writeContig_RepeatAnnotation(&repeatAnnotations,contigId,&repeatAnnotation,&path);
			//m_contig_paths.push_back(path);
			cout<<"Contig"<<contigId<<endl;
			contigId++;
		}
		sources=newSources;
	}
	amosFile.close();
	repeatAnnotation.close();
	contigsFileStream.close();
	coverageStream.close();
	cout<<endl;

/*
	(cout)<<"********** Sources discovery..."<<endl;
	(cout)<<"Iteration Sources Cumulative"<<endl;
	for(int i=0;i<(int)The_Discovery_Of_Sources.size();i++)
		(cout)<<i+1<<" "<<The_Discovery_Of_Sources[i]<<" "<<VisitsOfSources[i]<<endl;

	(cout)<<endl;
*/
	return;
}

void DeBruijnAssembler::Algorithm_Assembler_20090121(){
	Walk_In_GRAPH();
}



void DeBruijnAssembler::version2_Walker(uint64_t  a,vector<uint64_t>*path,
			vector<map<int,map<char,int> > >* currentReadPositions,
			vector<int>*repeatAnnotations){
	//cout<<"Starting from "<<a<<endl;
	vector<uint64_t>contig;
	vector<uint64_t>  children;
	children.push_back(a);

	hash_set<int> usedReads;
	hash_set<int>  readsInRange;

	// read, readposition
	hash_map<int,int> readsReadPosition;

	// read, position in contig
	hash_map<int,int> readsContigPositions;

	// read, read strand
	hash_map<int,char> readsReadStrands;

	while(children.size()==1){
/*
		if(contig.size()%1000==0)
			cout<<"CONTIG LENGTH "<<contig.size()<<endl;
*/
		uint64_t  currentVertex=children[0];
		VertexData*aData=m_data.get(currentVertex);
		if(aData->IsAssembled()&&contig.size()<200/*&&fromAPureParent*/){
			cout<<"Skipping spurious vertex"<<endl;
			return;
		}
		contig.push_back(currentVertex);
		//repeatAnnotations->push_back(-1);
		//cout<<idToWord(currentVertex,m_wordSize)<<endl;
		aData->assemble();
		
		/*
		map<int,map<char,int> > readPositions;
		currentReadPositions->push_back(readPositions);
		*/
		addAnnotations(aData,&usedReads,&readsInRange,
			&readsReadPosition,&readsContigPositions,&readsReadStrands,
			&contig,currentReadPositions);


		// process annotations
		children=aData->getChildren(currentVertex,m_wordSize);

		map<uint64_t,int> sumScores;
		map<uint64_t,vector<int> > annotationsForEach;



		bool skipThoroughtCheck=true;

		ThreadReads(aData,&usedReads,&readsInRange,
			&readsReadPosition,&readsContigPositions,&readsReadStrands,
			&contig,
			&sumScores,&annotationsForEach,&children,skipThoroughtCheck,
			currentReadPositions);

		// annotate children
	/*
		if(children.size()>1){
			if(getDebug()){
				cout<<endl;
			}
		}
*/


		double alpha=1.1;
		if(m_Solexa_detected)
			alpha=1.3;
		if(children.size()==2){
			int scoreA=sumScores[children[0]];
			int scoreB=sumScores[children[1]];

			if(scoreA >= alpha*scoreB){
				children.clear();
				children.push_back(children[0]);
			}else if(scoreB >= alpha*scoreA){
				children.clear();
				children.push_back(children[1]);
			}
		}

		if(children.size()>1){  // reduce it...
			//cout<<"MORE THAN 1"<<endl;
			for(map<uint64_t,vector<int> >::iterator i=annotationsForEach.begin();i!=annotationsForEach.end();i++){
				uint64_t currentVertex=i->first;
				int currentScore=sumScores[currentVertex];

				// try with sum
				bool isBest=true;
				for(map<uint64_t,vector<int> >::iterator j=annotationsForEach.begin();j!=annotationsForEach.end();j++){
					uint64_t otherVertex=j->first;
					int otherScore=sumScores[otherVertex];
					if(currentVertex!=otherVertex && currentScore < alpha*otherScore){
						isBest=false;
						break;
					}
				}
				if(isBest){
					children.clear();
					children.push_back(currentVertex);
					if(getDebug())
						cout<<"GOT BEST USING SUM "<<idToWord(currentVertex,m_wordSize)<<endl;
					break;
				}

				// try with number
				isBest=true;
				currentScore=annotationsForEach[currentVertex].size();
				for(map<uint64_t,vector<int> >::iterator j=annotationsForEach.begin();j!=annotationsForEach.end();j++){
					uint64_t otherVertex=j->first;
					int otherScore=annotationsForEach[otherVertex].size();
					if(currentVertex!=otherVertex && currentScore < 1.5*otherScore){
						isBest=false;
						break;
					}
				}
				if(isBest){
					children.clear();
					children.push_back(currentVertex);
		
					if(getDebug())
						cout<<"GOT BEST USING NUMBER"<<idToWord(currentVertex,m_wordSize)<<endl;
					break;
				}

			}
		}
		
		// detect homopolymer?
		if(children.size()==2){
			string first=idToWord(children[0],m_wordSize);
			string second=idToWord(children[1],m_wordSize);
			int iii=0;
			char theNuc='E';
			if(first[m_wordSize-2]==first[m_wordSize-1]){
				theNuc=first[m_wordSize-1];
			}else if(second[m_wordSize-2]==second[m_wordSize-1]){
				theNuc=second[m_wordSize-1];	
			}
			while(m_wordSize-1-iii>=0&&first[m_wordSize-2-iii]==second[m_wordSize-2-iii]&&
			first[m_wordSize-2-iii]==theNuc){
				if(theNuc=='E')
					break;
				iii++;
			}
			if(iii>=6){
				if(first[m_wordSize-1]==first[m_wordSize-2]||
				second[m_wordSize-1]==second[m_wordSize-2]){
					if(getDebug())
						cout<<"HOMOPOLYMER, n="<<iii<<endl;
					if(annotationsForEach[children[0]]>annotationsForEach[children[1]]){
						children.clear();
						children.push_back(children[0]);
					}else{
						children.clear();
						children.push_back(children[1]);
					}
				}
			}
		}

		if(children.size()==1&&(aData->Is_1_1()==false||skipThoroughtCheck==false)){ // check if it is ok
			uint64_t currentVertex=children[0];
			if(annotationsForEach[currentVertex].size()==0&&contig.size()>400){
				//if(getDebug())
				cout<<"No annotation found for "<<idToWord(currentVertex,m_wordSize)<<"."<<endl;
				break;
			}
		}


	}
	cout<<"STOP <<CHILDREN="<<children.size()<<endl;
	//cout<<"CONTIG LENGTH [FINAL] "<<contig.size()<<endl;
	*path=contig;
}







DeBruijnAssembler::~DeBruijnAssembler(){
}




void DeBruijnAssembler::setMinimumCoverage(string coverage){
	m_minimumCoverageParameter=coverage;
}

void DeBruijnAssembler::indexReadStrand(int readId,char strand,SequenceDataFull*sequenceData){
	Read*read=sequenceData->at(readId);
	string sequence=read->getSeq();

	if(strand=='R')
		sequence=reverseComplement(sequence);
	bool foundGoodHit=false;
	for(int readPosition=0;readPosition<(int)sequence.length();readPosition++){
		string myWord=sequence.substr(readPosition,m_wordSize);
		if(foundGoodHit)
			break;
		if(myWord.length()!=m_wordSize)
			continue;
		VERTEX_TYPE wordInBits=wordId(myWord.c_str());
		if(isValidDNA(myWord.c_str())
		&&m_data.hasNode(wordInBits)){
			bool thisIsTheFirst=false;
			if(foundGoodHit==false){
				foundGoodHit=true;
				thisIsTheFirst=true;
				//cout<<"Starting in: "<<readPosition<<endl;
				if(strand=='F'){
					sequenceData->at(readId)->setStartForward(readPosition);
				}else if(strand=='R'){
					sequenceData->at(readId)->setStartReverse(readPosition);
				}
			}

			// the position is the first, or
			// it is on a non-trivial node...
			// TODO: remove the true
			// TODO only annotate the  first
			VertexData*aVertexData=m_data.get(wordInBits);
			bool hasMixedParents=false;
			vector<VERTEX_TYPE> theParents=m_data.get(wordInBits)->getParents(wordInBits,m_wordSize);
			for(vector<VERTEX_TYPE>::iterator ii=theParents.begin();ii!=theParents.end();ii++){
				vector<VERTEX_TYPE> theChildren=m_data.get(*ii)->getChildren(*ii,m_wordSize);
				if(theChildren.size()>1)
					hasMixedParents=true;
			}
			if(/*m_data.get(wordInBits)->NotTrivial(wordInBits,m_wordSize)||*/(strand=='F'&&readPosition==sequenceData->at(readId)->getStartForward())||
			(strand=='R'&&readPosition==sequenceData->at(readId)->getStartForward())||false)
				m_data.get(wordInBits)->addAnnotation(readId,readPosition,strand);
		}
	}
}


void DeBruijnAssembler::writeContig_RepeatAnnotation(vector<int>*repeatAnnotations,int i,ofstream*file,vector<VERTEX_TYPE>*path){
	(*file)<<"Contig"<<i<<endl;
	int p=1;
	for(int vertexOffset=0;vertexOffset<path->size();vertexOffset++){
		if(vertexOffset==0){
			while(p<=m_wordSize){
				int edgeOffset=vertexOffset;
				(*file)<<p<<" "<<repeatAnnotations->at(edgeOffset)<<endl;
				p++;
			}
		}else{
			int edgeOffset=vertexOffset-1;
			(*file)<<p<<" "<<repeatAnnotations->at(edgeOffset)<<endl;
			p++;
		}
	}
}

void DeBruijnAssembler::writeContig_Coverage(vector<map<int,map<char,int> > >*currentReadPositions,vector<VERTEX_TYPE>*path,ofstream*file,int i){
	(*file)<<"Contig"<<i<<endl;
	int p=1;
	for(int vertexOffset=0;vertexOffset<path->size();vertexOffset++){
		if(vertexOffset==0){
			while(p<=m_wordSize){
				int edgeOffset=vertexOffset;
				(*file)<<p<<" "<<currentReadPositions->at(edgeOffset).size()<<endl;
				p++;
			}
		}else{
			int edgeOffset=vertexOffset;
			(*file)<<p<<" "<<currentReadPositions->at(edgeOffset).size()<<endl;
			p++;
		}
	}
}

/*

http://www.cbcb.umd.edu/research/contig_representation.shtml#AMOS

{CTG
iid:1
eid:1
seq:
CCTCTCCTGTAGAGTTCAACCGA-GCCGGTAGAGTTTTATCA
.
qlt:
DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
.
{TLE
src:1027
off:0
clr:618,0
gap:
250 612
.
}
}

*/
void DeBruijnAssembler::writeContig_Amos(vector<map<int,map<char,int> > >*currentReadPositions,vector<VERTEX_TYPE>*path,ofstream*file,int i){
	//(cout)<<"writeContig_Amos"<<endl;
	string sequenceDNA=pathToDNA(path);
	(*file)<<"{CTG"<<endl;
	//(*file)<<"com:generated by DNA"<<endl;
	(*file)<<"iid:"<<i<<endl;
	(*file)<<"eid:Contig"<<i<<endl;
	(*file)<<"seq:"<<endl;
	(*file)<<sequenceDNA<<endl;
	(*file)<<"."<<endl;
	(*file)<<"qlt:"<<endl;
	for(int i=0;i<sequenceDNA.length();i++)
		(*file)<<"D";
	(*file)<<endl;
	(*file)<<"."<<endl;

	map<int,int> readOffset;
	map<int,int> readStart;
	map<int,int> readEnd;
	map<int,char> readStrand;

	//(cout)<<"amos generator "<<currentReadPositions->size()<<endl;
	int contigPosition=0;
	for(vector<map<int,map<char,int> > >::iterator i=currentReadPositions->begin();
		i!=currentReadPositions->end();i++){
		//int contigPosition=i->first;
		for(map<int,map<char,int> >::iterator j=(*i).begin();
			j!=(*i).end();j++){
			int readId=j->first;
			char strand=j->second.begin()->first;
			int readPosition=j->second.begin()->second;
			if(readOffset.count(readId)==0){
				readOffset[readId]=contigPosition;
				readStart[readId]=0;
				readEnd[readId]=0;
				readStrand[readId]=strand;
			}
			if(readPosition>readEnd[readId])
				readEnd[readId]=readPosition;
		}
		contigPosition++;
	}

	for(map<int,int>::iterator i=readOffset.begin();i!=readOffset.end();i++){
		int readId=i->first;
		int readOffset=i->second;
		(*file)<<"{TLE"<<endl;
		//(*file)<<"src:"<<m_sequenceData->at(readId)->getId()<<endl;
		(*file)<<"src:"<<readId+1<<endl;
		int readLength=strlen(m_sequenceData->at(readId)->getSeq());
		(*file)<<"com:"<<m_sequenceData->at(readId)->getId()<<" "<<readLength<<" "<<readStrand[readId]<<" "<<readStart[readId]<<" "<<readEnd[readId]<<endl;

/*
 

---------------------------------------------
          ------------------------>



-------------------------------------------
            <--------------------------




*/
		char strand=readStrand[readId];
		if(strand=='F')
			(*file)<<"off:"<<readOffset+0<<endl;
		else
			(*file)<<"off:"<<readOffset+0<<endl;
		(*file)<<"clr:";
		if(strand=='F')
			(*file)<<readStart[readId]<<","<<readEnd[readId]+m_wordSize+1<<endl;
		else
			(*file)<<readLength-readStart[readId]<<","<<readLength-readEnd[readId]-m_wordSize-1<<endl;
			//(*file)<<readEnd[readId]+m_wordSize<<","<<readStart[readId];
		//(*file)<<endl;
		(*file)<<"gap:"<<endl;
		(*file)<<"."<<endl;
		(*file)<<"}"<<endl;
		
	}
	(*file)<<"}"<<endl;
}

void DeBruijnAssembler::writeContig_fasta(vector<VERTEX_TYPE>*path,ofstream*file,int i){
	int columns=60;
	string sequenceDNA=pathToDNA(path);
	(*file)<<">Contig"<<i<<"   "<<sequenceDNA.length()<<endl;
	int j=0;
	while(j<(int)sequenceDNA.length()){
		(*file)<<sequenceDNA.substr(j,columns);
		(*file)<<endl;
		j+=columns;
	}
}

void DeBruijnAssembler::debug(){
	m_DEBUG=true;
}




void DeBruijnAssembler::setSequenceData(SequenceDataFull*sequenceData){
	m_sequenceData=sequenceData;
}

void DeBruijnAssembler::setStrandUsage(string a){
	m_onlyOneStrand=a;
}

void DeBruijnAssembler::setMerUsage(string a){
	m_onlyFirstMer=a;
}

void DeBruijnAssembler::loadPairedInformation(){
	string pairedFile=m_assemblyDirectory+"/PairedReads.txt";
	int groups;

	ifstream f( pairedFile.c_str());

	if(!f){
		f.close();
		return;
	}
	f>>groups;
	string file1;
	string file2;
	int distance;
	for(int i=0;i<groups;i++){
		f>>file1>>file2>>distance;
		int file_1_start=m_sequenceData->getFirst(file1);
		int file_1_end=m_sequenceData->getLast(file1);
		int file_2_start=m_sequenceData->getFirst(file2);
		int file_2_end=m_sequenceData->getLast(file2);
		int amount=file_1_end-file_1_start+1;
		for(int j=0;j<amount;j++){
			int read1=file_1_start+j;
			int read2=file_2_start+j;
			PairedRead aPairedRead1;
			aPairedRead1.m_readNumber=read1;
			aPairedRead1.m_distance=distance;
			PairedRead aPairedRead2;
			aPairedRead2.m_readNumber=read2;
			aPairedRead2.m_distance=distance;
			m_paired_reads[read1]=aPairedRead2;
			m_paired_reads[read2]=aPairedRead1;
		}
	}
	cout<<m_paired_reads.size()<<" paired reads"<<endl;
	f.close();
}

void DeBruijnAssembler::addDebug(){
	m_DEBUG=true;
}

void DeBruijnAssembler::removeDebug(){
	m_DEBUG=false;
}

bool DeBruijnAssembler::getDebug(){
	return m_DEBUG;
}

void DeBruijnAssembler::updateDebug(){
	string file=m_assemblyDirectory+"/DEBUG";
	ifstream f(file.c_str());
	if(!f){
		removeDebug();
	}else{
		addDebug();
	}
	f.close();
}

void DeBruijnAssembler::addAnnotations(VertexData*aData,hash_set<int>*usedReads,hash_set<int>*readsInRange,
	hash_map<int,int>*readsReadPosition,hash_map<int,int>*readsContigPositions,hash_map<int,char>*readsReadStrands,
	vector<uint64_t>*contig,
		vector<map<int,map<char,int> > >* currentReadPositions){
	vector<AnnotationElement>*annotations=aData->getAnnotations();
	if(annotations->size()<m_REPEAT_DETECTION){
		for(int i=0;i<annotations->size();i++){
			int readId=annotations->at(i).readId;
			char readStrand=annotations->at(i).readStrand;
			int readPosition=annotations->at(i).readPosition;
			//int readFirstPosition=m_sequenceData->at(readId)->getStartForward();
			//if(readStrand=='R')
				//readFirstPosition=m_sequenceData->at(readId)->getStartReverse();
			if(/*readPosition==readFirstPosition&&*/usedReads->count(readId)==0){
				usedReads->insert(readId);	
				readsInRange->insert(readId);
				(*readsReadPosition)[readId]=readPosition;
				(*readsContigPositions)[readId]=contig->size()-1;
				(*readsReadStrands)[readId]=readStrand;
				//(*currentReadPositions)[contig->size()-1][readId][readStrand]=readPosition;
			}
		}
	}
}

void DeBruijnAssembler::ThreadReads(VertexData*aData,hash_set<int>*usedReads,hash_set<int>*readsInRange,
	hash_map<int,int>*readsReadPosition,hash_map<int,int>*readsContigPositions,hash_map<int,char>*readsReadStrands,
	vector<uint64_t>*contig,
	map<uint64_t,int>*sumScores,map<uint64_t,vector<int> >*annotationsForEach,vector<uint64_t>*children,bool skipThoroughtCheck,
		vector<map<int,map<char,int> > >* currentReadPositions){
	for(vector<uint64_t>::iterator i=children->begin();i!=children->end();i++){
		if(aData->Is_1_1()&&skipThoroughtCheck)
			break;
		uint64_t childVertex=*i;
		//string childSequence=idToWord(childVertex,m_wordSize);
		vector<int>  readNotInRangeAnymore;
		//cout<<readsInRange.size()<<" reads in range"<<endl;
		for(hash_set<int>::iterator i=readsInRange->begin();i!=readsInRange->end();i++){
			int readId=*i;
			int currentContigPosition=contig->size()-1+1;
			int lastContigPosition=(*readsContigPositions)[readId];
			int lastReadPosition=(*readsReadPosition)[readId];
			char readStrand=(*readsReadStrands)[readId];
			int distanceInContig=currentContigPosition-lastContigPosition;
			int inferedReadPosition=lastReadPosition+distanceInContig;
			int nucleotidePositionInRead=inferedReadPosition;
			if(nucleotidePositionInRead>=m_sequenceData->at(readId)->length()){
				//cout<<"OUT OF RANGE"<<endl;
				readNotInRangeAnymore.push_back(readId);
				continue;
			}
			//ccut<<endl;
			//cout<<"CHILD IS "<<childSequence<<endl;
			/*
			string readSequence=m_sequenceData->at(readId)->getSeq();
			if(readStrand=='R')
				readSequence=reverseComplement(readSequence);
			string aWord=readSequence.substr(inferedReadPosition,m_wordSize);
			*/
			uint64_t readVertex=m_sequenceData->at(readId)->Vertex(inferedReadPosition,m_wordSize,readStrand);
			//cout<<"READ IS  "<<readSequence.substr(inferedReadPosition,m_wordSize)<<endl;
			//if(aWord==childSequence){
			if(childVertex==readVertex){
				//cout<<"In range, threading"<<endl;
				(*annotationsForEach)[childVertex].push_back(nucleotidePositionInRead);
				(*sumScores)[childVertex]+=nucleotidePositionInRead;
				// check paired information
				if(m_paired_reads.count(readId)>0){
					PairedRead pairedInformation=m_paired_reads[readId];
					int otherReadNumber=pairedInformation.m_readNumber;
					if((*readsContigPositions).count(otherReadNumber)>0){
						int lastContigPositionForPairedMate=(*readsContigPositions)[otherReadNumber];
						int distance=pairedInformation.m_distance;
						int windowSemiSize=0.20*distance;
						int distanceInContig=contig->size()-lastContigPositionForPairedMate-nucleotidePositionInRead;
						if(distance-windowSemiSize<=distanceInContig&&
							distanceInContig<=distance+windowSemiSize){
							//cout<<"distance for mates "<<distanceInContig<<endl;
							//cout<<childSequence<<endl;
							(*annotationsForEach)[childVertex].push_back(distanceInContig);
							(*sumScores)[childVertex]+=distanceInContig;
						}
					}
				}
			}
		}
		for(int i=0;i<readNotInRangeAnymore.size();i++)
			(*readsInRange).erase(readNotInRangeAnymore[i]);
		int theScore=(*annotationsForEach)[childVertex].size();

		if(children->size()>1&&getDebug()){
			cout<<"Vertex: "<<idToWord(childVertex,m_wordSize)<<endl;
			cout<<"SUM="<<(*sumScores)[childVertex]<<endl;
			cout<<"n="<<(*annotationsForEach)[childVertex].size()<<endl;
			cout<<"LIST ";

			for(int i=0;i<(*annotationsForEach)[childVertex].size();i++)
				cout<<" "<<(*annotationsForEach)[childVertex][i];
			cout<<endl;
		}
	}
}
