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
#include<stack>
#include<queue>
#include"BinarySearch.h"
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


DeBruijnAssembler::DeBruijnAssembler(ostream*m_cout){
	m_DEBUG=false;
	this->m_cout=m_cout;
}

void DeBruijnAssembler::setWordSize(int k){
	m_wordSize=k;
	m_pairedAvailable=false;
	m_minimumCoverage=5;
}






void DeBruijnAssembler::build_From_Scratch(SequenceDataFull*sequenceData){
	ostream&m_cout=*(this->m_cout);
	m_cout<<endl;
	m_cout<<"********** Collecting mers from reads..."<<endl;
	cout<<endl;
	m_cout<<"k+1 = "<<m_wordSize+1<<endl;
	SortedList myList;
	int last_vertices_size=-1;

	for(int i=0;i<(int)sequenceData->size();i++){
		if(i%100000==0){
			m_cout<<"Reads: "<<i<<" / "<<sequenceData->size()<<endl;
		}
		vector<VERTEX_TYPE> highQualityMers=sequenceData->at(i)->getHighQualityMers(m_wordSize);
		for(vector<VERTEX_TYPE>::iterator iteratorMer=highQualityMers.begin();iteratorMer!=highQualityMers.end();iteratorMer++){
			myList.add(*iteratorMer);
/*
			if(!words.find(*iteratorMer))
				words.add(*iteratorMer,0);

			words.set(*iteratorMer,words.get(*iteratorMer)+1);
			if((int)words.size()%1000000==0&&(int)words.size()!=last_vertices_size){
				m_cout<<"Mers: "<<words.size()<<endl;
				last_vertices_size=words.size();
			}
*/
		}
	}	

	m_cout<<"Reads: "<<sequenceData->size()<<" / "<<sequenceData->size()<<endl;
	myList.sort();
	//m_cout<<"Mers: "<<words.size()<<endl;

/*
	cout<<"********** Buckets analysis"<<endl;
	cout<<"Buckets: "<<words.buckets()<<endl;
	cout<<"Elements: "<<words.size()<<endl;
	map<int,int> bucketsDistribution;
	for(int i=0;i<words.buckets();i++){
		Entry<int>*pointer=words.bucketAt(i);
		int count=0;
		while(pointer!=NULL){
			count++;
			pointer=pointer->m_next;
		}
		bucketsDistribution[count]++;
	}

	cout<<"Distribution of buckets, count, density"<<endl;
	for(map<int,int>::iterator i=bucketsDistribution.begin();i!=bucketsDistribution.end();i++){
		cout<<i->first<<" "<<i->second<<endl;
	}
*/
	int processed=0;
	int solid=0;

	m_cout<<endl;

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
	
	vector<VERTEX_TYPE> solidMers=myList.elementsWithALeastCCoverage(m_minimumCoverage);
	myList.clear();
	cout<<"c-confident mers: "<<solidMers.size()<<", c="<<m_minimumCoverage<<endl;
	if(solidMers.size()==0){
		m_cout<<"Error: mers are depleted..."<<endl;
		exit(0);
	}
	m_cout<<endl;
	//words.clear();


	m_cout<<"********** Creating vertices..."<<endl;
	cout<<endl;
	SortedList graphNodesList;
	for(vector<VERTEX_TYPE>::iterator i=solidMers.begin();i!=solidMers.end();i++){
		VERTEX_TYPE node=*i;
		string wordString=idToWord(node,m_wordSize+1);
		VERTEX_TYPE prefix=wordId(wordString.substr(0,m_wordSize).c_str());
		VERTEX_TYPE suffix=wordId(wordString.substr(1,m_wordSize).c_str());
		graphNodesList.add(prefix);
		graphNodesList.add(suffix);
	}
	graphNodesList.sort();
	vector<VERTEX_TYPE> nodes=graphNodesList.elementsWithALeastCCoverage(1);
	graphNodesList.clear();
	for(vector<VERTEX_TYPE>::iterator i=nodes.begin();i!=nodes.end();i++){
		m_data.add(*i);
	}
	cout<<nodes.size()<<" vertices"<<endl;
	nodes.clear();
	cout<<endl;
	m_data.makeMemory();

	m_cout<<"********** Creating edges..."<<endl;
	cout<<endl;
	for(vector<VERTEX_TYPE>::iterator i=solidMers.begin();i!=solidMers.end();i++){
		VERTEX_TYPE node=*i;
		string wordString=idToWord(node,m_wordSize+1);
		VERTEX_TYPE prefix=wordId(wordString.substr(0,m_wordSize).c_str());
		VERTEX_TYPE suffix=wordId(wordString.substr(1,m_wordSize).c_str());
		m_data.get(prefix)->addChild(suffix,m_wordSize);
		m_data.get(suffix)->addParent(prefix,m_wordSize);
	}
	cout<<solidMers.size()<<" edges"<<endl;
	solidMers.clear();
	cout<<endl;


	m_cout<<"********** Indexing solid mers in reads..."<<endl; // <-------
	cout<<endl;
/*
	VERTEX_TYPE*solidMersPTR=new VERTEX_TYPE[solidMers.size()];
	for(int i=0;i<solidMers.size();i++)
		solidMersPTR[i]=solidMers[i];
*/
	for(int readId=0;readId<(int)sequenceData->size();readId++){
		if(readId%10000==0)
			m_cout<<"Reads: "<<readId<<" / "<<sequenceData->size()<<endl;
		indexReadStrand(readId,'F',sequenceData);
		indexReadStrand(readId,'R',sequenceData);
	
	}
	m_cout<<"Reads: "<<sequenceData->size()<<" / "<<sequenceData->size()<<endl;
	m_cout<<endl;

	// paired information
	if(m_pairedInfoFile=="none")
		return;
	ifstream f(m_pairedInfoFile.c_str());
	if(!f){
		f.close();
		return;
	}
}



void DeBruijnAssembler::buildGraph(SequenceDataFull*sequenceData){
	m_sequenceData=sequenceData;
	ostream&m_cout=*(this->m_cout);
	bool debug=m_DEBUG;
	//debug=false;
	bool useCache=false;
	useCache=true;
	bool writeGraphFile=useCache;
	if(debug){
		useCache=true;
	}

	ifstream f(m_graphFile.c_str());
	if(!f)
		useCache=false;
	f.close();

	if(useCache){
		cout<<endl;
		load_graphFrom_file();
	}else{
		build_From_Scratch(sequenceData);
		if(writeGraphFile){
			writeGraph();
		}
	}

}


void DeBruijnAssembler::load_graphFrom_file(){
	cout<<"********** Reading graph file."<<endl;
	cout<<endl;
	ifstream f(m_graphFile.c_str());
	string version;
	string buffer;
	f>>version>>buffer>>m_minimumCoverage>>buffer>>m_coverage_mean>>buffer>>m_REPEAT_DETECTION>>buffer;
	int n;
	cout<<"Version: "<<version<<endl;
	cout<<"MinimumCoverage: "<<m_minimumCoverage<<endl;
	cout<<"PeakCoverage: "<<m_coverage_mean<<endl;
	cout<<"RepeatDetectionCoverage: "<<m_REPEAT_DETECTION<<endl;
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
	m_data.makeMemory();
	VertexData*dataPointer=m_data.getNodeData();
	f>>buffer>>buffer;
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
		}
	}

	cout<<"Loading edges: "<<n<<" / "<<n<<endl;
	f.close();
}

void DeBruijnAssembler::writeGraph(){
	cout<<"********** Writing graph file."<<endl;

	ofstream f(m_graphFile.c_str());
	f<<"GraphFormatVersion1"<<"\n";
	f<<"MinimumCoverage "<<m_minimumCoverage<<"\n";
	f<<"PeakCoverage "<<m_coverage_mean<<"\n";
	f<<"RepeatDetectionCoverage "<<m_REPEAT_DETECTION<<"\n";
	
	vector<VERTEX_TYPE>*nodes=m_data.getNodes();
	f<<"Vertices: "<<m_data.size()<<"\n";
	int k=0;
	//for(vector<VERTEX_TYPE>::iterator i=nodes->begin();i!=nodes->end();i++){
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


int DeBruijnAssembler::DFS_watch(VERTEX_TYPE a,int color){
	stack<VERTEX_TYPE> theQueue;
	theQueue.push(a);
	m_data.get(a)->setColor(color);
	int numberOfVertices=1;
	while(theQueue.size()>0){
		//cout<<theQueue.size()<<endl;
		VERTEX_TYPE b=theQueue.top();
		theQueue.pop();
		vector<VERTEX_TYPE> parents=m_data.get(b)->getParents(b,NULL,m_wordSize);
		if(parents.size()==1){
			VERTEX_TYPE parent=parents[0];
			VertexData*vData=m_data.get(parent);
			if(vData->getColor()==-1&&vData->getParents(parent,NULL,m_wordSize).size()<=1&&vData->getChildren(parent,m_wordSize).size()<=1){
				vData->setColor(color);
				numberOfVertices++;
				theQueue.push(parent);
			}
		}


		vector<VERTEX_TYPE> children=m_data.get(b)->getChildren(b,m_wordSize);
		if(children.size()==1){
			VERTEX_TYPE child=children[0];
			VertexData*vData=m_data.get(child);
			if(vData->getColor()==-1&&vData->getParents(child,NULL,m_wordSize).size()<=1&&vData->getChildren(child,m_wordSize).size()<=1){
				vData->setColor(color);
				numberOfVertices++;
				theQueue.push(child);
			}
		}
	}
	return numberOfVertices;
}


void DeBruijnAssembler::Walk_In_GRAPH(){
	cout<<endl;
	vector<VERTEX_TYPE>*theNodes=m_data.getNodes();
	(*m_cout)<<"********* Simplifying the graph"<<endl;
	cout<<endl;
	int color=0;
	for(int i=0;i<m_data.size();i++){
		if(i%100000==0){
			cout<<"Simplifying: "<<i<<" / "<<m_data.size()<<endl;
		}
		if(m_data.get((*theNodes)[i])->getColor()==-1){
			color++;
			int count=DFS_watch((*theNodes)[i],color);
			if(count>=100){
				//cout<<color<<" : "<<count<<endl;
			}
		}
	}
	cout<<"Simplifying: "<<m_data.size()<<" / "<<m_data.size()<<endl;

	cout<<endl;
	vector<VERTEX_TYPE> withoutParents;
	(*m_cout)<<"********* Inspecting the graph"<<endl;
	cout<<endl;
	map<int,map<int,int> > stats_parents_children;
	for(int i=0;i<m_data.size();i++){
		if(i%1000000==0){
			cout<<"Inspecting: "<<i<<" / "<<m_data.size()<<endl;
		}
		VERTEX_TYPE vertex=(*theNodes)[i];
		VertexData*dataNode=m_data.get(vertex);
		int parents=dataNode->getParents(vertex,NULL,m_wordSize).size();
		int children=dataNode->getChildren(vertex,m_wordSize).size();
		stats_parents_children[parents][children]++;
		if(parents==1&&children==1)
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
	(*m_cout)<<endl;
	(*m_cout)<<"********** Removing spurious edges"<<endl;
	bool removing=true;
	int spuriousRemoval=0;
	int id=0;
	cout<<endl;
	for(int myDataIterator=0;myDataIterator<m_data.size();myDataIterator++){
		break;
	//for(CustomMap<VertexData>::iterator i=m_data->begin();i!=m_data->end();i++){
		VERTEX_TYPE prefix=(*nodes)[(myDataIterator)];
		VertexData*nodeDataInstance=&(nodeData[myDataIterator]);
		if(nodeDataInstance->IsEliminated())
			continue;
		if(id%100000==0){
			cout<<id+1<<" / "<<m_data.size()<<endl;
		}
		vector<VERTEX_TYPE> theParents=nodeDataInstance->getParents(prefix,NULL,m_wordSize);
		//(*m_cout)<<theParents.size()<<endl;
		if(theParents.size()>1){
			//(*m_cout)<<theParents.size()<<endl;
			for(vector<VERTEX_TYPE>::iterator j=theParents.begin();j!=theParents.end();j++){
				if(m_data.get(*j)->IsEliminated())
					continue;
				int MaxDepth=30; // changed from 50 to 30 here
				VERTEX_TYPE currentNode=*j;
				set<VERTEX_TYPE> stuffVisited;
				int reachedDepth=visitVertices(currentNode,&stuffVisited,MaxDepth,true);
				if(reachedDepth<MaxDepth){
					//cout<<"Depth: "<<reachedDepth<<", "<<stuffVisited.size()<<" nodes"<<endl;
					m_data.get(*j)->eliminateNow();
					removing=true;
					spuriousRemoval++;
				}
			}
		}
		id++;
	}
	cout<<id+1<<" / "<<m_data.size()<<endl;
	(*m_cout)<<"Removed "<<spuriousRemoval<<" edges"<<endl;
	cout<<endl;

/*
	(*m_cout)<<"********** Collecting sources..."<<endl;
	cout<<endl;
	//for(MAP_TYPE<VERTEX_TYPE,MAP_TYPE<VERTEX_TYPE,vector<int> > >::iterator i=m_graph.begin();i!=m_graph.end();i++){
	for(int myDataIterator=0;myDataIterator<m_data.size();myDataIterator++){
		VERTEX_TYPE prefix=nodes[(myDataIterator)];
		
		if(nodeData[(myDataIterator)].IsEliminated())
			continue;
		vector<VERTEX_TYPE> theParents=nodeData[(myDataIterator)].getParents(prefix,&m_data);
		if(theParents.size()==0){
			withoutParents.push_back(prefix);
		}
	}
*/
	(*m_cout)<<"Done..., "<<withoutParents.size()<<" sources."<<endl;
	
	(*m_cout)<<endl;
	ostream&m_cout=*(this->m_cout);	
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
		m_cout<<endl;
		vector<VERTEX_TYPE> newSources;
		for(int i=0;i<(int)sources.size();i++){
			VERTEX_TYPE prefix=sources[i];
			if(m_data.get(prefix)->IsAssembled())
				continue;
			m_cout<<endl;
			m_cout<<"Source "<<i+1<<" / "<<sources.size()<<endl;//" REPEAT MODE"<<endl;
			m_cout<<"From: "<<idToWord(prefix,m_wordSize)<<endl;
		//m_cout<<m_contig_paths.size()<<" contigs"<<endl;
			vector<VERTEX_TYPE>path;
			vector<map<int,map<char,int> > > currentReadPositions;
			vector<int> repeatAnnotations;
			vector<VERTEX_TYPE> localNewSources;
			contig_From_SINGLE(&currentReadPositions,&path,&localNewSources,&repeatAnnotations,prefix);
			m_cout<<path.size()<<" vertices"<<endl;
/*
			int validSources=0;
			for(vector<VERTEX_TYPE>::iterator k=localNewSources.begin();k!=localNewSources.end();k++){
				if(m_data.get(*k)->IsAssembled())
					continue;
				// assemble vertices cannot be sources..
				validSources++;
				newSources.push_back(*k);
			}
			cout<<localNewSources.size()<<" new sources, "<<validSources<<" valid."<<endl;a
*/
			if(path.size()<=2)
				continue;
			writeContig_Amos(&currentReadPositions,&path,&amosFile,contigId);
			writeContig_fasta(&path,&contigsFileStream,contigId);
			writeContig_Coverage(&currentReadPositions,&path,&coverageStream,contigId);
			writeContig_RepeatAnnotation(&repeatAnnotations,contigId,&repeatAnnotation,&path);
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
	m_cout<<endl;

/*
	(m_cout)<<"********** Sources discovery..."<<endl;
	(m_cout)<<"Iteration Sources Cumulative"<<endl;
	for(int i=0;i<(int)The_Discovery_Of_Sources.size();i++)
		(m_cout)<<i+1<<" "<<The_Discovery_Of_Sources[i]<<" "<<VisitsOfSources[i]<<endl;

	(m_cout)<<endl;
*/
	return;
}

void DeBruijnAssembler::Algorithm_Assembler_20090121(){
	Walk_In_GRAPH();
}





bool DeBruijnAssembler::DETECT_BUBBLE(vector<VERTEX_TYPE>*path,VERTEX_TYPE a,VERTEX_TYPE b){
	return false;
	set<VERTEX_TYPE> stuffFromA;
	set<VERTEX_TYPE> stuffFromB;
	int maxDepth=m_wordSize+1;
	visitVertices(a,&stuffFromA,maxDepth,false);
	visitVertices(b,&stuffFromB,maxDepth,false);
	for(set<VERTEX_TYPE>::iterator i=stuffFromA.begin();i!=stuffFromA.end();i++){
		if(stuffFromB.count(*i)>0){
			cout<<"Bubble!"<<endl;
			return true;
		}
	}
	return false;
}



void DeBruijnAssembler::contig_From_SINGLE(vector<map<int,map<char,int> > >*currentReadPositions,vector<VERTEX_TYPE>*path,vector<VERTEX_TYPE>*newSources,vector<int>*repeatAnnotations,VERTEX_TYPE source){
	map<int,int> usedReads;
	map<int,int> firstReadOccurance;
	bool debug_print=m_DEBUG;
	//debug_print=true;
	//(*m_cout)<<"Depth: "<<path->size()<<endl;
	int lowCoverageLength=0;
	VERTEX_TYPE prefix;
	vector<VERTEX_TYPE> prefixNextVertices;
	prefixNextVertices.push_back(source);
	bool isSpurious=false;
	VERTEX_TYPE otherSource;
	set<VERTEX_TYPE> notSpurious;
	int otherPosition;
	set<int> reads_that_can_start_After_repeat;
	while(prefixNextVertices.size()==1){
		prefix=prefixNextVertices[0];
		if(m_data.get(prefix)->IsEliminated()==true){
			cout<<"Skipping spurious edge"<<endl;
			return;
		}
		VertexData*prefixVertexData=m_data.get(prefix);

		if(m_data.get(prefix)->IsAssembled()&&path->size()<50){
			cout<<"Skipping spurious source"<<endl;
			return;
		}
/*
		if(isSpurious==false&&prefixVertexData->IsAssembled()){
			map<VERTEX_TYPE,vector<int> >*thePositions=prefixVertexData->getPositions();
			isSpurious=true;
			otherSource=thePositions->begin()->first;
			otherPosition=thePositions->begin()->second[0];
			if(otherPosition==0){
				notSpurious.insert(otherSource);
				isSpurious=false;
			}
		}else if(isSpurious==true&&path->size()<=1000&&notSpurious.count(otherSource)==0){
			map<VERTEX_TYPE,vector<int> >*thePositions=prefixVertexData->getPositions();
			//cout<<"Might be spurious"<<endl;
			
			if((*thePositions).count(otherSource)>0&&(*thePositions)[otherSource][0]-800==otherPosition){
				cout<<"Halting, spurious vertex detected."<<endl;
				break;
			}
		}
*/
		prefixVertexData->addPositionInContig(source,path->size());
		//cout<<DeBruijnAssembler::idToWord(prefix,DeBruijnAssembler::m_wordSize)<<endl;
		//prefixVertexData->printPositions();
		//cout<<"Position:" <<path->size()<<" "<<idToWord(prefix,m_wordSize)<<endl;
		path->push_back(prefix);
		map<int,map<char,int> > a;
		(*currentReadPositions).push_back(a);
		//(*m_cout)<<"Pushing "<<idToWord(prefix,m_wordSize)<<endl;
		int added=0;

		vector<AnnotationElement>*annotations=m_data.get(prefix)->getAnnotations();
		if(debug_print&&annotations!=NULL){
			(*m_cout)<<idToWord(path->at(path->size()-1),m_wordSize)<<" ";
			(*m_cout)<<annotations->size();
			for(int j=0;j<annotations->size();j++){
				(*m_cout)<<" "<<m_sequenceData->at(annotations->at(j).readId)->getId()<<" "<<annotations->at(j).readPosition<<" "<<annotations->at(j).readStrand;
			}
			(*m_cout)<<endl;
		}



		//(*m_cout)<<"Path position "<<path->size()-1<<endl;
		if(path->size()%1000==0)
			(*m_cout)<<path->size()<<" progress."<<endl;
		//(*m_cout)<<"Threading.. reads "<<endl;
		//(*m_cout)<<"Coverage mean: "<<coverageMean<<" "<<annotations->size()<<endl; 
		if(annotations->size()>=m_REPEAT_DETECTION&&m_DEBUG){
			(*m_cout)<<"Coverage: "<<annotations->size()<<", refusing to start threading reads!"<<endl;
		}
		if(annotations->size()>=m_REPEAT_DETECTION){
			repeatAnnotations->push_back(1);
		}else{
			repeatAnnotations->push_back(0);
		}
		for(int h=0;h<(int)annotations->size();h++){
			bool addedInContig=false;
			//(*m_cout)<<m_sequenceData->at(annotations->at(h).readId)->getId()<<" "<<annotations->at(h).readPosition<<" "<<annotations->at(h).readStrand;
			if(usedReads.count(annotations->at(h).readId)>0){
				//cout<<" LastSeen="<<(*currentReadPositions)[usedReads[annotations->at(h).readId]][annotations->at(h).readId][annotations->at(h).readStrand]<<" "<<usedReads[annotations->at(h).readId];
			}
			if(firstReadOccurance.count(annotations->at(h).readId)==0){
				firstReadOccurance[annotations->at(h).readId]=annotations->at(h).readPosition;
			}
			if(usedReads.count(annotations->at(h).readId)==0){
					// add a read when it starts at its beginning...
				if(reads_that_can_start_After_repeat.count(annotations->at(h).readId)>0||
 			((annotations->at(h).readStrand=='F'&&annotations->at(h).readPosition==m_sequenceData->at(annotations->at(h).readId)->getStartForward())||
			(annotations->at(h).readStrand=='R'&&annotations->at(h).readPosition==m_sequenceData->at(annotations->at(h).readId)->getStartReverse())||
			path->size()<50
		//annotations->at(h).readPosition<5 // ||
		//(firstReadOccurance.count(annotations->at(h).readId)>0&&annotations->at(h).readPosition-firstReadOccurance[annotations->at(h).readId]>50)
			
				)){ // add at most a given amount of "new reads" to avoid depletion
					if(annotations->size()<m_REPEAT_DETECTION){
						(*currentReadPositions)[path->size()-1][annotations->at(h).readId][annotations->at(h).readStrand]=annotations->at(h).readPosition; // = 0
						//(*m_cout)<<path->size()<<" "<<idToWord(path->at(path->size()-2),m_wordSize)<<" -> "<<idToWord(path->at(path->size()-1),m_wordSize)<<endl;
						//(*m_cout)<<"Adding read "<<m_sequenceData->at(annotations->at(h).readId)->getId()<<" "<<annotations->at(h).readStrand<<" "<<annotations->at(h).readPosition<<endl;	
						//cout<<" Added ";
						added++;
						usedReads[(annotations->at(h).readId)]=path->size()-1;
						addedInContig=true;
					}else{
						reads_that_can_start_After_repeat.insert(annotations->at(h).readId);
					}
				}
			}else if(is_d_Threading(&(annotations->at(h)),currentReadPositions,path,&usedReads,false)){
				added++;
				if(m_carry_forward_offset==0||true){
					//cout<<" Threaded ";
					usedReads[annotations->at(h).readId]=path->size()-1;
					addedInContig=true;
					(*currentReadPositions)[path->size()-1][annotations->at(h).readId][annotations->at(h).readStrand]=annotations->at(h).readPosition;
					//cout<<"Threading "<<m_sequenceData->at(annotations->at(h).readId)->getId()<<endl;
				}
			}
			if(addedInContig){
				//cout<<" [IN CONTIG]";
			}
			//cout<<endl;
		}

		if(debug_print){
			(*m_cout)<<(*currentReadPositions)[path->size()-1].size();
			for(map<int,map<char,int> >::iterator i=(*currentReadPositions)[path->size()-1].begin();i!=(*currentReadPositions)[path->size()-1].end();i++){
				for(map<char,int>::iterator j=i->second.begin();j!=i->second.end();j++){
					(*m_cout)<<" "<<m_sequenceData->at(i->first)->getId()<<" "<<j->first<<" "<<j->second;
				}
			}
			(*m_cout)<<endl;
		}


		prefixNextVertices=nextVertices(path,currentReadPositions,newSources,&usedReads);
		if(added==0&&m_data.get(prefix)->hasManyChildren(prefix,m_wordSize)){
			(*m_cout)<<"Stop!, reason: No read threaded, "<<endl;
			break;
		}
		if(added>1)
			lowCoverageLength=0;
		if(added==1&&path->size()>400){
			lowCoverageLength++;
		}
		if(lowCoverageLength>30&&m_data.get(prefix)->NotTrivial(prefix,m_wordSize)){
			cout<<"Stop!, reason: Extensive 1-Coverage."<<endl;
			break;
		}
		//(*m_cout)<<"adding "<<(*currentReadPositions).at(currentReadPositions->size()-1).size()<<endl;
	}

	vector<VERTEX_TYPE>children=m_data.get(prefix)->getChildren(prefix,m_wordSize);

	if(children.size()==0){
		cout<<"Stop!, reason: this is a sink, no data available beyond this very nucleotide."<<endl;
	}else if(children.size()>1){
		cout<<"Stop!, reason: no choice possible (repeated region?)."<<endl;
	}else if(children.size()==1){
		cout<<"1 children (?)"<<endl;
	}
	(*m_cout)<<"Prefix "<<idToWord(prefix,m_wordSize)<<" "<<children.size()<<endl;
	// add newSources
	/*
	for(int j=0;j<(int)children.size();j++){
		VERTEX_TYPE vertexValue=children[j];
		(*m_cout)<<"Adding "<<idToWord(children[j],m_wordSize)<<endl;
		newSources->push_back(children[j]);
	}
	*/
	if(!debug_print&&false){
		return;
	}
	//return;
	if(children.size()>0&&path->size()>=1){
	/*
		vector<AnnotationElement>*annotations=m_data->get(path->at(path->size()-2)).getAnnotations(path->at(path->size()-1));
		map<int,int>allowedReads;
		for(int h=0;h<annotations->size();h++)
			allowedReads[(annotations->at(h).readId)]=annotations->at(h).readPosition;
*/
		if(path->size()>0){
			//(*m_cout)<<pathToDNA(path)<<endl;
			cout<<"Contig annotations"<<endl;
			for(int position=0;position<currentReadPositions->size();position++) {
			
				if(position<path->size()-500){
					//continue;
				}
				VERTEX_TYPE aVertex=path->at(position);
				cout<<"Position "<<position<<" "<<idToWord(aVertex,m_wordSize)<<endl;

				vector<AnnotationElement>*annotations=m_data.get(aVertex)->getAnnotations();
				for(int j=0;j<annotations->size();j++){
					(*m_cout)<<m_sequenceData->at(annotations->at(j).readId)->getId()<<" "<<annotations->at(j).readPosition<<" "<<annotations->at(j).readStrand;
					if((*currentReadPositions)[position].count(annotations->at(j).readId)>0&&
					(*currentReadPositions)[position][annotations->at(j).readId].count(annotations->at(j).readStrand)>0&&
					(*currentReadPositions)[position][annotations->at(j).readId][annotations->at(j).readStrand]==annotations->at(j).readPosition){
						cout<<" [IN CONTIG]";
					}
					cout<<endl;
				}
			}
			//cout<<
			(*m_cout)<<children.size()<<" Choices"<<endl;
			for(int i=0;i<children.size();i++){
				//int score=recThreading(path->at(path->size()-1),children[i],&allowedReads);
				//(*m_cout)<<"Score: "<<score<<endl;
				(*m_cout)<<idToWord(children[i],m_wordSize)<<" ";
				vector<AnnotationElement>*annotations=m_data.get(children[i])->getAnnotations();
				(*m_cout)<<annotations->size();
				for(int j=0;j<annotations->size();j++){
					(*m_cout)<<" "<<m_sequenceData->at(annotations->at(j).readId)->getId()<<" "<<annotations->at(j).readPosition<<" "<<annotations->at(j).readStrand;
				}
				(*m_cout)<<endl;
			}
		}
	}
}




/**
 * \param simple  force passFilter to use passFilter_ShortRead
 */
vector<VERTEX_TYPE> DeBruijnAssembler::nextVertices(vector<VERTEX_TYPE>*path,vector<map<int,map<char,int> > >*currentReadPositions,vector<VERTEX_TYPE>*newSources,map<int,int>*usedReads){
	bool debugPrint=m_DEBUG;
	VERTEX_TYPE prefix=path->at(path->size()-1);
	vector<VERTEX_TYPE> children=m_data.get(prefix)->getChildren(prefix,m_wordSize);
	if(children.size()==1)
		return children;
	// start when nothing is done yet
	//(*m_cout)<<currentReadPositions->size()<<" "<<path->size()<<endl;
	if(currentReadPositions->size()==0){//||currentReadPositions->size()<path->size())
		cout<<"Nothing here"<<endl;
		if(children.size()==2&&DETECT_BUBBLE(path,children[0],children[1])){
			vector<VERTEX_TYPE> output;
			cout<<"Early Bubble."<<endl;
			output.push_back(children[0]);
			return output;
		}
		return children;
	}
	
	if(children.size()==1){
		VERTEX_TYPE first=children[0];
		int colorToHave=m_data.get(prefix)->getColor();
		if(m_data.get(first)->getColor()==colorToHave){
			return children;
		}
	}


	if(children.size()>1&&m_DEBUG){
		(*m_cout)<<"More than 1 children: "<<children.size()<<endl;
	}

/*
	for(map<int,map<int,map<char,int> > >::iterator k=currentReadPositions->begin();k!=currentReadPositions->end();k++){
		for(map<int,map<char,int> >::iterator i=k->second.begin();i!=k->second.end();i++){
			for(map<char,int>::iterator j=i->second.begin();j!=i->second.end();j++){
				//(*m_cout)<<k->first<<" "<<m_sequenceData->at(i->first)->getId()<<" "<<j->first<<" "<<j->second<<endl;
			}
		}
	}
*/
	//(*m_cout)<<"previous position to check "<<path->size()-2<<endl;

	map<VERTEX_TYPE,int > scoresSum;
	map<VERTEX_TYPE,int> numbers;
	map<VERTEX_TYPE,int> scoresMax;
	map<VERTEX_TYPE,int>  coverageOfEdges;
	map<VERTEX_TYPE,vector<int> > valuesForEdge;
	ostringstream debugBuffer;
	//(*m_cout)<<"Children"<<endl;
	for(vector<VERTEX_TYPE>::iterator i=children.begin();i!=children.end();i++){

		//(*m_cout)<<"Child "<<idToWord(children[i],m_wordSize)<<endl;
		vector<AnnotationElement>*thisEdgeData=(m_data.get(*i)->getAnnotations());
		coverageOfEdges[*i]=0;
		if(thisEdgeData!=NULL){
			coverageOfEdges[*i]=thisEdgeData->size();
			for(int j=0;j<(int)thisEdgeData->size();j++){
				uint32_t readId=thisEdgeData->at(j).readId;

				if(is_d_Threading(&(thisEdgeData->at(j)),currentReadPositions,path,usedReads,true)){
					if(thisEdgeData->at(j).readPosition>=scoresMax[*i]){
						scoresMax[*i]=thisEdgeData->at(j).readPosition+m_carry_forward_offset;
					}

					scoresSum[*i]+=thisEdgeData->at(j).readPosition+m_carry_forward_offset;
					numbers[*i]++;
					valuesForEdge[*i].push_back(thisEdgeData->at(j).readPosition+m_carry_forward_offset);
				}
			}
		}
	}

	if(m_DEBUG)
		cout<<"Debug scores"<<endl;
	for(map<VERTEX_TYPE,int>::iterator i=scoresSum.begin();i!=scoresSum.end();i++){
		if(m_DEBUG==false)
			break;
		vector<int> values=valuesForEdge[i->first];
		for(vector<int>::iterator j=values.begin();j!=values.end();j++){
			cout<<" "<<*j;
		}
		cout<<endl;
	}


	VERTEX_TYPE best=0;
	bool foundBest=false;

	// best Sum
	for(map<VERTEX_TYPE,int>::iterator i=scoresSum.begin();i!=scoresSum.end();i++){
		if(foundBest==true)
			break;
		bool isBest=true;
		for(map<VERTEX_TYPE,int>::iterator j=scoresSum.begin();j!=scoresSum.end();j++){
			if(i->first==j->first)
				continue;
			if(!(i->second > 1.1*j->second)){
				isBest=false;
			}
		}
		if(isBest){
			foundBest=true;
			best=i->first;
		}
	}

	// best max
	for(map<VERTEX_TYPE,int>::iterator i=scoresMax.begin();i!=scoresMax.end();i++){
		if(foundBest==true)
			break;
		bool isBest=true;
		for(map<VERTEX_TYPE,int>::iterator j=scoresMax.begin();j!=scoresMax.end();j++){
			if(i->first==j->first)
				continue;
			if(!(i->second > 1.01*j->second)){
				isBest=false;
			}
		}
		if(isBest){
			foundBest=true;
			best=i->first;
		}
	}

	// best number
	for(map<VERTEX_TYPE,int>::iterator i=numbers.begin();i!=numbers.end();i++){
		if(foundBest==true)
			break;
		bool isBest=true;
		for(map<VERTEX_TYPE,int>::iterator j=numbers.begin();j!=numbers.end();j++){
			if(i->first==j->first)
				continue;
			if(!(i->second > 1.6*j->second)){
				isBest=false;
			}
		}
		if(isBest){
			foundBest=true;
			best=i->first;
		}
	}

	if(children.size()==2&&
		foundBest==false){
		// attempt to detect homopolymer...
		cout<<"Searching for homopolymer"<<endl;
		string firstOne=idToWord(children[0],m_wordSize);
		string secondOne=idToWord(children[1],m_wordSize);
		int trailingHomoPolymerSize=0;
		while(
m_wordSize-2-trailingHomoPolymerSize>0 &&
firstOne[m_wordSize-2-trailingHomoPolymerSize]==
			secondOne[m_wordSize-2-trailingHomoPolymerSize]&&
			firstOne[m_wordSize-2]==firstOne[m_wordSize-2-trailingHomoPolymerSize]
			){
			trailingHomoPolymerSize++;
		}
		if(trailingHomoPolymerSize>1){
			cout<<"Trailing homopolymer W00t "<<trailingHomoPolymerSize<<endl;
			cout<<firstOne<<endl;
			cout<<secondOne<<endl;
			foundBest=true;
			if(firstOne[m_wordSize-1]==firstOne[m_wordSize-2]){ // take the shortest one.
				best=children[1];
			}else{
				best=children[0];
			}
		}
	}

/*
	if(foundBest==false){
		for(int i=0;i<children.size();i++){
			if(newSources!=NULL){
				VERTEX_TYPE dataVertex=children[i];
				int nParents=m_data.get(dataVertex)->getParents(dataVertex,NULL).size();
				if(nParents>1||
				(m_data.get(prefix)->getAnnotations(dataVertex)!=NULL&&m_data.get(prefix)->getAnnotations(dataVertex)->size()>=m_REPEAT_DETECTION))
					continue;
				newSources->push_back(dataVertex);
				(*m_cout)<<"Adding alternative source: "<<idToWord(children[i],m_wordSize)<<", "<<nParents<<" parents"<<endl;
			}
		}
	}else{
		for(int i=0;i<children.size();i++){
			if(children[i]!=best&&newSources!=NULL&&coverageOfEdges[children[i]]<m_REPEAT_DETECTION&&
				!DETECT_BUBBLE(path,children[i],best)){
				VERTEX_TYPE dataVertex=children[i];
				int nParents=m_data.get(dataVertex)->getParents(dataVertex,NULL).size();
				if(nParents>1||
				(m_data.get(prefix)->getAnnotations(dataVertex)!=NULL&&m_data.get(prefix)->getAnnotations(dataVertex)->size()>=m_REPEAT_DETECTION))
					continue;
				newSources->push_back(children[i]);
				(*m_cout)<<"Adding alternative source: "<<idToWord(children[i],m_wordSize)<<", "<<nParents<<" parents"<<endl;
			}
		}

	}
*/
	if(foundBest==true){
		vector<VERTEX_TYPE> output;
		output.push_back(best);
		return output;
	}
	//(*m_cout)<<debugBuffer.str()<<endl;
	vector<VERTEX_TYPE> output;

	if(newSources!=NULL&&m_DEBUG){
		(*m_cout)<<"No children scored. "<<endl;
		vector<VERTEX_TYPE> allChildren=m_data.get(path->at(path->size()-1))->getChildren(path->at(path->size()-1),m_wordSize);
		(*m_cout)<<"Before filtering "<<allChildren.size()<<endl;
		for(int i=0;i<allChildren.size();i++)
			(*m_cout)<<idToWord(allChildren[i],m_wordSize)<<endl;

		(*m_cout)<<"After filtering "<<children.size()<<endl;
		for(int i=0;i<children.size();i++)
			(*m_cout)<<idToWord(children[i],m_wordSize)<<endl;
	}
	return output;
}






DeBruijnAssembler::~DeBruijnAssembler(){
}



VERTEX_TYPE DeBruijnAssembler::reverseComplement_VERTEX(VERTEX_TYPE a){
	return wordId(reverseComplement(idToWord(a,m_wordSize)).c_str());
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
	int maxSize=sequence.length()-m_wordSize;
	char*myWord=(char*)malloc(m_wordSize+2);
	for(int readPosition=0;readPosition<(int)maxSize;readPosition++){
		for(int i=0;i<m_wordSize;i++){
			myWord[i]=sequence[readPosition+i];
		}
		myWord[m_wordSize]='\0';
		VERTEX_TYPE wordInBits=wordId(myWord);
		if(read->isValidDNA(myWord)
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
			vector<VERTEX_TYPE> theParents=m_data.get(wordInBits)->getParents(wordInBits,NULL,m_wordSize);
			for(vector<VERTEX_TYPE>::iterator ii=theParents.begin();ii!=theParents.end();ii++){
				vector<VERTEX_TYPE> theChildren=m_data.get(*ii)->getChildren(*ii,m_wordSize);
				if(theChildren.size()>1)
					hasMixedParents=true;
			}
			if((strand=='F'&&readPosition==sequenceData->at(readId)->getStartForward())||
			(strand=='R'&&readPosition==sequenceData->at(readId)->getStartForward())||false)
				m_data.get(wordInBits)->addAnnotation(readId,readPosition,strand);
		}
	}
	free(myWord);
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
			int edgeOffset=vertexOffset-1;
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
	//(*m_cout)<<"writeContig_Amos"<<endl;
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

	//(*m_cout)<<"amos generator "<<currentReadPositions->size()<<endl;
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

bool DeBruijnAssembler::is_d_Threading(AnnotationElement*annotation,vector<map<int,map<char,int> > >*currentReadPositions,vector<VERTEX_TYPE>*path,map<int,int>*usedReads,bool beforeAdding){

	int readId=annotation->readId;
	if(usedReads->count(readId)==0){
		return false;
	}
	int lastPosition=(*usedReads)[readId];
	char readStrand=annotation->readStrand;
	if((*currentReadPositions)[lastPosition][readId].count(readStrand)==0){
		return false;
	}
	int lastPositionInRead=(*currentReadPositions)[lastPosition][readId][readStrand];
	int distanceInRead=annotation->readPosition-lastPositionInRead;
	int distanceInPath=path->size()-1-lastPosition;
/*
	string sequence=m_sequenceData->at(readId)->getSeq();
	if(readStrand=='R')
		sequence=reverseComplement(sequence);

	VERTEX_TYPE aNodeInRead=wordId(sequence.substr(lastPositionInRead+distanceInRead,m_wordSize).c_str());
	VERTEX_TYPE lastNodeInPath=path->at(path->size()-1);
	cout<<idToWord(aNodeInRead,m_wordSize)<<" "<<idToWord(lastNodeInPath,m_wordSize)<<endl;
	if(aNodeInRead==lastNodeInPath)
		return true;
*/

	if(beforeAdding) //&&lastPositionInRead==annotation->readPosition-1)
		distanceInPath++;
	/*
	if(lastPositionInRead!=annotation->readPosition-1)
		distanceInRead--;
	*/
	//(*m_cout)<<"R "<<distanceInRead<<" C "<<distanceInPath<<endl;
	m_carry_forward_offset=0;
	//return distanceInPath==distanceInRead;
	//cout<<"DistanceInPath "<<distanceInPath<<" DistanceInRead "<<distanceInRead;
	if(distanceInPath==0)
		return false;
	if(distanceInPath==distanceInRead){
		return true;
	}
	for(int i=0;i<=0;i++){
		m_carry_forward_offset=i;
		if(distanceInRead+i==distanceInPath)
			return true;
	}
	return false;
}






// DFS (?) accumulate 
int DeBruijnAssembler::visitVertices(VERTEX_TYPE a,set<VERTEX_TYPE>*nodes,int maxDepth,bool parents){
	queue<VERTEX_TYPE> dataQueue;
	queue<VERTEX_TYPE> depthQueue;
	dataQueue.push(a);
	depthQueue.push(0);
	int bestDepth=0;
	while(dataQueue.size()>0){
		VERTEX_TYPE aNode=dataQueue.front();
		dataQueue.pop();
		int aNodeDepth=depthQueue.front();
		depthQueue.pop();

		if(aNodeDepth>bestDepth)
			bestDepth=aNodeDepth;

		nodes->insert(aNode);
		vector<VERTEX_TYPE> elements;
		if(parents==true)
			elements=m_data.get(aNode)->getParents(aNode,NULL,m_wordSize);
		else
			elements=m_data.get(aNode)->getChildren(aNode,m_wordSize);

		for(vector<VERTEX_TYPE>::iterator i=elements.begin();i!=elements.end();i++){
			VERTEX_TYPE otherNode=*i;
			if(aNodeDepth+1>maxDepth)
				continue;
			if(nodes->count(otherNode)>0)
				continue;
			dataQueue.push(otherNode);
			depthQueue.push(aNodeDepth+1);
		}
	}
	return bestDepth;
}
