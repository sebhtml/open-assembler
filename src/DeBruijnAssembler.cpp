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


#include<cmath>
#include<cstdlib>
#include<map>
#include<iostream>
#include"DeBruijnAssembler.h"
#include<fstream>
#include<vector>
#include<string>
#include<vector>
#include"CustomMap.hpp"
using namespace std;




int min(int a,int b){
	if(a<b)
		return a;
	return b;
}

int abs(int a){
	if(a<0)
		return -a;
	return a;
}


int DeBruijnAssembler::m_WordSize=5;
int DeBruijnAssembler::m_NUCLEOTIDE_A=0;
int DeBruijnAssembler::m_NUCLEOTIDE_T=1;
int DeBruijnAssembler::m_NUCLEOTIDE_C=2;
int DeBruijnAssembler::m_NUCLEOTIDE_G=3;

DeBruijnAssembler::DeBruijnAssembler(ostream*m_cout){
	m_DEBUG=false;
	this->m_cout=m_cout;
}

void DeBruijnAssembler::setWordSize(int k){
	m_Coverage_From_DepletionCurve=5;
	m_wordSize=k;
	m_threshold=0.015;
	m_pairedAvailable=false;
	DeBruijnAssembler::m_WordSize=k;
	m_minimumCoverage=5;
}



char DeBruijnAssembler::complement(char a){
	if(a=='A')
		return 'T';
	if(a=='T')
		return 'A';
	if(a=='C')
		return 'G';
	if(a=='G')
		return 'C';
	return a;
}

string DeBruijnAssembler::reverseComplement(string a){
	ostringstream i;
	for(int p=a.length()-1;p>=0;p--){
		char b=complement(a[p]);
		i<< b;
	}
	return i.str();
}


void DeBruijnAssembler::build_From_Scratch(SequenceDataFull*sequenceData){
	ostream&m_cout=*(this->m_cout);
	m_cout<<endl;
	m_cout<<endl;
	m_cout<<"Collecting mers from reads"<<endl;
	m_cout<<"k+1 = "<<m_wordSize+1<<endl;
	CustomMap<int> words(m_buckets);
	int last_vertices_size=-1;
	for(int i=0;i<(int)sequenceData->size();i++){
		if(i%100000==0){
			m_cout<<"Reads: "<<i<<" / "<<sequenceData->size()<<endl;
		}
		vector<VERTEX_TYPE> highQualityMers=sequenceData->at(i)->getHighQualityMers(m_wordSize);
		for(vector<VERTEX_TYPE>::iterator iteratorMer=highQualityMers.begin();iteratorMer!=highQualityMers.end();iteratorMer++){
			if(!words.find(*iteratorMer))
				words.add(*iteratorMer,0);

			words.set(*iteratorMer,words.get(*iteratorMer)+1);
			if((int)words.size()%1000000==0&&(int)words.size()!=last_vertices_size){
				m_cout<<"Mers: "<<words.size()<<endl;
				last_vertices_size=words.size();
			}
		}
	}	

	m_cout<<"Reads: "<<sequenceData->size()<<" / "<<sequenceData->size()<<endl;
	m_cout<<"Mers: "<<words.size()<<endl;


	int processed=0;
	int solid=0;

	m_cout<<endl;

	(m_cout)<<"Quality analysis (Depletion curve), k+1 = "<<m_wordSize+1<<endl;
	(m_cout)<<"MinimumCoverage NumberOfSolidMers NormalizedNumberOfSolidMers Change"<<endl;
	map<int,int> solidCounts;
	for(int minimumCoverage=1;minimumCoverage<=32;minimumCoverage++){
		solidCounts[minimumCoverage]=0;
	}
	for(CustomMap<int>::iterator i=words.begin();i!=words.end();i++){
		for(int minimumCoverage=1;minimumCoverage<=32;minimumCoverage++){
			if(i.second()>=minimumCoverage){
				solidCounts[minimumCoverage]++;
			}
		}
	}

	map<int,double> solidRatio;
	
	for(int minimumCoverage=1;minimumCoverage<=32;minimumCoverage++){
		solidRatio[minimumCoverage]=solidCounts[minimumCoverage]/(0.0+solidCounts[1]);
	}
	map<int,double> solidDiff;
	for(int minimumCoverage=2;minimumCoverage<=32;minimumCoverage++){
		solidDiff[minimumCoverage]=solidRatio[minimumCoverage-1]-solidRatio[minimumCoverage];
	}

	m_minimumCoverage=5;
	double minimum=1;
	for(int minimumCoverage=1;minimumCoverage<=32;minimumCoverage++){
		m_cout<<minimumCoverage<<" "<<solidCounts[minimumCoverage]<<" "<<solidRatio[minimumCoverage];
		if(minimumCoverage!=1)
			m_cout<<" "<<solidDiff[minimumCoverage];
		m_cout<<endl;
	}
	m_cout<<endl;
	double cutOff=m_threshold;
	if(m_minimumCoverageParameter=="auto")
		m_cout<<"Cutoff: "<<cutOff<<endl;
	for(int minimumCoverage=2;minimumCoverage<=32;minimumCoverage++){
		if(solidDiff[minimumCoverage]<minimum){
			m_minimumCoverage=minimumCoverage;
			minimum=solidDiff[minimumCoverage];
		}
		if(solidDiff[minimumCoverage]>minimum)
			break;
		if(solidDiff[minimumCoverage]<cutOff){
			m_minimumCoverage=minimumCoverage;
			if(m_minimumCoverageParameter=="auto")
				m_cout<<"Best Coverage <- "<<m_minimumCoverage<<endl;
			break;
		}
	}

	m_Coverage_From_DepletionCurve=m_minimumCoverage;
	if(m_minimumCoverageParameter=="auto"){
		m_cout<<"Using depletion curve (-minimumCoverage auto)"<<endl;
		m_cout<<"m_minimumCoverage <- "<<m_minimumCoverage<<endl;
	}else{
		//m_cout<<""<<endl;
		m_cout<<"m_minimumCoverage <- "<<atoi(m_minimumCoverageParameter.c_str())<<endl;
		m_minimumCoverage=atoi(m_minimumCoverageParameter.c_str());
	}
	m_cout<<endl;

	uint64_t total_bases=0;
	uint64_t solid_bases=0;
	
	for(CustomMap<int>::iterator i=words.begin();i!=words.end();i++){
		if(i.second()>=m_minimumCoverage)
			solid_bases+=i.second()*(m_wordSize+1);
		total_bases+=i.second()*(m_wordSize+1);
	}

	m_cout<<"Solid mers: (weighted) "<<(solid_bases+0.0)/total_bases*100<<"%"<<endl;
	solid=0;
	CustomMap<int> solidMers(m_buckets);
	for(CustomMap<int>::iterator i=words.begin();i!=words.end();i++){
		processed++;
		if(i.second()>=m_minimumCoverage){
			solid++;
			VERTEX_TYPE w=i.first();
			solidMers.add(w,1);
		}
	}
	m_solidMers=solidMers.size();
	m_cout<<"Solid mers: "<<solid<<" / "<<words.size()<<" ---> "<<((solid+0.0)/words.size()+0.0)*100.0<<"%"<<endl;
	//m_cout<<" (this should be roughly twice the genome size)"<<endl;
	m_cout<<"Not-solid mers: "<<processed-solid<<" / "<<words.size()<<" ---> "<<(processed-solid+0.0)/words.size()*100.0<<"%"<<endl;

	if(m_solidMers==0){
		m_cout<<"Error: mers are depleted..."<<endl;
		exit(0);
	}
	m_cout<<endl;
	words.clear();

	// Don't load too much reads in edges, MAX: 15?
	m_cout<<"Indexing solid mers in reads, building graph with solid mers."<<endl; // <-------
	for(int readId=0;readId<(int)sequenceData->size();readId++){
		if(readId%10000==0)
			m_cout<<"Reads: "<<readId<<" / "<<sequenceData->size()<<endl;
		indexReadStrand(readId,'F',sequenceData,&solidMers);
		indexReadStrand(readId,'R',sequenceData,&solidMers);
	
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
	m_cout<<"Loading paired information"<<endl;
	m_pairedAvailable=true;
	int entries;
	f>>entries;
	for(int i=0;i<entries;i++){
		string file1;
		string file2;
		int distance;
		f>>file1>>file2>>distance;
		if(!sequenceData->hasFile(file1)||!sequenceData->hasFile(file2))
			continue;
		int start1=sequenceData->getFirst(file1);
		int start2=sequenceData->getFirst(file2);
		int last1=sequenceData->getLast(file1);
		while(start1<=last1){
			string read1=sequenceData->at(start1)->getSeq();
			string read2=sequenceData->at(start2)->getSeq();
			start1++;
			start2++;
			string word1=read1.substr(0,m_wordSize+1);
			string word2=read2.substr(0,m_wordSize+1);
			bool good=true;
			//no N
			for(int k=0;k<m_wordSize+1;k++){
				if(word1[k]=='N'||word2[k]=='N'||word1[k]=='.'||word2[k]=='.'){
					good=false;
					break;
				}
			}
			if(good==false)
				continue;
			VERTEX_TYPE word1_foward=wordId(word1.substr(0,m_wordSize).c_str());
			VERTEX_TYPE word2_foward=wordId(word2.substr(0,m_wordSize).c_str());
			if(solidMers.find(word1_foward)&&solidMers.find(word2_foward)){
				VERTEX_TYPE word1_reverse=wordId(reverseComplement(word1.substr(0,m_wordSize)).c_str());
				VERTEX_TYPE word2_reverse=wordId(reverseComplement(word2.substr(0,m_wordSize)).c_str());
				m_data->get(word2_foward).addPaired(word1_foward,distance);
				m_data->get(word1_reverse).addPaired(word2_reverse,distance);
			}
		}
	}
	f.close();


}

void DeBruijnAssembler::writeGraph(){
	//(*m_cout)<<"Writing graph."<<endl;
	ofstream graph(m_graphFile.c_str());
	string humanReadable=m_assemblyDirectory+"/Graph.txt";
	ofstream graph2(humanReadable.c_str());
	graph<<m_solidMers<<" "<<endl;
	graph2<<m_solidMers<<" edges (solid mers)"<<endl;
	//for(MAP_TYPE<VERTEX_TYPE,MAP_TYPE<VERTEX_TYPE,vector<int> > >::iterator i=m_graph.begin();i!=m_graph.end();i++){
	for(CustomMap<VertexData>::iterator i=m_data->begin();i!=m_data->end();i++){
		VERTEX_TYPE prefix=i.first();
		VertexData dataStructure= (i.second());
		vector<VERTEX_TYPE>children=dataStructure.getChildren(prefix);
		//for(MAP_TYPE<VERTEX_TYPE,vector<int> >::iterator j=i->second.begin();j!=i->second.end();j++){
		for(int j=0;j<(int)children.size();j++){
			VERTEX_TYPE suffix=children[j];
			vector<AnnotationElement>*reads=dataStructure.getAnnotations(suffix);
			graph<<prefix<<" "<<suffix<<" "<<reads->size();
			graph2<<idToWord(prefix,m_wordSize)<<" -> "<<idToWord(suffix,m_wordSize)<<" "<<reads->size();
			for(int k=0;k<(int)reads->size();k++){
				graph<<" "<<reads->at(k).readId<<" "<<reads->at(k).readPosition<<" "<<reads->at(k).readStrand;
				graph2<<" "<<reads->at(k).readId<<" "<<reads->at(k).readPosition<<" "<<reads->at(k).readStrand;
			}
			graph<<endl;
			graph2<<endl;
		}
	}
	graph.close();
	graph2.close();
}

// fill m_graph, and m_vertex_parents
void DeBruijnAssembler::load_graphFrom_file(){
	(*m_cout)<<"[load_graphFrom_file]"<<endl;
	
	ifstream graph(m_graphFile.c_str());
	int n;
	graph>>n;
	(*m_cout)<<n<<endl;
	for(int i=0;i<n;i++){
		if(i%100000==0){
			(*m_cout)<<i<<" / "<<n<<endl;
		}
		VERTEX_TYPE prefix;
		VERTEX_TYPE suffix;
		int reads;
		graph>>prefix>>suffix>>reads;
		if(!m_data->find(prefix)){
			VertexData data;
			m_data->add(prefix,data);
		}
		if(!m_data->find(suffix)){
			VertexData data;
			m_data->add(suffix,data);
		}
		//m_vertex_parents[suffix].insert(prefix);
		m_data->get(suffix).addParent(prefix);
		//vector<int>*theReads=&(m_graph[prefix][suffix]);
		VertexData*vertexData=&(m_data->get(prefix));
		for(int j=0;j<reads;j++){
			int read;
			int position;
			char strand;
			graph>>read>>position>>strand;
			//theReads->push_back(read);
			vertexData->addAnnotation(suffix,read,position,strand);
		}
	}
	(*m_cout)<<n<<" / "<<n<<endl;
	graph.close();
}

void DeBruijnAssembler::buildGraph(SequenceDataFull*sequenceData){
	m_sequenceData=sequenceData;
	ostream&m_cout=*(this->m_cout);
	bool debug=m_DEBUG;
	//debug=false;
	bool useCache=false;


	if(debug){
		useCache=true;
	}

	ifstream f(m_graphFile.c_str());
	if(!f)
		useCache=false;
	f.close();

	if(useCache){
		m_cout<<"Using cache information from "<<m_graphFile<<endl;
		load_graphFrom_file();
	}else{
		build_From_Scratch(sequenceData);
		if(debug==true)
			writeGraph();
	}

}



string DeBruijnAssembler::pathToDNA(vector<VERTEX_TYPE>*path){
	ostringstream contigSequence;
	for(vector<VERTEX_TYPE>::iterator i=path->begin();i!=path->end();i++){
		if(i==path->begin()){
			contigSequence<< idToWord(*i,m_wordSize);
		}else{
			contigSequence<<getLastSymbol(*i);
		}
	}
	return contigSequence.str();
}


// convert k-mer to VERTEX_TYPE
VERTEX_TYPE DeBruijnAssembler::wordId(const char*a){
	VERTEX_TYPE i=0;
	for(int j=0;j<(int)strlen(a);j++){
		VERTEX_TYPE k=0;
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
	return i;
}

string DeBruijnAssembler::idToWord(VERTEX_TYPE i,int wordSize){
	string a="";
	int maxSize=sizeof(VERTEX_TYPE)*8/2; // 32
	for(int p=0;p<wordSize;p++){
		VERTEX_TYPE j=i;
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
		}
		

	}
	return a;
}

char DeBruijnAssembler::getLastSymbol(VERTEX_TYPE i){
        VERTEX_TYPE j=i;
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
 
void DeBruijnAssembler::setAssemblyDirectory(string assemblyDirectory){
	m_assemblyDirectory=assemblyDirectory;
	ostringstream name;
	name<<m_assemblyDirectory<<"/graph";
	m_graphFile=name.str();
}





void DeBruijnAssembler::Walk_In_GRAPH(){

	(*m_cout)<<endl;
	(*m_cout)<<"Collecting sources"<<endl;
	vector<VERTEX_TYPE> withoutParents;
	//for(MAP_TYPE<VERTEX_TYPE,MAP_TYPE<VERTEX_TYPE,vector<int> > >::iterator i=m_graph.begin();i!=m_graph.end();i++){
	for(CustomMap<VertexData>::iterator i=m_data->begin();i!=m_data->end();i++){
		VERTEX_TYPE prefix=i.first();
		if(i.second().getParents(prefix).size()==0){
			withoutParents.push_back(prefix);
		}
	}
	(*m_cout)<<"Done..., "<<withoutParents.size()<<" sources."<<endl;
	
	m_coverage_mean=1000;
	m_coverage_stddev=100;

	(*m_cout)<<endl;
	(*m_cout)<<"(k+1)-mers coverage statistics"<<endl;
	int NUMBER_OF_ROUND_FOR_SHRINKING=5;
	for(int round=1;round<=NUMBER_OF_ROUND_FOR_SHRINKING;round++){
		ostringstream distributionFile;
		distributionFile<<m_assemblyDirectory<<"/Coverage-Round-"<<round<<".txt";
		ofstream coverageBuffer(distributionFile.str().c_str());
		int threshold=1*m_coverage_mean+3*m_coverage_stddev;
		uint64_t sumCoverage=0;
		int edges=0;
		for(CustomMap<VertexData>::iterator i=m_data->begin();i!=m_data->end();i++){
			VERTEX_TYPE prefix=i.first();
			vector<VERTEX_TYPE>children=(i).second().getChildren(prefix);
			for(int j=0;j<children.size();j++){
				int value=(i).second().getAnnotations(children[j])->size();
				if(value>=threshold)
					continue;
				sumCoverage+=value;
				coverageBuffer<<value<<endl;
				edges++;
			}
		}
		coverageBuffer.close();
		int meanCoverage=sumCoverage/edges;
		uint64_t sum_ofSquaredDiffs=0;
		for(CustomMap<VertexData>::iterator i=m_data->begin();i!=m_data->end();i++){
			VERTEX_TYPE prefix=i.first();
			vector<VERTEX_TYPE>children=(i).second().getChildren(prefix);
			for(int j=0;j<children.size();j++){
				int value=(i).second().getAnnotations(children[j])->size();
				if(value>=threshold)
					continue;
				int diff=meanCoverage-value;
				sum_ofSquaredDiffs+=diff*diff;
			}
		}

		int  stddev=sqrt(sum_ofSquaredDiffs/(edges-1));

	
		m_coverage_mean=meanCoverage;
		m_coverage_stddev=stddev;
		(*m_cout)<<"Round "<<round<<": Mean="<<m_coverage_mean<<", Standard deviation="<<m_coverage_stddev<<endl;
	}


	ostream&m_cout=*(this->m_cout);	
	vector<VERTEX_TYPE> sources=withoutParents;
	set<VERTEX_TYPE> sourcesVisited;
	vector<int> The_Discovery_Of_Sources;
	vector<int>VisitsOfSources;
	string assemblyAmos=m_assemblyDirectory+"/Assembly.afg";
	string contigsFile=m_assemblyDirectory+"/fasta.contigs";
	
	ofstream amosFile(assemblyAmos.c_str());
	ofstream contigsFileStream(contigsFile.c_str());
	int contigId=1;
	int round=1;
	while(sources.size()>0){
		m_cout<<endl;
		m_cout<<"Round: "<<round<<", "<<sources.size()<<" sources."<<endl;
		round++;
		The_Discovery_Of_Sources.push_back(sources.size());
		vector<VERTEX_TYPE> newSources;
		for(int i=0;i<(int)sources.size();i++){
			map<VERTEX_TYPE,int> visits;
			VERTEX_TYPE prefix=sources[i];
			if(sourcesVisited.count(prefix)>0)
				continue;
			m_cout<<endl;
			sourcesVisited.insert(prefix);
			{
				m_cout<<"Source "<<i+1<<" / "<<sources.size()<<endl;//" REPEAT MODE"<<endl;
				m_cout<<"From: "<<idToWord(prefix,m_wordSize)<<endl;
			//m_cout<<m_contig_paths.size()<<" contigs"<<endl;
				vector<VERTEX_TYPE>path;
				path.push_back(prefix);
				vector<map<int,map<char,int> > > currentReadPositions;
				contig_From_SINGLE(&currentReadPositions,&path,&newSources);
				m_cout<<path.size()<<" vertices"<<endl;
				if(path.size()<=2)
					continue;
				writeContig_Amos(&currentReadPositions,&path,&amosFile,contigId);
				writeContig_fasta(&path,&contigsFileStream,contigId);
				//m_contig_paths.push_back(path);
				contigId++;
			}
		}
		sources=newSources;
		VisitsOfSources.push_back(sourcesVisited.size());
	}
	amosFile.close();
	contigsFileStream.close();
	m_cout<<endl;

	(m_cout)<<"Source discovery"<<endl;
	(m_cout)<<"Iteration Sources Cumulative"<<endl;
	for(int i=0;i<(int)The_Discovery_Of_Sources.size();i++)
		(m_cout)<<i+1<<" "<<The_Discovery_Of_Sources[i]<<" "<<VisitsOfSources[i]<<endl;

	(m_cout)<<endl;
	return;
}

void DeBruijnAssembler::Algorithm_Assembler_20090121(){
	Walk_In_GRAPH();
}



// get a walk from a vertex, with a path, up to maxSize,
// the path is not added in the walk
vector<VERTEX_TYPE> DeBruijnAssembler::getWalk(VERTEX_TYPE prefix,vector<VERTEX_TYPE>*path,int length,vector<map<int,map<char,int> > >*currentReadPositions){
	vector<VERTEX_TYPE> subPath;
	vector<VERTEX_TYPE> path1=*path;
	path1.push_back(prefix);
	subPath.push_back(prefix);
	vector<VERTEX_TYPE> nextVertices1=m_data->get(path1.at(path1.size()-1)).getChildren(path1.at(path1.size()-1));
	while((int)nextVertices1.size()==1&&(int)subPath.size()<=length){
		prefix=nextVertices1[0];
		path1.push_back(prefix);
		nextVertices1=m_data->get(path1.at(path1.size()-1)).getChildren(path1.at(path1.size()-1));
		subPath.push_back(prefix);
	}
	return subPath;
}

bool DeBruijnAssembler::DETECT_BUBBLE(vector<VERTEX_TYPE>*path,VERTEX_TYPE a,VERTEX_TYPE b){
	int maxSize=200;
	vector<VERTEX_TYPE> n1=getWalk(a,path,maxSize,NULL);
	vector<VERTEX_TYPE> n2=getWalk(b,path,maxSize,NULL);
	set<VERTEX_TYPE> n1_Table;
	for(int i=0;i<(int)n1.size();i++){
		n1_Table.insert(n1[i]);
	}
	for(int i=0;i<(int)n2.size();i++){
		if(n1_Table.count(n2[i])>0){
			(*m_cout)<<"BUBBLE Length: "<<endl;
			return true;
		}
	}
	return false;
}



void DeBruijnAssembler::contig_From_SINGLE(vector<map<int,map<char,int> > >*currentReadPositions,vector<VERTEX_TYPE>*path,vector<VERTEX_TYPE>*newSources){
	VERTEX_TYPE prefix=path->at(path->size()-1);
	map<int,int> usedReads;
	bool debug_print=m_DEBUG;
	//debug_print=true;
	//(*m_cout)<<"Depth: "<<path->size()<<endl;
	vector<VERTEX_TYPE> prefixNextVertices=nextVertices(path,currentReadPositions,newSources,&usedReads);
	while(prefixNextVertices.size()==1){
		prefix=prefixNextVertices[0];
		path->push_back(prefix);

		map<int,map<char,int> > a;
		(*currentReadPositions).push_back(a);
		//(*m_cout)<<"Pushing "<<idToWord(prefix,m_wordSize)<<endl;
		int cumulativeCoverage=0;
		int added=0;

		vector<AnnotationElement>*annotations=m_data->get(path->at(path->size()-2)).getAnnotations(path->at(path->size()-1));
		if(debug_print){
			(*m_cout)<<idToWord(path->at(path->size()-2),m_wordSize)<<" -> "<<idToWord(path->at(path->size()-1),m_wordSize)<<" ";
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
		int HIGHCOVERAGETHRESHOLD=1*m_coverage_mean+3*m_coverage_stddev;
		if(annotations->size()>=HIGHCOVERAGETHRESHOLD&&m_DEBUG){
			(*m_cout)<<"Coverage: "<<annotations->size()<<", refusing to start threading reads!"<<endl;
		}
		for(int h=0;h<(int)annotations->size();h++){
			if(usedReads.count(annotations->at(h).readId)==0){
					// add a read when it starts at its beginning...
				if(annotations->size()<HIGHCOVERAGETHRESHOLD&&
			annotations->at(h).readPosition==0){ // add at most a given amount of "new reads" to avoid depletion
					(*currentReadPositions)[path->size()-2][annotations->at(h).readId][annotations->at(h).readStrand]=annotations->at(h).readPosition; // = 0
					//(*m_cout)<<path->size()<<" "<<idToWord(path->at(path->size()-2),m_wordSize)<<" -> "<<idToWord(path->at(path->size()-1),m_wordSize)<<endl;
					//(*m_cout)<<"Adding read "<<m_sequenceData->at(annotations->at(h).readId)->getId()<<" "<<annotations->at(h).readStrand<<" "<<annotations->at(h).readPosition<<endl;
					cumulativeCoverage++;
					added++;
					usedReads[(annotations->at(h).readId)]=path->size()-2;
				}
			}else if(is_d_Threading(&(annotations->at(h)),currentReadPositions,path,&usedReads,false)){
				added++;
				usedReads[annotations->at(h).readId]=path->size()-2;

				(*currentReadPositions)[path->size()-2][annotations->at(h).readId][annotations->at(h).readStrand]=annotations->at(h).readPosition;
			}
		}

		if(debug_print){
			(*m_cout)<<(*currentReadPositions)[path->size()-2].size();
			for(map<int,map<char,int> >::iterator i=(*currentReadPositions)[path->size()-2].begin();i!=(*currentReadPositions)[path->size()-2].end();i++){
				for(map<char,int>::iterator j=i->second.begin();j!=i->second.end();j++){
					(*m_cout)<<" "<<m_sequenceData->at(i->first)->getId()<<" "<<j->first<<" "<<j->second;
				}
			}
			(*m_cout)<<endl;
		}

		prefixNextVertices=nextVertices(path,currentReadPositions,newSources,&usedReads);

		if(added==0){
			(*m_cout)<<"Nothing threaded"<<endl;
			break;
		}

		//(*m_cout)<<"adding "<<(*currentReadPositions).at(currentReadPositions->size()-1).size()<<endl;
	}

	VertexData dataStructure= (m_data->get(prefix));
	vector<VERTEX_TYPE>children=dataStructure.getChildren(prefix);


	(*m_cout)<<"Prefix "<<idToWord(prefix,m_wordSize)<<" "<<children.size()<<endl;
	// add newSources
	for(int j=0;j<(int)children.size();j++){
		(*m_cout)<<"Adding "<<idToWord(children[j],m_wordSize)<<endl;
		newSources->push_back(children[j]);
	}
	if(!debug_print)
		return;
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
			int position=0;
			for(vector<map<int,map<char,int> > >::iterator k=currentReadPositions->begin();k!=currentReadPositions->end();k++){
				break;
				//position=k->first;
				if(position<path->size()-500)
					continue;
				for(map<int,map<char,int> >::iterator i=(*k).begin();i!=(*k).end();i++){
					for(map<char,int>::iterator j=i->second.begin();j!=i->second.end();j++){
						(*m_cout)<<position<<" "<<m_sequenceData->at(i->first)->getId()<<" "<<j->first<<" "<<j->second<<endl;
					}
				}
				position++;
			}
			for(int i=399;i>=2;i--){
				break;
				(*m_cout)<<idToWord(path->at(path->size()-i),m_wordSize)<<" -> "<<idToWord(path->at(path->size()-i+1),m_wordSize)<<" ";
				vector<AnnotationElement>*annotations=m_data->get(path->at(path->size()-i)).getAnnotations(path->at(path->size()-i+1));
				(*m_cout)<<annotations->size();
				for(int j=0;j<annotations->size();j++){
					(*m_cout)<<" "<<annotations->at(j).readId<<" "<<annotations->at(j).readPosition<<" "<<annotations->at(j).readStrand;
				}
				(*m_cout)<<endl;
			}
			(*m_cout)<<children.size()<<" Choices"<<endl;
			for(int i=0;i<children.size();i++){
				//int score=recThreading(path->at(path->size()-1),children[i],&allowedReads);
				//(*m_cout)<<"Score: "<<score<<endl;
				(*m_cout)<<idToWord(path->at(path->size()-1),m_wordSize)<<" -> "<<idToWord(children[i],m_wordSize)<<" ";
				vector<AnnotationElement>*annotations=m_data->get(path->at(path->size()-1)).getAnnotations(children[i]);
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
	vector<VERTEX_TYPE> children=m_data->get(path->at(path->size()-1)).getChildren(path->at(path->size()-1));
	// start when nothing is done yet
	//(*m_cout)<<currentReadPositions->size()<<" "<<path->size()<<endl;
	if(currentReadPositions->size()==0){//||currentReadPositions->size()<path->size())
		return children;
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

	map<int,int > scoresSum;
	map<int,int> scoresMax;
	ostringstream debugBuffer;
	//(*m_cout)<<"Children"<<endl;
	for(int i=0;i<(int)children.size();i++){

		//(*m_cout)<<"Child "<<idToWord(children[i],m_wordSize)<<endl;
		vector<AnnotationElement>*thisEdgeData=(m_data->get(path->at(path->size()-1)).getAnnotations(children[i]));
		for(int j=0;j<(int)thisEdgeData->size();j++){
			uint32_t readId=thisEdgeData->at(j).readId;

			//(*m_cout)<<readId<<" "<<thisEdgeData->at(j).readStrand<<" "<<thisEdgeData->at(j).readPosition<<endl;
			if(is_d_Threading(&(thisEdgeData->at(j)),currentReadPositions,path,usedReads,true)){
				if(thisEdgeData->at(j).readPosition>=scoresMax[i]){
					scoresMax[i]=thisEdgeData->at(j).readPosition;
				}

				scoresSum[i]+=thisEdgeData->at(j).readPosition;
			}
		}
	}

	int best=-1;
	double factor=1.1; // magic number
	for(map<int,int>::iterator i=scoresSum.begin();i!=scoresSum.end();i++){
		//(*m_cout)<<i->second<<endl;
		bool isBest=true;
		for(map<int,int>::iterator j=scoresSum.begin();j!=scoresSum.end();j++){
			if(i->first==j->first)
				continue;
			if(i->second> factor*  j->second){
			}else{
				isBest=false;
			}
		}
		if(isBest)
			best=i->first;
	}
	
/*
	if(best==-1){
		for(map<int,int>::iterator i=scoresMax.begin();i!=scoresMax.end();i++){
			//(*m_cout)<<i->second<<endl;
			bool isBest=true;
			for(map<int,int>::iterator j=scoresMax.begin();j!=scoresMax.end();j++){
				if(i->first==j->first)
					continue;
				if(i->second> factor*  j->second){
				}else{
					isBest=false;
				}
			}
			if(isBest)
				best=i->first;
		}

	}

*/

	if(best==-1){
		for(int i=0;i<children.size();i++){
			if(newSources!=NULL){
				newSources->push_back(children[i]);
				(*m_cout)<<"Adding alternative source: "<<idToWord(children[i],m_wordSize)<<endl;
			}
		}
	}else{
		for(int i=0;i<children.size();i++){
			if(children[i]!=children[best]&&newSources!=NULL&&!DETECT_BUBBLE(path,children[i],children[best])){
				newSources->push_back(children[i]);
				(*m_cout)<<"Adding alternative source: "<<idToWord(children[i],m_wordSize)<<endl;
			}
		}

	}

/*
	if(newSources!=NULL){
		children=removeBubblesAndTips(children,path,currentReadPositions);
	}
*/
	if(best!=-1){
		vector<VERTEX_TYPE> output;
		output.push_back(children[best]);
		return output;
	}
	//(*m_cout)<<debugBuffer.str()<<endl;
	vector<VERTEX_TYPE> output;
	if(newSources!=NULL&&m_DEBUG){
		(*m_cout)<<"No children scored. "<<endl;
		vector<VERTEX_TYPE> allChildren=m_data->get(path->at(path->size()-1)).getChildren(path->at(path->size()-1));
		(*m_cout)<<"Before filtering "<<allChildren.size()<<endl;
		for(int i=0;i<allChildren.size();i++)
			(*m_cout)<<idToWord(allChildren[i],m_wordSize)<<endl;

		(*m_cout)<<"After filtering "<<children.size()<<endl;
		for(int i=0;i<children.size();i++)
			(*m_cout)<<idToWord(children[i],m_wordSize)<<endl;
	}
	return output;
}




void DeBruijnAssembler::setPairedInfo(string a){
	m_pairedInfoFile=a;
}


void DeBruijnAssembler::setBuckets(uint64_t buckets){
	m_buckets=buckets;
	m_data=new CustomMap<VertexData>(m_buckets);
}

DeBruijnAssembler::~DeBruijnAssembler(){
	if(m_data==NULL)
		return;
	delete m_data;
	m_data=NULL;
}



VERTEX_TYPE DeBruijnAssembler::reverseComplement_VERTEX(VERTEX_TYPE a){
	return wordId(reverseComplement(idToWord(a,m_wordSize)).c_str());
}




void DeBruijnAssembler::setMinimumCoverage(string coverage){
	m_minimumCoverageParameter=coverage;
}

void DeBruijnAssembler::indexReadStrand(int readId,char strand,SequenceDataFull*sequenceData,CustomMap<int>*solidMers){
	Read*read=sequenceData->at(readId);
	string sequence=read->getSeq();

	if(strand=='R')
		sequence=reverseComplement(sequence);

	for(int readPosition=0;readPosition<(int)sequence.length();readPosition++){
		string wholeWord=sequence.substr(readPosition,m_wordSize+1);
		if(readPosition>10000){
			//(*m_cout)<<"WTF "<<sequence.length()<<" "<<strlen(read->getSeq())<<strand<<endl;
		}
		if((int)wholeWord.length()==m_wordSize+1&&read->isValidDNA(&wholeWord)
		&&solidMers->find(wordId(wholeWord.c_str()))){
			VERTEX_TYPE prefix=wordId(wholeWord.substr(0,m_wordSize).c_str());
			VERTEX_TYPE suffix=wordId(wholeWord.substr(1,m_wordSize).c_str());
			if(!m_data->find(prefix)){
				VertexData vertexData;
				m_data->add(prefix,vertexData);
			}
			if(!m_data->find(suffix)){
				VertexData vertexData;
				m_data->add(suffix,vertexData);
			}
			m_data->get(prefix).addAnnotation(suffix,readId,readPosition,strand);
			m_data->get(suffix).addParent(prefix);
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
	(*file)<<">Contig"<<i+1<<"   "<<sequenceDNA.length()<<endl;
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
	if(usedReads->count(readId)==0)
		return false;
	int lastPosition=(*usedReads)[readId];
	char readStrand=annotation->readStrand;
	if((*currentReadPositions)[lastPosition][readId].count(readStrand)==0)
		return false;

	int distanceInRead=annotation->readPosition-(*currentReadPositions)[lastPosition][readId][readStrand];
	int distanceInPath=path->size()-2-lastPosition;
	if(beforeAdding)
		distanceInPath++;
	//(*m_cout)<<"R "<<distanceInRead<<" C "<<distanceInPath<<endl;
	return distanceInRead==distanceInPath; // allow error in read threading, if d !=1
}
