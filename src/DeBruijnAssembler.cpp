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
	this->m_cout=m_cout;
}

void DeBruijnAssembler::setWordSize(int k){
	m_default_window=6;
	m_minimumCoverage_for_walk=1;
	m_wordSize=k;
	m_threshold=0.015;
	m_pairedAvailable=false;
	m_longReadAvailable=false;
	DeBruijnAssembler::m_WordSize=k;
	m_longReadMode_threshold=100;
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
		if(!m_longReadAvailable&&(int)strlen(sequenceData->at(i)->getSeq())>m_longReadMode_threshold){
			m_longReadAvailable=true;
		}
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

	if(m_longReadAvailable)
		m_cout<<"[longRead available]"<<endl;

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
			m_cout<<"Best Coverage <- "<<m_minimumCoverage<<endl;
			break;
		}
	}

	if(m_minimumCoverageParameter=="auto"){
		m_cout<<"Using depletion curve (-minimumCoverage auto)"<<endl;
		m_cout<<"m_minimumCoverage <- "<<m_minimumCoverage<<endl;
	}else{
		m_cout<<"Not using the depletion curve"<<endl;
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
	m_cout<<" (this should be roughly twice the genome size)"<<endl;
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
	bool debug=true;
	bool useCache=false;


	if(debug){
		useCache=true;
		m_longReadAvailable=true;
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
		writeGraph();
	}

}

void DeBruijnAssembler::setMinimumContigSize(int minimumContigSize){
	m_minimumContigSize=minimumContigSize;
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

void DeBruijnAssembler::outputContigs(){
	ostream&m_cout=*(this->m_cout);
	m_cout<<endl;
	m_cout<<"Writing contigs"<<endl;
	string contigsFile=m_assemblyDirectory+"/Contigs.fa";
	string paths=m_assemblyDirectory+"/Walks.txt";
	string coverageFile=m_assemblyDirectory+"/Coverage.txt";
	ofstream outputWalk(paths.c_str());
	ofstream output(contigsFile.c_str());
	ofstream coverageStream(coverageFile.c_str());
	int columns=60;
	if(m_contig_paths.size()==0){ 
		m_cout<<"Error: 0 contigs"<<endl;
		return;
	}
	m_cout<<endl;
	m_cout<<m_contig_paths.size()<<" contigs"<<endl;	
	int totalLength=0;
	multiset<int> contigSizes;
	int stat_minimumContigSize=9999999;
	int stat_maximumContigSize=0;
	for(int i=0;i<(int)m_contig_paths.size();i++){
		string sequenceDNA=pathToDNA(&m_contig_paths[i]);
		totalLength+=sequenceDNA.length();
		contigSizes.insert(sequenceDNA.length());
		if((int)sequenceDNA.length()<stat_minimumContigSize)
			stat_minimumContigSize=sequenceDNA.length();
		if((int)sequenceDNA.length()>stat_maximumContigSize)
			stat_maximumContigSize=sequenceDNA.length();
		
		m_cout<<"Contig"<<i+1<<" "<<sequenceDNA.length()<<" nucleotides"<<endl;
		output<<">Contig"<<i+1<<"   "<<sequenceDNA.length()<<endl;
		outputWalk<<"Contig"<<i+1<<" "<<m_contig_paths[i].size()-1<<" edges"<<endl;
		coverageStream<<"Contig"<<i+1<<" "<<sequenceDNA.length()<<" nucleotides"<<endl;
		for(int j=0;j<(int)sequenceDNA.length();j++){
			int p=max(0,j-25);
			//(m_cout)<<j<<" "<<p<<endl;
			VERTEX_TYPE prefix=m_contig_paths[i][p];
			VERTEX_TYPE suffix=m_contig_paths[i][p+1];
			vector<AnnotationElement>*edgeReads=(m_data->get(prefix).getAnnotations(suffix));
			coverageStream<<j+1<<" "<<edgeReads->size()<<endl;
		}

		for(int j=0;j<=(int)m_contig_paths[i].size()-2;j++){
			VERTEX_TYPE prefix=m_contig_paths[i][j];
			VERTEX_TYPE suffix=m_contig_paths[i][j+1];
			vector<AnnotationElement>*edgeReads=(m_data->get(prefix).getAnnotations(suffix));
			outputWalk<<j+1<<endl;
			outputWalk<<idToWord(prefix,m_wordSize)<<" -> "<<idToWord(suffix,m_wordSize)<<endl;
			outputWalk<<edgeReads->size()<<" reads"<<endl;
			for(int k=0;k<(int)edgeReads->size();k++)
				outputWalk<<edgeReads->at(k).readId<<" ";
			outputWalk<<endl;
		}
		int j=0;
		while(j<(int)sequenceDNA.length()){
			output<<sequenceDNA.substr(j,columns);
			output<<endl;
			j+=columns;
		}

	}
	output.close();
	coverageStream.close();
	outputWalk.close();
	m_cout<<endl;
/*
	int averageContigLength=totalLength/(m_contig_paths.size());
	m_cout<<"Asssembly statistics"<<endl;
	m_cout<<"Estimated genome size (total contig bases): "<<totalLength<<endl;
	m_cout<<"Average contig size: "<<averageContigLength<<endl;
	m_cout<<"Smallest contig size: "<<stat_minimumContigSize<<endl;
	m_cout<<"Largest contig size: "<<stat_maximumContigSize<<endl;
	m_cout<<"Contigs: "<<m_contig_paths.size()<<endl;
*/
	int cumulativeSize=0;
	set<int> stat_nx_done;
	for(multiset<int>::reverse_iterator i=contigSizes.rbegin();i!=contigSizes.rend();i++){
		cumulativeSize+=*i;
		for(int j=50;j<=50;j+=10){
			if(cumulativeSize>=totalLength*j/100.0&&stat_nx_done.count(j)==0){
				stat_nx_done.insert(j);
				//m_cout<<"N"<<j<<" size: "<<*i<<endl;
				break;
			}
		}
	}
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
	int maxSize=sizeof(long int)*8/2; // 32
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

/**
 * check that the last l edges of path have at least C reads in common
 */
bool DeBruijnAssembler::passFilter_ShortRead(vector<VERTEX_TYPE>*path,int l,int C){
	//(*m_cout)<<"Using short reads "<<path->size()<<endl;
	MAP_TYPE<int,int> votes;
	l=min(l,path->size()-1); // size of verification
/*

         0     1      2      3     4     5    6    7     8      9     10
					<---------  l: 5 edges  ------->
*/
	for(int i=(int)path->size()-l-1;i<=(int)path->size()-2;i++){
		//vector<int> reads=m_graph[path[i]][path[i+1]];
		//cout<<"case 1"<<endl;
		vector<AnnotationElement>*reads=m_data->get((*path)[i]).getAnnotations((*path)[i+1]);
		for(vector<AnnotationElement>::iterator j=reads->begin();j!=reads->end();j++){
			votes[(*j).readId]++;
		}
	}
	int ok=0;
	for(MAP_TYPE<int,int>::iterator i=votes.begin();i!=votes.end();i++){
		if(i->second>=l)
			ok++;
	}
	//bool res=ok>=C;
	//(*m_cout)<<"Short "<<res<<endl;
	return ok>=C;
}

vector<vector<VERTEX_TYPE> >DeBruijnAssembler::Filter_Remove_Smaller_Duplicates(vector<vector<VERTEX_TYPE> > contigs){
	(*m_cout)<<"[Filter_Remove_Smaller_Duplicates] "<<contigs.size()<<endl;
	if(contigs.size()==1)
		return contigs;
	vector<vector<VERTEX_TYPE> >filteredContigs;
	MAP_TYPE<int,set<VERTEX_TYPE> > dictionnary;
	CustomMap<vector<int > > walksIndex(m_buckets);
	set<int> eliminatedContigs;
	// fill dictionnary
	(*m_cout)<<"Indexing.."<<endl;
	for(int id=0;id<(int)contigs.size();id++){
		vector<VERTEX_TYPE> contig=contigs[id];
		for(vector<VERTEX_TYPE>::iterator k=contig.begin();k!=contig.end();k++){
			dictionnary[id].insert(*k);
			dictionnary[id].insert(wordId(reverseComplement(idToWord(*k,m_wordSize)).c_str()));
			vector<int> theIndex;
			VERTEX_TYPE reverseNode=wordId(reverseComplement(idToWord(*k,m_wordSize)).c_str());
			if(!walksIndex.find(*k))
				walksIndex.add(*k,theIndex);
			if(!walksIndex.find(reverseNode))
				walksIndex.add(reverseNode,theIndex);

			walksIndex.get(*k).push_back(id);
			walksIndex.get(reverseNode).push_back(id);
		}
	}

	(*m_cout)<<"Done.."<<endl;
	for(int progress=0;progress<(int)contigs.size();progress++){
		vector<VERTEX_TYPE>contig=contigs[progress];
		if(progress%400==0)
			(*m_cout)<<progress<<" / "<<contigs.size()<<endl;
		bool isDuplicate=false;
		VERTEX_TYPE toCheck=contig[contig.size()/2]; // the middle one
		vector<int> contigsToCheck=walksIndex.get(toCheck);
		for(vector<int> ::iterator j=contigsToCheck.begin();j!=contigsToCheck.end();j++){
			int otherContig=*j;
			set<VERTEX_TYPE>*otherContigIndex=&dictionnary[otherContig];
			if(otherContig!=progress&&eliminatedContigs.count(otherContig)==0
			&&contigs[progress].size()<=contigs[otherContig].size()){
				int notFound=0;
				for(vector<VERTEX_TYPE>::iterator k=contig.begin();k!=contig.end();k++){
					if(otherContigIndex->count(*k)==0){
						notFound++;
					}
				}
				isDuplicate=(notFound<=m_longReadMode_threshold);
				if(isDuplicate){
					(*m_cout)<<"Vertices not found: "<<notFound<<" / "<<contig.size()<<progress<<" is the same as "<<otherContig<<endl;
					break;
				}
			}
		}
		if(isDuplicate){
			eliminatedContigs.insert(progress);
		}else{
			filteredContigs.push_back(contig);
		}
	}
	(*m_cout)<<contigs.size()<<" / "<<contigs.size()<<endl;
	(*m_cout)<<contigs.size()<<" -> "<<filteredContigs.size()<<endl;
	return filteredContigs;
}

void DeBruijnAssembler::Walk_In_GRAPH(){
	string walksRawFile=m_assemblyDirectory+"/RawWalks";
	ifstream f(walksRawFile.c_str());

	if(f){
		(*m_cout)<<"Reading cached data from "<<walksRawFile<<endl;
		int numberOfWalks;
		f>>numberOfWalks;
		for(int i=0;i<numberOfWalks;i++){
			int numberOfVertices;
			vector<VERTEX_TYPE> path;
			f>> numberOfVertices;
			for(int j=0;j<numberOfVertices;j++){
				VERTEX_TYPE a;
				f>>a;
				path.push_back(a);
			}
			m_contig_paths.push_back(path);
		}
		f.close();
		return;
	}else{
		f.close();
	}
	(*m_cout)<<"Collecting sources"<<endl;
	vector<VERTEX_TYPE> withoutParents;
	//for(MAP_TYPE<VERTEX_TYPE,MAP_TYPE<VERTEX_TYPE,vector<int> > >::iterator i=m_graph.begin();i!=m_graph.end();i++){
	for(CustomMap<VertexData>::iterator i=m_data->begin();i!=m_data->end();i++){
		VERTEX_TYPE prefix=i.first();
		if(i.second().getParents(prefix).size()==0){
			withoutParents.push_back(prefix);
		}
	}

	ostream&m_cout=*(this->m_cout);	
	vector<VERTEX_TYPE> sources=withoutParents;
	set<VERTEX_TYPE> sourcesVisited;
	vector<int> The_Discovery_Of_Sources;
	vector<int>VisitsOfSources;
	while(sources.size()>0){
		m_cout<<endl;
		m_cout<<"[sources]: "<<sources.size()<<endl;
		The_Discovery_Of_Sources.push_back(sources.size());
		vector<VERTEX_TYPE> newSources;
		for(int i=0;i<(int)sources.size();i++){
			map<VERTEX_TYPE,int> visits;
			VERTEX_TYPE prefix=sources[i];
			if(sourcesVisited.count(prefix)>0)
				continue;
			m_cout<<endl;
			sourcesVisited.insert(prefix);
			m_cout<<"Source "<<i+1<<" / "<<sources.size()<<endl;
			m_cout<<"From: "<<idToWord(prefix,m_wordSize)<<endl;
			//m_cout<<m_contig_paths.size()<<" contigs"<<endl;
			vector<VERTEX_TYPE>path;
			path.push_back(prefix);
			map<int,map<char,int> >currentReadPositions;
			contig_From_SINGLE(&currentReadPositions,&path,&newSources);
			m_cout<<"Vertices: "<<path.size()<<""<<endl;
			m_contig_paths.push_back(path);
		}
		sources=newSources;
		VisitsOfSources.push_back(sourcesVisited.size());
	}
	m_cout<<endl;

	(m_cout)<<"Source discovery"<<endl;
	(m_cout)<<"Iteration Sources Cumulative"<<endl;
	for(int i=0;i<(int)The_Discovery_Of_Sources.size();i++)
		(m_cout)<<i+1<<" "<<The_Discovery_Of_Sources[i]<<" "<<VisitsOfSources[i]<<endl;

	(m_cout)<<endl;
	ofstream streamBuffer(walksRawFile.c_str());
	streamBuffer<<m_contig_paths.size()<<endl;
	for(int i=0;i<(int)m_contig_paths.size();i++){
		streamBuffer<<m_contig_paths[i].size()<<endl;
		for(int j=0;j<(int)m_contig_paths[i].size();j++)
			streamBuffer<<m_contig_paths[i][j]<<" ";
		streamBuffer<<endl;
	}
	streamBuffer.close();
}

void DeBruijnAssembler::Algorithm_Assembler_20090121(){
	Walk_In_GRAPH();
	//vector<vector<VERTEX_TYPE> > largeContigs=Remove_Small_Contigs(m_contig_paths);
	
	//vector<vector<VERTEX_TYPE> > mergedContigs=Filter_Remove_Smaller_Duplicates_Cached(largeContigs);

	//vector<vector<VERTEX_TYPE> > finalContigs=ExtendReverseComplements(mergedContigs);
	//finalContigs=mergedContigs;
	//m_contig_paths=finalContigs;
}

vector<vector<VERTEX_TYPE> > DeBruijnAssembler::Filter_Remove_Smaller_Duplicates_Cached(vector<vector<VERTEX_TYPE > > largeContigs){
	string mergedContigsFile=m_assemblyDirectory+"/MergedContigs";
	ifstream f(mergedContigsFile.c_str());
	if(f){
		vector<vector<VERTEX_TYPE> > mergedContigs;
		(*m_cout)<<"Reading cached data from "<<mergedContigsFile<<endl;
		int numberOfWalks;
		f>>numberOfWalks;
		for(int i=0;i<numberOfWalks;i++){
			int numberOfVertices;
			vector<VERTEX_TYPE> path;
			f>> numberOfVertices;
			for(int j=0;j<numberOfVertices;j++){
				VERTEX_TYPE a;
				f>>a;
				path.push_back(a);
			}
			mergedContigs.push_back(path);
		}
		f.close();
		return mergedContigs;
	}else{
		f.close();
	}

	vector<vector<VERTEX_TYPE> > mergedContigs=Filter_Remove_Smaller_Duplicates(largeContigs);
	while(mergedContigs.size()<largeContigs.size()){
		largeContigs=mergedContigs;
		mergedContigs=Filter_Remove_Smaller_Duplicates(largeContigs);
	}
	ofstream streamBuffer(mergedContigsFile.c_str());
	streamBuffer<<mergedContigs.size()<<endl;
	for(int i=0;i<(int)mergedContigs.size();i++){
		streamBuffer<<mergedContigs[i].size()<<endl;
		for(int j=0;j<(int)mergedContigs[i].size();j++)
			streamBuffer<<mergedContigs[i][j]<<" ";
		streamBuffer<<endl;
	}
	streamBuffer.close();

	return mergedContigs;
}

// get a walk from a vertex, with a path, up to maxSize,
// the path is not added in the walk
vector<VERTEX_TYPE> DeBruijnAssembler::getWalk(VERTEX_TYPE prefix,vector<VERTEX_TYPE>*path,int length,map<int,map<char,int> >*currentReadPositions){
	vector<VERTEX_TYPE> subPath;
	vector<VERTEX_TYPE> path1=*path;
	path1.push_back(prefix);
	subPath.push_back(prefix);
	vector<VERTEX_TYPE> nextVertices1=nextVertices(&path1,currentReadPositions);
	while((int)nextVertices1.size()==1&&(int)subPath.size()<=length){
		prefix=nextVertices1[0];
		path1.push_back(prefix);
		nextVertices1=nextVertices(&path1,currentReadPositions);
		subPath.push_back(prefix);
	}
	return subPath;
}

// remove bubble, if any
// remote tips also..
vector<VERTEX_TYPE> DeBruijnAssembler::removeBubblesAndTips(vector<VERTEX_TYPE> vertices,vector<VERTEX_TYPE>*path,map<int,map<char,int> >*currentReadPositions){
	if(vertices.size()==1)
		return vertices;

	int maxSize=200;
/*
	if(vertices.size()==2){ // 454 HOMOPOLYMER  DETECTION
		string word1=idToWord(vertices[0],m_wordSize);
		string word2=idToWord(vertices[1],m_wordSize);
		if(word1[m_wordSize-1]==word1[m_wordSize-2]||word2[m_wordSize-1]==word2[m_wordSize-2]){
			int homopolymerLength=0;
			while(m_wordSize-2-homopolymerLength>=0
			&&word1[m_wordSize-2-homopolymerLength]==word1[m_wordSize-2]
			&&word1[m_wordSize-2-homopolymerLength]==word2[m_wordSize-2-homopolymerLength])
				homopolymerLength++;

			if(homopolymerLength>=2){
				int votes1=m_data->get(path[path.size()-1]).getReads(vertices[0]).size();
				int votes2=m_data->get(path[path.size()-1]).getReads(vertices[1]).size();
				vector<VERTEX_TYPE> withoutHomopolymer;
				if(votes1>votes2)
					withoutHomopolymer.push_back(vertices[0]);
				else
					withoutHomopolymer.push_back(vertices[1]);
				(*m_cout)<<word1<<" "<<votes1<<endl;
				(*m_cout)<<word2<<" "<<votes2<<endl;
				(*m_cout)<<"454 Homopolymer "<<homopolymerLength<<endl;
				return withoutHomopolymer;
			}
		}
	}*/

	if(vertices.size()==2){
		vector<VERTEX_TYPE> n1=getWalk(vertices[0],path,maxSize,currentReadPositions);
		vector<VERTEX_TYPE> n2=getWalk(vertices[1],path,maxSize,currentReadPositions);
		set<VERTEX_TYPE> n1_Table;
		for(int i=0;i<(int)n1.size();i++){
			n1_Table.insert(n1[i]);
		}
		bool bubble=false;
		for(int i=0;i<(int)n2.size();i++){
			if(n1_Table.count(n2[i])>0){
				(*m_cout)<<"BUBBLE Length: "<<i<<endl;
				bubble=true;
				break;
			}
		}
		if(bubble==true){
			vector<VERTEX_TYPE> withoutBubbles;
			withoutBubbles.push_back(vertices[0]);
			return withoutBubbles;
		}
	}
	if(vertices.size()>1){
		vector<VERTEX_TYPE> withoutTips;
		for(int i=0;i<(int)vertices.size();i++){
			vector<VERTEX_TYPE> subPath=getWalk(vertices[i],path,maxSize,currentReadPositions);
			if((int)subPath.size()<2*m_wordSize&&nextVertices(&subPath,currentReadPositions).size()==0){
				(*m_cout)<<"TIP Length: "<<subPath.size()<<endl;
				//" From: "<<idToWord(vertices[i],m_wordSize)<<endl;
			}else{
				withoutTips.push_back(vertices[i]);
			}
		}
		return withoutTips;
	}
	return vertices;
}

void DeBruijnAssembler::contig_From_SINGLE(map<int,map<char,int> >*currentReadPositions,vector<VERTEX_TYPE>*path,vector<VERTEX_TYPE>*newSources){
	VERTEX_TYPE prefix=path->at(path->size()-1);


	//(*m_cout)<<"Depth: "<<path->size()<<endl;
	vector<VERTEX_TYPE> prefixNextVertices=removeBubblesAndTips(nextVertices(path,currentReadPositions),path,currentReadPositions);
	while(prefixNextVertices.size()==1){
		vector<VERTEX_TYPE> children=removeBubblesAndTips(m_data->get(prefix).getChildren(prefix),path,currentReadPositions);
		for(int i=0;i<(int)children.size();i++){
			if(children[i]!=prefixNextVertices[0]){
				(*m_cout)<<"Adding "<<idToWord(children[i],m_wordSize)<<endl;
				newSources->push_back(children[i]);
			}
		}
		prefix=prefixNextVertices[0];
		path->push_back(prefix);
		//(*m_cout)<<idToWord(prefix,m_wordSize)<<endl;
		int cumulativeCoverage=0;
		vector<AnnotationElement>*annotations=m_data->get(path->at(path->size()-2)).getAnnotations(path->at(path->size()-1));
		for(int h=0;h<(int)annotations->size();h++){
			if((*currentReadPositions).count(annotations->at(h).readId)==0){// add a read when it starts at its beginning...
				if(cumulativeCoverage<m_minimumCoverage&&annotations->at(h).readPosition==0){ // add at most a given amount of "new reads" to avoid depletion
					(*currentReadPositions)[annotations->at(h).readId][annotations->at(h).readStrand]=annotations->at(h).readPosition;
					cumulativeCoverage++;
				}
			}else if(
			(*currentReadPositions)[annotations->at(h).readId][annotations->at(h).readStrand] +1 ==  annotations->at(h).readPosition){
				(*currentReadPositions)[annotations->at(h).readId][annotations->at(h).readStrand]=annotations->at(h).readPosition;
			}
		}
		if(path->size()%1000==0){
			//(*m_cout)<<"Vertices: "<<path->size()<<endl;
		}
		prefixNextVertices=nextVertices(path,currentReadPositions);

		if(prefixNextVertices.size()>1)
			prefixNextVertices=removeBubblesAndTips(prefixNextVertices,path,currentReadPositions);
	}

	VertexData dataStructure= (m_data->get(prefix));
	vector<VERTEX_TYPE>children=dataStructure.getChildren(prefix);


	//(*m_cout)<<"Prefix "<<idToWord(prefix,m_wordSize)<<" "<<children.size()<<endl;
	// add newSources
	for(int j=0;j<(int)children.size();j++){
		(*m_cout)<<"Adding "<<idToWord(children[j],m_wordSize)<<endl;
		newSources->push_back(children[j]);
	}



	if(children.size()==0){
		//(*m_cout)<<"DEAD END"<<endl;
		return ;
	}




/*
	// show path
	for(int i=path->size()-20;i<=(int)path->size()-2;i++){
		if(i<0){
			continue;
		}
		int aIndex=i;
		int bIndex=aIndex+1;
		VERTEX_TYPE a=path->at(aIndex);
		VERTEX_TYPE b=path->at(bIndex);
		(*m_cout)<<idToWord(a,m_wordSize)<<" -> "<<idToWord(b,m_wordSize)<<" ";
		//vector<int>*theReads=&(m_graph[a][b]);
		//cout<<"case 2"<<endl;
		vector<int> theReads= (m_data->get(a).getReads(b));
		(*m_cout)<<theReads.size();
		for(int j=0;j<(int)theReads.size();j++){
			(*m_cout)<<" "<<(theReads)[j];
		}
		(*m_cout)<<endl;
	
	}
	(*m_cout)<<children.size()<<" CHOICES"<<endl;
	//MAP_TYPE<VERTEX_TYPE,vector<int> >* dataStructure=&(m_graph[prefix]);
	//for(MAP_TYPE<VERTEX_TYPE,vector<int> >::iterator i=dataStructure->begin();i!=dataStructure->end();i++){
	for(int i=0;i<(int)children.size();i++){
		VERTEX_TYPE suffix=children[i];
		//(*m_cout)<<idToWord(prefix,m_wordSize)<<" -> "<<idToWord(i->first,m_wordSize)<<" ";
		(*m_cout)<<idToWord(prefix,m_wordSize)<<" -> "<<idToWord(suffix,m_wordSize)<<" ";
		vector<int> theReads= (dataStructure.getReads(suffix));
		(*m_cout)<<theReads.size();
		for(int j=0;j<(int)theReads.size();j++){
			(*m_cout)<<" "<<( theReads)[j];
		}
		(*m_cout)<<endl;
		
	}
*/
	return ;
}

/**
 * \param simple  force passFilter to use passFilter_ShortRead
 */
vector<VERTEX_TYPE> DeBruijnAssembler::nextVertices(vector<VERTEX_TYPE>*path,map<int,map<char,int> >*currentReadPositions){
	vector<VERTEX_TYPE> children=m_data->get(path->at(path->size()-1)).getChildren(path->at(path->size()-1));

	// start when nothing is done yet
	if(currentReadPositions->size()==0)
		return children;


	map<int,int > scores;
	ostringstream debugBuffer;
	//(*m_cout)<<"Children"<<endl;
	for(int i=0;i<(int)children.size();i++){
		debugBuffer<<"Child "<<idToWord(children[i],m_wordSize)<<endl;
		vector<AnnotationElement>*thisEdgeData=m_data->get(path->at(path->size()-1)).getAnnotations(children[i]);
		for(int j=0;j<(int)thisEdgeData->size();j++){
			uint32_t readId=thisEdgeData->at(j).readId;
			if(
				// the read is there in the path already
				currentReadPositions->count(readId)>0
				// the strand is available
		&&		(*currentReadPositions)[readId].count(thisEdgeData->at(j).readStrand)>0
				// the position is greater than the one in the database
		&& 		(*currentReadPositions)[readId][thisEdgeData->at(j).readStrand] +1==thisEdgeData->at(j).readPosition
			){
					if(thisEdgeData->at(j).readPosition>scores[i])
						scores[i]=thisEdgeData->at(j).readPosition;
					debugBuffer<<m_sequenceData->at(readId)->getId()<<" "<<thisEdgeData->at(j).readStrand<<endl;

					if(thisEdgeData->at(j).readStrand=='F')
						debugBuffer<<m_sequenceData->at(readId)->getSeq()<<endl;
					else
						debugBuffer<<reverseComplement(m_sequenceData->at(readId)->getSeq())<<endl;

					for(int h=0;h<thisEdgeData->at(j).readPosition;h++)
						debugBuffer<<" ";
					debugBuffer<<"*"<<endl;
			}
		}
	}

	if(scores.size()==0){
		vector<VERTEX_TYPE> output;
		return output;
	}

	//(*m_cout)<<"Children"<<endl;

	if(scores.size()==1){
		vector<VERTEX_TYPE> output;
		output.push_back(children[scores.begin()->first]);
		return output;
	}

	int best=-1;
	for(map<int,int>::iterator i=scores.begin();i!=scores.end();i++){
		//(*m_cout)<<i->second<<endl;
		bool isBest=true;
		for(map<int,int>::iterator j=scores.begin();j!=scores.end();j++){
			if(i->first==j->first)
				continue;
			if(i->second> 1.2* j->second){
			}else{
				isBest=false;
			}
		}
		if(isBest)
			best=i->first;
	}
	if(best!=-1){
		vector<VERTEX_TYPE> output;
		output.push_back(children[best]);
		return output;
	}
	(*m_cout)<<debugBuffer.str()<<endl;
	return children;
}

bool DeBruijnAssembler::passFilter(vector<VERTEX_TYPE>*path,int l,int C){
	//(*m_cout)<<"simple "<<simple<<endl;
	if(l!=m_default_window){
		//(*m_cout)<<"Forcing passFilter_ShortRead"<<endl;
		return passFilter_ShortRead(path,l,C);
	}
	if((int)path->size()>m_longReadMode_threshold&&m_pairedAvailable){
		if(passFilter_Paired(path,l,C))
			return true;
	}

	if((int)path->size()>m_longReadMode_threshold&&m_longReadAvailable){
		if(passFilter_LongRead(path,l,C))
			return true;
	}
	return passFilter_ShortRead(path,l,C);
}

// TODO
bool DeBruijnAssembler::passFilter_Paired(vector<VERTEX_TYPE>*path,int l,int C){
	return false;
}

bool DeBruijnAssembler::passFilter_LongRead(vector<VERTEX_TYPE>*path,int l,int C){
	//(*m_cout)<<"Using long reads"<<endl;
/*
	for(int i=0;i<(int)path->size()-1;i++){
		//(*m_cout)<<i+1<<" "<<idToWord(path->at(i),m_wordSize)<<" -> "<<idToWord(path->at(i+1),m_wordSize)<<" ";
		vector<int> reads=m_data->get((*path)[i]).getReads((*path)[i+1]);
		//(*m_cout)<<reads.size();
		for(int j=0;j<(int)reads.size();j++){
			//(*m_cout)<<" "<<reads[j];
		}
		//(*m_cout)<<endl;
	}
*/
	int lastIndex=path->size()-1;
	map<int,int> votes;
	//cout<<"case 3"<<endl;
	vector<AnnotationElement>*reads=m_data->get((*path)[lastIndex-1]).getAnnotations((*path)[lastIndex]);
	for(int i=0;i<(int)reads->size();i++){
		//(*m_cout)<<" "<<reads[i];
		votes[reads->at(i).readId]++;
	}
	//(*m_cout)<<endl;

	int MAX_DISTANCE=500;
	//for(int firstIndex=m_wordSize+1;firstIndex<min(path->size(),400);firstIndex+=m_wordSize+1){
	for(int firstIndex=(int)path->size()-MAX_DISTANCE;firstIndex<=(int)path->size()-m_longReadMode_threshold;firstIndex+=m_wordSize+1){
		if(firstIndex<0)
			continue;
		map<int,int> votes2=votes;
		int prefixIndex=firstIndex;
		//(*m_cout)<<prefixIndex+1<<endl;
		int suffixIndex=prefixIndex+1;
		//cout<<"case 4"<<endl;
		reads=m_data->get((*path)[prefixIndex]).getAnnotations((*path)[suffixIndex]);
		for(int i=0;i<(int)reads->size();i++)
			votes2[reads->at(i).readId]++;
		int coverage=0;
		for(map<int,int>::iterator i=votes2.begin();i!=votes2.end();i++){
			if(i->second==2){
				coverage++;
			}
			//(*m_cout)<<i->first<<" "<<i->second<<endl;
		}

		//(*m_cout)<<C<<" "<<prefixIndex+1<<" "<<coverage<<endl;
		if(coverage>=C){
			//(*m_cout)<<"match!"<<endl;
			return true;
		}
	}
	//(*m_cout)<<"No match found"<<endl;
	return false;
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



// join
vector<vector<VERTEX_TYPE> > DeBruijnAssembler::ExtendReverseComplements(vector<vector<VERTEX_TYPE> > contigs){
	(*m_cout)<<endl;
	(*m_cout)<<"[DeBruijnAssembler::ExtendReverseComplements>>"<<endl;
	vector<vector<VERTEX_TYPE> >enhancedContigs=contigs;
	bool joiningOccured=true;
	while(joiningOccured==true){
		joiningOccured=false;

		vector<vector<VERTEX_TYPE> >newContigs;
		set<int>contigsProcessed;

		// -------------------------------------->
		// ----->                          ---->
		// forward_start                       forward_end
		//
		// <-----                          <-----
		// reverse_start                  reverse_end
		//

		map<VERTEX_TYPE,vector<int> > forward_start;
		map<VERTEX_TYPE,vector<int> > forward_end;
		map<VERTEX_TYPE,vector<int> > reverse_start;
		map<VERTEX_TYPE,vector<int> > reverse_end;
		
		// index
		for(int i=0;i<(int)enhancedContigs.size();i++){
			vector<VERTEX_TYPE>contig=enhancedContigs[i];
			forward_start[contig[0]].push_back(i);
			forward_end[contig[contig.size()-1]].push_back(i);
			reverse_start[reverseComplement_VERTEX(contig[0])].push_back(i);
			reverse_end[reverseComplement_VERTEX(contig[contig.size()-1])].push_back(i);
		}

		// Join ends
		for(int i=0;i<(int)enhancedContigs.size();i++){

			if(contigsProcessed.count(i)>0)
				continue;
			if(i%100==0)
				(*m_cout)<<i<<" / "<<enhancedContigs.size()<<endl;
			vector<VERTEX_TYPE>contig=enhancedContigs[i];
			

			vector<int> start_Matches_Foward;
			vector<int> start_Matches_Reverse;
			vector<int> end_Matches_Reverse;
			
			for(int l=0;l<(int)contig.size();l++){
				VERTEX_TYPE vertex=contig[l];
				if(forward_start.count(vertex)>0){
					vector<int> toAdd=not_Processed(forward_start[vertex],contigsProcessed,i);
					for(int o=0;o<(int)toAdd.size();o++)
						start_Matches_Foward.push_back(toAdd[o]);
				}
				if(reverse_start.count(vertex)>0){
					vector<int> toAdd=not_Processed(reverse_start[vertex],contigsProcessed,i);
					for(int o=0;o<(int)toAdd.size();o++)
						start_Matches_Reverse.push_back(toAdd[o]);
				}
				if(reverse_end.count(vertex)>0){
					vector<int> toAdd=not_Processed(reverse_end[vertex],contigsProcessed,i);
					for(int o=0;o<(int)toAdd.size();o++)
						end_Matches_Reverse.push_back(toAdd[o]);
				}
			}
/*
* 		start_Matches_Foward
 *
 *			---------------------------------------------->       Contig A
 *			<----------------------------------------------
 *                                                            -------------------------------------->  Contig B
 *                                                            <--------------------------------------
 */

			int WINDOW_SIZE=500;
			if(start_Matches_Foward.size()>0){
				//(*m_cout)<<"start_Matches_Foward"<<endl;
				if(start_Matches_Foward.size()==1){
					//(*m_cout)<<"Exactly 1"<<endl;
					int otherContigId=start_Matches_Foward[0];
					vector<VERTEX_TYPE>otherContig=contigs[otherContigId];
					VERTEX_TYPE vertex=otherContig[0];
					vector<VERTEX_TYPE> newContig;
					for(int u=0;u<(int)contig.size();u++){
						VERTEX_TYPE otherVertex=contig[u];
						if(otherVertex==vertex)
							break;
						newContig.push_back(otherVertex);
					}
					bool correct=true;
					int windowForJoin=0;
					for(int u=0;u<(int)otherContig.size();u++){
						if(newContig.size()==0){
							newContig.push_back(otherContig[u]);
							continue;
						}
						(*m_cout)<<"newContig.size "<<newContig.size()<<endl;
						if((newContig.size()-1)>=0&&!(m_data->find(newContig[newContig.size()-1]))){
							correct=false;
							break;
						}
						if((newContig.size()-1)>=0&&!(m_data->get(
							newContig[newContig.size()-1])
							.hasChild(otherContig[u]))){
							correct=false;
							break;
						}
						newContig.push_back(otherContig[u]);
						if(windowForJoin<=WINDOW_SIZE
						&&!passFilter(&newContig,m_default_window,m_minimumCoverage_for_walk)){
							correct=false;
							break;
						}
						windowForJoin++;
					}


					if(correct&&addNewContig(&newContigs,&newContig,i,otherContigId,&contigsProcessed))
						joiningOccured=true;
				}
/* *  		start_Matches_Reverse
*              	---------------------------------------------->    Contig A
*		<----------------------------------------------
*   ------------------->
*   <-------------------		Contig B
*/
			}else if(start_Matches_Reverse.size()>0){ 
				//(*m_cout)<<"start_Matches_Reverse"<<endl;
				if(start_Matches_Reverse.size()==1){
					//(*m_cout)<<"Exactly 1"<<endl;
					int otherContigId=start_Matches_Reverse[0];
					vector<VERTEX_TYPE>otherContig=contigs[otherContigId];
					VERTEX_TYPE vertex=reverseComplement_VERTEX(otherContig[0]);
					vector<VERTEX_TYPE> newContig;
	
					for(int u=otherContig.size()-1;u>=0;u--){
						newContig.push_back(reverseComplement_VERTEX(otherContig[u]));
					}
					bool ok=false;
					bool correct=true;
					int windowForJoin=0;
					for(int u=0;u<(int)contig.size();u++){
						if(ok){
							if((newContig.size()-1)>=0&&!(m_data->get(newContig[newContig.size()-1]).hasChild(contig[u]))){
								correct=false;
								break;
							}
							newContig.push_back(contig[u]);
							windowForJoin++;
						}
						if(contig[u]==vertex)
							ok=true;
						if(ok&&windowForJoin<=WINDOW_SIZE&&!passFilter(&newContig,m_default_window,m_minimumCoverage_for_walk)){
							correct=false;
							break;
						}
					}


					if(correct&&addNewContig(&newContigs,&newContig,i,otherContigId,&contigsProcessed))
						joiningOccured=true;
				}
/*
*  		end_Matches_Reverse
*  			
*  		--------------------------------------------------> Contig A
*  		<--------------------------------------------------
*
*  							------------------------------------>
*  							<------------------------------------ Contig B
*/

			}else if(end_Matches_Reverse.size()>0){ 
				//(*m_cout)<<"end_Matches_Reverse"<<endl;
				if(end_Matches_Reverse.size()==1){
					//(*m_cout)<<"Exactly 1"<<endl;
					int otherContigId=end_Matches_Reverse[0];
					vector<VERTEX_TYPE>otherContig=contigs[otherContigId];
					VERTEX_TYPE vertex=reverseComplement_VERTEX(otherContig[otherContig.size()-1]);
					vector<VERTEX_TYPE> newContig;

					for(int u=0;u<(int)contig.size();u++){
						VERTEX_TYPE otherVertex=contig[u];
						if(otherVertex==vertex)
							break;
						newContig.push_back(otherVertex);
						//(*m_cout)<<u<<" "<<idToWord(otherVertex,m_wordSize)<<endl;
					}
					bool correct=true;
					int windowForJoin=0;
					for(int u=otherContig.size()-1;u>=0;u--){
						if(newContig.size()==0){
							newContig.push_back(reverseComplement_VERTEX(otherContig[u]));
							continue;
						}
						if((newContig.size()-1)>=0&&
			!(m_data->get(newContig[newContig.size()-1]).hasChild(reverseComplement_VERTEX(otherContig[u])))){
							correct=false;
							break;
						}

						newContig.push_back(reverseComplement_VERTEX(otherContig[u]));
						//(*m_cout)<<u<<" "<<idToWord(reverseComplement_VERTEX(otherContig[u]),m_wordSize)<<endl;
						if(windowForJoin<=WINDOW_SIZE
					&&!passFilter(&newContig,m_default_window,m_minimumCoverage_for_walk)){
							correct=false;
							break;
						}
						windowForJoin++;
					}

					if(correct&&addNewContig(&newContigs,&newContig,i,otherContigId,&contigsProcessed))
						joiningOccured=true;
				}
			}

			if(contigsProcessed.count(i)==0){
				newContigs.push_back(contig);
				contigsProcessed.insert(i);
			}
		}
		(*m_cout)<<enhancedContigs.size()<<" / "<<enhancedContigs.size()<<endl;
		(*m_cout)<<enhancedContigs.size()<<" -> "<<newContigs.size()<<endl;
		enhancedContigs=newContigs;
	}
	return enhancedContigs;
}

VERTEX_TYPE DeBruijnAssembler::reverseComplement_VERTEX(VERTEX_TYPE a){
	return wordId(reverseComplement(idToWord(a,m_wordSize)).c_str());
}

vector<int> DeBruijnAssembler::not_Processed(vector<int>contigs,set<int>processed,int self){
	vector<int> out;
	for(int i=0;i<(int)contigs.size();i++)
		if(processed.count(contigs[i])==0&& contigs[i]!=self)
			out.push_back(contigs[i]);
	return out;
}

vector<vector<VERTEX_TYPE> >DeBruijnAssembler::Remove_Small_Contigs(vector<vector<VERTEX_TYPE> > contigs){
	vector<vector<VERTEX_TYPE> > largeContigs;
	for(int i=0;i<(int)contigs.size();i++){
		int nucleotides=m_contig_paths[i].size()+m_wordSize-1;
		if(nucleotides<m_minimumContigSize)
			continue;
		largeContigs.push_back(m_contig_paths[i]);
	}
	(*m_cout)<<"Remove_Small_Contigs] "<<contigs.size()<<" -> "<<largeContigs.size()<<endl;
	(*m_cout)<<endl;
	return largeContigs;
}

bool DeBruijnAssembler::addNewContig(vector<vector<VERTEX_TYPE> >*newContigs,vector<VERTEX_TYPE>*newContig,int currentContigId,int otherContigId,
	set<int>*contigsProcessed){
	bool isValid=true;
	for(int u=0;u<(int)newContig->size()-1;u++){
		if(!m_data->get(newContig->at(u)).hasChild(newContig->at(u+1))){
			isValid=false;
			//(*m_cout)<<"incorrect path"<<endl;
			for(int y=0;y<(int)newContig->size()-1;y++){
				//(*m_cout)<<idToWord(newContig->at(y),m_wordSize)<<" -> "<<idToWord(newContig->at(y+1),m_wordSize)<<endl;
			}
			break;
		}
	}
	if(isValid){
		newContigs->push_back(*newContig);
		contigsProcessed->insert(currentContigId);
		contigsProcessed->insert(otherContigId);
		//(*m_cout)<<"Joining"<<endl;
		return true;
	}
	return false;
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
