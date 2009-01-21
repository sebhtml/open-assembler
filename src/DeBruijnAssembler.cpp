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



#include<cstdlib>
#include"Logger.hpp"
#include<iostream>
#include"DeBruijnAssembler.h"
#include<fstream>
#include<vector>



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

// TODO: add paired end reads ?
// with 454 it is 3k

int DeBruijnAssembler::m_NUCLEOTIDE_A=0;
int DeBruijnAssembler::m_NUCLEOTIDE_T=1;
int DeBruijnAssembler::m_NUCLEOTIDE_C=2;
int DeBruijnAssembler::m_NUCLEOTIDE_G=3;

DeBruijnAssembler::DeBruijnAssembler(Logger*m_cout){
	this->m_cout=m_cout;
}




void DeBruijnAssembler::setWordSize(int k){
	m_wordSize=k;
}

void DeBruijnAssembler::setMinimumCoverage(int m){
	m_minimumCoverage=m;
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
	return 'E';
}

string DeBruijnAssembler::reverseComplement(string a){
	ostringstream i;
	for(int p=a.length()-1;p>=0;p--){
		char b=complement(a[p]);
		if(b=='E'){
			cout<<b<<endl;
			exit(0);
		}
		i<< complement(a[p]);
	}
	return i.str();
}



void DeBruijnAssembler::build_From_Scratch(){
	Logger&m_cout=*(this->m_cout);
	m_cout<<endl;
	m_cout<<"Collecting high-quality mers from reads"<<endl;
	map<VERTEX_TYPE,int> words;
	int last_vertices_size=-1;
	for(int i=0;i<(int)m_reads->size();i++){
		if(i%10000==0){
			m_cout<<"Reads: "<<i<<" / "<<m_reads->size()<<endl;
		}
		vector<VERTEX_TYPE>*highQualityMers=m_reads->at(i)->getHighQualityMers(m_wordSize);
		for(vector<VERTEX_TYPE>::iterator iteratorMer=highQualityMers->begin();iteratorMer!=highQualityMers->end();iteratorMer++){
			//m_cout<<endl;
			//m_cout<<wordFoward<<endl;
			//m_cout<<wordReverse<<endl;
			//m_cout<<idToWord(wordId(wordFoward))<<endl;
			words[*iteratorMer]++;
			if((int)words.size()%1000000==0&&(int)words.size()!=last_vertices_size){
				m_cout<<"High-quality mers: "<<words.size()<<endl;
				last_vertices_size=words.size();
			}
		}
	}	


	m_cout<<"High-quality mers: "<<words.size()<<endl;
	m_cout<<"Reads: "<<m_reads->size()<<" / "<<m_reads->size()<<endl;
	int processed=0;
	int solid=0;




	m_cout<<endl;
	//m_cout<<"Collecting solid mers"<<endl;
	for(map<VERTEX_TYPE,int>::iterator i=words.begin();i!=words.end();i++){
		processed++;
		if(i->second>=m_minimumCoverage){
			solid++;
			VERTEX_TYPE w=i->first;
			m_solidMers.insert(w);
		}
		if(processed%1000==0){
			//m_cout<<"Processed words: "<<processed<<" / "<<words.size()<<endl;
		}
		if(solid%1000==0){
			//m_cout<<"Solid words: "<<solid<<" / "<<words.size()<<endl;
		}
	}
	//m_cout<<"Processed mers: "<<processed<<" / "<<words.size()<<endl;
	m_cout<<"Solid mers: "<<solid<<" / "<<words.size()<<" ---> "<<((solid+0.0)/words.size()+0.0)*100.0<<"%"<<endl;
	m_cout<<" (this should be roughly twice the genome size)"<<endl;
	m_cout<<"Not-solid mers: "<<processed-solid<<" / "<<words.size()<<" ---> "<<(processed-solid+0.0)/words.size()*100.0<<"%"<<endl;

	m_cout<<endl;
	words.clear();
	//m_cout<<"Clearing mers"<<endl;
	//m_cout<<endl;



	// TODO: check if this can be done only once instead of twice
	// i.e. building the graph below is almost the same.
	m_cout<<"Indexing solid mers in reads, building graph with solid mers."<<endl; // <-------
	for(int readId=0;readId<(int)m_reads->size();readId++){
		if(readId%10000==0)
			m_cout<<readId<<" / "<<m_reads->size()<<endl;
		vector<VERTEX_TYPE>*highQualityMers=m_reads->at(readId)->getHighQualityMers(m_wordSize);
		for(vector<VERTEX_TYPE>::iterator merIterator=highQualityMers->begin();
			merIterator!=highQualityMers->end();merIterator++){
			if(m_solidMers.count(*merIterator)>0){
				VERTEX_TYPE solidMer=*merIterator;
				string wholeWord=idToWord(solidMer,m_wordSize+1);
				VERTEX_TYPE prefix=wordId(wholeWord.substr(0,m_wordSize).c_str());
				VERTEX_TYPE suffix=wordId(wholeWord.substr(1,m_wordSize).c_str());
				m_graph[prefix][suffix].push_back(readId);
				m_vertex_parents[suffix].insert(prefix);
			}
		}
	}
	m_cout<<m_reads->size()<<" / "<<m_reads->size()<<endl;
	m_cout<<endl;
	
	//m_cout<<"the graph is ready"<<endl;
	//m_cout<<"Clearing solid mers"<<endl;
	//solidMers.clear();
}

void DeBruijnAssembler::writeGraph(){
	ofstream graph(m_graphFile.c_str());
	string humanReadable=m_assemblyDirectory+"/Graph.txt";
	ofstream graph2(humanReadable.c_str());
	graph<<m_solidMers.size()<<" "<<endl;
	graph2<<m_solidMers.size()<<" edges (solid mers)"<<endl;
	for(map<VERTEX_TYPE,map<VERTEX_TYPE,vector<int> > >::iterator i=m_graph.begin();i!=m_graph.end();i++){
		VERTEX_TYPE prefix=i->first;
		for(map<VERTEX_TYPE,vector<int> >::iterator j=i->second.begin();j!=i->second.end();j++){
			VERTEX_TYPE suffix=j->first;
			vector<int> reads=j->second;
			graph<<prefix<<" "<<suffix<<" "<<reads.size();
			graph2<<idToWord(prefix,m_wordSize)<<" -> "<<idToWord(suffix,m_wordSize)<<" "<<reads.size();
			for(int k=0;k<(int)reads.size();k++){
				graph<<" "<<reads[k];
				graph2<<" "<<(*m_reads)[reads[k]]->getId();
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
		m_vertex_parents[suffix].insert(prefix);
		vector<int>*theReads=&(m_graph[prefix][suffix]);
		for(int j=0;j<reads;j++){
			int read;
			graph>>read;
			theReads->push_back(read);
		}
		//(*m_cout)<<m_graph[prefix][suffix].size()<<endl;
	}
	(*m_cout)<<n<<" / "<<n<<endl;
	graph.close();
}

void DeBruijnAssembler::buildGraph(vector<Read*>*reads){
	Logger&m_cout=*(this->m_cout);
	m_reads=reads;

	bool useCache=m_useCache=="yes";
	ifstream f(m_graphFile.c_str());
	if(!f)
		useCache=false;
	f.close();

	if(useCache){
		m_cout<<"Using cache information from "<<m_graphFile<<endl;
		load_graphFrom_file();
	}else{
		build_From_Scratch();
		writeGraph();
		m_solidMers.clear();
	}

}

void DeBruijnAssembler::setMinimumContigSize(int minimumContigSize){
	m_minimumContigSize=minimumContigSize;
}

string DeBruijnAssembler::pathToDNA(vector<VERTEX_TYPE> path){
	ostringstream contigSequence;
	for(vector<VERTEX_TYPE>::iterator i=path.begin();		i!=path.end();i++){
		if(i==path.begin()){
			contigSequence<< idToWord(*i,m_wordSize);
		}else{
			contigSequence<<getLastSymbol(*i);
		}
	}
	return contigSequence.str();
}

void DeBruijnAssembler::outputContigs(){
	Logger&m_cout=*(this->m_cout);
	// write contigs
	// TODO: don't duplicate each strand..
	m_cout<<endl;
	m_cout<<"Writing contigs"<<endl;
	string contigsFile=m_assemblyDirectory+"/Contigs.fa";
	ofstream output(contigsFile.c_str());
	int columns=60;
	vector<string> contigs;
	vector<int> contigNumbers;
	m_cout<<m_contig_paths.size()<<" contig paths"<<endl;
	for(int i=0;i<(int)m_contig_paths.size();i++){
		m_cout<<"Contig #"<<i+1<<": "<<m_contig_paths[i].size()+m_wordSize-1<<" nucleotides"<<endl;
		string contigSequence=pathToDNA(m_contig_paths[i]);
		contigs.push_back(contigSequence);
		contigNumbers.push_back(i+1);
	}


	if(contigs.size()==0){ 
		m_cout<<"Error: 0 contigs"<<endl;
		return;
	}
	m_cout<<endl;
	m_cout<<contigs.size()<<" contigs"<<endl;	
	//m_cout<<contigs.size()-duplicatedContig.size()<<" after mergin of complement contigs (the graph contains both strands)"<<endl;
	int totalLength=0;
	int contigId=1;
	multiset<int> contigSizes;
	int stat_minimumContigSize=9999999;
	int stat_maximumContigSize=0;
	for(int i=0;i<(int)contigs.size();i++){
		string sequenceDNA=contigs[i];
		totalLength+=sequenceDNA.length();
		contigSizes.insert(sequenceDNA.length());
		if((int)sequenceDNA.length()<stat_minimumContigSize)
			stat_minimumContigSize=sequenceDNA.length();
		if((int)sequenceDNA.length()>stat_maximumContigSize)
			stat_maximumContigSize=sequenceDNA.length();

		m_cout<<sequenceDNA.length()<<" nucleotides"<<endl;
		//m_cout<<"contig "<<iteratorContigFirstVertex->first<<endl;
		int j=0;
		output<<">Contig"<<contigId<<"   "<<sequenceDNA.length()<<endl;
		contigId++;
		while(j<(int)sequenceDNA.length()){
			output<<sequenceDNA.substr(j,columns);
			output<<endl;
			j+=columns;
		}

	}
	output.close();
	m_cout<<endl;
	int averageContigLength=totalLength/(contigs.size());
	m_cout<<"Asssembly statistics"<<endl;
	m_cout<<"Estimated genome size (total contig bases): "<<totalLength<<endl;
	m_cout<<"Average contig size: "<<averageContigLength<<endl;
	m_cout<<"Smallest contig size: "<<stat_minimumContigSize<<endl;
	m_cout<<"Largest contig size: "<<stat_maximumContigSize<<endl;
	m_cout<<"Contigs: "<<contigs.size()<<endl;
	int cumulativeSize=0;
	set<int> stat_nx_done;
	for(multiset<int>::reverse_iterator i=contigSizes.rbegin();i!=contigSizes.rend();i++){
		cumulativeSize+=*i;
		for(int j=50;j<=50;j+=10){
			if(cumulativeSize>=totalLength*j/100.0&&stat_nx_done.count(j)==0){
				stat_nx_done.insert(j);
				m_cout<<"N"<<j<<" size: "<<*i<<endl;
				break;
			}
		}
	}
}





// concert k-mer to VERTEX_TYPE
// algo works with (k+1)-mer
// TODO: replace *= by <<
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
	//m_cout<<"[DeBruijnAssembler,wordId] "<<a<<" -> "<<i<<endl;
	//string word=idToWord(i);
	//m_cout<<"Conversion <"<<a<<"> <"<<word<<">"<<endl;
	//m_cout<<"getLastSymbol "<<getLastSymbol(i)<<endl;
	return i;
}






// TODO: replace *= by << and /= by >>
string DeBruijnAssembler::idToWord(VERTEX_TYPE i,int wordSize){
	Logger&m_cout=*(this->m_cout);
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
			m_cout<<"j "<<j<<endl;
		}
		

	}
	return a;
}




// get the last symbol
// TODO: replace /= by >>1
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
	name<<m_assemblyDirectory<<"/assembly.graph,minimumCoverage="<<m_minimumCoverage<<",wordSize="<<m_wordSize<<".txt";
	m_graphFile=name.str();
}









bool DeBruijnAssembler::passFilterCoverage(vector<VERTEX_TYPE> path,int l,int C){
	map<int,int> votes;
	l=min(l,path.size()-1);
	for(int i=(int)path.size()-2-l+1;i<(int)path.size()-1;i++){
		vector<int> reads=m_graph[path[i]][path[i+1]];
		for(vector<int>::iterator j=reads.begin();j!=reads.end();j++){
			votes[*j]++;
		}
	}
	int ok=0;
	for(map<int,int>::iterator i=votes.begin();i!=votes.end();i++){
		if(i->second>=l)
			ok++;
	}
	//(*m_cout)<<ok<<" are ok"<<endl;
	return ok>=C;
}


vector<vector<VERTEX_TYPE> >DeBruijnAssembler::Filter_Remove_Smaller_Duplicates(vector<vector<VERTEX_TYPE> > contigs){
	//(*m_cout)<<"list "<<contigs.size()<<endl;
	(*m_cout)<<"[Filter_Remove_Smaller_Duplicates] "<<contigs.size()<<endl;
	if(contigs.size()==1)
		return contigs;
	vector<vector<VERTEX_TYPE> >filteredContigs;
	map<int,set<VERTEX_TYPE> > dictionnary;
	map<VERTEX_TYPE,vector<int> > walksIndex;
	set<int> eliminatedContigs;
	int id=0;
	// fill dictionnary
	// and index the last edge of contigs
	(*m_cout)<<"Indexing.."<<endl;
	for(vector<vector<VERTEX_TYPE> >::iterator i=contigs.begin();i!=contigs.end();i++){
		for(vector<VERTEX_TYPE>::iterator k=(*i).begin();k!=(*i).end();k++){
			dictionnary[id].insert(*k);
			dictionnary[id].insert(wordId(reverseComplement(idToWord(*k,m_wordSize)).c_str()));
			walksIndex[*k].push_back(id);
			walksIndex[ wordId(reverseComplement(idToWord(*k,m_wordSize)).c_str())].push_back(id);
		}
		id++;
	}

	(*m_cout)<<"Done.."<<endl;
	int progress=0;
	for(vector<vector<VERTEX_TYPE> >::iterator i=contigs.begin();i!=contigs.end();i++){
		//(*m_cout)<<(*i).size()<<endl;
		if(progress%100==0)
			(*m_cout)<<progress<<" / "<<contigs.size()<<endl;
		bool isDuplicate=false;
		vector<int> contigsToCheck=walksIndex[(*i)[(*i).size()/2]];
		for(vector<int> ::iterator j=contigsToCheck.begin();j!=contigsToCheck.end();j++){
			int otherContig=*j;
			if(otherContig!=progress&&eliminatedContigs.count(otherContig)==0){
				int notFound=0;
				for(vector<VERTEX_TYPE>::iterator k=(*i).begin();k!=(*i).end();k++){
					if(dictionnary[otherContig].count(*k)==0){
						notFound++;
					}
				}
				isDuplicate=(notFound<=2*m_wordSize);
				if(isDuplicate)
					break;
			}
		}
		if(!isDuplicate){
			filteredContigs.push_back(*i);
		}else{
			eliminatedContigs.insert(progress);
		}
		progress++;
	}
	(*m_cout)<<progress<<" / "<<contigs.size()<<endl;
	(*m_cout)<<contigs.size()<<" -> "<<filteredContigs.size()<<endl;
	return filteredContigs;
}



// new version of the assembler
void DeBruijnAssembler::run_New_Algorithm_Assembler_20090102(){

	vector<vector<VERTEX_TYPE> > assemblyContigs;
	
	vector<VERTEX_TYPE> withoutParents;
	for(map<VERTEX_TYPE,map<VERTEX_TYPE,vector<int> > >::iterator i=m_graph.begin();i!=m_graph.end();i++){
		VERTEX_TYPE prefix=i->first;
		if(m_vertex_parents[prefix].size()==0){
			withoutParents.push_back(prefix);
		}
	}

	int progress=0;
	Logger&m_cout=*(this->m_cout);	
	for(vector<VERTEX_TYPE>::iterator i=withoutParents.begin();
		i!=withoutParents.end();i++){
		map<VERTEX_TYPE,int> visits;
		progress++;
		m_cout<<endl;
		VERTEX_TYPE prefix=*i;
		m_cout<<"Starting to walk, source "<<progress<<" / "<<withoutParents.size()<<endl;
		m_cout<<assemblyContigs.size()<<" contigs"<<endl;
		vector<VERTEX_TYPE>path;
		path.push_back(prefix);
		vector< VERTEX_TYPE> localContig=contig_From_SINGLE(prefix,visits,path);
		m_cout<<localContig.size()<<" vertices "<<endl;
		int nucleotides=localContig.size()-m_wordSize+1;
		if(nucleotides<m_minimumContigSize)
			continue;
		assemblyContigs.push_back(localContig);
	}
	m_contig_paths=Filter_Remove_Smaller_Duplicates(assemblyContigs);

}

// get a walk from a vertex, with a path, up to maxSize,
// the path is not added in the walk
vector<VERTEX_TYPE> DeBruijnAssembler::getWalk(VERTEX_TYPE prefix,vector<VERTEX_TYPE>path,int length){
	vector<VERTEX_TYPE> subPath;
	path.push_back(prefix);
	subPath.push_back(prefix);
	while((int)nextVertices(prefix,path,m_windowSize,m_minimumCoverage).size()==1&&(int)subPath.size()<=length){
		prefix=nextVertices(prefix,path,m_windowSize,m_minimumCoverage)[0];
		path.push_back(prefix);
		subPath.push_back(prefix);
	}
	return subPath;
}

// remove bubble, if any
// remote tips also..
vector<VERTEX_TYPE> DeBruijnAssembler::removeBubblesAndTips(vector<VERTEX_TYPE> vertices,vector<VERTEX_TYPE>path){
	int maxDifference=5;
	int maxSize=2*m_wordSize+maxDifference;
	if(vertices.size()==2){
		vector<VERTEX_TYPE> n1=getWalk(vertices[0],path,maxSize);
		vector<VERTEX_TYPE> n2=getWalk(vertices[1],path,maxSize);
		set<VERTEX_TYPE> n1_Table;
		for(int i=0;i<(int)n1.size();i++){
			n1_Table.insert(n1[i]);
		}
		bool bubble=false;
		for(int i=0;i<(int)n2.size();i++){
			if(n1_Table.count(n2[i])>0){
				(*m_cout)<<i<<" BUBBLE"<<endl;
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
			vector<VERTEX_TYPE> subPath=getWalk(vertices[i],path,maxSize);
			if((int)subPath.size()<2*m_wordSize&&nextVertices(subPath[subPath.size()-1],subPath,m_windowSize,m_minimumCoverage).size()==0){
				(*m_cout)<<"TIP "<<getWalk(vertices[i],path,maxSize).size()<<" "<<idToWord(vertices[i],m_wordSize)<<endl;
			}else{
				withoutTips.push_back(vertices[i]);
			}
		}
		return withoutTips;
	}
	return vertices;
}

vector<VERTEX_TYPE> DeBruijnAssembler::contig_From_SINGLE(VERTEX_TYPE prefix,map<VERTEX_TYPE,int> visits,vector<VERTEX_TYPE> path){
	int MAX_VISITS=5;
	int l=m_windowSize;
	int C=m_minimumCoverage;
	if(visits[prefix]>MAX_VISITS){
		(*m_cout)<<"REPEAT? "<<idToWord(prefix,m_wordSize)<<endl;
		return path;
	}
	visits[prefix]++;

	(*m_cout)<<"Depth: "<<path.size()<<endl;
	vector<VERTEX_TYPE> prefixNextVertices=removeBubblesAndTips(nextVertices(prefix,path,l,C),path);
	//(*m_cout)<<"Only 1"<<endl;
	while(prefixNextVertices.size()==1){
		prefix=prefixNextVertices[0];
		if(visits[prefix]>MAX_VISITS){
			(*m_cout)<<"REPEAT? "<<idToWord(prefix,m_wordSize)<<endl;
			return path;
		}
		visits[prefix]++;

		path.push_back(prefix);
		prefixNextVertices=removeBubblesAndTips(nextVertices(prefix,path,l,C),path);
	}


	if(prefixNextVertices.size()>1){
		(*m_cout)<<"more than 1!"<<endl;

		// try to let go a branch
		int optimizedL=l;
		int optimizedC=C;
		vector<VERTEX_TYPE> nextOptimizedVertices=nextVertices(prefix,path,optimizedL,optimizedC);
		while(nextOptimizedVertices.size()>1){
			if(optimizedL>200)
				break;
			optimizedL++;
			optimizedC++;
			int maxC=7;
			if(optimizedC>maxC)
				optimizedC=maxC;
			nextOptimizedVertices=nextVertices(prefix,path,optimizedL,optimizedC);
		}
		(*m_cout)<<"Optimization completed"<<endl;
		if(nextOptimizedVertices.size()==0){
			(*m_cout)<<"Nothing found suitable"<<endl;
			return path;
		}
		if(nextOptimizedVertices.size()>1){
			(*m_cout)<<"Optimization failed"<<endl;
			return path;
		}
		(*m_cout)<<"Optimization successful!"<<endl;
		vector<VERTEX_TYPE> sPath=path;
		VERTEX_TYPE suffix=nextOptimizedVertices[0];
		sPath.push_back(suffix);
		vector<VERTEX_TYPE> sRecPath=contig_From_SINGLE(suffix,visits,sPath);
		return sRecPath;
	}
	(*m_cout)<<"DEAD END "<<endl;
	// show path
	for(int i=min(path.size(),20);i>=1;i--){
		if(path.size()<2){
			break;
		}
		VERTEX_TYPE a=path[path.size()-i-1];
		VERTEX_TYPE b=path[path.size()-i+1-1];
		(*m_cout)<<idToWord(a,m_wordSize)<<" -> "<<idToWord(b,m_wordSize)<<" ";
		vector<int>*theReads=&(m_graph[a][b]);
		(*m_cout)<<theReads->size();
		for(int j=0;j<(int)theReads->size();j++){
			(*m_cout)<<" "<<(*theReads)[j];
		}
		(*m_cout)<<endl;
	
	}
	(*m_cout)<<prefixNextVertices.size()<<" CHOICES"<<endl;
	map<VERTEX_TYPE,vector<int> >* dataStructure=&(m_graph[prefix]);
	for(map<VERTEX_TYPE,vector<int> >::iterator i=dataStructure->begin();i!=dataStructure->end();i++){
		(*m_cout)<<idToWord(prefix,m_wordSize)<<" -> "<<idToWord(i->first,m_wordSize)<<" ";
		(*m_cout)<<i->second.size();
		for(int j=0;j<(int)i->second.size();j++){
			(*m_cout)<<" "<<i->second[j];
		}
		(*m_cout)<<endl;
	}
	return path;
}

vector<VERTEX_TYPE> DeBruijnAssembler::nextVertices(VERTEX_TYPE prefix,vector<VERTEX_TYPE> path,int l,int C){
	vector<VERTEX_TYPE> next;

	map<VERTEX_TYPE,vector<int> >*prefixMap=&(m_graph[prefix]);
        for(map<VERTEX_TYPE,vector<int> >::iterator i=prefixMap->begin();i!=prefixMap->end();i++){
		VERTEX_TYPE suffix=i->first;
		vector<VERTEX_TYPE> sPath=path;
		if(m_graph[prefix][suffix].size()>500)
			continue;
		sPath.push_back(suffix);
		if(!passFilterCoverage(sPath,l,C))
			continue;
		next.push_back(suffix);
	}
	return next;
}

void DeBruijnAssembler::setWindowSize(int windowSize){
	m_windowSize=windowSize;
}

void DeBruijnAssembler::setUseCache(string useCache){
	m_useCache=useCache;
}

// THE END IS HERE  ------------> .
