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

DeBruijnAssembler::DeBruijnAssembler(Logger*m_cout){
	this->m_cout=m_cout;
	m_COLOR_REPEAT=-87;
	m_COLOR_NOT_ASSEMBLED=-2;
	m_COLOR_IN_PROGRESS_ASSEMBLER=-3;
	m_COLOR_IN_PROGRESS_SOLVER=-4;
	m_COLOR_DISCARDED=-99;
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

void DeBruijnAssembler::build_From_Scratch(){
	Logger&m_cout=*(this->m_cout);
	m_cout<<endl;
	m_cout<<"Collecting high-quality mers from reads"<<endl;
	map<unsigned long int,int> words;
	int last_vertices_size=-1;
	for(int i=0;i<(int)m_reads->size();i++){
		if(i%10000==0){
			m_cout<<"Reads: "<<i<<" / "<<m_reads->size()<<endl;
		}
		vector<unsigned long int>highQualityMers=getHighQualityMers(i);
		for(vector<unsigned long int>::iterator iteratorMer=highQualityMers.begin();iteratorMer!=highQualityMers.end();iteratorMer++){
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
	for(map<unsigned long int,int>::iterator i=words.begin();i!=words.end();i++){
		processed++;
		if(i->second>=m_minimumCoverage){
			solid++;
			unsigned long int w=i->first;
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
	m_cout<<"Indexing solid mers in reads."<<endl; // <-------
	for(int readId=0;readId<(int)m_reads->size();readId++){
		if(readId%10000==0)
			m_cout<<readId<<" / "<<m_reads->size()<<endl;
		vector<unsigned long int> highQualityMers=getHighQualityMers(readId);
		for(vector<unsigned long int>::iterator merIterator=highQualityMers.begin();
			merIterator!=highQualityMers.end();merIterator++){
			if(m_solidMers.count(*merIterator)>0){
				m_mer_to_read_table[*merIterator].push_back(readId);
			}
		}
	}

	m_cout<<m_reads->size()<<" / "<<m_reads->size()<<endl;

	m_cout<<endl;
	m_cout<<"Building de Bruijn graph with solid mers"<<endl;
	for(set<unsigned long int>::iterator i=m_solidMers.begin();i!=m_solidMers.end();i++){
		//m_cout<<endl;
		if(m_graph_mers.size()%100000==0){
			m_cout<<"Building graph: "<<m_graph_mers.size()<<" / "<<m_solidMers.size()<<endl;
		}
		string wholeWord=idToWord(*i,m_wordSize+1);
		unsigned long int prefix=wordId(wholeWord.substr(0,m_wordSize).c_str());
		unsigned long int suffix=wordId(wholeWord.substr(1,m_wordSize).c_str());
		m_graph_mers.insert(prefix);
		m_graph_mers.insert(suffix);
		for(vector<int>::iterator j=m_mer_to_read_table[*i].begin();
			j!=m_mer_to_read_table[*i].end();j++)
			m_graph[prefix][suffix].push_back(*j);
		m_vertex_parents[suffix].push_back(prefix);
	}

	m_cout<<"Building graph: "<<m_solidMers.size()<<" / "<<m_solidMers.size()<<endl;
	//m_cout<<"the graph is ready"<<endl;
	//m_cout<<"Clearing solid mers"<<endl;
	//solidMers.clear();
}

void DeBruijnAssembler::writeGraph(){
	ofstream graph(m_graphFile.c_str());
	graph<<m_solidMers.size()<<" "<<endl;
	for(map<unsigned long int,map<unsigned long int,vector<int> > >::iterator i=m_graph.begin();i!=m_graph.end();i++){
		unsigned long int prefix=i->first;
		for(map<unsigned long int,vector<int> >::iterator j=i->second.begin();j!=i->second.end();j++){
			unsigned long int suffix=j->first;
			vector<int> reads=j->second;
			graph<<prefix<<" "<<suffix<<" "<<reads.size();
			for(int k=0;k<(int)reads.size();k++){
				graph<<" "<<reads[k];
			}
			graph<<endl;
		}
	}
	graph.close();
}

// fill m_graph, m_graph_mers, and m_vertex_parents
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
		unsigned long int prefix;
		unsigned long int suffix;
		int reads;
		graph>>prefix>>suffix>>reads;
		m_graph_mers.insert(prefix);
		m_graph_mers.insert(suffix);
		m_vertex_parents[suffix].push_back(prefix);
		for(int j=0;j<reads;j++){
			int read;
			graph>>read;
			m_graph[prefix][suffix].push_back(read);
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
	}
}

void DeBruijnAssembler::setMinimumContigSize(int minimumContigSize){
	m_minimumContigSize=minimumContigSize;
}

string DeBruijnAssembler::pathToDNA(vector<unsigned long int> path){
	ostringstream contigSequence;
	for(vector<unsigned long int>::iterator i=path.begin();		i!=path.end();i++){
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
	unsigned long int stat_minimumContigSize=9999999999;
	int stat_maximumContigSize=0;
	for(int i=0;i<(int)contigs.size();i++){
		string sequenceDNA=contigs[i];
		totalLength+=sequenceDNA.length();
		contigSizes.insert(sequenceDNA.length());
		if(sequenceDNA.length()<stat_minimumContigSize)
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
		for(int j=50;j<=90;j+=10){
			if(cumulativeSize>=totalLength*j/100.0&&stat_nx_done.count(j)==0){
				stat_nx_done.insert(j);
				m_cout<<"N"<<j<<" size: "<<*i<<endl;
				break;
			}
		}
	}
	m_cout<<"See http://www.cbcb.umd.edu/research/castats.shtml#N50 for the definition of N50."<<endl;
}






// TODO: remove this, it is deprecated.
// returns a trivial for a vertex
string DeBruijnAssembler::getTrivialPath(unsigned long pathVertex,unsigned long int vertex){
	ostringstream path;
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
void DeBruijnAssembler::solveMultiPath(unsigned long int vertex,int depth){
	//Logger&m_cout=*(this->m_cout);

	if(depth<0)
		return;
	vector<unsigned long int> neighboursVector=getNeighbours(vertex);
	//m_cout<<endl;
	for(vector<unsigned long int>::iterator solverIterator=neighboursVector.begin();
		solverIterator!=neighboursVector.end();solverIterator++){
		//m_cout<<vertex<<" -> "<<*solverIterator<<endl;
		if(*solverIterator==vertex){
			m_edge_states[*solverIterator][*solverIterator]=m_COLOR_DISCARDED;
			//m_assembled_edges[*solverIterator].insert(*solverIterator);
		}else if(m_edge_states[vertex][*solverIterator]==m_COLOR_IN_PROGRESS_SOLVER){

		}else if(m_edge_states[vertex][*solverIterator]==m_COLOR_DISCARDED){
		}else{
			m_edge_states[vertex][*solverIterator]=m_COLOR_IN_PROGRESS_SOLVER;
		//if(m_assembled_edges
			solveMultiPath(*solverIterator,depth-1);
		}
	}

	// : remove tips

	neighboursVector=getNeighbours(vertex);
	if(neighboursVector.size()>1){
		// remove all with length < 2k (Like Velvet)
		for(int i=0;i<(int)neighboursVector.size();i++){
			string path=getTrivialPath(neighboursVector[i],vertex);
			if((int)path.length()< 2 * m_wordSize)
				m_edge_states[vertex][neighboursVector[i]]=m_COLOR_DISCARDED;
		}
	}

	// remove bubbles
	neighboursVector=getNeighbours(vertex);

	int maximumDifferenceInBubbles=3;
	if(neighboursVector.size()==2){
		(*m_solver_log)<<"2 alternative paths"<<endl;
		string path1=getTrivialPath(neighboursVector[0],vertex);
		string path2=getTrivialPath(neighboursVector[1],vertex);
		string lastKMerPath1=path1.substr(path1.length()-m_wordSize,m_wordSize);
		string lastKMerPath2=path2.substr(path2.length()-m_wordSize,m_wordSize);
		if(lastKMerPath1==lastKMerPath2){// bubble
			int diff=path1.length()-path2.length();
			if(diff<0)
				diff=-diff;
			if(diff<=maximumDifferenceInBubbles){
				(*m_solver_log)<<endl;
				// TODO: maybe chop it all
				(*m_solver_log)<<"Bubble found!, removing it... "<<path1.length()<<" "<<path2.length()<<endl;
				// TODO SNP detection goes here ?? or it can be sequencing errors also..
				m_edge_states[vertex][neighboursVector[1]]=m_COLOR_DISCARDED;
				(*m_solver_log)<<path1<<endl;
				(*m_solver_log)<<path2<<endl;
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

vector<unsigned long int> DeBruijnAssembler::getNeighbours(unsigned long int vertex){
        vector<unsigned long int> neighboursVector;


        for(map<unsigned long int,vector<int> >::iterator i=m_graph[vertex].begin();i!=m_graph[vertex].end();i++){
                unsigned long int currentNeighbour=i->first;
                if((m_edge_states[vertex][currentNeighbour]==m_COLOR_NOT_ASSEMBLED||
		m_edge_states[vertex][currentNeighbour]==m_COLOR_IN_PROGRESS_SOLVER||
		m_edge_states[vertex][currentNeighbour]==m_COLOR_IN_PROGRESS_ASSEMBLER)&&currentNeighbour!=vertex){
                //if(m_colors[currentNeighbour]==m_colors[vertex]&&vertex!=currentNeighbour){
                        neighboursVector.push_back(currentNeighbour);
                }
        }


        return neighboursVector;
}


 
void DeBruijnAssembler::setAssemblyDirectory(string assemblyDirectory){
	m_assemblyDirectory=assemblyDirectory;
	ostringstream name;
	name<<m_assemblyDirectory<<"/assembly.graph,minimumCoverage="<<m_minimumCoverage<<",wordSize="<<m_wordSize<<",minimumQuality="<<m_minimumQuality<<".txt";
	m_graphFile=name.str();
}



void DeBruijnAssembler::run_New_Algorithm_Assembler(){
	Logger&m_cout=*(this->m_cout);	
	string solveLogFile=m_assemblyDirectory+"/Solver.txt";
	m_solver_log=new ofstream(solveLogFile.c_str());

	int total_Edges=0;
	// set initial color
	for(map<unsigned long int,map<unsigned long int,vector<int> > >::iterator i=m_graph.begin();i!=m_graph.end();i++){

		for(map<unsigned long int,vector<int> >::iterator j=i->second.begin();j!=i->second.end();j++){
			m_edge_states[i->first][j->first]=m_COLOR_NOT_ASSEMBLED;
			total_Edges++;
		}
	}
	m_cout<<total_Edges<<" edges to assemble."<<endl;

	/*
	m_cout<<"fixing graph"<<endl;

	//solve strange structures
	int progress=0;
	int maximumDepth=1000;
	for(set<unsigned long int>::iterator i=m_graph_mers.begin();i!=m_graph_mers.end();i++){

		if(progress%10000==0){
			m_cout<<progress<<" / "<<total_Edges<<endl;
		}
		// remove bubbles and tips
		solveMultiPath(*i,maximumDepth);
		progress++;
	}

	m_cout<<progress<<" / "<<progress<<endl;

	*/
	vector<vector<unsigned long int> > assemblyContigs;
	map<unsigned long int,int> visits;
	
	vector<unsigned long int> withoutParents;
	for(set<unsigned long int>::iterator i=m_graph_mers.begin();
		i!=m_graph_mers.end();i++){
		bool ok=true;
		for(vector<unsigned long int>::iterator j=
			m_vertex_parents[*i].begin();
			j!=m_vertex_parents[*i].end();j++){
			if(m_edge_states[*j][*i]!=m_COLOR_DISCARDED){
				ok=false;
			}
		}
		if(ok==true){
			withoutParents.push_back(*i);
		}
	}

	string debugFileDna=m_assemblyDirectory+"/Walks.seq";
	string debugFile=m_assemblyDirectory+"/Walks.debug";
	ofstream debugContigs(debugFile.c_str());
	ofstream debugDna(debugFileDna.c_str());
	int progress=0;
	for(vector<unsigned long int>::iterator i=withoutParents.begin();
		i!=withoutParents.end();i++){
		progress++;
		m_cout<<endl;
		m_cout<<"Starting to walk."<<endl;
		m_cout<<progress<<" / "<<withoutParents.size()<<endl;
		m_cout<<assemblyContigs.size()<<" contigs"<<endl;
		vector<unsigned long int> previousEdges;
		//m_cout<<"From "<<idToWord(*i,m_wordSize)<<endl;
		map<int,map<unsigned long int,int> > readUsages;
		vector<vector< unsigned long int> > localContigs=Filter_Remove_Smaller_Duplicates(removeSmallContigs(contigs_From(*i,visits,0,previousEdges,readUsages)));

		for(vector<vector< unsigned long int> >::iterator j = localContigs.begin();j!=localContigs.end();j++){
			assemblyContigs.push_back(*j);
			debugContigs<<"Walk"<<assemblyContigs.size()<<" "<<(*j).size()<<" edges"<<endl;
			ostringstream sequence;
			
			for(vector<unsigned long int>::iterator k=(*j).begin();k!=(*j).end();k++){
				unsigned long int prefix=toPrefix(*k);
				unsigned long int suffix=toSuffix(*k);
				if(k==(*j).begin()){
					sequence<< idToWord(prefix,m_wordSize);
				}
				sequence<< getLastSymbol(suffix);
				debugContigs<<" "<<*k;
			}
			debugDna<<">Walk"<<assemblyContigs.size()<<" "<<sequence.str().length()<<" nucleotides"<<endl;
			int columns=70;
			string seq=sequence.str();
			int r=0;
			while((int)r<(int)seq.length()){
				debugDna<< seq.substr(r,columns)<<endl;
				r+=columns;
			}
			debugContigs<<endl;
			m_cout<<(*j).size()<<" edges"<<endl;
		}
	}
	debugContigs.close();
	debugDna.close();
	m_cout<<progress<<" / "<<withoutParents.size()<<endl;
	// push contigs in m_contig_paths
	
	int contigId=0;
	assemblyContigs=Filter_Remove_Smaller_Duplicates(assemblyContigs);
	(m_cout)<<"Now "<<assemblyContigs.size()<<endl;
	for(vector<vector<unsigned long int> >::iterator i=assemblyContigs.begin();i!=assemblyContigs.end();i++){
		m_contig_paths[contigId].push_back(toPrefix(*(*i).begin()));
		for(vector<unsigned long int>::iterator j=(*i).begin();
			j!=i->end();j++){
			m_contig_paths[contigId].push_back(toSuffix((*j)));
		}
		contigId++;
		m_cout<<(*i).size()<<" edges, "<<(*i).size()+m_wordSize-0<<endl;
	}
	m_solver_log->close();
}

unsigned long int DeBruijnAssembler::toKPlusOneMer(unsigned long int prefix,
	unsigned long int suffix){
	string token=idToWord(prefix,m_wordSize)+getLastSymbol(suffix);
	return wordId(token.c_str());
}

// return a list of edges,
// each edge is encoded in an unsigned long int
vector<vector<unsigned long int> > DeBruijnAssembler::contigs_From(unsigned long int prefix,map<unsigned long int,int> visits,int depth,vector<unsigned long int>previousEdges,map<int,map<unsigned long int,int> > readUsages){
	//cout<<"contigs from "<<idToWord(prefix,m_wordSize)<<endl;
	vector<vector<unsigned long int> > contigs;

	(*m_cout)<<"Depth: "<<depth<<endl;
	int maxVisits=5;
	int l=150;
	int C=m_minimumCoverage;
	int maxZ=2;
	// accelerator

	vector<unsigned long int> acceleratedPath=previousEdges;


	#define ACCELERATOR NOW   // comment this line to remove the accelerator
					// the accelerator avoids recursive calls when possible.

	#ifdef ACCELERATOR
	map<int,map<unsigned long int,int> > localReadUsages=readUsages;

	while(getNeighbours(prefix).size()==1){
		unsigned long int suffix=getNeighbours(prefix)[0];
		vector<unsigned long int> pathToTest=acceleratedPath;
		pathToTest.push_back(toKPlusOneMer(prefix,suffix));
		
		if(visits[toKPlusOneMer(prefix,suffix)]>maxVisits||!passFilterZ(pathToTest,l,C,maxZ)){
			contigs.push_back(acceleratedPath);
			return contigs;
		}
		//(*m_cout)<<"Extending "<<prefix<<" -> "<<suffix<<" visits: "<< visits[toKPlusOneMer(prefix,suffix)] <<endl;
		depth++;
		visits[toKPlusOneMer(prefix,suffix)]++;
		acceleratedPath=pathToTest;
		prefix=suffix;
	}
	#endif

	//(*m_cout)<<"depth (after trivial run): "<<depth<<endl;
	vector<unsigned long int> children=getNeighbours(prefix);
	//(*m_cout)<<"Next: "<<children.size()<<endl;
	for(vector<unsigned long int>::iterator i=children.begin();
		i!=children.end();i++){
		unsigned long int suffix=*i;
		//(*m_cout)<<"Checking: "<<prefix<<" -> "<<suffix<<endl;
		if(visits[toKPlusOneMer(prefix,suffix)]>maxVisits)
			continue;
		visits[toKPlusOneMer(prefix,*i)]++;
		vector<unsigned long int>path=acceleratedPath;
		path.push_back(toKPlusOneMer(prefix,suffix));
		map<int,map<unsigned long int,int> > localReadUsages2=localReadUsages;
		if(passFilterZ(path,l,C,maxZ)){
			//(*m_cout)<<"Going in,.."<<endl;
			vector<vector<unsigned long int> > sContigs=(contigs_From(suffix,visits,depth+1,path,localReadUsages2));
			for(vector<vector<unsigned long int> >::iterator j=sContigs.begin();j!=sContigs.end();j++){
				contigs.push_back(*j);
			}
		}
	}
	if(contigs.size()==0){
		//(*m_cout)<<"LEAF: no children."<<endl;
		contigs.push_back(acceleratedPath);
	}
	vector<vector<unsigned long int> >output=contigs ;//Filter_Remove_Smaller_Duplicates(contigs);
	//cout<<output.size()<<endl;
	return output;
}

bool DeBruijnAssembler::passFilterCoverage(vector<unsigned long int> path,int l,int C){
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

bool DeBruijnAssembler::passFilterZ(vector<unsigned long int>contig,int l,int C,int maxZ){
	l=min(l,contig.size());
	maxZ=min(maxZ,contig.size());
	for(int Z=1;Z<=maxZ;Z++){
		if(passFilter(contig,l,C,Z)){
			//(*m_cout)<<"ok..."<<endl;
			return true;
		}
	}
	//(*m_cout)<<"Failed..."<<endl;
	return false;
}

bool DeBruijnAssembler::passFilter(vector<unsigned long int>contig,int l,int C,int Z){
	// only check the last part, l edges
	//for(int i=contig.length()-1;i<verificationSize;i++){
	map<int,int> votesLeft;
	map<int,int> votesRight;
	vector<int> indexesToCheckLeft;
	vector<int> indexesToCheckRight;
	for(int toAdd=0;toAdd<Z;toAdd++){
		indexesToCheckLeft.push_back(contig.size()-l+toAdd);
		indexesToCheckRight.push_back(contig.size()-2-toAdd);
	}
	for(vector<int>::iterator i=indexesToCheckLeft.begin();i!=indexesToCheckLeft.end();i++){
		unsigned long int prefix=contig[*i];
		unsigned long int suffix=contig[*i+1];
		for(vector<int>::iterator k=m_graph[prefix][suffix].begin();k!=m_graph[prefix][suffix].end();k++){
			votesLeft[*k]++;
		}
	}
	for(vector<int>::iterator i=indexesToCheckRight.begin();i!=indexesToCheckRight.end();i++){
		unsigned long int prefix=contig[*i];
		unsigned long int suffix=contig[*i+1];
		for(vector<int>::iterator k=m_graph[prefix][suffix].begin();k!=m_graph[prefix][suffix].end();k++){
			votesRight[*k]++;
		}
	}
	int total=0;
	for(map<int,int>::iterator i=votesLeft.begin();i!=votesLeft.end();i++){
		if(votesRight.count(i->first)>0)
			total++;
	}
	//(*m_cout)<<"TOTAL: "<<total<<endl;
	return total>=C;
}


vector<vector<unsigned long int> >DeBruijnAssembler::removeSmallContigs(vector<vector<unsigned long int> > contigs){
	vector<vector<unsigned long int> >filteredContigs;
	for(vector<vector<unsigned long int> >::iterator i=contigs.begin();i!=contigs.end();i++){
		int nucleotides=(*i).size()-m_wordSize+1;
		if(nucleotides>=m_minimumContigSize){
			(*m_cout)<<nucleotides<<endl;
			filteredContigs.push_back(*i);
		}
	}
	(*m_cout)<<"[removeSmallContigs] "<<contigs.size()<<" -> "<<filteredContigs.size()<<endl;
	return filteredContigs;
}


// TODO: filter...
vector<vector<unsigned long int> >DeBruijnAssembler::Filter_Remove_Smaller_Duplicates(vector<vector<unsigned long int> > contigs){
	//(*m_cout)<<"list "<<contigs.size()<<endl;
	(*m_cout)<<"[Filter_Remove_Smaller_Duplicates] "<<contigs.size()<<endl;
	if(contigs.size()==1)
		return contigs;
	vector<vector<unsigned long int> >filteredContigs;
	map<int,set<unsigned long int> > dictionnary;
	map<unsigned long int,vector<int> > walksIndex;
	set<int> eliminatedContigs;
	int id=0;
	// fill dictionnary
	// and index the last edge of contigs
	(*m_cout)<<"Indexing.."<<endl;
	for(vector<vector<unsigned long int> >::iterator i=contigs.begin();i!=contigs.end();i++){
		for(vector<unsigned long int>::iterator k=(*i).begin();k!=(*i).end();k++){
			dictionnary[id].insert(*k);
			dictionnary[id].insert(wordId(reverseComplement(idToWord(*k,m_wordSize)).c_str()));
			walksIndex[*k].push_back(id);
			walksIndex[ wordId(reverseComplement(idToWord(*k,m_wordSize)).c_str())].push_back(id);
		}
		id++;
	}

	(*m_cout)<<"Done.."<<endl;
	int progress=0;
	for(vector<vector<unsigned long int> >::iterator i=contigs.begin();i!=contigs.end();i++){
		//(*m_cout)<<(*i).size()<<endl;
		if(progress%100==0)
			(*m_cout)<<progress<<" / "<<contigs.size()<<endl;
		bool isDuplicate=false;
		vector<int> contigsToCheck=walksIndex[(*i)[(*i).size()/2]];
		for(vector<int> ::iterator j=contigsToCheck.begin();j!=contigsToCheck.end();j++){
			int otherContig=*j;
			if(otherContig!=progress&&eliminatedContigs.count(otherContig)==0){
				int notFound=0;
				for(vector<unsigned long int>::iterator k=(*i).begin();k!=(*i).end();k++){
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

unsigned long int DeBruijnAssembler::toPrefix(unsigned long int edge){
	return wordId(idToWord(edge,m_wordSize+1).substr(0,m_wordSize).c_str());
}

unsigned long int DeBruijnAssembler::toSuffix(unsigned long int edge){
	return wordId(idToWord(edge,m_wordSize+1).substr(1,m_wordSize).c_str());
}

// new version of the assembler
void DeBruijnAssembler::run_New_Algorithm_Assembler_20090102(){

	string solveLogFile=m_assemblyDirectory+"/Solver.txt";
	m_solver_log=new ofstream(solveLogFile.c_str());
	vector<vector<unsigned long int> > assemblyContigs;
	
	vector<unsigned long int> withoutParents;
	for(set<unsigned long int>::iterator i=m_graph_mers.begin();
		i!=m_graph_mers.end();i++){
		bool ok=true;
		if(m_vertex_parents[*i].size()>0)
			ok=false;
		if(ok==true){
			withoutParents.push_back(*i);
		}
	}

	int progress=0;
	Logger&m_cout=*(this->m_cout);	
	for(vector<unsigned long int>::iterator i=withoutParents.begin();
		i!=withoutParents.end();i++){
		map<unsigned long int,int> visits;
		progress++;
		m_cout<<endl;
		unsigned long int prefix=*i;
		m_cout<<"Starting to walk."<<endl;
		m_cout<<progress<<" / "<<withoutParents.size()<<endl;
		m_cout<<assemblyContigs.size()<<" contigs"<<endl;
		vector<unsigned long int>path;
		path.push_back(prefix);
		vector< unsigned long int> localContig=contig_From_SINGLE(prefix,visits,path);
		m_cout<<localContig.size()<<" vertices "<<endl;
		assemblyContigs.push_back(localContig);
	}
	m_contig_paths=Filter_Remove_Smaller_Duplicates(removeSmallContigs(assemblyContigs));

	m_solver_log->close();
}

// get a walk from a vertex, with a path, up to maxSize,
// the path is not added in the walk
vector<unsigned long int> DeBruijnAssembler::getWalk(unsigned long int prefix,vector<unsigned long int>path,int length){
	vector<unsigned long int> subPath;
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
vector<unsigned long int> DeBruijnAssembler::removeBubbles(vector<unsigned long int> vertices,vector<unsigned long int>path){
	int maxDifference=5;
	int maxSize=2*m_wordSize+maxDifference;
	if(vertices.size()==2){
		vector<unsigned long int> n1=getWalk(vertices[0],path,maxSize);
		vector<unsigned long int> n2=getWalk(vertices[1],path,maxSize);
		set<unsigned long int> n1_Table;
		for(int i=0;i<(int)n1.size();i++){
			n1_Table.insert(n1[i]);
		}
		bool bubble=false;
		for(int i=0;i<(int)n2.size();i++){
			if(n1_Table.count(n2[i])>0){
				(*m_cout)<<i<<endl;
				bubble=true;
				break;
			}
		}
		if(bubble==true){
			(*m_cout)<<"BUBBLE"<<endl;
			vector<unsigned long int> withoutBubbles;
			withoutBubbles.push_back(vertices[0]);
			return withoutBubbles;
		}
	}
	if(vertices.size()>1){
		vector<unsigned long int> withoutTips;
		for(int i=0;i<(int)vertices.size();i++){
			if((int)getWalk(vertices[i],path,maxSize).size()<2*m_wordSize){
				(*m_cout)<<"TIP "<<getWalk(vertices[i],path,maxSize).size()<<endl;
			}else{
				withoutTips.push_back(vertices[i]);
			}
		}
		return withoutTips;
	}
	return vertices;
}

vector<unsigned long int> DeBruijnAssembler::contig_From_SINGLE(unsigned long int prefix,map<unsigned long int,int> visits,vector<unsigned long int> path){
	int MAX_VISITS=5;
	int l=m_windowSize;
	int C=m_minimumCoverage;
	if(visits[prefix]>MAX_VISITS){
		(*m_cout)<<"REPEAT? "<<idToWord(prefix,m_wordSize)<<endl;
		return path;
	}
	visits[prefix]++;

	(*m_cout)<<"Depth: "<<path.size()<<endl;
	vector<unsigned long int> prefixNextVertices=removeBubbles(nextVertices(prefix,path,l,C),path);
	//(*m_cout)<<"Only 1"<<endl;
	while(prefixNextVertices.size()==1){
		prefix=prefixNextVertices[0];
		if(visits[prefix]>MAX_VISITS){
			(*m_cout)<<"REPEAT? "<<idToWord(prefix,m_wordSize)<<endl;
			return path;
		}
		visits[prefix]++;

		path.push_back(prefix);
		prefixNextVertices=removeBubbles(nextVertices(prefix,path,l,C),path);
	}

	//(*m_cout)<<"now not 1, taking the longuest"<<endl;
/*
	vector<vector<unsigned long int> > paths;
	for(vector<unsigned long int>::iterator i=prefixNextVertices.begin();i!=prefixNextVertices.end();i++){
		// if the execution reaches this point, it means there are at least 2 choices
		vector<unsigned long int> sPath=path;
		unsigned long int suffix=*i;
		sPath.push_back(suffix);
		vector<unsigned long int> sRecPath=contig_From_SINGLE(suffix,visits,sPath);

		// TIPS HERE

		// BUBBLES

		paths.push_back(sRecPath);
	}
*/
	if(prefixNextVertices.size()>1){
		(*m_cout)<<"more than 1!"<<endl;
		(*m_solver_log)<<prefixNextVertices.size();
/*
		for(vector<vector<unsigned long int> >::iterator i=paths.begin();i!=paths.end();i++)
			(*m_solver_log)<<" "<<(*i).size();
	*/
		(*m_solver_log)<<endl;
		(*m_solver_log)<<"View different options"<<endl;
		/*
		for(vector<vector<unsigned long int> >::iterator i=paths.begin();i!=paths.end();i++)
			(*m_solver_log)<<(*i).size()-path.size()<<endl;
*/
/*
		for(int i=0;i<(int)paths.size();i++){
			(*m_solver_log)<<"Option "<<i+1<<endl;
			vector<unsigned long int> pathToPrint=paths[i];
			for(int j=0;j<(int)pathToPrint.size()-1;j++){
				(*m_solver_log)<<idToWord(pathToPrint[j],m_wordSize)<<" -> "<<idToWord(pathToPrint[j+1],m_wordSize)<<" ";
				vector<int> supportReads=m_graph[pathToPrint[j]][pathToPrint[j+1]];
				(*m_solver_log)<<supportReads.size();
				for(int k=0;k<(int)supportReads.size();k++){
					(*m_solver_log)<<" "<<supportReads[k];
				}
				if(j==(int)path.size()-1){
					(*m_solver_log)<<" <<<<<<<------- BREAK POINT ";
				}
				(*m_solver_log)<<endl;
			}
		}
*/
		// try to let go a branch
		int optimizedL=l;
		int optimizedC=C;
		vector<unsigned long int> nextOptimizedVertices=nextVertices(prefix,path,optimizedL,optimizedC);
		while(nextOptimizedVertices.size()>1){
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
		vector<unsigned long int> sPath=path;
		unsigned long int suffix=nextOptimizedVertices[0];
		sPath.push_back(suffix);
		vector<unsigned long int> sRecPath=contig_From_SINGLE(suffix,visits,sPath);
		return sRecPath;
	}
	(*m_cout)<<"DEAD END "<<endl;
	return path;
}

vector<unsigned long int> DeBruijnAssembler::nextVertices(unsigned long int prefix,vector<unsigned long int> path,int l,int C){
	vector<unsigned long int> next;

        for(map<unsigned long int,vector<int> >::iterator i=m_graph[prefix].begin();i!=m_graph[prefix].end();i++){
		unsigned long int suffix=i->first;
		vector<unsigned long int> sPath=path;
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
