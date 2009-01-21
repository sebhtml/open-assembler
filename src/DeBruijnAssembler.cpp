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
#include"Logger.hpp"
#include<iostream>
#include"DeBruijnAssembler.h"
#include<fstream>
#include<vector>



using namespace std;



// TODO: add paired end reads ?
// with 454 it is 3k

// TODO expand repeats using reads

DeBruijnAssembler::DeBruijnAssembler(Logger*m_cout){
	this->m_cout=m_cout;
	m_COLOR_REPEAT=-87;
	m_COLOR_NOT_ASSEMBLED=-2;
	m_COLOR_IN_PROGRESS_ASSEMBLER=-3;
	m_COLOR_IN_PROGRESS_SOLVER=-4;
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
		//m_cout<<endl;
		if(m_graph_mers.size()%100000==0){
			m_cout<<"Building graph: "<<m_graph_mers.size()<<" / "<<solidMers.size()<<endl;
		}
		string wholeWord=idToWord(*i,m_wordSize+1);
		unsigned long int prefix=wordId(wholeWord.substr(0,m_wordSize).c_str());
		unsigned long int suffix=wordId(wholeWord.substr(1,m_wordSize).c_str());
		//unsigned long int prefixRev=wordId(reverseComplement(wholeWord.substr(0,m_wordSize)).c_str());
		//unsigned long int suffixRev=wordId(reverseComplement(wholeWord.substr(1,m_wordSize)).c_str());
		m_graph_mers.insert(prefix);
		m_graph_mers.insert(suffix);
		for(vector<int>::iterator j=m_mer_to_read_table[*i].begin();
			j!=m_mer_to_read_table[*i].end();j++)
			m_graph[prefix][suffix].push_back(*j);
		m_vertex_parents[suffix].push_back(prefix);
		m_vertex_children[prefix].push_back(suffix);
	}

	m_cout<<"Building graph: "<<solidMers.size()<<" / "<<solidMers.size()<<endl;
	m_cout<<"the graph is ready"<<endl;
	m_cout<<"Clearing solid mers"<<endl;
	solidMers.clear();
}



int DeBruijnAssembler::number_of_reads_that_Agree_on_path(vector<unsigned long int> path){
	map<int,int> votes;
	for(int i=0;i<(int)path.size()-1;i++){
		if(m_graph.count(path[i])==0||m_graph[path[i]].count(path[i+1])==0)
			continue;
		for(vector<int>::iterator j=m_graph[path[i]][path[i+1]].begin();j!=m_graph[path[i]][path[i+1]].end();j++){
			votes[*j]++;
		}
	}
	int i=0;
	for(map<int,int>::iterator j=votes.begin();j!=votes.end();j++)
		if(j->second==(int)path.size()-1)
			i++;
	return i;
}

bool DeBruijnAssembler::areRelatives(vector<unsigned long int> path1,vector<unsigned long int> path2){
	if(path2.size()>path1.size())
		return areRelatives(path2,path1);
	// path2 is the smallest
	//cout<<"?"<<endl;
	//cout<<pathToDNA(path1)<<endl;
	//cout<<pathToDNA(path2)<<endl;
	
	for(int i=0;i<(int)path1.size();i++){
		//cout<<path1[i]<<" ";
	}
	//cout<<endl;
	for(int i=0;i<(int)path2.size();i++){
		//cout<<path2[i]<<" ";
	}
	//cout<<endl;
	for(int i=0;i<(int)path1.size();i++){
		bool found=true;
		for(int j=0;j<(int)path2.size();j++){
			if(path1[i+j]!=path2[j]){
				found=false;
				break;
			}
		}
		if(found){
			//cout<<"relatives"<<endl;
			return true;
		}
	}
	//cout<<"not relatives"<<endl;
	return false;
}

void DeBruijnAssembler::showPath(vector<unsigned long int> path){
	Logger&m_cout=*(this->m_cout);
	for(int i=0;i<(int)path.size();i++){
		m_cout<<idToWord(path[i],m_wordSize)<<" ";
		if(i!=(int)path.size()-1)
			m_cout<<m_edge_states[path[i]][path[i+1]]<<" ";
	}
	m_cout<<endl;
}



void DeBruijnAssembler::join_vertices(unsigned long int prefix,unsigned long int suffix){
	
	Logger&m_cout=*(this->m_cout);
	//m_cout<<"Joining "<<prefix<<" -> "<<suffix<<" state: "<<m_edge_states[prefix][suffix]<<endl;

	if(prefix==suffix){
		m_edge_states[prefix][suffix]=m_COLOR_DISCARDED;
	}

	if(m_edge_states[prefix][suffix]==m_COLOR_REPEAT)
		return;


	if(prefix==suffix){
		m_edge_states[prefix][suffix]=m_COLOR_DISCARDED;
		return;
	}

	if(m_edge_states[prefix][suffix]==m_COLOR_DISCARDED)
		return;

	if(m_edge_states[prefix][suffix]==m_COLOR_IN_PROGRESS_SOLVER)
		m_edge_states[prefix][suffix]=m_COLOR_IN_PROGRESS_ASSEMBLER;

	if(m_edge_states[prefix][suffix]==m_COLOR_NOT_ASSEMBLED)
		m_edge_states[prefix][suffix]=m_COLOR_IN_PROGRESS_ASSEMBLER;

	// state here should be m_COLOR_IN_PROGRESS_ASSEMBLER or a contig

	//cout<<"Joining vertices "<<idToWord(prefix,m_wordSize)<<" -> "<<idToWord(suffix,m_wordSize)<<endl;
	vector<unsigned long int> beforePrefix=m_vertex_parents[prefix];
	vector<unsigned long int> afterPrefix;
	for(map<unsigned long int,vector<int> >::iterator i=m_graph[prefix].begin();i!=m_graph[prefix].end();i++){
		afterPrefix.push_back(i->first);
	}

	vector<unsigned long int> beforeSuffix=m_vertex_parents[suffix];
	vector<unsigned long int> afterSuffix;
	for(map<unsigned long int,vector<int> >::iterator i=m_graph[suffix].begin();i!=m_graph[suffix].end();i++){
		afterSuffix.push_back(i->first);
	}

	int before_prefix_contig_id=-1;
	unsigned long int beforeKMer=0;
	int after_suffix_contig_id=-1;
	unsigned long int afterKMer=0;
	for(int i=0;i<(int)beforePrefix.size();i++){
		if(m_edge_states[beforePrefix[i]][prefix]>=1){
			if(m_edge_states[beforePrefix[i]][prefix]==m_COLOR_REPEAT){
			}else if(m_edge_states[beforePrefix[i]][prefix]==m_COLOR_DISCARDED){
			}else if(before_prefix_contig_id!=-1){
				m_cout<<"Warning: REPEAT found "<<idToWord(beforePrefix[i],m_wordSize)<<" -> "<<idToWord(prefix,m_wordSize)<<endl;
				//trim_REPEAT(beforePrefix[i],prefix);
				//m_edge_states[beforePrefix[i]][prefix]=m_COLOR_REPEAT;
				//m_edge_states[beforeKMer][prefix]=m_COLOR_DISCARDED;
				//m_edge_states[prefix][suffix]=m_COLOR_DISCARDED;
				//m_edge_states[beforePrefix[i]][prefix]=m_COLOR_DISCARDED;
				before_prefix_contig_id=-1;
				break;
			}
			before_prefix_contig_id=m_edge_states[beforePrefix[i]][prefix];
			//cout<<"Edge "<<idToWord(beforePrefix[i],m_wordSize)<<" -> "<<idToWord(prefix,m_wordSize)<<" is in contig "<<before_prefix_contig_id<<endl;
			beforeKMer=beforePrefix[i];
		}
	}
	for(int i=0;i<(int)afterSuffix.size();i++){
		if(m_edge_states[suffix][afterSuffix[i]]>=1){
			if(m_edge_states[suffix][afterSuffix[i]]==m_COLOR_REPEAT){
			}else if(m_edge_states[suffix][afterSuffix[i]]==m_COLOR_DISCARDED){
			}else if(after_suffix_contig_id!=-1){
				after_suffix_contig_id=-1;
				//m_edge_states[suffix][afterSuffix[i]]=m_COLOR_REPEAT;
				//m_edge_states[prefix][suffix]=m_COLOR_DISCARDED;
				//m_edge_states[suffix][afterKMer]=m_COLOR_DISCARDED;
				//m_edge_states[suffix][afterSuffix[i]]=m_COLOR_DISCARDED;
				//trim_REPEAT(suffix,afterSuffix[i]);
				m_cout<<"Warning: REPEAT found "<<idToWord(suffix,m_wordSize)<<" -> "<<idToWord(afterSuffix[i],m_wordSize)<<endl;
				//cout<<suffix<<" "<<afterSuffix[i]<<" "<<m_edge_states[suffix][afterSuffix[i]]<<endl;
				//cout<<suffix<<" "<<afterKMer<<" "<<m_edge_states[suffix][afterKMer]<<endl;
				break;
			}
			after_suffix_contig_id=m_edge_states[suffix][afterSuffix[i]];
			//cout<<"Edge "<<idToWord(suffix,m_wordSize)<<" -> "<<idToWord(afterSuffix[i],m_wordSize)<<" is in contig "<<after_suffix_contig_id<<endl;
			afterKMer=afterSuffix[i];
		}
	}
	//m_cout<<"Around: "<<before_prefix_contig_id<<" ("<<beforePrefix.size()<<") "<<after_suffix_contig_id<<" ("<<afterSuffix.size()<<")"<<endl;
	
	// count parent, and show an Alert if 1 or more
	int i=0;
	for(vector<unsigned long int>::iterator j=m_vertex_parents[suffix].begin();j!=m_vertex_parents[suffix].end();j++){
		if(m_edge_states[*j][suffix]>=1){ // already have a parent
			i++;	
			//cout<<m_edge_states[*j][suffix]<<endl;
		}
	}
	if(i>0){
		m_edge_states[prefix][suffix]=m_COLOR_DISCARDED;
		cout<<"Alert: "<<i<<" parents"<<" "<<idToWord(prefix,m_wordSize)<<" -> "<<idToWord(suffix,m_wordSize)<<endl;
		return;
	}
	/*
                  before_prefix_contig_id      prefix       ?      suffix        after_suffix_contig_id
     (         )-------------------------  (           )---------(          )--------------------------------(            )

*/
	if(m_edge_states[prefix][suffix]==m_COLOR_IN_PROGRESS_ASSEMBLER){
		m_edge_states[prefix][suffix]=m_contig_id;
		m_contig_edges[m_contig_id][prefix]=suffix;
		m_contig_id++;
	}

	if(m_edge_states[prefix][suffix]<=0){// not in a contig
		return;
	}

	if(before_prefix_contig_id!=-1&&m_edge_states[prefix][suffix]!=before_prefix_contig_id){ // before
		//m_cout<<"Merging with contig before ("<<before_prefix_contig_id<<") "<<idToWord(beforeKMer,m_wordSize)<<endl;
		changeAllEdgeColor(m_edge_states[prefix][suffix],before_prefix_contig_id);
	}else{
		//cout<<"nothing to do before"<<endl;
	}

	if(after_suffix_contig_id!=-1&&m_edge_states[prefix][suffix]!=after_suffix_contig_id){ // after
		//m_cout<<"Merging with contig after ("<<after_suffix_contig_id<<") "<<idToWord(afterKMer,m_wordSize)<<endl;
		changeAllEdgeColor(m_edge_states[prefix][suffix],after_suffix_contig_id);
	}else{
		//m_cout<<"nothing to do after"<<endl;
	}
}

void DeBruijnAssembler::changeAllEdgeColor(int contigToReplace,int replacementContig){
	//cout<<"Merging: contig "<<contigToReplace<<" ("<<m_contig_edges[contigToReplace].size()<<") ->  contig "<<replacementContig<<" ("<<m_contig_edges[replacementContig].size()<<")"<<endl;
	for(map<unsigned long int,unsigned long int>::iterator i=m_contig_edges[contigToReplace].begin();
		i!=m_contig_edges[contigToReplace].end();i++){
			//cout<<"Adding "<<idToWord(i->first,m_wordSize)<<" -> "<<idToWord(i->second,m_wordSize)<<" in contig "<<replacementContig<<endl;
			m_edge_states[i->first][i->second]=replacementContig;
			m_contig_edges[replacementContig][i->first]=i->second;
	}
	m_contig_edges.erase(m_contig_edges.find(contigToReplace));
	//cout<<"Merging: contig "<<replacementContig<<" ("<<m_contig_edges[replacementContig].size()<<")"<<endl;
}

void DeBruijnAssembler::assemble_edge(unsigned long int prefix, unsigned long int suffix){

	Logger&m_cout=*(this->m_cout);
	if(m_edge_states[prefix][suffix]!=m_COLOR_NOT_ASSEMBLED&&m_edge_states[prefix][suffix]!=m_COLOR_IN_PROGRESS_SOLVER)
		return;
	m_edge_states[prefix][suffix]=m_COLOR_IN_PROGRESS_ASSEMBLER;
	list<unsigned long int> trivialPath;
	trivialPath.push_front(prefix);
	vector<unsigned long int> beforePrefix=m_vertex_parents[prefix];
	while(beforePrefix.size()==1){
		trivialPath.push_front(beforePrefix[0]);
		beforePrefix=m_vertex_parents[beforePrefix[0]];
	}
	map<unsigned long int,vector<int> > afterSuffix=m_graph[suffix];
	trivialPath.push_back(suffix);
	while(afterSuffix.size()==1){
		trivialPath.push_back(afterSuffix.begin()->first);
		afterSuffix=m_graph[afterSuffix.begin()->first];
	}

	vector<unsigned long int> path;
	for(list<unsigned long int>::iterator i=trivialPath.begin();i!=trivialPath.end();i++)
		path.push_back(*i);
	
	m_cout<<"[assemble_edge] "<<idToWord(prefix,m_wordSize)<<" -> "<<idToWord(suffix,m_wordSize)<<" :: "<<path.size()+m_wordSize<<endl;
	for(int i=0;i<(int)path.size()-1;i++){
		join_vertices(path[i],path[i+1]);
	}
	return ;
	//Logger&m_cout=*(this->m_cout);
	//m_cout<<endl;
	//m_cout<<"assemble_edge "<<idToWord(prefix,m_wordSize)<<" -> "<<idToWord(suffix,m_wordSize)<<endl;

#ifdef TOTO
	if(prefix==suffix){
		m_edge_states[prefix][suffix]=m_COLOR_REPEAT;
	}
	if(m_edge_states[prefix][suffix]==m_COLOR_REPEAT)
		return;

	if(m_edge_states[prefix][suffix]==m_COLOR_NOT_ASSEMBLED)
		m_edge_states[prefix][suffix]=m_COLOR_IN_PROGRESS_ASSEMBLER;

	vector<unsigned long int> beforePrefix=m_vertex_parents[prefix];
	vector<unsigned long int> afterPrefix;
	for(map<unsigned long int,vector<int> >::iterator i=m_graph[prefix].begin();i!=m_graph[prefix].end();i++){
		afterPrefix.push_back(i->first);
	}

	vector<unsigned long int> beforeSuffix=m_vertex_parents[suffix];
	vector<unsigned long int> afterSuffix;
	for(map<unsigned long int,vector<int> >::iterator i=m_graph[suffix].begin();i!=m_graph[suffix].end();i++){
		afterSuffix.push_back(i->first);
	}

	// discard if  near repeat.
	for(vector<unsigned long int>::iterator i=beforePrefix.begin();i!=beforePrefix.end();i++){
		if(m_edge_states[*i][prefix]==m_COLOR_REPEAT)
			m_edge_states[prefix][suffix]=m_COLOR_DISCARDED;
	}
	for(vector<unsigned long int>::iterator i=afterPrefix.begin();i!=afterPrefix.end();i++){
		if(m_edge_states[prefix][*i]==m_COLOR_REPEAT)
			m_edge_states[prefix][suffix]=m_COLOR_DISCARDED;
	}
	for(vector<unsigned long int>::iterator i=beforeSuffix.begin();i!=beforeSuffix.end();i++){
		if(m_edge_states[*i][suffix]==m_COLOR_REPEAT)
			m_edge_states[prefix][suffix]=m_COLOR_DISCARDED;
	}
	for(vector<unsigned long int>::iterator i=afterSuffix.begin();i!=afterSuffix.end();i++){
		if(m_edge_states[suffix][*i]==m_COLOR_REPEAT)
			m_edge_states[prefix][suffix]=m_COLOR_DISCARDED;
	}

	if(m_edge_states[prefix][suffix]==m_COLOR_DISCARDED)
		return;



	vector<vector<unsigned long int> > pathsToTry;

	// (prefix)-(suffix)
	vector<unsigned long int> path1;
	path1.push_back(prefix);
	path1.push_back(suffix);
	pathsToTry.push_back(path1);

	// (before prefix)-(prefix)-(after prefix)
	for(int i=0;i<(int)beforePrefix.size();i++){
		for(int j=0;j<(int)afterPrefix.size();j++){
			vector<unsigned long int> path2;
			path2.push_back(beforePrefix[i]);
			path2.push_back(prefix);
			path2.push_back(afterPrefix[j]);
			pathsToTry.push_back(path2);
		}
	}

	// (before suffix)-(suffix)-(after suffix)
	for(int i=0;i<(int)beforeSuffix.size();i++){
		for(int j=0;j<(int)afterSuffix.size();j++){
			vector<unsigned long int> path2;
			path2.push_back(beforeSuffix[i]);
			path2.push_back(suffix);
			path2.push_back(afterSuffix[j]);
			pathsToTry.push_back(path2);
		}
	}

	// (before prefix)-(prefix)-(suffix)-(after suffix)
	for(int i=0;i<(int)beforePrefix.size();i++){
		for(int j=0;j<(int)afterSuffix.size();j++){
			vector<unsigned long int> path2;
			path2.push_back(beforePrefix[i]);
			path2.push_back(prefix);
			path2.push_back(suffix);
			path2.push_back(afterSuffix[j]);
			pathsToTry.push_back(path2);
		}
	}
	
	int bestPath=-1;
	//m_cout<<endl;
	for(int i=0;i<(int)pathsToTry.size();i++){
		//cout<<pathToDNA(pathsToTry[i])<<endl;
		//showPath(pathsToTry[i]);
		//cout<<number_of_reads_that_Agree_on_path(pathsToTry[i])<<endl;
		if(number_of_reads_that_Agree_on_path(pathsToTry[i])>0){
			if(bestPath==-1){
				bestPath=i;
			}else if(bestPath!=-1||areRelatives(pathsToTry[bestPath],pathsToTry[i])){
				if(pathsToTry[bestPath].size()<pathsToTry[i].size())
					bestPath=i;
			}else{
				bestPath=-1;
				break;
			}
		}
	}
	// invariant: only one path 
	if(bestPath!=-1){
		//m_cout<<"Winner "<<endl;
		//showPath(pathsToTry[bestPath]);
		for(int i=0;i<(int)pathsToTry[bestPath].size()-1;i++){
			join_vertices(pathsToTry[bestPath][i],pathsToTry[bestPath][i+1]);
		}
		//m_cout<<"after "<<endl;
		//showPath(pathsToTry[bestPath]);
		//m_cout<<(m_contig_edges[m_edge_states[prefix][suffix]].size())<<" edges "<<endl;
	}else{
		m_cout<<"DIscarded"<<endl;
		// TODO: search more to choose between two things
		m_edge_states[prefix][suffix]=m_COLOR_DISCARDED;
	}
#endif
}

void DeBruijnAssembler::run_Assembler(){
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
	m_cout<<"fixing graph"<<endl;

	//solve strange structures
	int progress=0;
	int maximumDepth=1000;
	for(set<unsigned long int>::iterator i=m_graph_mers.begin();i!=m_graph_mers.end();i++){

		if(progress%10000==0){
			m_cout<<progress<<" / "<<total_Edges<<endl;
		}
		solveMultiPath(*i,maximumDepth);
		progress++;
	}

	m_cout<<progress<<" / "<<progress<<endl;
	m_cout<<"Assembling contigs"<<endl;

	// remove when 2 parents.
	for(set<unsigned long int>::iterator i=m_graph_mers.begin();i!=m_graph_mers.end();i++){
		if(m_vertex_parents[*i].size()>1){
			for(vector<unsigned long int>::iterator j=m_vertex_parents[*i].begin();j!=m_vertex_parents[*i].end();j++){
				if(m_edge_states[*j][*i]==m_COLOR_NOT_ASSEMBLED){
					m_cout<<"Tagging REPEAT.. "<<idToWord(*j,m_wordSize)<<" -> "<<idToWord(*i,m_wordSize)<<endl;
					m_edge_states[*j][*i]=m_COLOR_REPEAT;
				}
			}
		}
	}


	// run assembler on edges
	progress=0;
	for(map<unsigned long int,map<unsigned long int,vector<int> > >::iterator i=m_graph.begin();i!=m_graph.end();i++){
		// sanity check, should be <= 4 (A T C G)
		if(i->second.size()>4){
			m_cout<<endl;
			m_cout<<"1 Error: more than 4, "<<i->second.size()<<endl;
			m_cout<<idToWord(i->first,m_wordSize)<<endl;
			for(map<unsigned long int,vector<int> >::iterator j=i->second.begin();j!=i->second.end();j++){
				m_cout<<" "<<idToWord(j->first,m_wordSize)<<" "<<j->second.size()<<endl;
			}
		}

		for(map<unsigned long int,vector<int> >::iterator j=i->second.begin();j!=i->second.end();j++){
			if(progress%10000==0){
				m_cout<<progress<<" / "<<total_Edges<<endl;
			}
			assemble_edge(i->first,j->first);
			if(m_edge_states[i->first][j->first]==m_COLOR_NOT_ASSEMBLED||
				m_edge_states[i->first][j->first]==m_COLOR_IN_PROGRESS_ASSEMBLER){
				m_cout<<"ERROR: "<<m_edge_states[i->first][j->first]<<endl;
			}
		
			progress++;
		}
	}
	m_cout<<progress<<" / "<<total_Edges<<endl;

	// sanity check, should be ok always
	for(set<unsigned long int>::iterator i=m_graph_mers.begin();i!=m_graph_mers.end();i++){
		map<int,int> parentsForContig;
		for(vector<unsigned long int>::iterator j=m_vertex_parents[*i].begin();j!=m_vertex_parents[*i].end();j++){
			if(m_edge_states[*j][*i]>=1)
				parentsForContig[m_edge_states[*j][*i]]++;
		}
		for(map<int,int>::iterator j=parentsForContig.begin();j!=parentsForContig.end();j++){
			if(j->second>1){
				cout<<"Error: "<<idToWord(*i,m_wordSize)<<" have more than 1 parent in contig, (REPEAT region)"<<j->first<<endl;
			}
		}
	}

	// write edges
	string edgesFile=m_assemblyDirectory+"/Edges.txt";
	m_cout<<"Writing edges file: "<<edgesFile<<endl;
	ofstream edgesStream(edgesFile.c_str());
	for(map<unsigned long int,map<unsigned long int,vector<int> > >::iterator i=m_graph.begin();i!=m_graph.end();i++){
		for(map<unsigned long int,vector<int> >::iterator j=i->second.begin();j!=i->second.end();j++){
			edgesStream<<idToWord(i->first,m_wordSize)<<" -> "<<idToWord(j->first,m_wordSize)<<" "<<m_edge_states[i->first][j->first]<<endl;
			//edgesStream<<(i->first)<<" -> "<<(j->first)<<" "<<m_edge_states[i->first][j->first]<<endl;
		}
	}
	edgesStream.close();

	// make sure that all contigs have at least a startddj
	m_cout<<endl;
	m_cout<<"Walking on contigs ("<<m_contig_edges.size()<<")"<<endl;
	map<int,set<unsigned long int> > contigVertices;
	for(map<int,map<unsigned long int,unsigned long int> >::iterator i=m_contig_edges.begin();
		i!=m_contig_edges.end();i++){
		for(map<unsigned long int, unsigned long int>::iterator j=i->second.begin();
			j!=i->second.end();j++){
			contigVertices[i->first].insert(j->first);
			contigVertices[i->first].insert(j->second);
		}
	}
	for(map<int,set<unsigned long int> >::iterator contigIt=contigVertices.begin();
		contigIt!=contigVertices.end();contigIt++){
		int contigIndex=contigIt->first;
		unsigned long int startPoint=0;
		m_cout<<contigIndex<<endl;
		bool foundStart=true;
		int maxDistance=0;
		map<unsigned long int,int> cache;
		// iterate and find a vertex without parent, it will be the starting point of the contig
		for(set<unsigned long int>::iterator i=contigIt->second.begin();i!=contigIt->second.end();i++){
			//cout<<*i<<endl;
			unsigned long int toTest=*i;
			int  distance=0;
			bool expanding=true;
			while(expanding==true){
				expanding=false;
				if(cache.count(toTest)>0){
					distance+=cache[toTest];
					break;
				}
				for(map<unsigned long int,vector<int> >::iterator j=m_graph[toTest].begin();j!=m_graph[toTest].end();j++){
					if(m_edge_states[toTest][j->first]==contigIndex){
						distance++;
						//cout<<"Expanding "<<distance<<endl;
						toTest=j->first;
						expanding=true;
						break;
					}
				}
			}
			//cout<<"d= "<<distance<<endl;
			cache[*i]=distance;
			if(distance>maxDistance){
				maxDistance=distance;
				startPoint=*i;
			}
		}
		m_cout<<contigVertices[contigIndex].size()+m_wordSize<<" nucleotides"<<endl;
		
		if(foundStart==false){
			m_cout<<"contig has no start..."<<endl;
		}else{
			m_cout<<"Making the contig path de novo, start= "<<idToWord(startPoint,m_wordSize)<<endl;
			// sanity check
			if(m_vertex_parents[startPoint].size()>0){
				//m_cout<<"start has parent: state is "<<m_edge_states[m_vertex_parents[startPoint][0]][startPoint]<<endl;
			}
			makeContigPath(startPoint,contigIndex);
		}
	}
	
	
	//m_cout<<m_graph_mers.size()<<" / "<<m_graph_mers.size()<<endl;
	m_solver_log->close();
}

void DeBruijnAssembler::makeContigPath(unsigned long int prefix,int contig){
	m_contig_paths[contig].push_back(prefix);
	bool becomingLarger=true;
	map<unsigned long int,set<unsigned long int> > doneSoFar;
	while(becomingLarger){
		becomingLarger=false;
		for(map<unsigned long int,vector<int> >::iterator j=m_graph[prefix].begin();j!=m_graph[prefix].end();j++){
			unsigned long int suffix=j->first;
			if(doneSoFar[prefix].count(suffix)>0){
				cout<<"Warning(2): REPEAT found, "<<idToWord(prefix,m_wordSize)<<" -> "<<idToWord(suffix,m_wordSize)<<endl;
			}
			if(m_edge_states[prefix][suffix]==contig){
				doneSoFar[prefix].insert(suffix);
				m_contig_paths[contig].push_back(suffix);
				//cout<<idToWord(suffix,m_wordSize)<<endl;
				prefix=suffix;
				becomingLarger=true;
				break;
			}
		}
		//cout<<i<<endl;
		//i++;
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
	for(map<int,vector<unsigned long int> >::iterator iteratorContigFirstVertex=m_contig_paths.begin();
		iteratorContigFirstVertex!=m_contig_paths.end();iteratorContigFirstVertex++){
		m_cout<<"Contig #"<<iteratorContigFirstVertex->first<<": "<<iteratorContigFirstVertex->second.size()+m_wordSize<<" nucleotides"<<endl;
		if((int)iteratorContigFirstVertex->second.size()+m_wordSize<m_minimumContigSize)
			continue;
		string contigSequence=pathToDNA(iteratorContigFirstVertex->second);
		contigs.push_back(contigSequence);
		contigNumbers.push_back(iteratorContigFirstVertex->first);
	}


	// remove duplicates
	set<int> duplicatedContig;

	for(int contigIndex=0;contigIndex<(int)contigs.size();contigIndex++){
		//cout<<contigIndex<<endl;
		if(duplicatedContig.count(contigIndex)>0)// if contig is not already tagged as duplicated
			continue;
		for(int otherContigIndex=0;otherContigIndex<(int)contigs.size();otherContigIndex++){
			//cout<<"  "<<otherContigIndex<<endl;
			if(otherContigIndex==contigIndex)// is the same
				continue;
			if(duplicatedContig.count(otherContigIndex)>0)// already tagged as duplicated
				continue;
			if(contigs[contigIndex].length()!=contigs[otherContigIndex].length()){
				int diff=contigs[contigIndex].length()-contigs[otherContigIndex].length();
				if(diff<0)
					diff=-diff;
				int smallContig=contigIndex;
				int otherContig=otherContigIndex;
				if(contigs[otherContigIndex].length()<contigs[smallContig].length()){
					smallContig=otherContigIndex;
					otherContig=contigIndex;
				}
				if(diff==1){
					string revComp=reverseComplement(contigs[otherContig]).substr(1,contigs[smallContig].length());
					if(revComp==contigs[contigIndex]){
						//m_cout<<"duplicate"<<endl;
						duplicatedContig.insert(contigIndex);
						m_cout<<"Duplicate -+1, contigs "<<contigNumbers[otherContig]<<" and "<<contigNumbers[smallContig]<<endl;
					}
					revComp=reverseComplement(contigs[otherContig].substr(1,contigs[smallContig].length()));
					if(revComp==contigs[contigIndex]){
						//m_cout<<"duplicate"<<endl;
						duplicatedContig.insert(contigIndex);
						m_cout<<"Duplicate -+1, contigs "<<contigNumbers[otherContig]<<" and "<<contigNumbers[smallContig]<<endl;
					}



				}else{
					//m_cout<<revComp.substr(0,70)<<endl;
					//m_cout<<contigs[otherContigIndex].substr(0,70)<<endl;
				}
			}
			if(contigs[contigIndex].length()==contigs[otherContigIndex].length()){
				string revComp=reverseComplement(contigs[contigIndex]);
				if(revComp==contigs[otherContigIndex]){
					m_cout<<"Duplicate +-0, contigs "<<contigNumbers[otherContigIndex]<<" and "<<contigNumbers[contigIndex]<<endl;
					duplicatedContig.insert(otherContigIndex);
				}else{
					//m_cout<<revComp.substr(0,70)<<endl;
					//m_cout<<contigs[otherContigIndex].substr(0,70)<<endl;
				}
			}
		}
	}

	if(contigs.size()==0){ 
		m_cout<<"Error: 0 contigs"<<endl;
		return;
	}
	m_cout<<endl;
	m_cout<<endl;
	m_cout<<contigs.size()-duplicatedContig.size()<<" contigs ("<<contigs.size()<<" in the graph)"<<endl;	
	//m_cout<<contigs.size()-duplicatedContig.size()<<" after mergin of complement contigs (the graph contains both strands)"<<endl;
	int totalLength=0;
	int contigId=1;
	multiset<int> contigSizes;
	unsigned long int stat_minimumContigSize=9999999999;
	int stat_maximumContigSize=0;
	for(int i=0;i<(int)contigs.size();i++){
		if(duplicatedContig.count(i)>0){
			continue;
		}

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
		output<<">Contig";
		output<<contigId;
		output<<endl;
		contigId++;
		while(j<(int)sequenceDNA.length()){
			output<<sequenceDNA.substr(j,columns);
			output<<endl;
			j+=columns;
		}

	}
	output.close();
	m_cout<<endl;
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
				m_cout<<"N"<<j<<" size: "<<*i<<endl;
				break;
			}
		}
	}
	m_cout<<"See http://www.cbcb.umd.edu/research/castats.shtml#N50 for the definition of N50."<<endl;
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
	// majority vote
	neighboursVector=getNeighbours(vertex);
	if(neighboursVector.size()>1){
		// TODO: put debug here only.
		// new algo here
		(*m_solver_log)<<"new algo"<<endl;
		map<int,map<int,int> > votes;
		int maximumThreshold=20;
		for(int i=0;i<(int)neighboursVector.size();i++){
			for(int threshold=0;threshold<=maximumThreshold;threshold++)
				votes[threshold][i]=0;
			// TODO: do majority vote in binary space here.
			// this will give the same results faster
			string path=getTrivialPath(neighboursVector[i],vertex);
			(*m_solver_log)<<"Path "<<path<<endl;
			// count the reads that agree..
			//map<unsigned long int,int> treeForPath;
			// build hash
			map<int,int> readVotes;
			for(int iteratorPath=0;iteratorPath<(int)path.length()-m_wordSize;iteratorPath++){
				unsigned long int kMerInPath=wordId(path.substr(iteratorPath,m_wordSize+1).c_str());
				//treeForPath[kMerInPath]=iteratorPath;
				//(*m_solver_log)<<m_mer_to_read_table[kMerInPath].size()<<" reads for "<<idToWord(kMerInPath,m_wordSize+1)<<endl;
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
				(*m_solver_log)<<"Found a winner ;P, threshold="<<threshold<<endl;
				break; // found a winner
			}
		}
		if(winner!=-1){
			(*m_solver_log)<<winner<<" wins, has "<<maxVotes<<" other 0."<<endl;
			for(int iteratorOnNext=0;iteratorOnNext<(int)neighboursVector.size();iteratorOnNext++){
				if(iteratorOnNext!=winner){
					m_edge_states[vertex][neighboursVector[iteratorOnNext]]=m_COLOR_DISCARDED;
				}
			}
		}else{
			// TODO: add this as a parameter?
			// bubbles are caused by 454's homopolymers that are cursed by wizards
			int maximumDifferenceInBubbles=3;
			(*m_solver_log)<<"No winner..."<<endl;
			bool allTheSameLength=true;
			for(int u=0;u<(int)neighboursVector.size();u++){
				if(getTrivialPath(neighboursVector[u],vertex).length()!=getTrivialPath(neighboursVector[0],vertex).length()){
					allTheSameLength=false;
				}
			}
			if((int)(allTheSameLength&&getTrivialPath(neighboursVector[0],vertex).length())==m_wordSize+1){
				(*m_solver_log)<<endl;
				(*m_solver_log)<<"All too short, skipping."<<endl;
				for(int u=0;u<(int)neighboursVector.size();u++){
					(*m_solver_log)<<getTrivialPath(neighboursVector[u],vertex)<<endl;
					//m_edge_states[vertex][neighboursVector[u]]=m_COLOR_DISCARDED;
				}
			}else if(neighboursVector.size()==2){
				(*m_solver_log)<<"2 alternative paths"<<endl;
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
						(*m_solver_log)<<endl;
						// TODO: maybe chop it all
						(*m_solver_log)<<"Bubble found!, removing it... "<<path1.length()<<" "<<path2.length()<<endl;
						// TODO SNP detection goes here ?? or it can be sequencing errors also..
						m_edge_states[vertex][neighboursVector[1]]=m_COLOR_DISCARDED;
						(*m_solver_log)<<path1<<endl;
						(*m_solver_log)<<path2<<endl;
					}else{
						(*m_solver_log)<<endl;
						(*m_solver_log)<<"Non-trivial ("<<diff<<") Bubble found!, skipping... "<<path1.length()<<" "<<path2.length()<<endl;
						//m_edge_states[vertex][neighboursVector[0]]=m_COLOR_DISCARDED;
						//m_edge_states[vertex][neighboursVector[1]]=m_COLOR_DISCARDED;
						(*m_solver_log)<<path1<<endl;
						(*m_solver_log)<<path2<<endl;
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
						(*m_solver_log)<<endl;
						(*m_solver_log)<<mismatches<<" mismatches in first "<<minimumPathSize<< " bases (2 things, not a bubble): taking the longuest"<<endl;
						(*m_solver_log)<<path1<<endl;
						(*m_solver_log)<<path2<<endl;
						if((int)minimumPathSize==(int)path1.length()){
							//m_colors[neighboursVector[0]]=m_COLOR_DISCARDED;
							m_edge_states[vertex][neighboursVector[0]]=m_COLOR_DISCARDED;
						}else{
							//m_colors[neighboursVector[1]]=m_COLOR_DISCARDED;
							m_edge_states[vertex][neighboursVector[1]]=m_COLOR_DISCARDED;

						}
					}else{
						(*m_solver_log)<<endl;
						(*m_solver_log)<<"Working on it (2 things, not a bubble), "<<mismatches<<" mismatches in the first "<<minimumPathSize<<" bases"<<endl;
						if(path1IsRepeat){
							(*m_solver_log)<<" Path 1 is a REPEAT ("<<getTrivialPath(getNeighbours(repeatPath1)[0],repeatPath1)<<")"<<endl;
							//m_edge_states[vertex][repeatPath1]=m_COLOR_REPEAT;
						}
						if(path2IsRepeat){
							(*m_solver_log)<<" Path 2 is a REPEAT ("<<getTrivialPath(getNeighbours(repeatPath2)[0],repeatPath2)<<")"<<endl;
							//m_edge_states[vertex][repeatPath2]=m_COLOR_REPEAT;
						}
						(*m_solver_log)<<path1<<endl;
						(*m_solver_log)<<path2<<endl;
					}
				}
			}else if(neighboursVector.size()==3){
				// TODO here
				(*m_solver_log)<<endl;
				(*m_solver_log)<<"3 alternative paths"<<endl;
				for(int u=0;u<(int)neighboursVector.size();u++){
					(*m_solver_log)<<getTrivialPath(neighboursVector[u],vertex)<<endl;
					//m_edge_states[vertex][neighboursVector[u]]=m_COLOR_DISCARDED;
				}
			}else if(neighboursVector.size()==3){
				// TODO here
				(*m_solver_log)<<endl;
				(*m_solver_log)<<"4 alternative paths"<<endl;
				for(int u=0;u<(int)neighboursVector.size();u++){
					(*m_solver_log)<<getTrivialPath(neighboursVector[u],vertex)<<endl;
					//m_edge_states[vertex][neighboursVector[u]]=m_COLOR_DISCARDED;
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
}



// THE END IS HERE  ------------> .
