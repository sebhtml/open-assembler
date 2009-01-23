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

#include"SequenceData.h"
#include"Loader.h"
#include<fstream>

SequenceData::SequenceData(vector<string>*files,ostream*logger){
	m_inputFiles=*files;
	m_cout=logger;
	m_reads=0;
	m_bases=0;
	for(int i=0;i<(int)m_inputFiles.size();i++){
		(*m_cout)<<"Files: "<<i+1<<" / "<<m_inputFiles.size()<<endl;
		vector<Read*>reads;
		Loader loader(m_cout);
		loader.load(m_inputFiles[i],&reads);
		m_bases+=loader.getBases();
		m_file_first[m_inputFiles[i]]=m_reads;
		//(*m_cout)<<reads.size()<<endl;
		m_reads+=reads.size();
		m_file_last[m_inputFiles[i]]=m_reads-1;
		for(int j=0;j<(int)reads.size();j++)
			delete reads[j];
		(*m_cout)<<endl;
	}

	m_file_1=0;
	m_file_2=1;
	loadFiles();
}

int SequenceData::size(){
	return m_reads;
}


Read*SequenceData::at(int i){
	if(m_start_1<=i&&i<=m_end_1)
		return m_reads_1[i-m_start_1];
	if(m_start_2<=i&&i<=m_end_2)
		return m_reads_2[i-m_start_2];
	
	m_file_1+=2;
	m_file_2+=2;
	loadFiles();
	return at(i);
}





uint64_t SequenceData::bases(){
	return m_bases;
}

void SequenceData::loadFiles(){
	//(*m_cout)<<"LOAD"<<endl;
	(*m_cout)<<endl;
	(*m_cout)<<"Loading "<<m_file_1+1<<" / "<<m_inputFiles.size()<<endl;
	if(m_reads_1.size()>0)
		for(int i=0;i<(int)m_reads_1.size();i++){
			delete m_reads_1[i];
			m_reads_1[i]=NULL;
		}
	if(m_reads_2.size()>0)
		for(int i=0;i<(int)m_reads_2.size();i++){
			delete m_reads_2[i];
			m_reads_2[i]=NULL;
		}
	m_reads_1.clear();
	m_reads_2.clear();
	if(m_file_1>=(int)m_inputFiles.size()){
		//(*m_cout)<<"Error file 2 not set"<<m_file_1<<endl;
		//exit(0);
		m_file_1=0;
		m_file_2=1;
		loadFiles();
		return;
	}
	string file1=m_inputFiles[m_file_1];

	m_start_1=m_file_first[file1];
	m_end_1=m_file_last[file1];
	Loader loader1(m_cout);
	(*m_cout)<<file1<<endl;
	loader1.load(file1,&m_reads_1);
	if(m_file_2>=(int)m_inputFiles.size())
		return;
	
	(*m_cout)<<endl;
	(*m_cout)<<"Loading "<<m_file_2+1<<" / "<<m_inputFiles.size()<<endl;
	Loader loader2(m_cout);
	string file2=m_inputFiles[m_file_2];
	loader2.load(file2,&m_reads_2);
	m_start_2=m_file_first[file2];
	m_end_2=m_file_last[file2];
}

int SequenceData::getFirst(string a){
	return m_file_first[a];
}

int SequenceData::getLast(string a){
	return m_file_last[a];
}

bool SequenceData::hasFile(string a){
	return m_file_last.count(a)>0;
}
