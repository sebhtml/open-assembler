/*
	dna: De Novo Assembler
    Copyright (C) 2008, 2009 Sébastien Boisvert
	$Id: SequenceData.cpp 11 2009-01-23 00:57:32Z boiseb01 $

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

#include"SequenceDataFull.h"
#include"Loader.h"
#include<fstream>

SequenceDataFull::SequenceDataFull(vector<string>*files,ostream*logger){
	m_inputFiles=*files;
	m_cout=logger;
	m_reads=0;
	m_bases=0;
	for(int i=0;i<(int)m_inputFiles.size();i++){
		(*m_cout)<<endl;
		(*m_cout)<<"Files: "<<i+1<<" / "<<m_inputFiles.size()<<endl;
		Loader loader(m_cout);
		m_file_first[m_inputFiles[i]]=m_reads;
		loader.load(m_inputFiles[i],&m_reads_vector);
		m_bases+=loader.getBases();
		//(*m_cout)<<reads.size()<<endl;
		m_reads=m_reads_vector.size();
		m_file_last[m_inputFiles[i]]=m_reads-1;
	}

}

int SequenceDataFull::size(){
	return m_reads;
}


Read*SequenceDataFull::at(int i){
	return (m_reads_vector[i]);
}





uint64_t SequenceDataFull::bases(){
	return m_bases;
}



int SequenceDataFull::getFirst(string a){
	return m_file_first[a];
}

int SequenceDataFull::getLast(string a){
	return m_file_last[a];
}

bool SequenceDataFull::hasFile(string a){
	return m_file_last.count(a)>0;
}
