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

#include<sstream>
#include<iostream>
#include<string>
#include<vector>
#include"Loader.h"
#include"SffLoader.h"
#include<stdlib.h>
#include"Read.h"
using namespace std;

Loader::Loader(){
	m_bases=0;
}

void Loader::load(string file,vector<Read*>*reads){
	(cout)<<"Loading "<<file<<endl;
	if(file.length()<4){
		(cout)<<"Error: "<<file<<endl;
		exit(0);
	}
	if(file.substr(file.length()-4,4)==".sff"){
		(cout)<<"Format: SFF"<<endl;
		SffLoader sffLoader;
		int num=reads->size();
		sffLoader.load(file,reads);
		int num2=reads->size()-num;
		m_bases=sffLoader.getBases();
		(cout)<<num2<<" reads"<<endl;
		(cout)<<m_bases<<" bases"<<endl;
		return;
	}

	ifstream f(file.c_str());
	int fasta=0;
	int fastq=1;
	int type=fasta;
	m_total=0;
	string buffer;
	f>>buffer;
	if(buffer[0]=='>'){
		type=fasta;
		(cout)<<"Format: fasta."<<endl;
	}else if(buffer[0]=='@'){
		type=fastq;
		(cout)<<"Format: fastq."<<endl;
	}else{
		(cout)<<"Format: unknown."<<endl;
	}
	//TODO: put this in FastaLoader.cpp
	f.seekg(0,ios_base::beg);
	if(type==fasta){
		string id;
		ostringstream sequence;
		string buffer;
		while(!f.eof()){
			buffer="";
			f>>buffer;
			if(buffer=="")
				continue;
			if(buffer[0]=='>'){
				char bufferForLine[1024];
				f.getline(bufferForLine,1024);
				if(id!=""){
					string sequenceStr=sequence.str();
					ostringstream quality;
					for(int i=0;i<(int)sequenceStr.length();i++){
						quality<< "F";
					}
					add(reads,&id,&sequence,&quality);
				}
				id=buffer;
				sequence.str("");
			}else{
				sequence<< buffer;
			}
		}
		string sequenceStr=sequence.str();
		ostringstream quality;
		for(int i=0;i<(int)sequenceStr.length();i++){
			quality<< "F";
		}
		add(reads,&id,&sequence,&quality);
		
	//TODO: put this in FastqLoader.cpp
	}else if(type==fastq){
		string id;
		ostringstream sequence;
		ostringstream quality;
		int seq=0;
		int qual=1;
		int mode=seq;
		string buffer;
		while(!f.eof()){
			buffer="";
			f>>buffer;
			if(buffer=="")
				continue;
			if(buffer[0]=='@'&&!(mode==qual&&quality.str().length()==0)){
				char bufferForLine[1024];
				f.getline(bufferForLine,1024);
				if(id!=""){
					add(reads,&id,&sequence,&quality);
				}
				id=buffer;
				sequence.str("");
				quality.str("");
				mode=seq;
			}else if(buffer[0]=='+'&&!(mode==qual&&quality.str().length()==0)){
				char bufferForLine[1024];
				f.getline(bufferForLine,1024);
				mode=qual;
			}else if(mode==qual){
				quality<< buffer;
			}else if(mode==seq){
				sequence<< buffer;
			}
		}
		add(reads,&id,&sequence,&quality);
	}
	f.close();
	(cout)<<"Reads: "<<m_total<<endl;
}

int Loader::getBases(){
	return m_bases;
}

void Loader::add(vector<Read*>*reads,string*id,ostringstream*sequence,ostringstream*quality){
	if(id->length()==0)
		return;
	if(sequence->str().length()==0){
		//(cout)<<"0 22"<<endl;
		(cout)<<"Empty sequence? "<<*id<<endl;
		return;
	}
	if(sequence->str().length()!=quality->str().length()){
		(cout)<<*id<<endl;
		(cout)<<sequence->str()<<endl;
		(cout)<<quality->str()<<endl;

		(cout)<<"ERROR length"<<endl;
		exit(0);
	}
	string sequenceStr=sequence->str();
	string theId=*id;
	if(theId[0]=='>')
		theId=theId.substr(1);

	Read*read=new Read(theId.c_str(),sequenceStr.c_str());
	m_bases+=sequenceStr.length();
	m_total++;
	reads->push_back(read);
}
