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
#include<iostream>
#include"Read.h"
#include"DeBruijnAssembler.h"
#include<vector>
#include"Logger.h"
#include<fstream>

#define endl '\n'

using namespace std;


int main(int argc,char*argv[]){
	Logger m_cout;
	m_cout<<"Command line:"<<endl;
	m_cout<<endl;
	for(int i=0;i<argc;i++){
		m_cout<<argv[i]<<" ";
	}
	m_cout<<endl;
	m_cout<<endl;
	// check processor
	m_cout<<(int)sizeof(unsigned long int)*8<<" bits available"<<endl;
	if(sizeof(unsigned long int)!=8){
		m_cout<<"Abording, you need a 64-bit machine."<<endl;
		return 0;
	}

	m_cout<<endl;


	// show usage
	m_cout<<"Welcome to 454dna, the 454 De Novo Assembler."<<endl;
	m_cout<<"http://DeNovoAssembler.sf.net"<<endl;
	m_cout<<"License: GPLv3"<<endl;
	m_cout<<endl;
	m_cout<<"The input file format is the fastq format."<<endl;
	m_cout<<endl;
	m_cout<<"To convert your sff files, use flower."<<endl;
	m_cout<<"See http://blog.malde.org/index.php/2008/11/14/454-sequencing-and-parsing-the-sff-binary-format/"<<endl;
	m_cout<<endl;
	m_cout<<"Usage:"<<endl;
	m_cout<<"454dna [options] <fastq files>"<<endl;



	// show options
	m_cout<<" OPTIONS"<<endl;
	vector<string> inputFiles;

	string assemblyDirectory="454Assembly";
	m_cout<<" -assemblyDirectory   default: "<<assemblyDirectory<<endl;
	int minimumCoverage=5;
	m_cout<<" -minimumCoverage     default: "<<minimumCoverage<<endl;
	int wordSize=25;
	m_cout<<" -wordSize            default: "<<wordSize<<" (maximum: 31)"<<endl;
	int minimumQuality=40;
	m_cout<<" -minimumQuality      default: "<<minimumQuality<<endl;
	int minimumContigSize=500;
	m_cout<<" -minimumContigSize   default: "<<minimumContigSize<<endl;

	// collect arguments
	for(int i=1;i<argc;i++){
		string option=argv[i];
		if(option=="-assemblyDirectory"){
			i++;
			assemblyDirectory=argv[i];
		}else if(option=="-minimumCoverage"){
			i++;
			minimumCoverage=atoi(argv[i]);
		}else if(option=="-minimumContigSize"){
			i++;
			minimumContigSize=atoi(argv[i]);
		}else if(option=="-wordSize"){
			i++;
			if(wordSize>31){
				m_cout<<"Wordsize cannot be greater than 31"<<endl;
			}else{
				wordSize=atoi(argv[i]);
			}
		}else if(option=="-minimumQuality"){
			i++;
			minimumQuality=atoi(argv[i]);
		}else{
			inputFiles.push_back(argv[i]);
		}
	}
	string command=" mkdir -p "+assemblyDirectory+" # [454dna] ";
	system(command.c_str());
	m_cout.setOutput(assemblyDirectory+"/Log.txt");
	
	m_cout<<endl;
	m_cout<<"Writing log to "<<assemblyDirectory+"/Log.txt"<<endl;

	m_cout<<endl;

	// check files
	if(inputFiles.size()==0){
		m_cout<<"You must provide fastq files"<<endl;
		return 0;
	}

	// show user options
	m_cout<<"[454dna] assemblyDirectory = "<<assemblyDirectory<<endl;
	m_cout<<" minimumCoverage = "<<minimumCoverage<<endl;
	m_cout<<" wordSize = "<<wordSize<<endl;
	m_cout<<" minimumQuality = "<<minimumQuality<<endl;
	m_cout<<" minimumContigSize = "<<minimumContigSize<<endl;
	m_cout<<" inputFiles = "<<endl;
	for(int i=0;i<(int)inputFiles.size();i++){
		m_cout<<"        "<<inputFiles[i]<<endl;
	}
	m_cout<<endl;


	// loading reads
	set<string> readNames;
	vector<Read*> reads;
	m_cout<<"Loading reads"<<endl;
	int _minimumQuality=999999;
	int _maximumQuality=-999999;
	unsigned long int bases=0;
	for(int i=0;i<(int)inputFiles.size();i++){
		unsigned long int file_quality=0;
		m_cout<<endl;
		m_cout<<"  Loading "<<inputFiles[i]<<endl;
		int file_reads=0;
		ifstream f(inputFiles[i].c_str());
		if(!f){
			m_cout<<"Error: invalid file "<<inputFiles[i]<<endl;
			return 0;
		}
		unsigned long int file_bases=0;
		while(!f.eof()){
			string id;
			string id2;
			string sequence;
			string quality;
			f>>id>>sequence>>id2>>quality;
			if(id=="")
				break;

			if(readNames.count(id)!=0){
				m_cout<<"Warning: read "<<id<<" is duplicated."<<endl;
				continue;
			}

			for(int u=0;u<(int)quality.length();u++){
				int qualityValue=(int)quality[u];
				file_quality+=qualityValue;
				if(qualityValue<_minimumQuality)
					_minimumQuality=qualityValue;
				if(qualityValue>_maximumQuality)
					_maximumQuality=qualityValue;
			}
			if(id!=""){
				readNames.insert(id);
				file_bases+=sequence.length();
				bases+=sequence.length();
				Read*read=new Read(id.c_str(),sequence.c_str(),quality.c_str());
				reads.push_back(read);
				file_reads++;
			}
		}
		f.close();
		if(file_reads==0)
			continue;
		m_cout<<"   reads: "<<file_reads<<" "<<endl;
		m_cout<<"   bases: "<<file_bases<<endl;
		m_cout<<"   average read length: "<<file_bases/file_reads<<endl;
		m_cout<<"   average quality: "<<file_quality/file_bases<<endl;
		m_cout<<"   minimum quality = "<<_minimumQuality<<endl;
		m_cout<<"   maximum quality = "<<_maximumQuality<<endl;
	}

	
	m_cout<<endl;
	m_cout<<" Total reads: "<<reads.size()<<endl;
	m_cout<<endl;



	// starting the assembler
	m_cout<<"Starting  De Bruijn assembler!."<<endl;
	DeBruijnAssembler assembler(&m_cout);
	// word size must be odd (2k+1)
	// to avoid palindromes
	assembler.setMinimumQuality(minimumQuality);
	assembler.setWordSize(wordSize);
	assembler.setMinimumCoverage(minimumCoverage);
	assembler.buildGraph(&reads);
	assembler.run_Assembler();


	m_cout<<"Writing contigs"<<endl;
	string contigsFile=assemblyDirectory+"/Contigs.fa";
	ofstream output(contigsFile.c_str());
	assembler.setMinimumContigSize(minimumContigSize);
	assembler.outputContigs(&output);
	output.close();


	m_cout<<"Done."<<endl;
	// finish
	return 0;
}
