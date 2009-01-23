/*
	dna: De Novo Assembler
    Copyright (C) 2008, 2009 Sébastien Boisvert
	$Id: De_Bruijn_De_Novo_Assembler_main.cpp 298 2009-01-17 01:52:03Z sebhtml $

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
#include"SequenceData.h"
#include"Loader.h"
#include"DeBruijnAssembler.h"
#include"CustomMap.hpp"
#include<fstream>
#include<stdint.h>
#include<stdlib.h>


using namespace std;


int main(int argc,char*argv[]){
	ostringstream buffer;

	// show usage
	buffer<<"Welcome to dna, De Novo Assembler."<<endl;
	buffer<<"Notice: You need a 64-bits processor"<<endl;
	buffer<<"Website: http://DeNovoAssembler.sf.net"<<endl;
	buffer<<"Version: "<<SOFTWARE_VERSION<<endl;
	buffer<<"License: GPLv3 (see http://fsf.org/)"<<endl;
	buffer<<endl;
	buffer<<"Usage:"<<endl;
	buffer<<"dna [options] <files (fasta, sff, or fastq)>"<<endl;

	// show options
	buffer<<" OPTIONS"<<endl;
	vector<string> inputFiles;

	string assemblyDirectory="Assembly";
	buffer<<" -assemblyDirectory   default: "<<assemblyDirectory<<endl;
	buffer<<"                      description: the directory where files will be written."<<endl;
	int wordSize=25;
	buffer<<" -wordSize            default: "<<wordSize<<" (maximum: 31)"<<endl;
	buffer<<"                      description: the length of strings inside vertices, edges will be defined for words of length <wordSize>+1."<<endl;
	string m_minimumCoverageParameter="auto";
	buffer<<" -minimumCoverage     default: auto (with depletion curve)"<<endl;
	
	int minimumContigSize=500;
	//buffer<<" -minimumContigSize   default: "<<minimumContigSize<<endl;
	//buffer<<"                      description: the minimum length of contigs generated with the graph."<<endl;
	uint64_t buckets=100000000;
	buffer<<" -buckets             default: "<<buckets<<endl;
	buffer<<"                      description: number of buckets"<<endl;
	string pairedInfo="none";
	buffer<<" -pairedInfo          default: none"<<endl;
	buffer<<"                      description: a file that contains paired information (For Solexa, and others)"<<endl;
	buffer<<"                                   Note that the files must be in the same order in the <pairedInfo> file and in the command line."<<endl;
	buffer<<endl;


	// collect arguments
	for(int i=1;i<argc;i++){
		string option=argv[i];
		if(option=="-assemblyDirectory"){
			i++;
			assemblyDirectory=argv[i];
		}else if(option=="-minimumContigSize"){
			i++;
			minimumContigSize=atoi(argv[i]);
		}else if(option=="-pairedInfo"){
			i++;
			pairedInfo=argv[i];
		}else if(option=="-minimumCoverage"){
			i++;
			m_minimumCoverageParameter=argv[i];
		}else if(option=="-buckets"){
			i++;
			buckets=atoi(argv[i]);
		}else if(option=="-wordSize"){
			i++;
			if(wordSize>31){
				buffer<<"Wordsize cannot be greater than 31"<<endl;
			}else{
				wordSize=atoi(argv[i]);
			}
		}else{
			inputFiles.push_back(argv[i]);
		}
	}

	buffer<<endl;
	buffer<<endl;
	buffer<<"dna -assemblyDirectory "<<assemblyDirectory<<" -minimumCoverage "<<m_minimumCoverageParameter<<" -buckets "<<buckets<<" -pairedInfo "<<pairedInfo<<  " -wordSize "<<wordSize <<" -minimumContigSize "<<minimumContigSize;
	for(int i=0;i<(int)inputFiles.size();i++)
		buffer<<" "<<inputFiles[i];
	buffer<<endl;



	string command=" mkdir -p "+assemblyDirectory+" # [dna] ";
	system(command.c_str());
	
	string logFile=assemblyDirectory+"/Log.txt";
	ofstream m_cout(logFile.c_str());
	cout<<buffer.str()<<endl;
	// switching to m_cout instead of cout
	cout<<endl;
	cout<<"See "<<logFile<<endl;


	m_cout<<endl;
	if(inputFiles.size()==0){
		m_cout<<"Error: no files provided."<<endl;
		exit(0);
	}
	m_cout<<buffer.str();
	m_cout<<endl;
	m_cout<<"Writing log to "<<assemblyDirectory+"/Log.txt"<<endl;

	m_cout<<endl;

	// TODO: a class SequenceData
	// to avoid loading all files in memory (only 1 at any moment)
	m_cout<<"Indexing files"<<endl;
	m_cout<<endl;
	SequenceData sequenceData(&inputFiles,&m_cout);

	if(sequenceData.size()==0){
		m_cout<<"Error: no reads provided."<<endl;
		exit(0);
	}


	m_cout<<endl;
	m_cout<<" Total reads: "<<sequenceData.size()<<endl;
	m_cout<<" Total bases: "<<sequenceData.bases()<<endl;
	m_cout<<endl;



	// starting the assembler
	m_cout<<"Starting  De Bruijn assembler!."<<endl;
	DeBruijnAssembler assembler(&m_cout);
	// word size must be odd (2k+1)
	// to avoid palindromes
	assembler.setBuckets(buckets);
	assembler.setWordSize(wordSize);
	//assembler.setMinimumContigSize(minimumContigSize);
	assembler.setMinimumCoverage(m_minimumCoverageParameter);
	assembler.setAssemblyDirectory(assemblyDirectory);
	assembler.setPairedInfo(pairedInfo);
	assembler.buildGraph(&sequenceData);
	assembler.Algorithm_Assembler_20090121();

	assembler.outputContigs();

	m_cout<<endl;
	m_cout<<"Files written in "+assemblyDirectory<<endl;
	m_cout<<"Done."<<endl;
	m_cout.close();
	return 0;
}
