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
#include<iostream>
#include"Read.h"
#include"SequenceDataFull.h"
#include"Loader.h"
#include"DeBruijnAssembler.h"
#include<fstream>
#include<stdint.h>
#include<stdlib.h>
#include"Merger.h"
#include"common_functions.h"
#include"VertexData.h"

using namespace std;


int main(int argc,char*argv[]){
	//cout<<"VertexDataSize <- "<<sizeof(VertexData)<<endl;
	CommonHeader(&cout);
	// show usage
	cout<<endl;
	cout<<"Usage:"<<endl;
	cout<<argv[0]<<" [options] <files (fasta, sff, or fastq)>"<<endl;

	// show options
	cout<<" OPTIONS"<<endl;
	vector<string> inputFiles;
	
	bool assemblyDirectoryWasSet=false;
	string assemblyDirectory="Assembly";
	cout<<" -directory   default: "<<assemblyDirectory<<endl;
	cout<<"                      description: the directory where files will be written."<<endl;
	int wordSize=21;
	cout<<" -wordSize            default: "<<wordSize<<" (maximum: 31)"<<endl;
	cout<<"                      description: the length of strings inside vertices, edges will be defined for words of length <wordSize>+1."<<endl;
	string m_minimumCoverageParameter="2";
	cout<<" -minimumCoverage     default: "<<m_minimumCoverageParameter<<""<<endl;
	cout<<"    auto if you want to use the distribution curve"<<endl;
	bool DEBUGMODE=false;
	//cout<<" [ -debug ]"<<endl;
	string pairedInfo="none";
	//cout<<" -pairedInfo          default: none"<<endl;
	//cout<<"                      description: a file that contains paired information (For Solexa, and others)"<<endl;
	//cout<<"                                   Note that the files must be in the same order in the <pairedInfo> file and in the command line."<<endl;
	cout<<endl;


	// collect arguments
	for(int i=1;i<argc;i++){
		string option=argv[i];
		if(option=="-directory"&&i!=argc-1){
			i++;
			assemblyDirectoryWasSet=true;
			assemblyDirectory=argv[i];
		}else if(option=="-minimumCoverage"&&i!=argc-1){
			i++;
			m_minimumCoverageParameter=argv[i];
		}else if(option=="-wordSize"&&i!=argc-1){
			i++;
			if(wordSize>31){
				cout<<"Wordsize cannot be greater than 31"<<endl;
			}else{
				wordSize=atoi(argv[i]);
			}
		}else{
			inputFiles.push_back(argv[i]);
		}
	}

	if(assemblyDirectoryWasSet==false){
		cout<<"Error: please provide a value for -directory"<<endl;
		return 0;
	}
	cout<<"  -directory="<<assemblyDirectory<<endl;
	cout<<"  -minimumCoverage="<<m_minimumCoverageParameter<<endl;

	cout<<"  -wordSize="<<wordSize<<endl;
	cout<<" <FILES>"<<endl;


	string command=" mkdir -p "+assemblyDirectory+" # [dna] ";
	system(command.c_str());

	string filesFile=assemblyDirectory+"/InputFiles.txt";
	ofstream aFilesStream(filesFile.c_str());
	for(int i=0;i<(int)inputFiles.size();i++){
		cout<<" "<<inputFiles[i]<<endl;
		aFilesStream<<inputFiles[i]<<endl;
	}
	aFilesStream.close();
	cout<<endl;

	cout<<"writing files to "<<filesFile<<endl;

	
	// switching to cout instead of cout
	if(inputFiles.size()==0){
		cout<<"Error: no files provided."<<endl;
		exit(0);
	}

	//cout<<"Indexing files"<<endl;
	cout<<"********** Loading sequence data..."<<endl;
	SequenceDataFull sequenceData(&inputFiles,&cout);

	if(sequenceData.size()==0){
		cout<<"Error: no reads provided."<<endl;
		exit(0);
	}


	cout<<endl;
	cout<<"Total reads: "<<sequenceData.size()<<endl;
	cout<<"Total bases: "<<sequenceData.bases()<<endl;



	// starting the assembler
	//cout<<"Starting  De Bruijn assembler!."<<endl;
	DeBruijnAssembler assembler(&cout);
	// word size must be odd (2k+1)
	// to avoid palindromes
	//assembler.setBuckets(buckets);
	assembler.setWordSize(wordSize);
	//assembler.setMinimumContigSize(minimumContigSize);
	assembler.setMinimumCoverage(m_minimumCoverageParameter);
	assembler.setAssemblyDirectory(assemblyDirectory);
	//assembler.setPairedInfo(pairedInfo);
	if(DEBUGMODE)
		assembler.debug();

	assembler.setSequenceData(&sequenceData);
	assembler.buildGraph();
	string applicationName=argv[0];
	cout<<"Graph is built."<<endl;	
	return 0;
	//assembler.Algorithm_Assembler_20090121();

	//assembler.outputContigs();

	//cout<<"Files written in "+assemblyDirectory<<endl;
	//cout<<"Done!"<<endl;

	//return 0;
}
