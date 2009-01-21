/*
	De Novo Assembler
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
#include"Loader.h"
#include"DeBruijnAssembler.h"
#include"Logger.hpp"
#include<fstream>

using namespace std;


int main(int argc,char*argv[]){
	Logger m_cout;
	// check processor
	m_cout<<(int)sizeof(unsigned long int)*8<<" bits available"<<endl;
	if(sizeof(unsigned long int)!=8){
		m_cout<<"Abording, you need a 64-bit machine."<<endl;
		return 0;
	}

	m_cout<<endl;


	// show usage
	m_cout<<"Welcome to dna, the De Novo Assembler."<<endl;
	m_cout<<"http://DeNovoAssembler.sf.net"<<endl;
	m_cout<<"License: GPLv3"<<endl;
	m_cout<<endl;
	m_cout<<"Usage:"<<endl;
	m_cout<<"dna [options] <files (fasta or fastq)>"<<endl;



	// show options
	m_cout<<" OPTIONS"<<endl;
	vector<string> inputFiles;

	string assemblyDirectory="Assembly";
	m_cout<<" -assemblyDirectory   default: "<<assemblyDirectory<<endl;
	m_cout<<"                      description: the directory where files will be written."<<endl;
	string useCache="no";
	m_cout<<" -useCache            default: "<<useCache<<endl;
	m_cout<<"                      description: use the graph file from a previous assembly in <assemblyDirectory>, if any."<<endl;
	int wordSize=25;
	m_cout<<" -wordSize            default: "<<wordSize<<" (maximum: 31)"<<endl;
	m_cout<<"                      description: the length of strings inside vertices, edges will be defined for words of length <wordSize>+1."<<endl;
	int minimumCoverage=3;
	m_cout<<" -minimumCoverage     default: "<<minimumCoverage<<""<<endl;
	m_cout<<"                      description: the number of occurances of a word needed to make it solid."<<endl;
	int minimumContigSize=500;
	m_cout<<" -minimumContigSize   default: "<<minimumContigSize<<endl;
	m_cout<<"                      description: the minimum length of contigs generated with the graph."<<endl;
	int windowSize=2;
	m_cout<<" -windowSize          default: "<<windowSize<<endl;
	m_cout<<"                      description: the window size to assess the quality of paths in the graph, it is the number of edges to be checked."<<endl;
	m_cout<<"                                   note that the algorithm removes tips, bubbles, and includes an optimization that automatically increases the cutoff to eliminate incorrect edges."<<endl;

	// collect arguments
	for(int i=1;i<argc;i++){
		string option=argv[i];
		if(option=="-assemblyDirectory"){
			i++;
			assemblyDirectory=argv[i];
		}else if(option=="-useCache"){
			i++;
			useCache=argv[i];
		}else if(option=="-minimumCoverage"){
			i++;
			minimumCoverage=atoi(argv[i]);
		}else if(option=="-minimumContigSize"){
			i++;
			minimumContigSize=atoi(argv[i]);
		}else if(option=="-windowSize"){
			i++;
			windowSize=atoi(argv[i]);
		}else if(option=="-wordSize"){
			i++;
			if(wordSize>31){
				m_cout<<"Wordsize cannot be greater than 31"<<endl;
			}else{
				wordSize=atoi(argv[i]);
			}
		}else{
			inputFiles.push_back(argv[i]);
		}
	}

	m_cout<<endl;
	m_cout<<" Your command was:"<<endl;
	m_cout<<endl;
	m_cout<<"dna -assemblyDirectory "<<assemblyDirectory<<" -useCache "<<useCache<<" -wordSize "<<wordSize<<" -minimumCoverage "<<minimumCoverage<<" -windowSize "<<windowSize<<" -minimumContigSize "<<minimumContigSize;
	for(int i=0;i<inputFiles.size();i++)
		m_cout<<" "<<inputFiles[i];
	m_cout<<endl;

	m_cout<<endl;
	if(inputFiles.size()==0){
		m_cout<<"Error: no files provided."<<endl;
		exit(0);
	}

	string command=" mkdir -p "+assemblyDirectory+" # [dna] ";
	system(command.c_str());
	m_cout.setOutput(assemblyDirectory+"/Log.txt");
	
	m_cout<<endl;
	m_cout<<"Writing log to "<<assemblyDirectory+"/Log.txt"<<endl;

	m_cout<<endl;




	// loading reads
	vector<Read*> reads;
	m_cout<<"Loading reads"<<endl;
	unsigned long int total_bases_count=0;
	for(int i=0;i<(int)inputFiles.size();i++){
		Loader loader(&m_cout);
		loader.load(inputFiles[i],&reads);
		total_bases_count+=loader.getBases();
	}

	if(reads.size()==0){
		m_cout<<"Error: no reads provided."<<endl;
		exit(0);
	}


	m_cout<<endl;
	m_cout<<" Total reads: "<<reads.size()<<endl;
	m_cout<<" Total bases (including low quality bases): "<<total_bases_count<<endl;
	m_cout<<endl;



	// starting the assembler
	m_cout<<"Starting  De Bruijn assembler!."<<endl;
	DeBruijnAssembler assembler(&m_cout);
	// word size must be odd (2k+1)
	// to avoid palindromes
	assembler.setWordSize(wordSize);
	assembler.setMinimumContigSize(minimumContigSize);
	assembler.setMinimumCoverage(minimumCoverage);
	assembler.setWindowSize(windowSize);
	assembler.setUseCache(useCache);
	assembler.setAssemblyDirectory(assemblyDirectory);
	assembler.buildGraph(&reads);
	//assembler.run_New_Algorithm_Assembler();
	assembler.run_New_Algorithm_Assembler_20090102();

	assembler.outputContigs();

	m_cout<<endl;
	m_cout<<"Files written in "+assemblyDirectory<<endl;
	m_cout<<"Done."<<endl;
	// finish
	return 0;
}
