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
#include"CustomMap.hpp"
#include<fstream>
#include<stdint.h>
#include<stdlib.h>


using namespace std;


int main(int argc,char*argv[]){
	DeBruijnAssembler::CommonHeader(&cout);
	// show usage
	cout<<endl;
	cout<<"Usage:"<<endl;
	cout<<"dna_DeBruijnAssembler [options] <files (fasta, sff, or fastq)>"<<endl;

	// show options
	cout<<" OPTIONS"<<endl;
	vector<string> inputFiles;

	string assemblyDirectory="Assembly";
	cout<<" -assemblyDirectory   default: "<<assemblyDirectory<<endl;
	cout<<"                      description: the directory where files will be written."<<endl;
	int wordSize=21;
	cout<<" -wordSize            default: "<<wordSize<<" (maximum: 31)"<<endl;
	cout<<"                      description: the length of strings inside vertices, edges will be defined for words of length <wordSize>+1."<<endl;
	string m_minimumCoverageParameter="2";
	cout<<" -minimumCoverage     default: "<<m_minimumCoverageParameter<<""<<endl;
	cout<<"    auto if you want to use the depletion curve"<<endl;
	int minimumContigSize=500;
	//cout<<" -minimumContigSize   default: "<<minimumContigSize<<endl;
	//cout<<"                      description: the minimum length of contigs generated with the graph."<<endl;
	uint64_t buckets=100000000;
	cout<<" -buckets             default: "<<buckets<<endl;
	cout<<"                      description: number of buckets"<<endl;
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
		if(option=="-assemblyDirectory"){
			i++;
			assemblyDirectory=argv[i];
		}else if(option=="-pairedInfo"){
			i++;
			pairedInfo=argv[i];
		}else if(option=="-minimumCoverage"){
			i++;
			m_minimumCoverageParameter=argv[i];
		}else if(option=="-debug"){
			DEBUGMODE=true;
		}else if(option=="-buckets"){
			i++;
			buckets=atoi(argv[i]);
		}else if(option=="-wordSize"){
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

	cout<<endl;
	cout<<endl;
	cout<<"  -assemblyDirectory="<<assemblyDirectory<<endl;
	cout<<"  -minimumCoverage="<<m_minimumCoverageParameter<<endl;
	cout<<"  -buckets="<<buckets<<endl;
	cout<<"  -wordSize="<<wordSize<<endl;
	cout<<" <FILES>"<<endl;
	//cout<<"dna_DeBruijnAssembler -assemblyDirectory "<<assemblyDirectory<<" -minimumCoverage "<<m_minimumCoverageParameter<<" -buckets "<<buckets<<" -pairedInfo "<<pairedInfo<<  " -wordSize "<<wordSize ;
	for(int i=0;i<(int)inputFiles.size();i++)
		cout<<" "<<inputFiles[i]<<endl;
	cout<<endl;



	string command=" mkdir -p "+assemblyDirectory+" # [dna] ";
	system(command.c_str());
	
	// switching to cout instead of cout
	cout<<endl;


	cout<<endl;
	if(inputFiles.size()==0){
		cout<<"Error: no files provided."<<endl;
		exit(0);
	}
	cout<<endl;

	cout<<endl;

	// TODO: a class SequenceData
	// to avoid loading all files in memory (only 1 at any moment)
	//cout<<"Indexing files"<<endl;
	cout<<endl;
	cout<<"********** Loading sequence data..."<<endl;
	SequenceDataFull sequenceData(&inputFiles,&cout);

	if(sequenceData.size()==0){
		cout<<"Error: no reads provided."<<endl;
		exit(0);
	}


	cout<<endl;
	cout<<"Total reads: "<<sequenceData.size()<<endl;
	cout<<"Total bases: "<<sequenceData.bases()<<endl;
	cout<<endl;



	// starting the assembler
	//cout<<"Starting  De Bruijn assembler!."<<endl;
	DeBruijnAssembler assembler(&cout);
	// word size must be odd (2k+1)
	// to avoid palindromes
	assembler.setBuckets(buckets);
	assembler.setWordSize(wordSize);
	//assembler.setMinimumContigSize(minimumContigSize);
	assembler.setMinimumCoverage(m_minimumCoverageParameter);
	assembler.setAssemblyDirectory(assemblyDirectory);
	assembler.setPairedInfo(pairedInfo);
	if(DEBUGMODE)
		assembler.debug();

	assembler.buildGraph(&sequenceData);
	assembler.Algorithm_Assembler_20090121();

	//assembler.outputContigs();

	cout<<endl;
	cout<<"Files written in "+assemblyDirectory<<endl;
	cout<<"Done."<<endl;
	string hawkeyeFile=assemblyDirectory+"/RunHawkeye.sh";
	string bankFile=assemblyDirectory+"/CreateBank.sh";
	ofstream fileStreamBank(bankFile.c_str());
	ofstream fileStream(hawkeyeFile.c_str());
	string currentDirectory;
	string pwdFile=assemblyDirectory+"/pwd.txt";
	string pwdCommand="pwd > "+pwdFile;
	system(pwdCommand.c_str());
	ifstream pwdFileStream(pwdFile.c_str());
	pwdFileStream>>currentDirectory;
	pwdFileStream.close();
	fileStream<<"if test -d bank"<<endl;
	fileStream<<"then"<<endl;
	fileStream<<"	true"<<endl;
	fileStream<<"else"<<endl;
	fileStream<<"	bash CreateBank.sh > bank.log"<<endl;
	fileStreamBank<<"bank-transact -m "<<AMOS_FILE_NAME<<" -b bank -c"<<endl;
	fileStreamBank<<"cat";
	for(int i=0;i<inputFiles.size();i++){
		if(inputFiles[i][0]!='/'){
			fileStreamBank<<" "<<currentDirectory<<"/"<<inputFiles[i];
		}else{
			fileStreamBank<<" "<<inputFiles[i];
		}
	}
			
	fileStreamBank<<" > reads.fasta"<<endl;
	fileStreamBank<<"dna_fastaToAMOS reads.fasta reads.afg"<<endl;
	fileStreamBank<<"bank-transact -m reads.afg -b bank"<<endl;
	fileStream<<"fi"<<endl;
	fileStream<<"hawkeye bank"<<endl;
	fileStreamBank.close();
	fileStream.close();
	string mergerFile=assemblyDirectory+"/Merge.sh";
	ofstream fileStreamMerger(mergerFile.c_str());
	fileStreamMerger<<"dna_Merger "<<FASTA_FILE_NAME<<" Merger.fasta > merger.log"<<endl;
	fileStreamMerger.close();
	return 0;
}
