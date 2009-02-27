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




#include"DeBruijnAssembler.h"
#include<fstream>
#include<string>
#include<vector>
#include<iostream>
using namespace std;

int main(int argc,char*argv[]){
	if(argc!=3){
		cout<<"Usage"<<endl;
		cout<<argv[0]<<" -directory <directory>"<<endl;
		return 0;
	}
	string directory=argv[2];
	string filesFile=directory+"/InputFiles.txt";
	string parametersFile=directory+"/Parameters.txt";
	DeBruijnAssembler assembler(&cout);
	assembler.setAssemblyDirectory(directory);
	assembler.loadParameters();
	assembler.load_graphFrom_file();
	vector<string> files;
	ifstream input(filesFile.c_str());
	cout<<"Loading files from "<<filesFile<<endl;
	string aFile;
	if(!input){
		cout<<"InputFiles.txt is not present.."<<endl;	
		return 0;
	}
	while(!input.eof()){
		aFile="";
		input>>aFile;
		cout<<aFile<<endl;
		if(aFile!="")
			files.push_back(aFile);
	}
	input.close();
	SequenceDataFull sequenceData(&files,&cout);
	assembler.setSequenceData(&sequenceData);
	assembler.Algorithm_Assembler_20090121();
	return 0;
}
