
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
	string aFile;
	while(input.eof()==false){
		aFile="";
		input>>aFile;
		if(aFile!="")
			files.push_back(aFile);
	}
	input.close();
	SequenceDataFull sequenceData(&files,&cout);
	assembler.setSequenceData(&sequenceData);
	assembler.Algorithm_Assembler_20090121();
	return 0;
}
