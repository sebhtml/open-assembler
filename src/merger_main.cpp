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

#include<iostream>
#include"Read.h"
#include<set>
#include<map>
#include<string>
#include"Loader.h"
#include"DeBruijnAssembler.h"
#include"Merger.h"

using namespace std;

int main(int argc,char*argv[]){
	DeBruijnAssembler::CommonHeader(&cout);
	cout<<"This is dna_Merger."<<endl;
	if(argc!=3){
		cout<<"usage"<<endl;
		cout<<"dna_Merger contigs.fa mergedContig.fa"<<endl;
		return 0;
	}
	int minimumContigSize=500;
	string contigsFile=argv[1];
	string outputFile=argv[2];
	vector<Read*> contigs;
	Loader loader;
	loader.load(contigsFile,&contigs);

	//for(int i=0;i<contigs.size();i++)
		//cout<<strlen(contigs[i]->getSeq())<<endl;
	Merger merger;
	vector<Read*> theBestContigs=merger.mergeContigs(contigs);

	ofstream f(outputFile.c_str());
	for(int j=0;j<theBestContigs.size();j++){
		string sequence=theBestContigs[j]->getSeq();
		f<<">"<<theBestContigs.at(j)->getId()<<" "<<sequence.length()<<" nucleotides"<<endl;
		cout<<">"<<theBestContigs.at(j)->getId()<<" "<<sequence.length()<<" nucleotides"<<endl;
		int columns=70;
		int k=0;
		while(k<sequence.length()){
			f<<sequence.substr(k,columns)<<endl;
			k+=columns;
		}
	}
	f.close();
	return 0;
}

