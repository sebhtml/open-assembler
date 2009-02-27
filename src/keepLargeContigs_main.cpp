/*
	dna: De Novo Assembler
    Copyright (C) 2008, 2009 Sébastien Boisvert
	$Id: DeBruijnAssembler.cpp 116 2009-02-16 21:19:41Z boiseb01 $

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




#include"Loader.h"
#include<fstream>
#include<string>
#include<iostream>
using namespace std;

int main(int argc,char*argv[]){
	if(argc!=4){
		cout<<"Usage"<<endl;
		cout<<"<program> <contigs> <minimumLength> <output>"<<endl;
		return 0;
	}
	
	string contigsFile=argv[1];
	int minimumLength=atoi(argv[2]);
	string output=argv[3];
	vector<Read*> contigs;
	Loader loader;
	loader.load(contigsFile,&contigs);
	int columns=60;
	ofstream f(output.c_str());
	for(vector<Read*>::iterator i=contigs.begin();i!=contigs.end();i++){
		string sequence=(*i)->getSeq();
		if(sequence.length()>=minimumLength){
			f<<">"<<(*i)->getId()<<endl;
			int j=0;
			while(j<sequence.length()){
				f<<sequence.substr(j,columns)<<endl;
				j+=columns;
			}
		}
	}
	f.close();
	return 0;
}
