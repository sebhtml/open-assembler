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


using namespace std;



int main(int argc,char*argv[]){
	DeBruijnAssembler::CommonHeader(&cout);
	if(argc!=3){
		cout<<"usage"<<endl;
		cout<<"dna_fastaToAMOS fastaFile afgFile"<<endl;
		return 0;
	}
	string contigsFile=argv[1];
	string outputFile=argv[2];
	vector<Read*> contigs;
	Loader loader;
	loader.load(contigsFile,&contigs);
	ofstream f(outputFile.c_str());
/*
 {RED
iid:1
eid:E7DK28409FQDJD
seq:
TAGTATTGGATAAGTTGCTATCGTGGTTGGGTTTGAAAGTGCACCACCAATCCCCAAAAG
TAGACCTGCTGCAGGAAGTATAGCTATAGGTAACATAAAAGCTTTTCACTTT
.
qlt:
DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
.
clr:0,112
}
*/
	for(int i=0;i<contigs.size();i++){
		f<<"{RED"<<endl;
		f<<"iid:"<<i+1<<endl;
		f<<"eid:"<<contigs[i]->getId()<<endl;
		f<<"seq:"<<endl;
		f<<contigs[i]->getSeq()<<endl;
		f<<"."<<endl;
		f<<"qlt:"<<endl;
		for(int j=0;j<strlen(contigs[i]->getSeq());j++)
			f<<"D";
		f<<endl;
		f<<"."<<endl;
		f<<"clr:0,"<<strlen(contigs[i]->getSeq())-1<<endl;
		f<<"}"<<endl;
	}
	f.close();
	return 0;
}

