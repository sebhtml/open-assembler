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

#include"Loader.h"
#include<iostream>
#include"Read.h"
#include"DeBruijnAssembler.h"
using namespace std;

int main(int argc,char*argv[]){
	DeBruijnAssembler::CommonHeader(&cout);
	if(argc!=3){
		cout<<"usage"<<endl;
		cout<<"dna_Scaffolder reads.fasta contigs.fasta scaffolds.fasta"<<endl;
		return 0;
	}
	string contigsFile=argv[2];
	string readsFile=argv[1];
	vector<Read*> contigs;
	Loader loader(&cout);
	loader.load(contigsFile,&contigs);
	vector<Read*> reads;
	Loader loaderReads(&cout);
	loaderReads.load(readsFile,&reads);
	return 0;
}

