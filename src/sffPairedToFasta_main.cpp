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

#include"Read.h"
#include<fstream>
#include<sstream>
#include"SffLoader.h"
#include"Loader.h"
#include<iostream>
#include<vector>
#include"DeBruijnAssembler.h"

using namespace std;

int main(int argc,char*argv[]){
	DeBruijnAssembler::CommonHeader(&cout);
	cout<<argv[0]<<" <sffFile>"<<endl;
	if(argc!=2){
		cout<<"Incorrect usage"<<endl;
		return 0;
	}
	vector<Read*> reads;
	string file=argv[1];
	cout<<"Loading file"<<endl;
	Loader loader(&cout);
	loader.load(file,&reads);
	int num2=reads.size();
	string standardReads=file+"_normal.fasta";
	string paired1=file+"_paired_2.fasta";
	string paired2=file+"_paired_1.fasta";
	ofstream*standard=NULL;
	ofstream*left=NULL;
	ofstream*right=NULL;
	string linker_FLX="GTTGGAACCGAAAGGGTTTGAATTCAAACCCTTTCGGTTCCAAC";
	string linker_TITANIUM="TCGTATAACTTCGTATAATGTATGCTATACGAAGTTATTACG";
	int normal=0;
	int flxLinker=0;
	int titaniumLinker=0;
	for(int i=0;i<(int)reads.size();i++){
		ostringstream name;
		string readName=reads[i]->getId();
		string sequence=reads[i]->getSeq();
		size_t positionFLX=sequence.find(linker_FLX);
		size_t positionTITANIUM=sequence.find(linker_TITANIUM);
		if(positionFLX!=string::npos){
			flxLinker++;
			if(left==NULL)
				left=new ofstream(paired1.c_str());
			if(right==NULL)
				right=new ofstream(paired2.c_str());

			string leftPart=sequence.substr(0,positionFLX);
			int multiplicity=1;
			string rightPart=sequence.substr(positionFLX+multiplicity*linker_FLX.length());
			while(rightPart.find(linker_FLX)==0){
				multiplicity++;
				rightPart=sequence.substr(positionFLX+multiplicity*linker_FLX.length());
			}
			if(leftPart.length()>0&&rightPart.length()>0){
				(*left)<<">"<<readName<<"_1 "<<leftPart.length()<<endl;
				(*left)<<leftPart<<endl;
				(*right)<<">"<<readName<<"_2 "<<rightPart.length()<<endl;
				(*right)<<rightPart<<endl;
			}
		}else if(positionTITANIUM!=string::npos){
			titaniumLinker++;
			if(left==NULL)
				left=new ofstream(paired1.c_str());
			if(right==NULL)
				right=new ofstream(paired2.c_str());
			string leftPart=sequence.substr(0,positionTITANIUM);
			string rightPart=sequence.substr(positionTITANIUM+linker_TITANIUM.length());
			(*left)<<">"<<readName<<"_1"<<endl;
			(*left)<<leftPart<<endl;
			(*right)<<">"<<readName<<"_2"<<endl;
			(*right)<<rightPart<<endl;
		}else{
			normal++;
			if(standard==NULL)
				standard=new ofstream(standardReads.c_str());
			(*standard)<<">"<<readName<<" "<<sequence.length()<<endl;
			(*standard)<<sequence<<endl;
		}
	}
	cout<<"Normal: "<<normal<<endl;
	cout<<"FLX Linker: "<<flxLinker<<endl;
	cout<<"Titanium Linker: "<<titaniumLinker<<endl;

	(cout)<<file<<"\t"<<num2<<"\t"<<loader.getBases()<<"\t"<<flxLinker<<endl;
	if(standard!=NULL){
		standard->close();
		delete standard;
		standard=NULL;
		cout<<"Wrote "<<standardReads<<endl;
	}
	if(left!=NULL){
		left->close();
		delete left;
		left=NULL;
		cout<<"Wrote "<<paired1<<endl;
	}
	if(right!=NULL){
		right->close();
		delete right;
		right=NULL;
		cout<<"Wrote "<<paired2<<endl;
	}
	cout<<endl;
	return 0;
}
