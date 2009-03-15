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



#include"common_functions.h"
#include<iostream>
#include<string>
#include<cstring>
#include<sstream>
using namespace std;

int m_NUCLEOTIDE_A=0;
int m_NUCLEOTIDE_T=1;
int m_NUCLEOTIDE_C=2;
int m_NUCLEOTIDE_G=3;


char complement(char a){
	if(a=='A')
		return 'T';
	if(a=='T')
		return 'A';
	if(a=='C')
		return 'G';
	if(a=='G')
		return 'C';
	return a;
}

string reverseComplement(string a){
	ostringstream i;
	for(int p=a.length()-1;p>=0;p--){
		char b=complement(a[p]);
		i<< b;
	}
	return i.str();
}

// convert k-mer to VERTEX_TYPE
VERTEX_TYPE wordId(const char*a){
	uint64_t i=0;
	for(int j=0;j<(int)strlen(a);j++){
		VERTEX_TYPE k=0;
		if(a[j]=='A'){
			k=0;
		}else if(a[j]=='T'){
			k=1;
		}else if(a[j]=='C'){
			k=2;
		}else if(a[j]=='G'){
			k=3;
		}
		i=(i|(k<<(j<<1)));
	}
	return i;
}

string idToWord(VERTEX_TYPE i,int wordSize){
	string a="";
	for(int p=0;p<wordSize;p++){
		VERTEX_TYPE j=(i<<(62-2*p))>>62;
		if(j==0){
			a+='A';
		}else if(j==1){
			a+='T';
		}else if(j==2){
			a+='C';
		}else if(j==3){
			a+='G';
		}else{
		}
	}
	return a;
}

void CommonHeader(ostream*out){
	*out<<"********** Starting..."<<endl;
	*out<<endl;
	*out<<"DNA: De Novo Assembler"<<endl;
	*out<<"Documentation: http://denovoassembler.sf.net/"<<endl;
	*out<<"License: http://www.gnu.org/licenses/gpl.html"<<endl;
	*out<<"Publication: in preparation"<<endl;
	*out<<"Version: "<<"$Id$"<<endl;
}

char getLastSymbol(VERTEX_TYPE i,int m_wordSize){
        VERTEX_TYPE j=i;
        for(int i=0;i<m_wordSize;i++){
                j/=2;
                j/=2;
        }

        if((int)j==m_NUCLEOTIDE_A)
                return 'A';
        if((int)j==m_NUCLEOTIDE_T)
                return 'T';
        if((int)j==m_NUCLEOTIDE_C)
                return 'C';
        if((int)j==m_NUCLEOTIDE_G)
                return 'G';
        return 'E';
}

bool isValidDNA(const char*x){
	int len=strlen(x);
	for(int i=0;i<len;i++){
		char a=x[i];
		if(!(a=='A'||a=='T'||a=='C'||a=='G'))
			return false;
	}
	return true;
}


/*
 * 
 *   63 62 ... 1 0
 *
 */

void coutBIN(uint64_t a){
	cout<<hex<<a<<dec<<endl;
	for(int i=63;i>=0;i--){
		cout<<(int)((a<<(63-i))>>63);
	}
	cout<<endl;
}
