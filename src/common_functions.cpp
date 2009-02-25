#include"common_functions.h"
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
	VERTEX_TYPE i=0;
	for(int j=0;j<(int)strlen(a);j++){
		VERTEX_TYPE k=0;
		if(a[j]=='A'){
			// binary=00
			// dec = 0
			k=m_NUCLEOTIDE_A;
		}else if(a[j]=='T'){
			// binary=01
			//  dec = 1
			k=m_NUCLEOTIDE_T;
		}else if(a[j]=='C'){
			// binary=10
			// dec = 2
			k=m_NUCLEOTIDE_C;
		}else if(a[j]=='G'){
			// binary=11
			// dec = 3
			k=m_NUCLEOTIDE_G;
		}
		for(int l=0;l<=j;l++){
			k*=4; // right shift two times two positions
		}
		i+=k;
	}
	return i;
}

string idToWord(VERTEX_TYPE i,int wordSize){
	string a="";
	int maxSize=sizeof(VERTEX_TYPE)*8/2; // 32
	for(int p=0;p<wordSize;p++){
		VERTEX_TYPE j=i;
		for(int k=0;k<(maxSize-p-2);k++){
			j*=4;
		}
		for(int k=0;k<(maxSize-1);k++){
			j/=4;
		}
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
	*out<<"Version: "<<"$Id: DeBruijnAssembler.cpp 180 2009-02-22 21:55:22Z boiseb01 $"<<endl;
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
