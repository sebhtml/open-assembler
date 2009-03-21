#include"common_functions.h"
#include<iostream>
#include<string>
using namespace std;

int main(){
	int k=21;
	string kPlus1Mer="GCTTGTTTATGTTGAGCTACAC";
	uint64_t kPlusOne=wordId(kPlus1Mer.c_str());
	uint64_t kPrefix=getKPrefix(kPlusOne,k);
	uint64_t kSuffix=getKSuffix(kPlusOne,k);
	cout<<" (k+1)-mer = "<<idToWord(kPlusOne,k+1)<<endl;
	cout<<" k-prefix =  "<<idToWord(kPrefix,k)<<endl;
	cout<<" k-suffix =   "<<idToWord(kSuffix,k)<<endl;
	
	return 0;
}
