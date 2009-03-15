#include"common_functions.h"
#include"VertexData.h"
#include<string>
#include<vector>
#include<iostream>
using namespace std;

int main(){
	vector<string> tests;
	tests.push_back("GGGAGGCGGAAGAAACAAAATT");
	tests.push_back("GGGAGGCGGAAGAAACAAAAAC");
	tests.push_back("CCCCCCCCCCCCCCCCCCCCCC");
	for(vector<string>::iterator i=tests.begin();i!=tests.end();i++){
		string a=*i;
		cout<<a<<" "<<idToWord(wordId(a.c_str()),22)<<endl;
		cout<<wordId(a.c_str())<<" "<<wordId(idToWord(wordId(a.c_str()),22).c_str())<<endl;
		cout<<"BIN ";
		uint64_t c=wordId(a.c_str());
		coutBIN(c);
		int k=22;
		uint64_t a2=(c<<(64-22*2));
		uint64_t a3=(c>>2);
		uint64_t A=3;
		A=(A<<22*2);
		uint64_t vv=a3|A;
		cout<<idToWord(vv,22)<<endl;
	}

	VertexData aaa;
	string parent="CGGGAGGCGGAAGAAACAAAAT";
	 string prefix="GGGAGGCGGAAGAAACAAAATT";
	 string suffix= "GGAGGCGGAAGAAACAAAATTG";
	cout<<"ADding child "<<suffix<<endl;
	aaa.addChild(wordId(suffix.c_str()),22);
	vector<uint64_t> children=aaa.getChildren(wordId(prefix.c_str()),22);
	cout<<"CHILD "<<idToWord(children[0],22)<< " REF "<<suffix<<endl;

	aaa.addParent(wordId(parent.c_str()),22);
	cout<<"PARENT  "<<idToWord(aaa.getParents(wordId(prefix.c_str()),22)[0],22)<<" "<<parent<<endl;
	
	{
	string prefix="ATAGACTATCGATCAGCTAGA";
	string suffix=  "TAGACTATCGATCAGCTAGAG";
	VertexData gg;
	gg.addChild(wordId(suffix.c_str()),21);
	uint64_t prefixInt=wordId(prefix.c_str());
	vector<uint64_t> children=gg.getChildren(prefixInt,21);
	cout<<idToWord(children[0],21)<<" "<<suffix<<endl;
	}
	return 0;
}

