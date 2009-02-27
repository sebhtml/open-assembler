
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
