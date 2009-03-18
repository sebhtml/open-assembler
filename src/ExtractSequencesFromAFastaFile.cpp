/*
--
SÃ©bastien Boisvert
<Sebastien.Boisvert.3@ulaval.ca>
http://agora.ulaval.ca/~seboi46/
Summer Intern in Bioinformatics
Ken Dewar - Comparative Genetics
McGill University and GÃ©nome QuÃ©bec Innovation Centre
740, Dr Penfield Avenue, room 7211
MontrÃ©al (QuÃ©bec)
Canada H3A 1A4
(514) 388-3311 x00471

*/


#include<iostream>
#include<fstream>
#include<cstring>
#include<vector>
#include<string>
using namespace std;

bool verifyBufferWithEntries(char*buffer,vector<string>entries){
	for(int i=0;i<(int)entries.size();i++){
		string entry=entries[i];
		if(strstr(buffer,entry.c_str())!=NULL){
			return true;
		}
	}
	return false;
}

int main(int argc,char*argv[]){
	if(argc!=4){
		cout<<"This program extracts each sequence from"<<endl;
		cout<<"<fastaFile> whose header includes en entry"<<endl;
		cout<<"from <entriesFile>. The resulting set is saved"<<endl;
		cout<<"in the <outputFile> file."<<endl;
		cout<<"Usage:"<<endl;
		cout<<"PoolExtractSequenceFromFile.cpp <fastaFile> <entriesFile> <outputFile>"<<endl;
		return 0;
	}
	char*fastaFile=argv[1];
	vector<string> entries;
	char buffer[1024];
	char*entriesFile=argv[2];
	char*outputFile=argv[3];
	ifstream f(entriesFile);
	int acceptedSequences=0;
	int rejectedSequences=0;

	while(!f.eof()){
		buffer[0]='\0';
		f>> buffer;
		//cout<<buffer<<endl;
		if(strlen(buffer)>0){//empty lines are bad
			entries.push_back(buffer);
			//cout<<"Entry "<<buffer<<endl;
		}
	}
	ofstream out(outputFile);
	bool theCurrentSequenceMustBeSaved=false;
	f.close();
	ifstream fastaInput(fastaFile);
	while(!fastaInput.eof()){
		fastaInput.getline(buffer,1024);
		if(buffer[0]=='>'){
			theCurrentSequenceMustBeSaved=verifyBufferWithEntries(buffer,entries);
			if(theCurrentSequenceMustBeSaved){
				acceptedSequences++;
			}else{
				rejectedSequences++;
			}
		}
		if(theCurrentSequenceMustBeSaved){
			out<<buffer<<endl;
		}
	}

	fastaInput.close();
	out.close();
	cout<<"fastaFile: "<<fastaFile<<endl;
	cout<<"sequences: "<<acceptedSequences+rejectedSequences<<endl;
	cout<<"entriesFile: "<<entriesFile<<endl;
	cout<<"entries: "<<entries.size()<<endl;
	cout<<"outputFile: "<<outputFile<<endl;
	cout<<"accepted sequences: "<<acceptedSequences<<endl;
	cout<<"rejected sequences: "<<rejectedSequences<<endl;
	return 0;
}
