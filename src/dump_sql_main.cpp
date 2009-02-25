#include<mysql.h>
#include"Read.h"
#include<fstream>
#include<sstream>
#include<string>
#include"DeBruijnAssembler.h"
#include<iostream>
#include<stdlib.h>
using namespace std;

bool validMer(string*a,int wordSize){
	if(a->length()!=wordSize)
		return false;
	for(int i=0;i<a->length();i++){
		char b=a->at(i);
		if(b=='A'||b=='T'||b=='C'||b=='G'){
		}else{
			return false;
		}
	}
	return true;
}

void process_Read(MYSQL*mysql,int readNumber,int wordSize,string*readSequence,char readStrand){
	string theSequence=*readSequence;
	if(readStrand=='R')
		theSequence=reverseComplement(*readSequence);
	for(int readPosition=0;readPosition<theSequence.length();readPosition++){
		string prefix=theSequence.substr(readPosition,wordSize);
		string suffix=theSequence.substr(readPosition+1,wordSize);
		bool prefixValid=validMer(&prefix,wordSize);
		bool suffixValid=validMer(&suffix,wordSize);
		VERTEX_TYPE prefixInt;
		VERTEX_TYPE suffixInt;
		if(prefixValid){
			prefixInt=wordId(prefix.c_str());
			ostringstream query;
			query<<"insert into vertex_annotations(vertex,readNumber,readStrand,readPosition) values ("<<prefixInt<<","<<readNumber<<",'"<<readStrand<<"',"<<readPosition<<")";
			//cout<<query.str()<<endl;
			mysql_query(mysql,query.str().c_str());
		}
		if(suffixValid){
			suffixInt=wordId(suffix.c_str());
			ostringstream query;
			query<<"insert into vertex_annotations(vertex,readNumber,readStrand,readPosition) values("<<suffixInt<<","<<readNumber<<",'"<<readStrand<<"',"<<readPosition+1<<")";
			mysql_query(mysql,query.str().c_str());
			//cout<<query.str()<<endl;
		}
		if(prefixValid&&suffixValid){
			ostringstream query;
			query<<"insert into edges (prefix,suffix) values("<<prefixInt<<","<<suffixInt<<")";
			mysql_query(mysql,query.str().c_str());
			//cout<<query.str()<<endl;
		}
	}

/*
	cout<<"ReadNumber="<<readNumber<<endl;
	cout<<"Nucleotides="<<readSequence->length()<<endl;
	cout<<*readSequence<<endl;
*/
	
}

int main(int argc,char*argv[]){
	MYSQL mysql;
	MYSQL_ROW row;

	cout<<"Usage:"<<endl;
	cout<<"dna_FastaToSQL <wordSize> <fasta files>"<<endl;
	if(argc<2)
		return 0;
	mysql_init(&mysql);
	mysql_options(&mysql,MYSQL_READ_DEFAULT_GROUP,"your_prog_name");

	if (!mysql_real_connect(&mysql,"localhost","root","root","assembly",0,NULL,0)){
    		fprintf(stderr, "Failed to connect to database: Error: %s\n",
          	mysql_error(&mysql));
	}
	mysql_query(&mysql,"show tables");
	unsigned int num_fields;
 	MYSQL_RES *result=mysql_use_result(&mysql)  ;
	num_fields = mysql_num_fields(result);
	
	cout<<"Checking tables"<<endl;
	while(row = mysql_fetch_row(result)){
		for(int i=0;i<num_fields;i++){
			cout<<" "<<row[i];
		}
		cout<<endl;
	}

	cout<<"Done"<<endl;

	int wordSize=atoi(argv[1]);
	cout<<"Wordsize: "<<wordSize<<endl;
	cout<<argc-2<<" files"<<endl;
	int readNumber=0;
	for(int fileNumber=2;fileNumber<argc;fileNumber++){
		if(readNumber%1000==0){
			cout<<"Progress: "<<readNumber<<endl;
		}
		string fileName=argv[fileNumber];
		ifstream f(fileName.c_str());
		ostringstream sequence;
		char buffer[2000];
		while(!f.eof()){
			f.getline(buffer,2000);
			string line=buffer;
			//cout<<"Line="<<line<<endl;
			if(line[0]=='>'){
				string readSequence=sequence.str();
				if(readSequence.length()>0){
					process_Read(&mysql,readNumber,wordSize,&readSequence,'F');
					process_Read(&mysql,readNumber,wordSize,&readSequence,'R');
					readNumber++;
				}
				sequence.str("");
			}else{
				sequence<<line;
			}
		}
		string readSequence=sequence.str();
		process_Read(&mysql,readNumber,wordSize,&readSequence,'F');
		process_Read(&mysql,readNumber,wordSize,&readSequence,'R');

		f.close();
		cout<<"Loading "<<fileName<<endl;
	}
	cout<<"Reads: "<<readNumber<<endl;
	cout<<"Closing MYSQL now."<<endl;
	mysql_close(&mysql);
	return 0;
}
