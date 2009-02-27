/*
	dna: De Novo Assembler
    Copyright (C) 2008, 2009 Sébastien Boisvert
	$Id: DeBruijnAssembler.cpp 116 2009-02-16 21:19:41Z boiseb01 $

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

void process_Read(MYSQL*mysql,int readNumber,int wordSize,string*readSequence,char readStrand,ostream*out){

	if(readNumber%1000==0&&readStrand=='F'){
		//cout<<"Progress: "<<readNumber<<endl;
	}
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
		prefixInt=wordId(prefix.c_str());
		suffixInt=wordId(suffix.c_str());
		if(prefixValid){
			ostringstream query;
			query<<"insert into vertex_annotations(vertex,readNumber,readStrand,readPosition) values ("<<prefixInt<<","<<readNumber<<",'"<<readStrand<<"',"<<readPosition<<") ;";
			//cout<<query.str()<<endl;
			//mysql_query(mysql,query.str().c_str());
			*out<<query.str()<<endl;
		}
/*
		if(suffixValid&&readPosition!=theSequence.length()-1){
			ostringstream query;
			query<<"insert into vertex_annotations(vertex,readNumber,readStrand,readPosition) values("<<suffixInt<<","<<readNumber<<",'"<<readStrand<<"',"<<readPosition+1<<") ;";
			//mysql_query(mysql,query.str().c_str());
			//cout<<query.str()<<endl;
			*out<<query.str()<<endl;
		}
*/
		if(prefixValid&&suffixValid){
			ostringstream query;
			query<<"insert into edges (prefix,suffix) values("<<prefixInt<<","<<suffixInt<<") ;";
			//mysql_query(mysql,query.str().c_str());
			//cout<<query.str()<<endl;
			*out<<query.str()<<endl;
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

	if(argc<2){

		cout<<"Usage:"<<endl;
		cout<<"dna_FastaToSQL <wordSize> <fasta files>"<<endl;
		return 0;
	}
	mysql_init(&mysql);
	mysql_options(&mysql,MYSQL_READ_DEFAULT_GROUP,"your_prog_name");
/*
	string mysql_hostname=argv[1];
	string mysql_user=argv[2];
	string mysql_password=argv[3];
	string mysql_database=argv[4];
	if (!mysql_real_connect(&mysql,mysql_hostname.c_str(),mysql_user.c_str(),mysql_password.c_str(),mysql_database.c_str(),0,NULL,0)){
    		fprintf(stderr, "Failed to connect to database: Error: %s\n",
          	mysql_error(&mysql));
	}
*/
	//mysql_query(&mysql,"show tables");
	unsigned int num_fields;
 	//MYSQL_RES *result=mysql_use_result(&mysql)  ;
	//num_fields = mysql_num_fields(result);
	/*
	//cout<<"Checking tables"<<endl;
	while(row = mysql_fetch_row(result)){
		for(int i=0;i<num_fields;i++){
			//cout<<" "<<row[i];
		}
		//cout<<endl;
	}
*/
	//cout<<"Done"<<endl;

	int wordSize=atoi(argv[1]);
	//cout<<"Wordsize: "<<wordSize<<endl;
	//cout<<argc-2<<" files"<<endl;
	int readNumber=0;
	for(int fileNumber=2;fileNumber<argc;fileNumber++){
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
					process_Read(&mysql,readNumber,wordSize,&readSequence,'F',&cout);
					process_Read(&mysql,readNumber,wordSize,&readSequence,'R',&cout);
					readNumber++;
				}
				sequence.str("");
			}else{
				sequence<<line;
			}
		}
		string readSequence=sequence.str();
		process_Read(&mysql,readNumber,wordSize,&readSequence,'F',&cout);
		process_Read(&mysql,readNumber,wordSize,&readSequence,'R',&cout);

		f.close();
		//cout<<"Loading "<<fileName<<endl;
	}
	//graphSql.close();
	cout<<"graph.sql written"<<endl;
	cout<<"Reads: "<<readNumber<<endl;
	cout<<"Closing MYSQL now."<<endl;
	mysql_close(&mysql);
	return 0;
}
