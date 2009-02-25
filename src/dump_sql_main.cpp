#include<mysql.h>
#include<fstream>
#include<sstream>
#include<string>
#include<iostream>
#include<stdlib.h>
using namespace std;

void process_Read(MYSQL*mysql,int readNumber,int wordSize,string*readSequence){
	cout<<"ReadNumber="<<readNumber<<endl;
	cout<<"Nucleotides="<<readSequence->length()<<endl;
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
		string fileName=argv[fileNumber];
		ifstream f(fileName.c_str());
		ostringstream sequence;
		char buffer[2000];
		while(!f.eof()){
			f.getline(buffer,2000);
			string line=buffer;
			if(line[0]=='>'){
				string readSequence=sequence.str();
				if(readSequence.length()>0){
					process_Read(&mysql,readNumber,wordSize,&readSequence);
					readNumber++;
				}
				sequence.clear();
			}else{
				sequence<<line;
			}
		}
		string readSequence=sequence.str();
		process_Read(&mysql,readNumber,wordSize,&readSequence);

		f.close();
		cout<<"Loading "<<fileName<<endl;
	}
	cout<<"Reads: "<<readNumber<<endl;
	cout<<"Closing MYSQL now."<<endl;
	mysql_close(&mysql);
	return 0;
}
