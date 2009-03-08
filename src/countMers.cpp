#include<string>
#include<iostream>
#include<map>
#include<sstream>
#include<stdint.h>
#include<fstream>
#include<hash_map>
using namespace std;
using namespace __gnu_cxx;

string revComp(string a){
	ostringstream c;
	for(int b=a.length()-1;b>=0;b--){
		char d=a[b];
		if(d=='A')
			c<<'T';
		else if(d=='T')
			c<<'A';
		else if(d=='C')
			c<<'G';
		else if(d=='G')
			c<<'C';
		else
			c<<d;
	}
	return c.str();
}

bool correct(string a){
	for(int i=0;i<a.length();i++){
		if(!(a[i]=='A'||a[i]=='T'||a[i]=='C'||a[i]=='G'))
			return false;
	}
	return true;
}

int main(){
	string  file="ETJITFZ02.sff_normal.fasta";
	ifstream f(file.c_str());
	int k=21;
	char buffer[10000];
	map<string,int> dataStore;
	int r=0;
	while(!f.eof()){
		if(r%10000==0)
			cout<<r<<endl;
		f.getline(buffer,10000);
		f.getline(buffer,10000);
		string dna=buffer;
		for(int p=0;p<dna.length();p++){
			string word=dna.substr(p,k+1);
			if(word.length()==k+1&&correct(word)){
				dataStore[word]++;
				dataStore[revComp(word)]++;
			}
		}
		r++;
	}
	f.close();

	int c=2;
	ofstream f2("Edges.txt");
	for(map<string,int>::iterator i=dataStore.begin();i!=dataStore.end();i++){
		if(i->second>=c)
			f2<<i->first<<endl;
	}
	f2.close();
	return 0;
}

