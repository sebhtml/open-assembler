#include<iostream>
#include<fstream>
#include<string>
#include<set>
using namespace std;

void load(string f,set<string>*a){
	ifstream c(f.c_str());

	int i=0;
	while(!c.eof()){
		if(i%100000==0)
			cout<<i<<endl;
		string buffer;
		c>>buffer;
		if(buffer!="")
			a->insert(buffer);
		i++;
	}


	c.close();
}

int main(){
	string r103="r103/Edges.txt";
	string last="new/Edges.txt";
	set<string> r103Set;
	set<string> lastSet;
	load(r103,&r103Set);
	load(last,&lastSet);
	for(set<string>::iterator i=r103Set.begin();i!=r103Set.end();i++){
		if(lastSet.count(*i)==0){
			cout<<*i<<" not in new?"<<endl;
		}
	}
	for(set<string>::iterator i=lastSet.begin();i!=lastSet.end();i++){
		if(r103Set.count(*i)==0){
			cout<<*i<<" not in r103Set?"<<endl;
		}
	}

	return 0;
}
