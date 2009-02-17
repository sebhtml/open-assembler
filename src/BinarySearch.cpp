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

#include"BinarySearch.h"

int BinarySearch(vector<VERTEX_TYPE>*a,VERTEX_TYPE b){
	//cout<<"BinarySearch"<<endl;
	int start=0;
	int end=a->size()-1;
	while(start<=end){
		int mid=(start+end)/2;
		//cout<<start<<" "<<middle<<" "<<end<<endl;
		if(b<a->at(mid)){
			//cout<<"lower"<<endl;
			end=mid-1;
		}else if(a->at(mid)<b){	
			//cout<<"greater"<<endl;
			start=mid+1;
		}else{
			return mid;
		}
	}
	return -1;
}


