/*
	dna: De Novo Assembler
    Copyright (C) 2008, 2009 Sébastien Boisvert
	$Id$

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

#ifndef _CustomMap
#define _CustomMap
#include<cmath>
#include<stdint.h>
#include<stdlib.h>
#include<iostream>
using namespace std;

template <class VALUE>
class Entry{
public:
	Entry(uint64_t key,VALUE value){
		m_key=key;
		m_value=value;
		m_next=NULL;
	}
	uint64_t m_key;
	VALUE m_value;
	Entry*m_next;
	~Entry(){
		m_next=NULL;
	}
};

template <class VALUE>
class CustomMap{
public:

	class iterator{
		uint32_t m_bucket;
		Entry<VALUE>*m_current;
		bool m_end;
		CustomMap*m_custom_map;
	public:
		bool operator!=(iterator b){
			if(m_end==true){
				//cout<<"THE END "<<endl;
				return m_end!=b.m_end;
			}
			return m_current!=b.m_current;
		}
		iterator(uint32_t bucket,Entry<VALUE>*current, bool isEnd,CustomMap<VALUE>*customMap){
			m_bucket=bucket;
			m_current=current;
			m_end=isEnd;
			m_custom_map=customMap;
		}
		/**
 			precondition: m_current is not NULL or m_end is true
			postcondition: m_current is not NULL or m_end is true
		*/

		void operator++(int i){
			m_current=m_current->m_next; // take the next in the current bucket
			if(m_current==NULL){ // if the bucket is done
				while(m_current==NULL){ // find another bucket containing data
					if(m_bucket==m_custom_map->m_table_size-1){ //if buckets are depleted, you are done
						m_end=true;
						return;
					}

					m_bucket++; // else, try the next bucket.
					m_current=m_custom_map->m_buckets[m_bucket];
				}
			}
		}
		uint64_t first(){
			return m_current->m_key;
		}
		VALUE second(){
			return m_current->m_value;
		}
	};

	Entry<VALUE>**m_buckets;
	uint64_t m_size;
	uint64_t m_table_size;
	VALUE m_defaultValue;
	// http://www.concentric.net/~Ttwang/tech/inthash.htm
	uint32_t hash6432shift(uint64_t key)
	{
  		key = (~key) + (key << 18); // key = (key << 18) - key - 1;
  		key = key ^ (key >> 31);
  		key = key * 21; // key = (key + (key << 2)) + (key << 4);
  		key = key ^ (key >> 11);
  		key = key + (key << 6);
  		key = key ^ (key >> 22);
		//cout<<t<<" "<<(uint32_t) key<<endl;
		uint32_t value=key%m_table_size;
  		return value;
	}


	uint32_t hash6432shift2(uint64_t key){
		return ((uint32_t)key)%m_table_size;
		uint32_t a=0;
		int nBitsOutput=log(m_table_size)/log(2);
		for(int i=0;i<64;i+=nBitsOutput){
			uint64_t key2=(key<<(64-i-nBitsOutput));
			key2=(key2>>(64-nBitsOutput));
			a=a | key2;
		}
		if(a>=m_table_size){
			cout<<"Error hash"<<endl;
			exit(0);
		}
		return a;
	}

	iterator begin(){
		uint32_t bucket=0;
		Entry<VALUE>*current=m_buckets[bucket];
		while(current==NULL){
			bucket++;
			current=m_buckets[bucket];
		}
		iterator i(bucket,current,false,this);
		//cout<<"Begin "<<current->m_key<<endl;
		return i;
	}
	iterator end(){
		//cout<<"end"<<endl;
		iterator i(0,NULL,true,this);
		return i;
	}
	CustomMap(uint64_t table_size){
		m_table_size=table_size;
		//cout<<m_table_size<<endl;
		m_buckets=new Entry<VALUE>*[m_table_size];
		//(*m_cout)<<"[m_table_size] "<<m_table_size<<endl;
		//(*m_cout)<<"Space: "<<m_table_size*sizeof(uint64_t)/1024/1024<<" MB"<<endl;
		for(uint64_t i=0;i<m_table_size;i++){
			if(i%1000==0){
			}
			m_buckets[i]=NULL;
		}
		m_size=0;
	}
	void add(uint64_t key,VALUE value){
		uint32_t hashKey=hash6432shift(key);
		//cout<<"Adding "<<key<<endl;
		if(m_buckets[hashKey]==NULL){
			m_buckets[hashKey]=new Entry<VALUE>(key,value);
			m_size++;
			return;
		}
		Entry<VALUE>*entry=m_buckets[hashKey];
		while(entry!=NULL){
			if(entry->m_key==key)
				return;
			entry=entry->m_next;
		}
		entry=m_buckets[hashKey];
		while(entry->m_next!=NULL)
			entry=entry->m_next;
		entry->m_next=new Entry<VALUE>(key,value);
		m_size++;
	}
	VALUE&get(uint64_t key){
		uint32_t hashKey=hash6432shift(key);
		Entry<VALUE>*entry=m_buckets[hashKey];
		while(entry!=NULL){
			if(entry->m_key==key)
				return entry->m_value;
			entry=entry->m_next;
		}
		
		return m_defaultValue;
	}
	void set(uint64_t key,VALUE value){
		uint32_t hashKey=hash6432shift(key);
		Entry<VALUE>*entry=m_buckets[hashKey];
		while(entry!=NULL){
			if(entry->m_key==key)
				entry->m_value=value;
			entry=entry->m_next;
		}
	}
	uint64_t size(){
		return m_size;
	}
	bool find(uint64_t key){
		uint32_t hashKey=hash6432shift(key);
		Entry<VALUE>*entry=m_buckets[hashKey];
		while(entry!=NULL){
			if(entry->m_key==key)
				return true;
			entry=entry->m_next;
		}
		return false;
	}
	void clear(){
		if(m_buckets==NULL)
			return;
		for(uint64_t i=0;i<m_table_size;i++){
			Entry<VALUE>*c=m_buckets[i];
			while(c!=NULL){
				Entry<VALUE>*d=c;
				c=c->m_next;
				delete d;
			}
		}
		delete[] m_buckets;
		m_buckets=NULL;
	}
	~CustomMap(){
		clear();
	}
	int buckets(){
		return m_table_size;
	}
	Entry<VALUE>*bucketAt(int i){
		return m_buckets[i];
	}
};

#endif
