/*
	454 De Novo Assembler
    Copyright (C) 2008 Sébastien Boisvert

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

#ifndef _Logger
#define _Logger

#include<sstream>
#include<cstdlib>
#include<stdlib.h>
#include<string>
#include<fstream>
#include<iostream>
using namespace std;

#define endl '\n'


class Logger{
	string m_outputDirectory;
	ofstream*m_stream;
	ostringstream m_buffer;
public:
	Logger(){
		m_stream=NULL;
	}
	void setOutput(string file){
		m_stream=new ofstream(file.c_str());
	}
	template<class A>
	Logger&operator<<(const A a){
		cout<<a;
		cout.flush();
		if(m_stream==NULL){
			m_buffer<<a;
		}else{
			if(m_buffer.str().length()>0){
				(*m_stream)<<m_buffer.str();
				m_buffer.str("");
			}
			(*m_stream)<<(a);
			m_stream->flush();
		}
		return *this;
	}
	~Logger(){
		m_stream->close();
		if(m_stream!=NULL){
			delete m_stream;
			m_stream=NULL;
		}

	}
};

#endif
