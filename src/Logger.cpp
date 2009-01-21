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


#include<cstdlib>
#include<stdlib.h>
#include"Logger.h"
#include<string>
#include<fstream>
#include<iostream>
using namespace std;
Logger::Logger(){
	m_stream=NULL;
}

void Logger::setOutput(const string output){
	m_stream=new ofstream(output.c_str());
}

Logger&Logger::operator<<(const string a){
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

// TODO: template operator<<
Logger&Logger::operator<<(const char a){
	ostringstream b;
	b<<a;
	return operator<<(b.str());
}

Logger&Logger::operator<<(const double a){
	ostringstream b;
	b<<a;
	return operator<<(b.str());
}

Logger&Logger::operator<<(const unsigned long int a){
	ostringstream b;
	b<<a;
	return operator<<(b.str());
}

Logger&Logger::operator<<(const int a){
	ostringstream b;
	b<<a;
	return operator<<(b.str());
}

Logger&Logger::operator<<(const char* a){
	ostringstream b;
	b<<a;
	return operator<<(b.str());
}

Logger::~Logger(){
	m_stream->close();
	if(m_stream!=NULL){
		delete m_stream;
		m_stream=NULL;
	}
}
