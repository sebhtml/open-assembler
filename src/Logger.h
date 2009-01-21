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
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
using namespace std;

class Logger{
	string m_outputDirectory;
	ofstream*m_stream;
	ostringstream m_buffer;
public:
	Logger();
	void setOutput(string file);
	Logger&operator<<(const string a);
	Logger&operator<<(const int a);
	Logger&operator<<(const double a);
	Logger&operator<<(const unsigned long int a);
	Logger&operator<<(const char*a);
	Logger&operator<<(const char a);
	~Logger();
};

#endif
