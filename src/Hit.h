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



#ifndef _Hit
#define _Hit

class Hit{
	int m_contigNumber;
	int m_contigPosition;
	char m_contigStrand;
	
	int m_readNumber;
	int m_readPosition;
	char m_readStrand;
	char m_L_or_R;
public:
	Hit(int contigNumber,int contigPosition,char contigStrand,
		int readNumber,int readPosition,char readStrand);
	Hit();
	void setL_or_R(char l);
	void show();
	int getContigNumber();
	int getContigPosition();
	char getContigStrand();
	int getReadPosition();
	int getReadNumber();
	char getReadStrand();
	char side();
};

#endif
