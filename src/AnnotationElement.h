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

#ifndef _AnnotationElement
#define _AnnotationElement

//#define POSITION_TYPE uint32_t
#define POSITION_TYPE uint16_t

#include<stdint.h>

class AnnotationElement{
	uint32_t m_readId;
	POSITION_TYPE m_readPosition;
	uint8_t m_readStrand;  // TODO: replace  this with 2 objects instead
public:
	AnnotationElement(uint32_t readId,uint16_t readPosition,uint8_t readStrand);
	uint32_t getReadId();
	uint16_t getReadPosition();
	uint8_t getReadStrand();
	void print();
};


#endif

