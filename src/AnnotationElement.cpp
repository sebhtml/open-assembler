#include"AnnotationElement.h"
#include<iostream>
using namespace std;

AnnotationElement::AnnotationElement(uint32_t readId,uint16_t readPosition,uint8_t readStrand){
	m_readId=readId;
	m_readPosition=readPosition;
	m_readStrand=readStrand;
}

uint32_t AnnotationElement::getReadId(){
	return m_readId;
}

uint16_t AnnotationElement::getReadPosition(){
	return m_readPosition;
}

uint8_t AnnotationElement::getReadStrand(){
	return m_readStrand;
}

void AnnotationElement::print(){
	cout<<" ("<<getReadId()<<" "<<getReadStrand()<<" "<<getReadPosition()<<")";
}
