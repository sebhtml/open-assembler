
#include"HitPair.h"
#include<iostream>
using namespace std;

Hit*HitPair::getLeft(){
	return &m_left;
}

Hit*HitPair::getRight(){
	return &m_right;
}

HitPair::HitPair(Hit*left,Hit*right){
	m_left=*left;
	m_right=*right;
}

HitPair::HitPair(){
}

bool HitPair::valid(){
	if(m_left.getContigPosition()==0&&m_right.getContigPosition()==0)
		return false;
	if(//m_left.getContigStrand()==m_right.getContigStrand()&&
		((m_left.getContigPosition()==0&&m_right.getContigPosition()==0)||
		(m_left.getContigPosition()!=0&&m_right.getContigPosition()!=0)))
		return false;

	return true;
}

void HitPair::show(){
	cout<<"Hitpair"<<endl;
	m_left.show();
	m_right.show();
}
