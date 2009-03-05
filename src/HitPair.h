
#ifndef _HitPair
#define _HitPair

#include"Hit.h"

class HitPair{
	Hit m_left;
	Hit m_right;
public:
	HitPair(Hit*left,Hit*right);
	HitPair();
	bool valid();
	void show();
	Hit*getLeft();
	Hit*getRight();
};



#endif
