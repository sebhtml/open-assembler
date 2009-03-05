
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
};

#endif
