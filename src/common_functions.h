
#ifndef _common_functions
#define _common_functions

/*
 *  unsigned long int
 *			for 64-bit machines, allows wordSize<=31
 *  unsigned int
 *  			for 32-bit machines, allows wordSize<=15
 */

//#define VERTEX_TYPE unsigned long int
//#define VERTEX_TYPE uint32_t
#define VERTEX_TYPE uint64_t

#include<string>
using namespace std;

string reverseComplement(string a);
VERTEX_TYPE wordId(const char*a);
string idToWord(VERTEX_TYPE i,int wordSize);
char complement(char a);
void CommonHeader(ostream*out);


char getLastSymbol(VERTEX_TYPE i,int w);

#endif
