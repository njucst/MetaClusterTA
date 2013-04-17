/*
 * last fixed: 2012.01.05.
 * by Wang Yi.
 * */
#ifndef MCH_HYBRID_ULLN_H_

#define MCH_HYBRID_ULLN_H_
#include <iostream>
#include "PARAMETER.h"
using namespace std;
class ULLN
{
public:
	static const int ULL_SIZE = (PARA_READ-1)/32+1;
//	int ULL_SIZE;
	unsigned long long U64[ULL_SIZE];

	ULLN();
	ULLN(const ULLN& obj);
	ULLN(const unsigned long long *t);
	ULLN(const string& str, bool isRev = false);
	void setzero();
	int non0base();//not non 0 bits

	ULLN& keeptopk(unsigned k);//keep top k bases(i.e. 2k bits), set tail to 0;

	bool operator < (const ULLN& obj) const;
	bool operator == (const ULLN& obj) const;
	bool operator != (const ULLN& obj) const;
	ULLN& operator <<= (const int obj);
	ULLN  operator << (const int obj) const;
	ULLN& operator >>= (const int obj);
	ULLN  operator >> (const int obj) const;
	ULLN& operator |= (const int obj);
	ULLN& operator |= (const ULLN& obj);
	ULLN operator | (const ULLN& obj)const;
	ULLN& operator &= (const ULLN& obj);
	ULLN operator & (const ULLN& obj)const;
	ULLN& operator ^= (const ULLN& obj);
	ULLN  operator ^ (const ULLN& obj) const;
	ULLN operator~() const;

	ULLN& operator=(const ULLN&obj);
 //   void setbase(int position,unsigned long long base);
	friend ostream& operator<<(ostream& os,const ULLN & obj);
};
#endif
