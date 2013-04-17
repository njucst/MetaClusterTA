/*
 * last fixed: 2012.11.11.
 * by Wang Yi.
 * */
#ifndef MCH_HYBRID_ShortKMER_H_

#define MCH_HYBRID_ShortKMER_H_
#include <iostream>
#include "PARAMETER.h"
using namespace std;
class ShortKMER
{
public:
	static const unsigned PARA_SKMER = 28;//short kmer length
	//1: static const int ULL_SIZE = (PARA_SKmer-1)/32+1;
	//0: static const int BROKEN_U_ID = PARA_SKmer%32==0?-1:(ULL_SIZE-1);
//	static const unsigned long long MASK = (1ULL<<(PARA_SKMER*2))-1;
	static const unsigned long long MASK = PARA_SKMER<32?((1ULL<<(PARA_SKMER*2))-1):(~(0ULL));
	unsigned long long U64;

	ShortKMER();
	ShortKMER(const ShortKMER& obj);
	ShortKMER(const unsigned long long t);
	ShortKMER(const string& str, bool isRev = false,bool is_K_minus_1_mer = false);

	void shiftInLow(unsigned base);
	void shiftInHigh(unsigned base);
	void cutBroken();
	void setzero();
	int non0base();//not non 0 bits
	unsigned hash() const;

	bool operator < (const ShortKMER& obj) const;
	bool operator == (const ShortKMER& obj) const;
	bool operator != (const ShortKMER& obj) const;
	ShortKMER& operator <<= (const int obj);
	ShortKMER  operator << (const int obj) const;
	ShortKMER& operator >>= (const int obj);
	ShortKMER  operator >> (const int obj) const;
	ShortKMER& operator |= (const unsigned long long obj);
	ShortKMER& operator |= (const ShortKMER& obj);
	ShortKMER operator | (const ShortKMER& obj)const;
	ShortKMER& operator &= (const ShortKMER& obj);
	ShortKMER operator & (const ShortKMER& obj)const;
	ShortKMER& operator ^= (const ShortKMER& obj);
	ShortKMER  operator ^ (const ShortKMER& obj) const;
	ShortKMER operator~() const;

	ShortKMER& operator=(const ShortKMER&obj);
	friend ostream& operator<<(ostream& os,const ShortKMER & obj);
};
#endif
