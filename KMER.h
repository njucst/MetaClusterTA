/*
 * last fixed: 2012.11.11.
 * by Wang Yi.
 * */
#ifndef MCH_HYBRID_KMER_H_

#define MCH_HYBRID_KMER_H_
#include <iostream>
#include "PARAMETER.h"
using namespace std;
class KMER
{
public:
	static const int ULL_SIZE = (PARA_KMER-1)/32+1;
	static const int BROKEN_U_ID = PARA_KMER%32==0?-1:(ULL_SIZE-1);
	static const unsigned long long MASK = (1ULL<<((PARA_KMER%32)*2))-1;
	unsigned long long U64[ULL_SIZE];

	KMER();
	KMER(int i);
	KMER(const KMER& obj);
	KMER(const unsigned long long *t);
	KMER(const string& str, bool isRev = false,bool is_K_minus_1_mer = false);

	void shiftInLow(unsigned base);
	void shiftInHigh(unsigned base);
	void cutBroken();
	void setzero();
	int non0base();//not non 0 bits
	unsigned hash() const;

	bool operator < (const KMER& obj) const;
	bool operator == (const KMER& obj) const;
	bool operator != (const KMER& obj) const;
	KMER& operator <<= (const int obj);
	KMER  operator << (const int obj) const;
	KMER& operator >>= (const int obj);
	KMER  operator >> (const int obj) const;
	KMER& operator |= (const int obj);
	KMER& operator |= (const KMER& obj);
	KMER operator | (const KMER& obj)const;
	KMER& operator &= (const KMER& obj);
	KMER operator & (const KMER& obj)const;
	KMER& operator ^= (const KMER& obj);
	KMER  operator ^ (const KMER& obj) const;
	KMER operator~() const;

	KMER& operator=(const KMER&obj);
 //   void setbase(int position,unsigned long long base);
	friend ostream& operator<<(ostream& os,const KMER & obj);
};
#endif
