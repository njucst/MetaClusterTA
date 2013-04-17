#include <cassert>
#include <iostream>

#include "ShortKMER.h"

using namespace std;

ShortKMER::ShortKMER()
{
	U64 = 0ULL;
}

ShortKMER::ShortKMER(const ShortKMER& obj)
{
	U64 = obj.U64;
}

ShortKMER::ShortKMER(const unsigned long long t)
{
	U64 = t;
}

ShortKMER::ShortKMER(const string& str, bool isRev, bool is_K_minus_1_mer)
{
	U64 = 0ULL;
	int uend;
	if(is_K_minus_1_mer)
		uend = str.length()<(PARA_SKMER-1)?str.length():PARA_SKMER-1;
	else
		uend = str.length()< PARA_SKMER ? str.length():PARA_SKMER;

	if(isRev)
	{
		for(int i=uend-1;i>=0;--i)
		{
			unsigned b = 0;
			if(str[i]=='C' || str[i]=='c')
				b = 2;
			else if(str[i]=='G' || str[i]=='g')
				b = 1;
			else if(str[i]=='A' || str[i]=='a')
				b = 3;
			U64 <<= 2; 
			U64 |= b;
		}
	}
	else
	{
		for(int i=0;i<uend;++i)
		{
			unsigned b = 0;
			if(str[i]=='C' || str[i]=='c')
				b = 1;
			else if(str[i]=='G' || str[i]=='g')
				b = 2;
			else if(str[i]=='T' || str[i]=='t')
				b = 3;
			U64 <<= 2; 
			U64 |= b;
		}
	}
}

void ShortKMER::shiftInLow(unsigned base)
{
	assert(base>=0 && base<4);
	U64 <<= 2;
	U64 |= base;
	cutBroken();
}
void ShortKMER::shiftInHigh(unsigned base)
{
	assert(base>=0 && base<4);
	U64 >>= 2;
	U64 |= ((unsigned long long)base)<<(((PARA_SKMER-1)%32)*2);
}

void ShortKMER::cutBroken()
{
	if(PARA_SKMER != 32)
		U64 &=  MASK;
}

void ShortKMER::setzero()
{
	U64 = 0ULL;
}

bool ShortKMER::operator < (const ShortKMER& obj) const
{
	return U64 < obj.U64;
}

bool ShortKMER::operator == (const ShortKMER& obj) const
{
	return U64 == obj.U64;
}

inline bool ShortKMER::operator != (const ShortKMER& obj) const
{
	return U64 != obj.U64;
}

ShortKMER ShortKMER::operator << (const int obj) const
{
	return ShortKMER(U64 << obj);
}

ShortKMER& ShortKMER::operator <<= (const int obj) 
{
	U64 <<= obj;
	return (*this);
}

int ShortKMER::non0base()
{
	int result = 0;
	unsigned long long v=((U64>>1)|(U64))&0x5555555555555555;
	for(;v;++result)
		v &= v-1;
	return result;
}

unsigned ShortKMER::hash()const
{
	//return U64 & DBHASHMASK;
	//return U64 % DBHASHSIZE;
	return ((U64*1299709+104729)%323780508946331ULL) % DBHASHSIZE;
}

ShortKMER ShortKMER::operator >> (const int obj) const
{
	return ShortKMER(U64>>obj);
}

ShortKMER& ShortKMER::operator >>= (const int obj) 
{
	U64 >>= obj;
	return (*this);
}

ShortKMER& ShortKMER::operator |= (const unsigned long long obj)
{
	U64 |= obj;
	return *this;
}

ShortKMER& ShortKMER::operator |= (const ShortKMER& obj)
{
	U64 |= obj.U64;
	return *this;
}

ShortKMER ShortKMER::operator | (const ShortKMER& obj)const
{
	return ShortKMER(U64 | obj.U64);
}

ShortKMER ShortKMER::operator ^ (const ShortKMER& obj)const
{
	return ShortKMER(U64 ^ obj.U64);
}

ShortKMER ShortKMER::operator & (const ShortKMER& obj)const
{
	return ShortKMER(U64 & obj.U64);
}

ShortKMER& ShortKMER::operator ^= (const ShortKMER& obj)
{
	U64 ^= obj.U64;
	return *this;
}
ShortKMER& ShortKMER::operator &= (const ShortKMER& obj)
{
	U64 &= obj.U64;
	return *this;
}

ShortKMER ShortKMER::operator~() const
{
	return ShortKMER(~U64);
}

ShortKMER& ShortKMER::operator=(const ShortKMER&obj)
{
	U64 = obj.U64;
}

ostream& operator<<(ostream& os,const ShortKMER &obj)
{
	os << hex << obj.U64;
	return os;
}
