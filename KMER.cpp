#include <cassert>
#include <iostream>

#include "KMER.h"
#include "GLOBAL.h"

using namespace std;

KMER::KMER()
{
	for(int i=0;i<ULL_SIZE;++i)
		U64[i] = 0ULL;
}
KMER::KMER(int i)
{
	U64[0] = i;
}
KMER::KMER(const KMER& obj)
{
	for(int i=0;i<ULL_SIZE;++i)
		U64[i] = obj.U64[i];
}

KMER::KMER(const unsigned long long *t)
{
	for(int i=0;i<ULL_SIZE;++i)
		U64[i] = t[i];
}

KMER::KMER(const string& str, bool isRev, bool is_K_minus_1_mer)
{
	for(int i=0;i<ULL_SIZE;++i)
		U64[i] = 0ULL;
	int uend;
	if(is_K_minus_1_mer)
		uend = str.length()<(PARA_KMER-1)?str.length():PARA_KMER-1;
	else
		uend = str.length()< PARA_KMER ? str.length():PARA_KMER;

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
			U64[i>>5] <<= 2; 
			U64[i>>5] |= b;
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
			U64[(uend-1-i)>>5] <<= 2; 
			U64[(uend-1-i)>>5] |= b;
		}
	}
}

void KMER::shiftInLow(unsigned base)
{
	assert(base>=0 && base<4);
	(*this) <<= 2;
	U64[0] |= base;
	cutBroken();
}
void KMER::shiftInHigh(unsigned base)
{
	assert(base>=0 && base<4);
	(*this) >>= 2;
	U64[ULL_SIZE-1] |= ((unsigned long long)base)<<(((PARA_KMER-1)%32)*2);
}

void KMER::cutBroken()
{
	if(BROKEN_U_ID >= 0)
		U64[BROKEN_U_ID] &=  MASK;
}

void KMER::setzero()
{
	for(int i=0;i<ULL_SIZE;++i)
		U64[i] = 0ULL;
}

bool KMER::operator < (const KMER& obj) const
{
	for(int i=ULL_SIZE-1;i>=0;--i)
	{
		if(U64[i] < obj.U64[i])
			return true;
		else if(U64[i] > obj.U64[i])
			return false;
	}
	return false;
}

bool KMER::operator == (const KMER& obj) const
{
	for(int i=0;i<ULL_SIZE;++i)
	{
		if(U64[i]!=obj.U64[i])
			return false;
	}
	return true;
}

inline bool KMER::operator != (const KMER& obj) const
{
	for(int i=0;i<ULL_SIZE;++i)
	{
		if(U64[i]!=obj.U64[i])
			return true;
	}
	return false;
}

KMER KMER::operator << (const int obj) const
{
	const int sid = obj/64;
	const int sft = obj%64;
	unsigned long long result[ULL_SIZE];
	for(int i=ULL_SIZE-1;i>=sid+1;--i)
	{
		result[i] = U64[i-sid] << sft;
		if(sft)
		result[i]|= U64[i-sid-1] >> (64-sft);
	}
	result[sid] = U64[0] << sft;
	for(int i=sid-1;i>=0;--i)
		result[i] = 0ULL;
	return KMER(result);
}

KMER& KMER::operator <<= (const int obj) 
{
	const int sid = obj/64;
	const int sft = obj%64;
	for(int i=ULL_SIZE-1;i>=sid+1;--i)
	{
		U64[i] = U64[i-sid] << sft;
		if(sft)
		U64[i]|= U64[i-sid-1] >> (64-sft);
	}
	U64[sid] = U64[0] << sft;
	for(int i=sid-1;i>=0;--i)
		U64[i] = 0ULL;
	return (*this);
}

int KMER::non0base()
{
	int result = 0;
	for(int i=0;i<ULL_SIZE;++i)
	{
		unsigned long long v=((U64[i]>>1)|(U64[i]))&0x5555555555555555;
		for(;v;++result)
			v &= v-1;
	}
	return result;
}

unsigned KMER::hash()const
{
	unsigned ans = U64[0] % HASHSIZE;
	for(int i=1;i<ULL_SIZE;++i)
	{
		ans += (U64[i]%HASHSIZE)*U64HASH[i];
		ans %= HASHSIZE;
	}
	return ans;
}

KMER KMER::operator >> (const int obj) const
{
	const int sid = obj/64;
	const int sft = obj%64;
	unsigned long long result[ULL_SIZE];
	for(int i=0;i<ULL_SIZE-1-sid;++i)
	{
		result[i] = U64[i+sid] >> sft;
		if(sft)
		result[i]|= U64[i+sid+1] << (64-sft);
	}
	result[ULL_SIZE-1-sid] = U64[ULL_SIZE-1] >> sft;
	for(int i=ULL_SIZE-sid;i<ULL_SIZE;++i)
		result[i] = 0ULL;
	return KMER(result);
}

KMER& KMER::operator >>= (const int obj) 
{
	const int sid = obj/64;
	const int sft = obj%64;
	for(int i=0;i<ULL_SIZE-1-sid;++i)
	{
		U64[i] = U64[i+sid] >> sft;
		if(sft)
		U64[i]|= U64[i+sid+1] << (64-sft);
	}
	U64[ULL_SIZE-1-sid] = U64[ULL_SIZE-1] >> sft;
	for(int i=ULL_SIZE-sid;i<ULL_SIZE;++i)
		U64[i] = 0ULL;
	return (*this);
}

KMER& KMER::operator |= (const int obj)
{
	U64[0] |= obj;
	return *this;
}

KMER& KMER::operator |= (const KMER& obj)
{
	for(int i=0;i<ULL_SIZE;++i)
		U64[i] |= obj.U64[i];
	return *this;
}

KMER KMER::operator | (const KMER& obj)const
{
	KMER result;
	for(int i=0;i<ULL_SIZE;++i)
		result.U64[i] = this->U64[i] | obj.U64[i];
	return result;
}

KMER KMER::operator ^ (const KMER& obj)const
{
	KMER result;
	for(int i=0;i<ULL_SIZE;++i)
		result.U64[i] = this->U64[i] ^ obj.U64[i];
	return result;
}

KMER KMER::operator & (const KMER& obj)const
{
	KMER result;
	for(int i=0;i<ULL_SIZE;++i)
		result.U64[i] = this->U64[i] & obj.U64[i];
	return result;
}

KMER& KMER::operator ^= (const KMER& obj)
{
	for(int i=0;i<ULL_SIZE;++i)
		U64[i] ^= obj.U64[i];
	return *this;
}
KMER& KMER::operator &= (const KMER& obj)
{
	for(int i=0;i<ULL_SIZE;++i)
		U64[i] &= obj.U64[i];
	return *this;
}

KMER KMER::operator~() const
{
	KMER result;
	for(int i=0;i<ULL_SIZE;++i)
		result.U64[i] = ~(this->U64[i]);
	return result;
}

KMER& KMER::operator=(const KMER&obj)
{
	for(int i=0;i<ULL_SIZE;++i)
		U64[i] = obj.U64[i];
}

ostream& operator<<(ostream& os,const KMER &obj)
{
	for(int i=obj.ULL_SIZE-1;i>0;--i)
		os << hex << obj.U64[i] << '-';
	os << hex << obj.U64[0];
	return os;
}
