#include <iostream>

#include "ULLN.h"
using namespace std;

ULLN::ULLN()
{
	for(int i=0;i<ULL_SIZE;++i)
		U64[i] = 0ULL;
}

ULLN::ULLN(const ULLN& obj)
{
	for(int i=0;i<ULL_SIZE;++i)
		U64[i] = obj.U64[i];
}

ULLN::ULLN(const unsigned long long *t)
{
	for(int i=0;i<ULL_SIZE;++i)
		U64[i] = t[i];
}

ULLN::ULLN(const string& str, bool isRev)
{
	for(int i=0;i<ULL_SIZE;++i)
		U64[i] = 0ULL;
	int lastbase = str.length()-1;
	if(isRev)//reverse compliment
	{
		for(int i=str.length()-1;i>=0;--i)
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
		for(int i=0;i<str.length();++i)
		{
			unsigned b = 0;
			if(str[i]=='C' || str[i]=='c')
				b = 1;
			else if(str[i]=='G' || str[i]=='g')
				b = 2;
			else if(str[i]=='T' || str[i]=='t')
				b = 3;
			U64[(lastbase-i)>>5] <<= 2; 
			U64[(lastbase-i)>>5] |= b;
		}
	}
}

void ULLN::setzero()
{
	for(int i=0;i<ULL_SIZE;++i)
		U64[i] = 0ULL;
}

ULLN& ULLN::keeptopk(unsigned k)//keep top k bases(i.e. 2k bits), set tail to 0;
{
	if(k >= ULL_SIZE*32)return(*this);
	if(k&0x1f)
		U64[k>>5] &= ((1ULL<<((k & 0x1f) << 1))-1);
	else
		U64[k>>5] = 0ULL;
	for(int i=(k/32+1);i<ULL_SIZE;++i)
		U64[i] = 0ULL;
	return (*this);
}

bool ULLN::operator < (const ULLN& obj) const
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

inline bool ULLN::operator == (const ULLN& obj) const
{
	for(int i=0;i<ULL_SIZE;++i)
	{
		if(U64[i]!=obj.U64[i])
			return false;
	}
	return true;
}

inline bool ULLN::operator != (const ULLN& obj) const
{
	for(int i=0;i<ULL_SIZE;++i)
	{
		if(U64[i]!=obj.U64[i])
			return true;
	}
	return false;
}

ULLN ULLN::operator << (const int obj) const
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
	return ULLN(result);
}

ULLN& ULLN::operator <<= (const int obj) 
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

inline unsigned diffsearch(unsigned long long a)
{
	unsigned long  pos = 0;
	__asm__("bsrq %1,%0\n\t"
					: "+r" (pos)
					: "rm" (a));
	return pos;
}

int ULLN::non0base()
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

ULLN ULLN::operator >> (const int obj) const
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
	return ULLN(result);
}

ULLN& ULLN::operator >>= (const int obj) 
{
	const int sid = obj/64;
	const int sft = obj%64;
	for(int i=0;i<ULL_SIZE-1-sid;++i)
	{
		U64[i] = U64[i+sid] >> sft;
		if(sft)	U64[i]|= U64[i+sid+1] << (64-sft);
	}
	U64[ULL_SIZE-1-sid] = U64[ULL_SIZE-1] >> sft;
	for(int i=ULL_SIZE-sid;i<ULL_SIZE;++i)
		U64[i] = 0ULL;
	return (*this);
}

ULLN& ULLN::operator |= (const int obj)
{
	U64[0] |= obj;
	return *this;
}

ULLN& ULLN::operator |= (const ULLN& obj)
{
	for(int i=0;i<ULL_SIZE;++i)
		U64[i] |= obj.U64[i];
	return *this;
}

ULLN ULLN::operator | (const ULLN& obj)const
{
	ULLN result;
	for(int i=0;i<ULL_SIZE;++i)
		result.U64[i] = this->U64[i] | obj.U64[i];
	return result;
}

ULLN ULLN::operator ^ (const ULLN& obj)const
{
	ULLN result;
	for(int i=0;i<ULL_SIZE;++i)
		result.U64[i] = this->U64[i] ^ obj.U64[i];
	return result;
}

ULLN ULLN::operator & (const ULLN& obj)const
{
	ULLN result;
	for(int i=0;i<ULL_SIZE;++i)
		result.U64[i] = this->U64[i] & obj.U64[i];
	return result;
}

ULLN& ULLN::operator ^= (const ULLN& obj)
{
	for(int i=0;i<ULL_SIZE;++i)
		U64[i] ^= obj.U64[i];
	return *this;
}
ULLN& ULLN::operator &= (const ULLN& obj)
{
	for(int i=0;i<ULL_SIZE;++i)
		U64[i] &= obj.U64[i];
	return *this;
}

ULLN ULLN::operator~() const
{
	ULLN result;
	for(int i=0;i<ULL_SIZE;++i)
		result.U64[i] = ~(this->U64[i]);
	return result;
}

ULLN& ULLN::operator=(const ULLN&obj)
{
	for(int i=0;i<ULL_SIZE;++i)
		U64[i] = obj.U64[i];
}

ostream& operator<<(ostream& os,const ULLN &obj)
{
	for(int i=obj.ULL_SIZE-1;i>0;--i)
		os << hex << obj.U64[i] << '-';
	os << hex << obj.U64[0];
	return os;
}

//position is [1,75];
/*void ULLN::setbase(int position,unsigned long long base)
{
	if(position > 97)
	{
		position -= 97;
		unsigned long long mask = 3ULL << (position*2);
		high &= (~mask);
		high |= (base << (position*2));
	}
	else if(position == 97)
	{
		high &= (~3ULL);
		high |= base;
	}
	else if(position > 65)
	{
		position -= 65;
		unsigned long long mask = 3ULL << (position*2);
		high &= (~mask);
		high |= (base << (position*2));
	}
	else if(position == 65)
	{
		high &= (~3ULL);
		high |= base;
	}
	else if(position > 33)
	{
		position -= 33;
		unsigned long long mask = 3ULL << (position*2);
		mid &= (~mask);
		mid |= (base << position*2);
	}
	else if(position == 33)
	{
		mid &= (~3ULL);
		mid |= base;
	}
	else if(position > 1)
	{
		--position;
		unsigned long long mask = 3ULL << (position*2);
		low &= (~mask);
		low |= (base << position*2);
	}
	else
	{
		low &= (~3ULL);
		low |= base;
	}
}*/
