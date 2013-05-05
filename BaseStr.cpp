#include <cassert>
#include <cmath>
#include <iostream>

#include "BaseStr.h"
using namespace std;

/*
BaseStr::BaseStr()
{
	ULL_SIZE = DEFAULT_S;
	size = 0
	U64 = new unsigned long long[ULL_SIZE];
	for(int i=0;i<ULL_SIZE;++i)
		U64[i] = 0ULL;
}
BaseStr::BaseStr(int capa)
{
	ULL_SIZE = capa;
	U64 = new unsigned long long[ULL_SIZE];
	for(int i=0;i<ULL_SIZE;++i)
		U64[i] = 0ULL;
}
*/

BaseStr::BaseStr(const unsigned long long *t, int usize, int size_)
{
	ULL_SIZE = usize;
	length = size_;
	U64 = new unsigned long long[ULL_SIZE];
	for(int i=0;i<ULL_SIZE;++i)
		U64[i] = t[i];
}

BaseStr::BaseStr(const BaseStr& obj)
{
	ULL_SIZE = obj.ULL_SIZE;
	length = obj.length;
	U64 = new unsigned long long[ULL_SIZE];
	for(int i=0;i<ULL_SIZE;++i)
		U64[i] = obj.U64[i];
}

//BaseStr::BaseStr(const string& str_)
BaseStr::BaseStr(const string& str_,bool isRev)
{
	ULL_SIZE = (str_.length()-1)/32+1;
	length = str_.length();
	U64 = new unsigned long long[ULL_SIZE];
	for(int i=0;i<ULL_SIZE;++i)
		U64[i] = 0ULL;
	if(!isRev)
	{
		for(int i=0;i<length;++i)
		{
			unsigned b = 0;
			if(str_[i]=='C' || str_[i]=='c')
				b = 1;
			else if(str_[i]=='G' || str_[i]=='g')
				b = 2;
			else if(str_[i]=='T' || str_[i]=='t')
				b = 3;
			U64[(length-1-i)>>5] <<= 2; 
			U64[(length-1-i)>>5] |= b;
		}
	}
	else
	{
		for(int i=length-1;i>=0;--i)
		{
			unsigned b = 3;
			if(str_[i]=='C' || str_[i]=='c')
				b = 2;
			else if(str_[i]=='G' || str_[i]=='g')
				b = 1;
			else if(str_[i]=='T' || str_[i]=='t')
				b = 0;
			U64[(length-1-i)>>5] <<= 2; 
			U64[(length-1-i)>>5] |= b;
		}
	}
	this->str = str_;
}

bool BaseStr::operator < (const BaseStr& obj) const
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

inline bool BaseStr::operator == (const BaseStr& obj) const
{
	if(length != obj.length)return false;
	for(int i=0;i<ULL_SIZE;++i)
		if(U64[i]!=obj.U64[i])
			return false;
	return true;
}

inline bool BaseStr::operator != (const BaseStr& obj) const
{
	if(length != obj.length)return true;
	for(int i=0;i<ULL_SIZE;++i)
		if(U64[i]!=obj.U64[i])
			return true;
	return false;
}

BaseStr BaseStr::operator << (const int obj) const
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
	return BaseStr(result,ULL_SIZE,min(length+obj/2,32*ULL_SIZE));
}

BaseStr& BaseStr::operator <<= (const int obj) 
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

int BaseStr::non0base()
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

/*void BaseStr::resize(int ulls)
{
	if(ULL_SIZE != ulls && ulls>0)
	{
		delete[]U64;
		U64 = new unsigned long long[ulls];
	}
}*/

BaseStr BaseStr::operator >> (const int obj) const
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
	return BaseStr(result,ULL_SIZE,max(length-obj/2,0));
}

BaseStr& BaseStr::operator >>= (const int obj) 
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

BaseStr& BaseStr::operator |= (const int obj)
{
	U64[0] |= obj;
	return *this;
}

BaseStr& BaseStr::operator |= (const BaseStr& obj)
{
	for(int i=0;i<ULL_SIZE;++i)
		U64[i] |= obj.U64[i];
	return *this;
}

BaseStr BaseStr::operator | (const BaseStr& obj)const
{
	BaseStr result(*this);
	for(int i=0;i<ULL_SIZE;++i)
		result.U64[i] = this->U64[i] | obj.U64[i];
	return result;
}

BaseStr BaseStr::operator ^ (const BaseStr& obj)const
{
	BaseStr result(*this);
	for(int i=0;i<ULL_SIZE;++i)
		result.U64[i] = this->U64[i] ^ obj.U64[i];
	return result;
}

BaseStr BaseStr::operator & (const BaseStr& obj)const
{
	BaseStr result(*this);
	for(int i=0;i<ULL_SIZE;++i)
		result.U64[i] = this->U64[i] & obj.U64[i];
	return result;
}

BaseStr& BaseStr::operator ^= (const BaseStr& obj)
{
	for(int i=0;i<ULL_SIZE;++i)
		U64[i] ^= obj.U64[i];
	return *this;
}
BaseStr& BaseStr::operator &= (const BaseStr& obj)
{
	for(int i=0;i<ULL_SIZE;++i)
		U64[i] &= obj.U64[i];
	return *this;
}

BaseStr BaseStr::operator~() const
{
	BaseStr result(*this);
	for(int i=0;i<ULL_SIZE;++i)
		result.U64[i] = ~(this->U64[i]);
	return result;
}

BaseStr& BaseStr::operator=(const BaseStr& obj)
{
	if(ULL_SIZE != obj.ULL_SIZE)
	{
		ULL_SIZE = obj.ULL_SIZE;
		if(U64!=NULL)
			delete[]U64;
		U64 = new unsigned long long[ULL_SIZE];
	}
	for(int i=0;i<ULL_SIZE;++i)
		U64[i] = obj.U64[i];
}
inline unsigned long long mask(int b)
{
	if(b==64 || b==0)
		return ~0ULL;
	return (1ULL<<b)-1;
}

ULLN BaseStr::subStr(int idx, int len)const// append 0 at the end of ULLN.
{
	len = (len <= (length-idx)) ? len : (length-idx);
	ULLN ans;
	assert(len < 32*ans.ULL_SIZE);
	unsigned long long T64[ans.ULL_SIZE+1];
	for(int i=0;i<ans.ULL_SIZE+1;++i)T64[i]=0ULL;
	int minU = (length-len-idx)>>5, maxU = (length-1-idx)>>5;
	for(int i=minU;i<=maxU;++i)
		T64[i-minU] = U64[i];
	int rft = (1+((length-1-idx)&0x1f))*2;
	maxU -= minU;
	if(rft<64) T64[maxU] &= ((1ULL<<rft)-1);
	int rsft = ((length-len-idx)&0x1f)*2,lsft = 64-rsft;
	if(rsft==0)
		for(int i=0;i<ans.ULL_SIZE;++i)
			ans.U64[i] = T64[i];
	else
		for(int i=0;i<ans.ULL_SIZE;++i)
			ans.U64[i] = (T64[i]>>rsft)|(T64[i+1]<<lsft);
	//////////////////////////////////////////////////////////////
	/*
	len = (len <= (length-idx)) ? len : (length-idx);
	ULLN ans;
	int j=0;
	for(int i=((length-len-idx)>>5);i<=(length-1-idx)>>5;++i,++j)
		ans.U64[j] = U64[i];
	int rft = (1+((length-1-idx)&0x1f))*2;
	if(rft<64) 
		ans.U64[--j] &= ((1ULL<<rft)-1);
	ans >>= (((length-len-idx)&0x1f)*2);
	*/
	//////////////////////////////////////////////////////////////
	return ans;
}

void BaseStr::checkRead(int strposi,const ULLN& read,int readposi,const int readLen, int& match, int& mismatch)const//check reads alignment result. return # of same bases
{
	int idx = strposi - readposi;
	if(idx < 0)
	{
		idx = -idx;
		ULLN tread = subStr(0, readLen - idx);
		mismatch = (tread^=(read>>(idx*2))).non0base();
		match = readLen-idx;
	}
	else if (idx > length-readLen)
	{
		ULLN tread = subStr(idx, length - idx);
		tread ^= read>>(2*(readLen-length+idx));
		mismatch = tread.non0base();
		match = length-idx-mismatch;
	}
	else
	{
		ULLN tread = subStr(idx,readLen);
		mismatch = (tread ^= read).non0base();
		match = readLen - mismatch;
	}
	return;
}

ostream& operator<<(ostream& os,const BaseStr& obj)
{
	for(int i=obj.ULL_SIZE-1;i>0;--i)
		os << hex << obj.U64[i] << '-';
	os << hex << obj.U64[0];
	return os;
}

//position is [1,75];
/*void BaseStr::setbase(int position,unsigned long long base)
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
