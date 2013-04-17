/*
 * last fixed: 2012.11.11.
 * by Wang Yi.
 * */
#ifndef MCH_HYBRID_BASESTR_H_

#define MCH_HYBRID_BASESTR_H_
#include <iostream>
#include "ULLN.h"
using namespace std;
class BaseStr
{
public:
	const static int DEFAULT_S = 10;
	int ULL_SIZE, length;
	unsigned long long* U64;
	string str;//to calculate kmer-distribution. 

//	BaseStr(int capa);
	BaseStr(const BaseStr& obj);
//	BaseStr(const string& str);
	BaseStr(const string& str,bool isRev = false);
	BaseStr(const unsigned long long *t, int usize, int size_);
	virtual ~BaseStr()
	{
		delete[]U64;
	}

	int non0base();//# of non-0 bases
	ULLN subStr(int idx, int len)const;//idx & len are both based on bps. append 0 at the end of ULLN.
	void checkRead(int strposi,const ULLN& read,int readposi,const int readLen, int& match, int& mismatch)const;//check reads alignment result. return # of same bases

//	void setzero();
//	BaseStr substr(int idx,int len);
//	void resize(int ulls);

	bool operator < (const BaseStr& obj) const;
	bool operator == (const BaseStr& obj) const;
	bool operator != (const BaseStr& obj) const;
	BaseStr& operator <<= (const int obj);
	BaseStr  operator << (const int obj) const;
	BaseStr& operator >>= (const int obj);
	BaseStr  operator >> (const int obj) const;
	BaseStr& operator |= (const int obj);
	BaseStr& operator |= (const BaseStr& obj);
	BaseStr operator | (const BaseStr& obj)const;
	BaseStr& operator &= (const BaseStr& obj);
	BaseStr operator & (const BaseStr& obj)const;
	BaseStr& operator ^= (const BaseStr& obj);
	BaseStr  operator ^ (const BaseStr& obj) const;
	BaseStr operator~() const;

	BaseStr& operator=(const BaseStr&obj);
	friend ostream& operator<<(ostream& os,const BaseStr & obj);
private:
	BaseStr(){};
	BaseStr(int size){};
//	const BaseStr &operator=(const BaseStr &t){return *this;}
};
#endif
