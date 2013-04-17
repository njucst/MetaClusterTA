#ifndef __USET_H_

#define __USET_H_
#include <iostream>
using namespace std;
class USet
{
public:
	explicit USet(unsigned size_)
	{
		parent = new unsigned[size_];
		rank = new unsigned[size_];
		readnum = new unsigned[size_];
		ctglen = new unsigned[size_];
		Size = size_;
		CtgNum = 0;
		for(unsigned i=0;i<size_;++i)
			parent[i]=i;
		for(unsigned i=0;i<size_;++i)
			rank[i]=0;
		for(unsigned i=0;i<size_;++i)
			readnum[i]=1;
		for(unsigned i=0;i<size_;++i)
			ctglen[i]=0;
	}
	USet()
	{
		parent = NULL;
		rank = NULL;
		readnum = NULL;
		ctglen = NULL;
		Size = 0;
		CtgNum = 0;
	}
	void init(const unsigned size_,const unsigned CtgNum_,const unsigned*CtgLen)
	{
		parent = new unsigned[size_];
		rank = new unsigned[size_];
		readnum = new unsigned[size_];
		ctglen = new unsigned[size_];
		Size = size_;
		CtgNum = CtgNum_;
		for(unsigned i=0;i<size_;++i)
			parent[i]=i;
		for(unsigned i=0;i<size_;++i)
			rank[i]=0;
		/////////////////////////////////
		for(unsigned i=0;i<CtgNum;++i)
			readnum[i] = 0;
		for(unsigned i=CtgNum;i<size_;++i)
			readnum[i]=1;
		for(unsigned i=0;i<CtgNum;++i)
			ctglen[i] = CtgLen[i];
		for(unsigned i=CtgNum;i<size_;++i)
			ctglen[i]=0;
	}

	virtual ~USet()
	{
		delete []parent;
		delete []rank;
		delete []readnum;
		delete []ctglen;
	}

//	void ReInitial();
//	void make_set(const unsigned &x);
//	void incctglen(unsigned i);
	unsigned length()const{return Size;}
	unsigned size()const{return Size;}
	unsigned ctgsize()const{return CtgNum;}

	unsigned getReadNum(const unsigned x);
	unsigned getCtgLen(const unsigned x);
	unsigned find(const unsigned x);
	void Union(unsigned x,unsigned y);

private:
	unsigned* parent,*rank;
	unsigned* readnum,*ctglen;
	unsigned CtgNum;
	unsigned Size;
	USet(const USet &uset){}
	const USet &operator=(const USet &uset){return *this;}
};
#endif
