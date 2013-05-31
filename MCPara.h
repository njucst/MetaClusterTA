/*
 * last fixed: 2012.11.11.
 * by Wang Yi.
 * */
#ifndef MCH_HYBRID_MCPARA_H_

#define MCH_HYBRID_MCPARA_H_
#include <map>
#include "MetaCluster.h"
#include "CtgsReads.h"
class MCPara
{
public:
	int KmerLen;
	int Size;
	int **ReverKmerDistri;
	bool* isoutlier;
	map<unsigned, unsigned>toNewId;
	map<unsigned, unsigned>toOldId;
	///////////////////////////
	int GenoNum;
	int **Component;
	///////////////////////////
	MCPara(int KmerLen_, const ContigsClass& Ctgs, USet& uset, const AccTester& acc_tester)
	{
		///////////////////////set outliers
		int usetsize = uset.size();
		bool* isoutlier = new bool[usetsize];
		for(unsigned i=0;i<usetsize;++i)
		{
			if(uset.getCtgLen(i)<1000)
				isoutlier[i] = true;
			else isoutlier[i] = false;
		}
		//last mark
		//////////////////////get new Id
		unsigned NewIdx = 0;
		for(unsigned i=0;i<usetsize;++i)
		{
			if(uset.find(i)!=i || isoutlier[i])continue;
			toNewId[i] = NewIdx;
			toOldId[NewIdx] = i;
			++NewIdx;
		}
		////////////////////////////////initialize kmer distribution
		Size = NewIdx;
		KmerLen = KmerLen_;
		ReverSize = tomReverSize(KmerLen);
		NonReverSize = 1<<(KmerLen<<1);
		KmerMask = (1<<(KmerLen<<1))-1;
		RCIdx = new int[NonReverSize];
		{
			int tidx = 0;
			for(int i=0;i<NonReverSize;++i)
			{
				int rev = tomReverComple(i,KmerLen);
				if(i <= rev)
				{
					RCIdx[i] = RCIdx[rev] = tidx;
					++tidx;
				}
			}
		}
		ReverKmerDistri = new int*[Size];
		for(int i=0;i<Size;++i)
		{
			ReverKmerDistri[i] = new int[ReverSize];
			for(int j=0;j<ReverSize;++j)
				ReverKmerDistri[i][j] = 0;
		}
		////////////////////////////////get k-mer distribution
		int CtgNum = Ctgs.CtgNum;
		int* DisBuf = new int[ReverSize];
		for(unsigned i=0;i<CtgNum;++i)
		{
			if(isoutlier[i])continue;
			int idx = uset.find(i);
			getRever((Ctgs.contigs[i]->str).c_str(),DisBuf);
			assert(toNewId.find(idx)!=toNewId.end());
			idx = toNewId[idx];
			for(int j=0;j<ReverSize;++j)
				ReverKmerDistri[idx][j] += DisBuf[j];
		}
		delete[] DisBuf;
		/////////////////////////////////////for testing
		GenoNum = acc_tester.GenomeNum;
		if(GenoNum > 0)
		{
			Component = new int*[Size];
			for(int i=0;i<Size;++i)
			{
				Component[i] = new int[GenoNum];
				for(int j=0;j<GenoNum;++j)
					Component[i][j] = 0;
			}
			for(unsigned i=0;i<CtgNum;++i)
			{
				if(isoutlier[i])continue;
				int idx = uset.find(i);
				assert(toNewId.find(idx)!=toNewId.end());
				idx = toNewId[idx];
				for(int j=0;j<GenoNum;++j)
					Component[idx][j] += acc_tester.CtgComp[i][j];
			}
		}
		else Component = NULL;
	}
	long long getOldId(unsigned cid)const
	{
		map<unsigned,unsigned>::const_iterator itr = toOldId.find(cid);
		if(itr==toOldId.end())
			return -1;
		return itr->second;
	}
	long long getNewId(unsigned cid)const
	{
		map<unsigned,unsigned>::const_iterator itr = toNewId.find(cid);
		if(itr==toNewId.end())
			return -1;
		return itr->second;
	}

	virtual ~MCPara()
	{
	/*	for(int i=0;i<Size;++i)
			delete[]ReverKmerDistri[i];
		delete[]ReverKmerDistri;
		if(Component!=NULL)
		{
			for(int i=0;i<Size;++i)
				delete[]Component[i];
			delete[]Component;
		}
		delete[]isoutlier;
		delete[]RCIdx;*/
	}
private:
	int ReverSize;
	int NonReverSize;
	unsigned KmerMask;
	int* RCIdx;
	void getRever(const char* str_,int *ReverDis)
	{
		vector<int>nonRever(NonReverSize,0);
		for(int i=0;i<ReverSize;++i)
			ReverDis[i] = 0;
		unsigned prek = 0;
		for(int i=0;i<KmerLen-1;++i)
		{
			prek <<= 2;
			prek |= str_[i];
		}
		const char* str = str_+(KmerLen-1);
		while(*str)
		{
			prek <<= 2;
			if(*str == 'C')
				prek |= 0x1; 
			else if(*str == 'G')
				prek |= 0x2; 
			else if(*str == 'T')
				prek |= 0x3; 
			else if(*str != 'A')
			{
				++str;
				continue;
			}
			++nonRever[prek & KmerMask];
			++str;
		}
		for(int i=0;i<NonReverSize;++i)
		{
			if(tomReverComple(i,KmerLen)==i)
				ReverDis[RCIdx[i]] = nonRever[i]*2;
			else
				ReverDis[RCIdx[i]] += nonRever[i];
		}
	}
	MCPara(){}
};
#endif
