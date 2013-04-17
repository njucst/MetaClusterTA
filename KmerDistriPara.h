/*
 * last fixed: 2012.11.11.
 * by Wang Yi.
 * */
#ifndef MCH_HYBRID_KMERDISTRIPARA_H_

#define MCH_HYBRID_KMERDISTRIPARA_H_
#include<vector>
#include"TomAlgorithm.h"
using namespace std;
class KmerDistriPara
{
public:
	int KmerLen;
	int ReverSize;
	int NonReverSize;
	unsigned KmerMask;
	vector<int>RCIdx;

	KmerDistriPara(int KmerLen_)
	{
		KmerLen = KmerLen_;
		if(KmerLen & 1U)
			ReverSize = ((1<<(KmerLen<<1))>>1);
		else
			ReverSize = (((1<<(KmerLen<<1))+(1<<KmerLen))>>1);
		NonReverSize = 1<<(KmerLen<<1);
		KmerMask = (1<<(KmerLen<<1)) - 1;
		RCIdx.resize(NonReverSize);
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
	}
private:
	KmerDistriPara(){}
};
#endif
