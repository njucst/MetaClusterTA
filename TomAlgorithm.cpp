#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <utility>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include "TomAlgorithm.h"

using namespace std;

double* tomNormalize(double* distri, int size,bool newresult)
{
	if(newresult)
	{
		double sum=0.;
		for(int i=0;i<size;++i)
			sum += distri[i];
		if(sum==0)return NULL;

		double *result = new double[size];
		for(int i=0;i<size;++i)
			result[i] = distri[i]/sum;
		return result;
	}
	else
	{
		double sum=0.;
		for(int i=0;i<size;++i)
			sum += distri[i];
		if(sum==0)return NULL;
	
		for(int i=0;i<size;++i)
			distri[i] /= sum;
		return NULL;
	}
}

double* tomNormalize(int *distri, int size,bool remove)
{
	double sum=0.;
	for(int i=0;i<size;++i)
		sum += distri[i];
	if(sum==0)return NULL;

	double *result = new double[size];
	for(int i=0;i<size;++i)
		result[i] = distri[i]/sum;
	if(remove)
		delete[]distri;
	return result;
}
/*
inline int tomReverSize(int kmerLen)
{
	if(kmerLen & 1U)
		return ((1<<(kmerLen<<1))>>1);
	else
		return (((1<<(kmerLen<<1))+(1<<kmerLen))>>1);
}*/

//kmer is a number represents a kmer like AGGT or TTAC, return the reverse complement.
int tomReverComple(int kmer,int kmerLen)
{
	int result=0;
	unsigned orig=kmer;
	for(int i=0;i<kmerLen;++i)
	{
		unsigned digit=orig&3U;
		orig >>= 2;
		result <<= 2;
		result |= (~digit)&3U;
	}
	return result;
}

//return the position of each kmer(rever complement is considered to be the same) in the old kmer-vector
int* getReverPosition(const int KmerLen)
{
	size_t orisize = 1<<(KmerLen<<1);
	size_t newsize;
	if(KmerLen & 1)
		newsize = orisize>>1;
	else
		newsize = (orisize+(1<<KmerLen))>>1;
	int *result = new int[newsize];
	int count = 0;
	for(int i=0;i<orisize;++i)
	{
		if(tomReverComple(i,KmerLen)<i)continue;
		result[count++]=i;
	}
	return result;
}

////added at 2011.11.30
void getRever(char* str_,int *ReverDis, const int KmerLen, const unsigned KmerMask, const int NonReverSize,const int ReverSize,const int* RCIdx)
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
	char* str = str_+(KmerLen-1);
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
