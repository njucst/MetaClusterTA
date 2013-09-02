#include "BWTDtStr.h"
#include <iostream>
#include <fstream>
#include <omp.h>
using namespace std;

inline void set_ith_bit_1(unsigned short* arr,unsigned i)
{
	arr[i/16] |= (1<<(i%16));
}
BWTDtStr::BWTDtStr()
{
	for(int i=0;i<4;++i)
	{
		bits[i] = NULL;
		AppearL1[i] = NULL;
		AppearL2[i] = NULL;
	}
	Sid_of_SA = NULL;
	NSeq = 0;
	Id2Taxon = NULL;
	for(unsigned i=0;i<(1U<<16);++i)
	{ 
		int oneN = 0;
		for(unsigned k = i;k;k &= (k-1), ++ oneN);
		No_Of_Ones_In[i] = oneN;
	}
}

void BWTDtStr::preprocess(char* bwt_in_DESIGN_format, int* Sid_of_SA_, long long nuc_length, int* Id2Taxon_, int NSeq_)
{
	Sid_of_SA = new int[nuc_length];
	memcpy(Sid_of_SA, Sid_of_SA_, nuc_length*4);
	NSeq = NSeq_;
	Id2Taxon = new int[NSeq];
	memcpy(Id2Taxon, Id2Taxon_, NSeq*4);

	char* bwt = bwt_in_DESIGN_format;
	////////get bwt length
	length = strlen(bwt);
	if(bwt[length-1]&0x7)
		length *=2;
	else
		length = length*2-1;
cerr << "bwtdtstr length: " << length << endl;
	///////get bwt
	SL1 = length/(L1)+5;
	SL2 = length/(L2)+5;
	bitL = length/(sizeof(unsigned short)*8)+5;
	for(int i=0;i<4;++i)
	{
		AppearL1[i] = new unsigned[SL1];
		AppearL2[i]  = new unsigned char[SL2];
		bits[i] = new unsigned short[bitL];
		for(int j=0;j<bitL;++j)
			bits[i][j] = 0;
		for(int j=0;j<SL1;++j)
			AppearL1[i][j] = 0;
		for(int j=0;j<SL2;++j)
			AppearL2[i][j] = 0;
	}
	//////get Count[]&appear[] for bwt
	{
		for(int i=0;i<9;++i)
			Count[i] = 0;

		int idx = 0;
		unsigned Sum[4],tsum[4];
		for(int i=0;i<4;++i)
			Sum[i] = tsum[i] = 0;
		for(char* str = bwt;*str;++str)
		{
			unsigned char c1=((*str))>>4;
			unsigned char c2=(*str)&0x7;
			if(idx%L1==0)
				for(int i=0;i<4;++i)
				{
					AppearL1[i][idx/L1] = Sum[i];
					tsum[i] = 0;
				}
			if(idx%L2==0)
				for(int i=0;i<4;++i)
					AppearL2[i][idx/L2] = tsum[i];
			if(!c1)break;
			++Count[c1];
			if(c1>=4 && c1<=7)
			{
				set_ith_bit_1(bits[c1-4],idx);
				++tsum[c1-4];
				++Sum[c1-4];
			}
			++idx;

			if(idx%L1==0)
				for(int i=0;i<4;++i)
				{
					AppearL1[i][idx/L1] = Sum[i];
					tsum[i] = 0;
				}
			if(idx%L2==0)
				for(int i=0;i<4;++i)
					AppearL2[i][idx/L2] = tsum[i];
			if(!c2)break;
			++Count[c2];
			if(c2>=4 && c2<=7)
			{
				set_ith_bit_1(bits[c2-4],idx);
				++tsum[c2-4];
				++Sum[c2-4];
			}
			++idx;
		}
		////////////////////modify Count[]
		{
			Count[0] = 0;
			for(int i=1;i<9;++i)
				Count[i] += Count[i-1];
			for(int i=0;i<6;++i)
				Count[i] = Count[i+3];
			for(int i=6;i<9;++i)
				Count[i] = 0;
		}
	}
	//////////////////////////////////////////////////
	//build index
	IdxStart.resize(IdxNumber);
	IdxEnd.resize(IdxNumber);
#pragma omp parallel for 
	for(unsigned i=0;i<IdxNumber;++i)
	{
		vector<int> query(IdxNuc,0);
		for(unsigned tidx = IdxNuc-1, x=i;x;--tidx)
		{
			query[tidx] = x&0x3;
			x >>= 2;
		}
		pair<long long, long long> ans = naiveSearch(query);
		if(ans.second < ans.first)
		{
			IdxStart[i] = 1;
			IdxEnd[i] = 0;
		}
		else
		{
			IdxStart[i] = ans.first;
			IdxEnd[i] = ans.second;
		}
	}
}
long long BWTDtStr::getAppear(int ch, long long posi)const
{
	if(posi<=0)return 0;
	long long ans = AppearL1[ch][posi/L1] + AppearL2[ch][posi/L2];
	for(long long i=posi/L2*(L2/16);i<posi/16;++i)
		ans += No_Of_Ones_In[bits[ch][i]];
	ans += No_Of_Ones_In[bits[ch][posi/16]&((1<<(posi&0xf))-1)];
//	if(posi%16 != 0)
//		ans += No_Of_Ones_In[(bits[ch][posi/16] << (16-posi%16))&0xffff];
	return ans;
}

inline int nucToInt(const char ch)
{
	switch(ch)
	{
		case 'A':return 0;
		case 'C':return 1;
		case 'G':return 2;
		case 'T':return 3;
		case 'a':return 0;
		case 'c':return 1;
		case 'g':return 2;
		case 't':return 3;
	}
	return -1;
}
pair<long long, long long> BWTDtStr::search(const string& P)const
{
	if(P.length()==0)
		return make_pair(1ULL,0ULL);
	for(unsigned i=0;i<P.length();++i)
		if(nucToInt(P[i])<0)return make_pair(1ULL,0ULL);
	if(P.length() <= IdxNuc)
	{
		vector<int> query(P.length());
		for(unsigned i=0;i<P.length();++i)
			query[i] = nucToInt(P[i]);
		return naiveSearch(query);
	}

	unsigned anchor = 0;
	for(int i=P.length()-IdxNuc;i<P.length();++i)
	{
		anchor <<= 2;
		anchor |= nucToInt(P[i]);
	}
	vector<int> query(P.length()-IdxNuc);
	for(unsigned i=0;i<(P.length()-IdxNuc);++i)
		query[i] = nucToInt(P[i]);
	return naiveSearchBeginWith(query,IdxStart[anchor],IdxEnd[anchor]);
}
pair<long long, long long> BWTDtStr::naiveSearch(const vector<int>& query)const
{
	vector<int>::const_reverse_iterator itr= query.rbegin();
	long long first = Count[*itr];
	long long last = Count[(*itr)+1]-1;
	for(++itr;itr!=query.rend() && first<=last;++itr)
	{
		first = Count[*itr] + getAppear(*itr, first);
		last  = Count[*itr] + getAppear(*itr, last+1)-1;
	}
	return make_pair(first,last);
}
pair<long long, long long> BWTDtStr::naiveSearchBeginWith(const vector<int>& query,long long s_start,long long s_end)const
{
	long long first = s_start;
	long long last = s_end;
	for(vector<int>::const_reverse_iterator itr= query.rbegin();itr!=query.rend() && first<=last;++itr)
	{
		first = Count[*itr] + getAppear(*itr, first);
		last  = Count[*itr] + getAppear(*itr, last+1)-1;
	}
	return make_pair(first,last);
}
BWTDtStr::~BWTDtStr()
{
/*	if(Sid_of_SA)
		delete[] Sid_of_SA;
	if(Id2Taxon)
		delete[] Id2Taxon;
	for(int i=0;i<4;++i)
	{
		if(bits[i])
			delete[] bits[i];
		if(AppearL1[i])
			delete[] AppearL1[i];
		if(AppearL2[i])
			delete[] AppearL2[i];
	}*/
}
void BWTDtStr::clear()
{
	if(Sid_of_SA)
		delete[] Sid_of_SA;
	Sid_of_SA = NULL;
	if(Id2Taxon)
		delete[] Id2Taxon;
	Id2Taxon = NULL;
	for(int i=0;i<4;++i)
	{
		if(bits[i])
			delete[] bits[i];
		if(AppearL1[i])
			delete[] AppearL1[i];
		if(AppearL2[i])
			delete[] AppearL2[i];
		bits[i] = NULL;
		AppearL1[i] = NULL;
		AppearL2[i] = NULL;
	}
}
