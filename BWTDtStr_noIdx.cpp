#include "BWTDtStr.h"
#include <iostream>
#include <fstream>
using namespace std;

void BWTDtStr::iniStatic()
{
}

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
	Taxid_of_SA = NULL;
	for(unsigned i=0;i<(1U<<16);++i)
	{ 
		int oneN = 0;
		for(unsigned k = i;k;k &= (k-1), ++ oneN);
		No_Of_Ones_In[i] = oneN;
	}
}
void BWTDtStr::preprocess(char* bwt_in_DESIGN_format, int* Taxid_of_SA_)
{
	Taxid_of_SA = Taxid_of_SA_;
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
pair<long long, long long> BWTDtStr::search(const string& P)const
{
	if(P.length()==0)
		return make_pair(1ULL,0ULL);
	vector<int> query(P.length());
	for(unsigned i=0;i<P.length();++i)
	{
		switch(P[i])
		{
			case 'A':query[i]=0;break;
			case 'C':query[i]=1;break;
			case 'G':query[i]=2;break;
			case 'T':query[i]=3;break;
			case 'a':query[i]=0;break;
			case 'c':query[i]=1;break;
			case 'g':query[i]=2;break;
			case 't':query[i]=3;break;
			default: return make_pair(1ULL,0ULL);
		}
	}
	return search(query);
}
pair<long long, long long> BWTDtStr::search(const vector<int>& query)const
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
BWTDtStr::~BWTDtStr()
{
/*	delete[] Taxid_of_SA;
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
	if(Taxid_of_SA)
		delete[] Taxid_of_SA;
	Taxid_of_SA = NULL;
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
