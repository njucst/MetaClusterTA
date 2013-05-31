/*
 * last fixed: 2013.01.21.
 * by Wang Yi.
 * */
#ifndef MCH_HYBRID_BWTDTSTR_H_

#define MCH_HYBRID_BWTDTSTR_H_
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include "BWTPARA.h"
using namespace std;
class BWTDtStr
{
public:
	pair<long long, long long> search(const string& P)const;

	void iniStatic();
	unsigned getlength(){return length;}

	int* Taxid_of_SA;

	void preprocess(char* bwt_in_DESIGN_format, int* Taxid_of_SA_);
	void clear();
	BWTDtStr();
	virtual ~BWTDtStr();
private:
	static const int L1 = 128;//in terms of bits
	static const int L2 = 32;//in terms of bits

	static const unsigned IdxNuc = 10;
	static const unsigned IdxNumber = 1<<(2*IdxNuc);
	vector<long long> IdxStart;
	vector<long long> IdxEnd;

	unsigned char No_Of_Ones_In[1U<<16];
	unsigned short* bits[4];
	unsigned* AppearL1[4];
	unsigned char* AppearL2[4];
	int Count[9];unsigned length;
	int SL1,SL2,bitL;

	long long getAppear(int ch, long long posi)const;
	pair<long long, long long> naiveSearch(const vector<int>& query)const;
	pair<long long, long long> naiveSearchBeginWith(const vector<int>& query,long long s_start,long long s_end)const;
};
#endif
