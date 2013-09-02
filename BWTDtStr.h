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

	pair<int,int> getSidTaxidOfPosition(int i)const{return make_pair(Sid_of_SA[i],Id2Taxon[Sid_of_SA[i]]);}
	unsigned getlength()const{return length;}
	int getNSeq()const{return NSeq;}
//	void preprocess(char* bwt_in_DESIGN_format, int* Sid_of_SA_);
	void preprocess(char* bwt_in_DESIGN_format, int* Sid_of_SA_,long long nuc_length, int* Id2Taxon, int NSeq);
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

	int* Sid_of_SA;
	///////////////////////////////////
	unsigned char No_Of_Ones_In[1U<<16];
	unsigned short* bits[4];
	unsigned* AppearL1[4];
	unsigned char* AppearL2[4];
	int Count[9];unsigned length;
	int SL1,SL2,bitL;
	///////////////////////////////////
	int NSeq;
	int*Id2Taxon;

	long long getAppear(int ch, long long posi)const;
	pair<long long, long long> naiveSearch(const vector<int>& query)const;
	pair<long long, long long> naiveSearchBeginWith(const vector<int>& query,long long s_start,long long s_end)const;
};
#endif
