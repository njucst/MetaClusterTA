/*
 * last fixed: 2012.01.05.
 * by Wang Yi.
 * */
#ifndef MCH_HYBRID_CTGSREADS_H_

#define MCH_HYBRID_CTGSREADS_H_

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include "Utils.h"
#include "KmerDistriPara.h"
#include "BaseStr.h"
#include "Structs.h"
#include "AccTester.h"
using namespace std;

class ContigsClass
{ 
public:
	string* info;

	BaseStr** contigs;
	int CtgNum;
	static const int MAXLINE = 4000000;

	ContigsClass()
	{
		contigs = NULL;
		CtgNum = 0;
	}
	~ContigsClass()
	{
		delete[]contigs;
		delete[]info;
	}
	ContigsClass(string path, int Thresh, KmerNode** KmerMap,KmerNodeAloc& AllNodes)
	{
		init(path,Thresh, KmerMap,AllNodes);
	}
/*	int** getSpear(const KmerDistriPara& K5mer)
	{
		int** ans = new int*[CtgNum];
		////////////////////////////////get k-mer distribution
		for(unsigned i=0;i<CtgNum;++i)
		{
			ans[i] = new int[K5mer.ReverSize];
			getRever((contigs[i]->str).c_str(),ans[i],K5mer);
		}
		return ans;
	}*/
	void insertHash(KmerNodeAloc& AllNodes, KmerNode** KmerMap, const KMER& kmer_, const unsigned hashidx, int id, int posi)
	{
		unsigned idx = hashidx;
		if(KmerMap[idx] == NULL)
			KmerMap[idx] = AllNodes.getNew(kmer_,id,posi);
		else
		{
			KmerNode* p = KmerMap[idx];
			for(;p->next!=NULL;p=p->next)
				if(p->kmer == kmer_)
				{
					p->push_back(id,posi);
					return;
				}
			if(p->kmer == kmer_)
				p->push_back(id,posi);
			else
				p->next = AllNodes.getNew(kmer_,id,posi);
		}
	}
	KmerNode* findKmer(KmerNode** KmerMap, const KMER & kmer_)
	{
		KmerNode* ans = KmerMap[kmer_.hash()];
		for(;ans!=NULL;ans = ans->next)
			if(ans->kmer == kmer_)
				return ans;
		return NULL;
	}
	void init(string path, int Thresh, KmerNode** KmerMap,KmerNodeAloc& AllNodes)
	{
		ifstream ifs(path.c_str());
		if(ifs.fail())
		{
			cerr << "File open failed:\t"<<path<<endl;
			exit(-1);
		}
		char* Buf = new char[MAXLINE];
		CtgNum = getFileLine(path.c_str())/2;
		BaseStr** contigs2 = new BaseStr*[CtgNum];
		int idx = 0;
		string* ctgstr= new string[CtgNum];
		info = new string[CtgNum];
		while(!ifs.eof())
		{
			ifs.getline(Buf,MAXLINE);
			info[idx] = Buf;
			ifs.getline(Buf,MAXLINE);
			if(strlen(Buf) >= Thresh)
			{
				ctgstr[idx] = Buf;
				contigs2[idx++] = new BaseStr(Buf);
			}
		}
		ifs.close();
		delete[]Buf;
		CtgNum = idx;
		contigs = new BaseStr*[CtgNum];
		for(int i=0;i<CtgNum;++i)
			contigs[i] = contigs2[i];
		delete[]contigs2;
/////////////////////////////////////////////////////////////
		omp_lock_t* hash_lock = new omp_lock_t[1U<<20];
		for(int i=0;i<(1U<<20);++i)
		omp_init_lock(&hash_lock[i]);
#pragma omp parallel for schedule(dynamic)
		for(int i=0;i<idx;++i)
			mapto(hash_lock, ctgstr[i], i, AllNodes, KmerMap);
		for(int i=0;i<(1U<<20);++i)
			omp_destroy_lock(&hash_lock[i]);
		delete[] hash_lock;
		delete[]ctgstr;
	}
	void mapto(omp_lock_t* hash_lock,const string& str, int idx, KmerNodeAloc& AllNodes, KmerNode** KmerMap)
	{
		assert(KmerMap !=NULL);
		KMER kmer0(str);
		KMER kmer1(str,true);
		for(int i=PARA_KMER;i<str.length();++i)
		{
			unsigned b=0;
			if(str[i]=='C' || str[i]=='c')
				b=1;
			else if(str[i]=='G' || str[i]=='g')
				b=2;
			else if(str[i]=='T' || str[i]=='t')
				b=3;
			kmer0.shiftInLow(b);
			kmer1.shiftInHigh(3-b);
			if(kmer0<kmer1)
			{
				unsigned hashidx = kmer0.hash();
				omp_set_lock(&hash_lock[hashidx&0xfffff]);
				insertHash(AllNodes, KmerMap, kmer0,hashidx, idx, (i<<2)|0x2);
				omp_unset_lock(&hash_lock[hashidx&0xfffff]);
			}
			else
			{
				unsigned hashidx = kmer1.hash();
				omp_set_lock(&hash_lock[hashidx&0xfffff]);
				insertHash(AllNodes, KmerMap, kmer1,hashidx, idx, (i<<2)|0x3);
				omp_unset_lock(&hash_lock[hashidx&0xfffff]);
			}
		}
	}
private:
	ContigsClass(const ContigsClass& t){}
	const ContigsClass& operator=(const ContigsClass& t){return *this;}

	void getRever(const char* str_,int *ReverDis,const KmerDistriPara& K5mer)
	{
		int NonReverSize = K5mer.NonReverSize;
		int ReverSize = K5mer.ReverSize;
		int KmerLen = K5mer.KmerLen;
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
			++nonRever[prek & K5mer.KmerMask];
			++str;
		}
		for(int i=0;i<NonReverSize;++i)
		{
			if(tomReverComple(i,K5mer.KmerLen)==i)
				ReverDis[K5mer.RCIdx[i]] = nonRever[i]*2;
			else
				ReverDis[K5mer.RCIdx[i]] += nonRever[i];
		}
	}
};

class ReadsClass
{
public:
	static const int MAXLINE = 1000000;
	///////////////////////////////////////////
	int ReadLen ;
	int AlignThresh ;
	///////////////////////////////////////////
	ULLN** R1;//reads
	ULLN** R2;//reverse complement of reads
	int* MatchId;//for unmapped reads, this value<0;
	int* Score;
	int* MPosi;//for unmapped reads, it's new Id.
	int* NewIdToOldId;
	unsigned long long ReadNum;
	unsigned long long TotalNum;
	vector<unsigned>ReadNumInCtg;
	ReadsClass()
	{
		R1 = R2 = NULL;
		ReadNum = 0;
	}
	ReadsClass(string path, int ReadLen_, int AlignThresh_, KmerNode** KmerMap,KmerNodeAloc& AllNodes,ContigsClass& Ctgs, AccTester& acc_tester)
	{
		init(path, ReadLen_, AlignThresh_, KmerMap,AllNodes,Ctgs,acc_tester);
	}
	void init(string path, int ReadLen_, int AlignThresh_, KmerNode** KmerMap,KmerNodeAloc& AllNodes,ContigsClass& Ctgs,AccTester& acc_tester)
	{
		ifstream ifs(path.c_str());
		if(ifs.fail())
		{
			cerr << "File open failed:\t"<<path<<endl;
			exit(-1);
		}
		ReadLen = ReadLen_;
		AlignThresh = AlignThresh_;
		ReadNumInCtg.resize(Ctgs.CtgNum,0);

		char* Buf = new char[MAXLINE];
		ReadNum = getFileLine(path.c_str())/2;
		TotalNum = ReadNum;
		if(INTEST)acc_tester.init_total(TotalNum);
		MatchId = new int[ReadNum];
		Score = new int[ReadNum];
		MPosi = new int[ReadNum];
		NewIdToOldId = new int[ReadNum];
		for(unsigned long long i=0;i<TotalNum;++i)
			MatchId[i] = -1;
		for(unsigned long long i=0;i<TotalNum;++i)
			Score[i] = 0;

		ULLN** T1 = new ULLN*[ReadNum];
		ULLN** T2 = new ULLN*[ReadNum];
		long long totalidx=-1;
		string* readstrs= new string[ReadNum];

		string* allreads = new string[ReadNum];
		while(!ifs.eof())
		{
			++totalidx;
			ifs.getline(Buf,MAXLINE);
			if(ifs.eof())break;
			if(INTEST)	acc_tester.getGId(totalidx,Buf);
			////////////////////////////////////////////////////
			ifs.getline(Buf,MAXLINE);
			Buf[ReadLen]=0;
			allreads[totalidx] = Buf;
		}
		ifs.close();
		delete[]Buf;
		totalidx = ((totalidx+1)/2)<<1;

		long long idx = -1;
		omp_lock_t input_lock;
		omp_init_lock(&input_lock);
#pragma omp parallel for schedule(dynamic)
		for(int i=0;i<totalidx;++i)
		{
			bool isAligned =  align(allreads[i],i,KmerMap,Ctgs);
			if(!isAligned)
			{
				omp_set_lock(&input_lock);
				long long tidx = ++idx;
				omp_unset_lock(&input_lock);
				readstrs[tidx] = allreads[i];
				T1[tidx] = new ULLN(allreads[i]);
				T2[tidx] = new ULLN(allreads[i],true);
				//////////////////////////
				MatchId[i] = -2-MatchId[i];
				MPosi[i] = tidx;
				NewIdToOldId[tidx] = i;
				///////////////////////////
			}
		}
		omp_destroy_lock(&input_lock);
		delete[]allreads;
		++idx;

	/*	while(!ifs.eof())
		{
			++totalidx;
			ifs.getline(Buf,MAXLINE);
			if(ifs.eof())break;
			////////////////////////////////////////////////////
			if(INTEST)	acc_tester.getGId(totalidx,Buf);
			////////////////////////////////////////////////////
			ifs.getline(Buf,MAXLINE);
			Buf[ReadLen]=0;
			////////////////////////////////////////////////////
			bool isAligned =  align(Buf,totalidx,KmerMap,Ctgs);
			if(!isAligned)
			{
				readstrs[idx] = Buf;
				T1[idx] = new ULLN(Buf);
				T2[idx] = new ULLN(Buf,true);
				//////////////////////////
				MatchId[totalidx] = -1-MatchId[totalidx];
				MPosi[totalidx] = idx;
				NewIdToOldId[idx] = totalidx;
				///////////////////////////
				++idx;
			}
			////////////////////////////////////////////////////
		}
		ifs.close();
		delete[]Buf;*/
		ReadNum = idx;
		R1 = new ULLN*[ReadNum];
		R2 = new ULLN*[ReadNum];
		for(int i=0;i<ReadNum;++i)
			R1[i] = T1[i];
		for(int i=0;i<ReadNum;++i)
			R2[i] = T2[i];
		delete[]T1;
		delete[]T2;
/////////////////////////////////////////////////////////////
		omp_lock_t* hash_lock = new omp_lock_t[1U<<20];
		for(int i=0;i<(1U<<20);++i)
		omp_init_lock(&hash_lock[i]);
#pragma omp parallel for schedule(dynamic)
		for(int i=0;i<ReadNum;++i)
			mapto(hash_lock, readstrs[i], i, AllNodes, KmerMap);
		for(int i=0;i<(1U<<20);++i)
			omp_destroy_lock(&hash_lock[i]);
		delete[] hash_lock;
		delete[]readstrs;
////////////////////////////////////////////////////////////////
	}
	int toNewId(int id)const
	{
		return MPosi[id];
	}
	virtual ~ReadsClass()
	{
		delete[]MatchId;
		delete[]Score;
		delete[]MPosi;
		delete[]NewIdToOldId;
		for(int i=0;i<ReadNum;++i)
			delete R1[i];
		delete[]R1;
		for(int i=0;i<ReadNum;++i)
			delete R2[i];
		delete[]R2;
	}
	void insertHash(KmerNodeAloc& AllNodes, KmerNode** KmerMap, const KMER& kmer_, const unsigned hashidx, int id, int posi)
	{
		unsigned idx = hashidx;
		if(KmerMap[idx] == NULL)
			KmerMap[idx] = AllNodes.getNew(kmer_,id,posi);
		else
		{
			KmerNode* p = KmerMap[idx];
			for(;p->next!=NULL;p=p->next)
				if(p->kmer == kmer_)
				{
					p->push_back(id,posi);
					return;
				}
			if(p->kmer == kmer_)
				p->push_back(id,posi);
			else
				p->next = AllNodes.getNew(kmer_,id,posi);
		}
	}
	KmerNode* findKmer(KmerNode** KmerMap, const KMER & kmer_)
	{
		KmerNode* ans = KmerMap[kmer_.hash()];
		for(;ans!=NULL;ans = ans->next)
			if(ans->kmer == kmer_)
				return ans;
		return NULL;
	}

	bool align(string& str, const int totalidx, KmerNode** KmerMap, ContigsClass& Ctgs)
	{
		KMER kmer0(str);
		KMER kmer1(str,true);
		ULLN tR1(str);
		ULLN tR2(str,true);
		assert(str.length()==ReadLen);
		for(int i=PARA_KMER;i<ReadLen;++i)
		{
			unsigned b=0;
			if(str[i]=='C' || str[i]=='c')
				b=1;
			else if(str[i]=='G' || str[i]=='g')
				b=2;
			else if(str[i]=='T' || str[i]=='t')
				b=3;
			kmer0.shiftInLow(b);
			kmer1.shiftInHigh(3-b);
			///////////////////////////////////////////////////////////////
			//align
			int match = 0, mismatch = ReadLen;
			if(kmer0<kmer1)
			{
				KmerNode* cur = findKmer(KmerMap,kmer0);
				if(cur == NULL)continue;
				NodeIDPosi* tVec = cur->myvector;
				for(int j=0;j<cur->VSize && (tVec[j].posi & 0x2);++j)
				{
					if(tVec[j].posi & 0x1)
						Ctgs.contigs[tVec[j].id]->checkRead(tVec[j].posi>>2,tR2,PARA_KMER-2+ReadLen-i,ReadLen,match,mismatch);
					else
						Ctgs.contigs[tVec[j].id]->checkRead(tVec[j].posi>>2,tR1,i,ReadLen,match,mismatch);
					if(match > Score[totalidx])
					{
						Score[totalidx] = match;
						MatchId[totalidx] = tVec[j].id;
						if(match==ReadLen)
						{
							++ReadNumInCtg[MatchId[totalidx]];
							return true;
						}
					}
				}
			}
			else
			{
				KmerNode* cur = findKmer(KmerMap,kmer1);
				if(cur == NULL)continue;
				NodeIDPosi* tVec = cur->myvector;
				for(int j=0;j<cur->VSize && (tVec[j].posi & 0x2);++j)
				{
					if(tVec[j].posi & 0x1)
						Ctgs.contigs[tVec[j].id]->checkRead(tVec[j].posi>>2,tR1,i,ReadLen,match,mismatch);
					else
						Ctgs.contigs[tVec[j].id]->checkRead(tVec[j].posi>>2,tR2,PARA_KMER-2+ReadLen-i,ReadLen,match,mismatch);
					if(match > Score[totalidx])
					{
						Score[totalidx] = match;
						MatchId[totalidx] = tVec[j].id;
						if(match==ReadLen)
						{
							++ReadNumInCtg[MatchId[totalidx]];
							return true;
						}
					}
				}
			}
		}
		if(Score[totalidx] >= AlignThresh)
		{
			++ReadNumInCtg[MatchId[totalidx]];
			return true;
		}
		return false;
	}
	void mapto(omp_lock_t* hash_lock,const string str, int idx, KmerNodeAloc& AllNodes, KmerNode** KmerMap)
	{
		assert(KmerMap !=NULL);
		assert(str.length()==ReadLen);
		KMER kmer0(str);
		KMER kmer1(str,true);
		int RL = ReadLen;
		for(int i=PARA_KMER;i< RL;++i)
		{
			unsigned b=0;
			if(str[i]=='C' || str[i]=='c')
				b=1;
			else if(str[i]=='G' || str[i]=='g')
				b=2;
			else if(str[i]=='T' || str[i]=='t')
				b=3;
			kmer0.shiftInLow(b);
			kmer1.shiftInHigh(3-b);
			if(kmer0<kmer1)
			{
				unsigned hashidx = kmer0.hash();
				omp_set_lock(&hash_lock[hashidx&0xfffff]);
				insertHash(AllNodes, KmerMap, kmer0,hashidx, idx, (1-PARA_KMER+i)<<2);
				omp_unset_lock(&hash_lock[hashidx&0xfffff]);
			}
			else
			{
				unsigned hashidx = kmer1.hash();
				omp_set_lock(&hash_lock[hashidx&0xfffff]);
				insertHash(AllNodes, KmerMap, kmer1, hashidx, idx, ((RL-i-1)<<2)|0x1);
				omp_unset_lock(&hash_lock[hashidx&0xfffff]);
			}
		}
	}
private:
	ReadsClass(const ReadsClass& t){}
	const ReadsClass& operator=(const ReadsClass& t){return *this;}
};
#endif
