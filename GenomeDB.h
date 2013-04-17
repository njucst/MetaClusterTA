/*
 * last fixed: 2012.11.30.
 * by Wang Yi.
 * */
#ifndef MCH_HYBRID_GENOMESCLASS_H_

#define MCH_HYBRID_GENOMESCLASS_H_

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <omp.h>
#include "BaseStr.h"
#include "Structs.h"
#include "Utils.h"
#include "SparseMatrix.h"

using namespace std;

class GenomesClass
{ 
private:
	const static int GMAX = 3000;
	string Path;
	int GenomeNum;
public:
	int **Taxo;
	BaseStr** G1;
	BaseStr** G2;
	string* description;
	static const int MAXLINE = 20000000;
	static const int MAXLENGTH = 20000000;

	int size(){return GenomeNum;}
	int length(){return GenomeNum;}
	GenomesClass()
	{
		Taxo = NULL;
		G1 = NULL;
		G2 = NULL;
		description = NULL;
		GenomeNum = 0;
	}
	GenomesClass(string path, KmerNode** KmerMap,KmerNodeAloc& AllNodes)
	{
		init(path, KmerMap,AllNodes);
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
	void init(string path, KmerNode** KmerMap,KmerNodeAloc& AllNodes)
	{
		ifstream ifs(path.c_str());
		if(ifs.fail())
		{
			std::cerr<<"File open failed: "<<path<<endl;
			GenomeNum = 0;
			return;
		}
		Path = path;
		char* gffbuf = new char[MAXLINE];
		char* genomebuf = new char[MAXLENGTH];
		const int TaxoBufMax = 10000;
		char Buf[TaxoBufMax];
		ifs.getline(Buf,TaxoBufMax);
		int gidx = 0;
		int taxo[7];

		GenomeNum = getFileLine(path.c_str())-1;
		G1 = new BaseStr*[GenomeNum];
		G2 = new BaseStr*[GenomeNum];
		description = new string[GenomeNum];
		string* ctgstr= new string[GenomeNum];
		Taxo = new int*[GenomeNum];
		Taxo[0] = new int[7];
		while(!ifs.eof())
		{
			if(gidx > 200)break;
			ifs.getline(Buf,TaxoBufMax);
			char* inp = Buf;
			int tabc = 0;
			while(tabc<4 && (*inp))
			{
				if(*inp == '\t')++tabc;
				++inp;
			}
			int tmp;
			sscanf(inp,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",&taxo[0],&taxo[1],&taxo[2],&taxo[3],&taxo[4],&taxo[5],&tmp,&taxo[6]);
			bool allTaxoExist = true;
			for(int i=0;i<7;++i)
				if(taxo[i]<=0)
				{
					allTaxoExist = false;
					break;
				}
			if(allTaxoExist)
			{
				for(int i=0;i<7;++i)
					Taxo[gidx][i] = taxo[i];
				if((int)inp[strlen(inp)-1]==13)
					inp[strlen(inp)-1] = 0;
				inp = Buf;
				while(*inp != '/')
					++inp;
				//////////////////////////////////////////////////////////////////
				FILE* fp = fopen(inp,"rt");
				int glength = 0;
				if(fp == NULL)
				{
					cerr << "can not open " << inp <<'\t'<<(int)inp[strlen(inp)-1]<< endl;
					continue;
				}
				fgets(gffbuf,MAXLINE,fp);
				if(gffbuf[0]!='>')
				{
					cerr << "First line is in wrong format:\t"<<gffbuf<<endl;
					throw exception();
				}
				description[gidx] = gffbuf;
				while(fgets(gffbuf, MAXLINE, fp) != NULL)
				{
					if(gffbuf[0] != '>')
					{
						if(gffbuf[strlen(gffbuf)-1]=='\n')
							gffbuf[strlen(gffbuf)-1]=0;
						strcpy(genomebuf+glength,gffbuf);
						glength += strlen(gffbuf);
					}
					else
					{
						cerr << "Wrong format:\t"<<gffbuf<<endl;
						throw exception();
					}
				}
				fclose(fp);
				//////////////////////////////////////////////////////////////////
				if(glength>0)
				{
					cerr << gidx << " has been processed. " << endl;
					ctgstr[gidx] = genomebuf;
					G1[gidx] = new BaseStr(genomebuf);
					G2[gidx] = new BaseStr(genomebuf,true);
					++gidx;
					Taxo[gidx] = new int[7];
				}
				else
					cerr << "It's not added into database: " << inp << endl;
			}
		}
		cerr<<gidx<<" out of "<<GenomeNum<<" genomes are loaded."<<endl;
		GenomeNum = gidx;
		ifs.close();
		delete[]gffbuf;
		delete[]genomebuf;
//////////////////////////////////////////////////////////////////////////////
		omp_lock_t* hash_lock = new omp_lock_t[1U<<20];
		for(int i=0;i<(1U<<20);++i)
		omp_init_lock(&hash_lock[i]);
#pragma omp parallel for
		for(int i=0;i<GenomeNum;++i)
			mapto(hash_lock, ctgstr[i], i, AllNodes, KmerMap);
		for(int i=0;i<(1U<<20);++i)
			omp_destroy_lock(&hash_lock[i]);
		delete[] hash_lock;
		delete[]ctgstr;
		cerr << "Genome DB Initialization is finished."<<endl;
	}
	void mapto(omp_lock_t* hash_lock,const string& str, int idx, KmerNodeAloc& AllNodes, KmerNode** KmerMap)
	{
		assert(KmerMap !=NULL);
		KMER kmer0(str);
		KMER kmer1(str,true);
		int tlength = str.length()-1;
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
				insertHash(AllNodes, KmerMap, kmer0,hashidx, idx, (1-PARA_KMER+i<<2)|0x2);
				omp_unset_lock(&hash_lock[hashidx&0xfffff]);
			}
			else
			{
				unsigned hashidx = kmer1.hash();
				omp_set_lock(&hash_lock[hashidx&0xfffff]);
				insertHash(AllNodes, KmerMap, kmer1,hashidx, idx, ((tlength-i)<<2)|0x3);
				omp_unset_lock(&hash_lock[hashidx&0xfffff]);
			}
	//		assert(!b || (kmer0.U64[0] && kmer1.U64[0]));
/*			if(b && kmer0.U64[0]==0 || kmer1.U64[0]==0)
			{
				cerr << "kmer error. " <<endl;
				cerr << str[i] << '\t' << kmer0 << '\t' << kmer1 << '\t' << b << endl;
				exit(-1);
			}*/
		}
	}
	void getStatis(const KmerNodeAloc& AllNodes)
	{
		cerr << "Start to get statistics. " << endl;
		cout << "vector size:\t" << endl;
		map<int,int>vecSize;
		for(int i=AllNodes.getNodeNum()-1;i>=0;--i)
		{
			int ct = AllNodes.getVSize(i);
			if(ct > 10)ct = (ct/10)*10;
			++vecSize[ct];
			if(ct > 1000)
				cout <<hex << AllNodes.getKmer(i).U64[0] << '\t' << dec << ct << endl;
		}
		for(map<int,int>::const_iterator itr = vecSize.begin();itr!=vecSize.end();++itr)
			cout << (itr->first) << '\t' << itr->second << endl;

		long long shareMtxMulti[3000][8];
		long long shareMtxSingle[3000][8];
		for(int i=0;i<3000;++i)
			for(int j=0;j<8;++j)
				shareMtxMulti[i][j]=shareMtxSingle[i][j]=0;

		SparseMatrix MTX1;
		SparseMatrix MTX2;
		for(int i=AllNodes.getNodeNum()-1;i>=0;--i)
		{
			int tsize = AllNodes.getVSize(i);
			NodeIDPosi* vec = AllNodes.getVector(i);
			map<int,int> comp;
			for(int j=0;j<tsize;++j)
				++comp[vec[j].id];
			if(comp.size()<=1)continue;
			for(map<int,int>::const_iterator itr1 = comp.begin();itr1!=comp.end();++itr1)
			{
				map<int,int>::const_iterator itr2 = itr1;
				for(++itr2;itr2!=comp.end();++itr2)
				{
					int sidx;
					for(sidx =6;sidx>=0;--sidx)
					{
						if(Taxo[itr1->first][sidx]==Taxo[itr2->first][sidx])
							break;
					}
					if(sidx<0)sidx = 7;
					MTX1.insert(itr1->first,itr2->first);
					MTX2.insert(itr1->first,itr2->first,min(itr1->second,itr2->second));
				}
			}
		}
		map<unsigned,map<unsigned,unsigned> > Neighbor;
		cout << "Using Multi-Occurence: "<<endl;
		MTX1.toNeighbor(Neighbor);
		for(map<unsigned,map<unsigned,unsigned> >::const_iterator itr1 = Neighbor.begin();itr1!=Neighbor.end();++itr1)
		{
			for(map<unsigned,unsigned>::const_iterator itr2 = itr1->second.begin();itr2!=itr1->second.end();++itr2)
			{
				int sidx;
				for(sidx =6;sidx>=0;--sidx)
				{
					if(Taxo[itr1->first][sidx]==Taxo[itr2->first][sidx])
						break;
				}
				if(sidx<0)sidx = 7;
				cout << (itr1->first) << '\t' << (itr2->first) << '\t' << sidx <<'\t'<<itr2->second << endl;
			}
		}

		cout << "Using Single-Occurence: "<<endl;
		Neighbor.clear();
		MTX2.toNeighbor(Neighbor);
		for(map<unsigned,map<unsigned,unsigned> >::const_iterator itr1 = Neighbor.begin();itr1!=Neighbor.end();++itr1)
		{
			for(map<unsigned,unsigned>::const_iterator itr2 = itr1->second.begin();itr2!=itr1->second.end();++itr2)
			{
				int sidx;
				for(sidx =6;sidx>=0;--sidx)
				{
					if(Taxo[itr1->first][sidx]==Taxo[itr2->first][sidx])
						break;
				}
				if(sidx<0)sidx = 7;
				cout << (itr1->first) << '\t' << (itr2->first) << '\t' << sidx <<'\t'<<itr2->second << endl;
			}
		}
		/*
		for(int i=AllNodes.getNodeNum()-1;i>=0;--i)
		{
			int tsize = AllNodes.getVSize(i);
			NodeIDPosi* vec = AllNodes.getVector(i);
			map<int,int> comp;
			for(int j=0;j<tsize;++j)
				++comp[vec[j].id];
			if(comp.size()<=1)continue;
			for(map<int,int>::const_iterator itr1 = comp.begin();itr1!=comp.end();++itr1)
			{
				map<int,int>::const_iterator itr2 = itr1;
				for(++itr2;itr2!=comp.end();++itr2)
				{
					int sidx;
					for(sidx =6;sidx>=0;--sidx)
					{
						if(Taxo[itr1->first][sidx]==Taxo[itr2->first][sidx])
							break;
					}
					if(sidx<0)sidx = 7;
					shareMtxMulti[itr1->first][sidx] += itr1->second;
					shareMtxMulti[itr2->first][sidx] += itr2->second;
					++shareMtxSingle[itr1->first][sidx];
					++shareMtxSingle[itr2->first][sidx];
				}
			}
		}
		cout << "Using Multi-Occurence: "<<endl;
		for(int i=0;i<3000;++i)
		{
			int tsum = 0;
			for(int j=0;j<8;++j)
				tsum += shareMtxMulti[i][j];
			if(tsum==0)continue;
			cout<<dec << i<<'\t';
			for(int j=0;j<8;++j)
				cout << shareMtxMulti[i][j] << '\t';
			cout << endl;
		}
		cout << "Using Single-Occurence: "<<endl;
		for(int i=0;i<3000;++i)
		{
			int tsum = 0;
			for(int j=0;j<8;++j)
				tsum += shareMtxSingle[i][j];
			if(tsum==0)continue;
			cout << i<<'\t';
			for(int j=0;j<8;++j)
				cout << shareMtxSingle[i][j] << '\t';
			cout << endl;
		}*/
	}
private:
	GenomesClass(const GenomesClass& t){}
	const GenomesClass& operator=(const GenomesClass& t){return *this;}
};

#endif
