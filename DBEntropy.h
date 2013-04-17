/*
 * last fixed: 2013.04.08.
 * by Wang Yi.
 * 1. make DBKmerType a class
 * */
#ifndef MCH_HYBRID_DBENTROPY_H_

#define MCH_HYBRID_DBENTROPY_H_
#include <cstring>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <omp.h>
#include "PARAMETER.h"
#include "ShortKMER.h"
#include "Utils.h"
#include "SparseMatrix.h"
using namespace std;
////////////////////////////////////////////////////
////////////////////////////////////////////////////
typedef ShortKMER DBKmerType;
////////for same NodeAloc
struct DBKmerNode
{
	unsigned VSize,Capacity;
	DBKmerType kmer;
	unsigned* GeID;
	DBKmerNode* next;
	DBKmerNode(){}
	void set(const DBKmerType& kmer_,unsigned VSize_,unsigned Capacity_)
	{
		kmer = kmer_;
		VSize = VSize_;
		Capacity = Capacity_;
	}
	void naiveset(DBKmerType kmer_)
	{
		kmer = kmer_;
		next = NULL;
		VSize = 1;
	}
	void incVSize()
	{
		++VSize;
	}
	void push_back(unsigned sqid)
	{
		GeID[VSize++] = sqid;
	}
	void mallocVector()
	{
		assert(VSize>0);
		GeID = new unsigned[VSize];
		Capacity = VSize;
		VSize = 0;
	}
	void sortVector()
	{
		sort(GeID,GeID+VSize);
	}
	void Capacity2UniqueIDNum()
	{
		assert(Capacity==VSize);
		unsigned* pend = GeID + VSize;
		for(unsigned* itr = GeID+1;itr<pend;++itr)
			if(*itr == *(itr-1))
				--Capacity;
	}
	void passData(DBKmerType &kmer_,unsigned &VSize_, unsigned &Capacity_)const
	{
		kmer_ = kmer;
		VSize_ = VSize;
		Capacity_ = Capacity;
	}
	void clear()
	{
		VSize = Capacity = 0;
		if(GeID != NULL)
		{
			delete[] GeID;
			GeID = NULL;
		}
	}
	void moveFrom(DBKmerNode& node)
	{
		VSize = node.VSize;
		Capacity = node.Capacity;
		GeID = node.GeID;
		next = node.next;
		node.VSize = 0;
		node.next = NULL;
	}
};

//have to process memory myself.
class DBNodeAloc
{
public:
	const static int RowSizeBit = 20;
	const static int RowSize = 1U<<(RowSizeBit);
	const static int RowNum = 1U<<19;
	const static unsigned MASK = 0xfffff;
	DBKmerNode** AllNodes;
	omp_lock_t getnew_lock;

	DBNodeAloc()
	{
		NodeNum = 0;
		AllNodes = new DBKmerNode*[RowNum];
		for(unsigned i=0;i<RowNum;++i)
		{
			AllNodes[i] = NULL;
		}
		omp_init_lock(&getnew_lock);
	}
	void clear()
	{
		for(unsigned i=0;i<NodeNum;++i)
		{
			AllNodes[i>>RowSizeBit][i&MASK].clear();
		}
		for(unsigned i=0;i<RowNum;++i)
		{
			if(AllNodes[i] != NULL)
				delete[] AllNodes[i];
			AllNodes[i] = NULL;
		}
		delete[] AllNodes;
		NodeNum = 0;
		AllNodes = NULL;
	}
	virtual ~DBNodeAloc()
	{
		for(int i=0;i<NodeNum;++i)
			AllNodes[i>>RowSizeBit][i&MASK].clear();
		for(unsigned i=0;i<RowNum;++i)
			if(AllNodes[i] != NULL)
				delete[] AllNodes[i];
		delete[]AllNodes;
		AllNodes = NULL;
		omp_destroy_lock(&getnew_lock);
	}

	DBKmerNode* getNew(DBKmerType kmer_)
	{
		omp_set_lock(&getnew_lock);
		unsigned rowNo = NodeNum>>RowSizeBit;
		if(AllNodes[rowNo]==NULL)
			AllNodes[rowNo] = new DBKmerNode[RowSize];
		unsigned curNum = NodeNum++;
		omp_unset_lock(&getnew_lock);
		AllNodes[rowNo][curNum&MASK].naiveset(kmer_);
		return &AllNodes[rowNo][curNum&MASK];
	}

	unsigned long long getNodeNum() const
	{
		return NodeNum;
	}

	DBKmerNode* getRef(unsigned id) const
	{
		return &AllNodes[id>>RowSizeBit][id&MASK];
	}

	void shrinkSize(DBKmerNode** KmerMap, const int Thresh=2)//shrinkSize1(INodeRef* KmerMap)
	{
		unsigned long long srcid = NodeNum-1, tgtid = 0;//move src node to tgt(target) node
		while(tgtid < srcid && AllNodes[tgtid>>RowSizeBit][tgtid&MASK].VSize >= Thresh)
			++tgtid;
		while(tgtid < srcid && AllNodes[srcid>>RowSizeBit][srcid&MASK].VSize < Thresh)
			--srcid;
		while(srcid > tgtid)
		{
			if(KmerMap != NULL)
			{
				DBKmerNode* result = updateNext(KmerMap,getRef(srcid), getKmer(srcid).hash(), getRef(tgtid));
				assert(result!=NULL);
//				if(result==INULL)std::cerr<<"ERROR in shrinksize 1:\t"<<tgtid<<"\t"<<std::hex<<getKmer(srcid)<<std::dec<<"\t"<<srcid<<std::endl;
			}
			AllNodes[tgtid>>RowSizeBit][tgtid&MASK].moveFrom(AllNodes[srcid>>RowSizeBit][srcid&MASK]);
			--srcid;
			++tgtid;
			while(tgtid < srcid && AllNodes[tgtid>>RowSizeBit][tgtid&MASK].VSize >= Thresh)
				++tgtid;
			while(tgtid < srcid && AllNodes[srcid>>RowSizeBit][srcid&MASK].VSize < Thresh)
				--srcid;
		}
		unsigned oriNodeNum = NodeNum;
		const unsigned oldnum = (NodeNum-1)>>RowSizeBit;
		NodeNum=0;
		while(NodeNum<oriNodeNum && getVSize(NodeNum)>=Thresh)
			++NodeNum;
		for(unsigned i=((NodeNum-1)>>RowSizeBit)+1;i<=oldnum;++i)
		{
			delete[]AllNodes[i];
			AllNodes[i] = NULL;
		}
		std::cerr<<"deleted rows:\t"<<oldnum-((NodeNum-1)>>RowSizeBit)<<std::endl;
	}
	unsigned getVSize(unsigned id) const
	{
		return AllNodes[id>>RowSizeBit][id&MASK].VSize;
	}
	unsigned getCapacity(unsigned id) const
	{
		return AllNodes[id>>RowSizeBit][id&MASK].Capacity;
	}
	unsigned* getVector(unsigned id)const
	{
		return AllNodes[id>>RowSizeBit][id&MASK].GeID;
	}
	DBKmerType getKmer(unsigned id) const
	{
		return AllNodes[id>>RowSizeBit][id&MASK].kmer;
	}
	void mallocVector(unsigned id)
	{
		AllNodes[id>>RowSizeBit][id&MASK].mallocVector();
	}
	void sortVector(unsigned id)
	{
		AllNodes[id>>RowSizeBit][id&MASK].sortVector();
	}
	void Capacity2UniqueIDNum(unsigned id)
	{
		AllNodes[id>>RowSizeBit][id&MASK].Capacity2UniqueIDNum();
	}
private:
	inline DBKmerNode* updateNext(DBKmerNode** KmerMap, const DBKmerNode* next_,unsigned hashid, DBKmerNode* nextnew)
	{
		if(KmerMap[hashid] == next_)
		{
			KmerMap[hashid] = nextnew;
			return KmerMap[hashid];
		}
		for(DBKmerNode* curr = KmerMap[hashid];curr != NULL;curr = curr -> next)
		{
			if(curr->next == next_)
			{
				curr->next = nextnew;
				return curr;
			}
		}
		return NULL;
	}
	DBKmerNode* getNext(unsigned id)
	{
		return AllNodes[id>>RowSizeBit][id&MASK].next;
	}
	void setVSize(unsigned id,unsigned VSize_)
	{
		AllNodes[id>>RowSizeBit][id&MASK].VSize = VSize_;
	}
	void incVSize(unsigned id)
	{
		++AllNodes[id>>RowSizeBit][id&MASK].VSize;
	}
	unsigned setNext(unsigned id,DBKmerNode* next_)
	{
		AllNodes[id>>RowSizeBit][id&MASK].next = next_;
	}
	unsigned setKmer(unsigned id,DBKmerType kmer_)
	{
		AllNodes[id>>RowSizeBit][id&MASK].kmer = kmer_;
	}
	void push_back(unsigned id,unsigned ele)
	{
		AllNodes[id>>RowSizeBit][id&MASK].push_back(ele);
	}
	void clearVect(unsigned id)
	{
		AllNodes[id>>RowSizeBit][id&MASK].clear();
	}

	unsigned long long NodeNum;
	DBNodeAloc(const DBNodeAloc &uset){}
	const DBNodeAloc &operator=(const DBNodeAloc &uset){return *this;}
};//INodePool;

struct SimpleTaxoInfoClass
{
	int taxid;
	double score;
	void set(int id_,double score_)
	{
		taxid = id_;
		score = score_;
	}
};
struct ClustTaxoInfoClass
{
	SimpleTaxoInfoClass taxo[40];
	double alignscore[40];
};

class DBEntropy
{
private:
	int GenomeNum;
	static const int MAXLINE = 20000000;
	static const int MAXLENGTH = 20000000;
public:
	int **Taxo;
	int *GenomeLength;

	DBNodeAloc NodePool;
	DBKmerNode** KmerMap;
	DBEntropy(){KmerMap = NULL;}
	DBEntropy(string file)
	{
		init(file);
	}
	void init(string path)
	{
		ifstream ifs(path.c_str());
		if(ifs.fail())
		{
			std::cerr<<"File open failed: "<<path<<std::endl;
			Taxo = NULL;
			GenomeNum = 0;
			return;
		}
////////////////////////////////////////////////////////////////////
		KmerMap = new DBKmerNode*[DBHASHSIZE];
		for(unsigned i=0;i<DBHASHSIZE;++i)
			KmerMap[i] = NULL;
////////////////////////////////////////////////////////////////////
		char* gffbuf = new char[MAXLINE];
		char* genomebuf = new char[MAXLENGTH];
		const int TaxoBufMax = 10000;
		char Buf[TaxoBufMax];
		ifs.getline(Buf,TaxoBufMax);
		int gidx = 0;
		int taxo[7];

		GenomeNum = getFileLine(path.c_str())-1;
	/*	G1 = new BaseStr*[GenomeNum];
		G2 = new BaseStr*[GenomeNum];
		description = new string[GenomeNum];*/
		string* ctgstr= new string[GenomeNum];
		Taxo = new int*[GenomeNum];
		GenomeLength = new int[GenomeNum];
		for(int i=0;i<GenomeNum;++i)
			Taxo[i] = NULL;
		Taxo[0] = new int[7];
		int ttcount = 0;
		cerr << "before loading genomes." << endl;
		while(!ifs.eof())
		{
//			if(gidx > 10)break;
			ifs.getline(Buf,TaxoBufMax);
			int CurGid = atoi(Buf);
			char* inp = Buf;
			int tabc = 0;
			while(tabc<4 && (*inp))
			{
				if(*inp == '\t')++tabc;
				++inp;
			}
			int tmp;
			sscanf(inp,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",&taxo[0],&taxo[1],&taxo[2],&taxo[3],&taxo[4],&taxo[5],&tmp,&taxo[6]);
			///////////////////////
		//	if(ttcount++ < 1500)continue;
		//	int Forbid[]={399742,548,61645,-1,79967,182337,215689,338565,563,564};
		//	int Forbid[] = {235,234,236,233,232,772,771,775,773,190,191,189,194,195,192,193,790,784,782,779,780,783,794,781,32,45,37,36,34,1436,1437,1434,1450,1449,1464,1465,1466,70,67,68,69,1507,1505,1506,1503,1504,1508,1670,1671,1673,1674,1675,1676,1672,382,383,284,285,272,273,282,280,1328,1535,1536,396,1547,1559,1534,1544,1546,1545,1548,1568,1569,1542,1549,1543,1550,105,107,106,108,675,674,672,673,2032,2031,2030,2029,2028,1299,1300,833,18,19,832,831,838};
		//	For 100 species of Level 5 
		//	int Forbid[] = {772,771,775,191,1398,1090,788,782,37,1440,1449,1445,1455,1460,955,1666,683,693,1615,1660,2056,667,1386,1492,1356,1355,1367,662,586,582,1307,1966,2038,1087,854,1653,761,938,1474,2020,1364,22,266,935,1504,280,1536,1300,838,1078,1972,1359,946,1468,665,1688,1689,1382,1314,12,1475,1997,75,1167,104,1374,620,631,1870,1872,967,878,1899,926,1892,1339,1905,1874,1091,1326,1748,24,664,1850,1939,1274,26,1698,206,968,1623,1627,412,1669,2023,1936,1146,1388,1125,1142};
		//	For 100 species of Level 4
		//	int Forbid[] = {191,1398,1997,75,955,1666,788,782,37,1440,1449,1445,1455,1460,683,693,1615,582,1660,2056,1653,2020,1364,22,1504,773,232,193,1404,781,1466,753,962,2037,1652,2002,851,266,201,1672,280,1158,311,1494,1284,823,642,884,347,1731,838,1078,1520,1958,946,1468,665,1688,1689,1382,1945,12,1475,631,620,964,1101,1899,342,128,1091,1326,1748,24,932,1090,854,79,879,1474,1125,1142,531,858,1679,1669,2023,1274,26,1698,206,968,1623,1627,412,56,411,1874,1339,926};
		//	For 100 species of Level 4_2
			int Forbid[] = {1573,1398,1997,75,1490,1666,585,659,37,1995,1121,1862,863,873,683,693,1615,582,1352,2056,1332,2012,1364,22,1504,773,232,193,1404,781,1466,753,962,2037,1652,2002,851,266,201,1672,280,1158,311,1494,1284,823,642,884,347,1731,838,1078,1520,1958,946,1468,665,1688,1689,1382,1945,12,1475,631,620,964,1101,1899,342,128,1091,1326,1748,24,932,1090,854,79,879,1474,1125,1142,531,858,1679,1669,2023,1274,26,1698,206,968,1623,1627,412,56,411,1874,1339,926};
			bool isForbidden = false;
			for(int i=0;i<100;++i)
				if(CurGid==Forbid[i])
				{
					isForbidden = true;
					break;
				}
			if(isForbidden)
				continue;
			
		//	if(taxo[6]==292 || taxo[6]==13373)continue;
		//	if(taxo[5]==32008)continue;
		//	if(taxo[5]==1386)continue;
			//////////////////////////
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
		//		description[gidx] = gffbuf;
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
//					cerr << gidx << " has been processed. " << endl;
					ctgstr[gidx] = genomebuf;
					GenomeLength[gidx] = glength;
			/*		G1[gidx] = new BaseStr(genomebuf);
					G2[gidx] = new BaseStr(genomebuf,true);
					*/
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
		//////count occurence & hash k-mers.
		omp_lock_t* hash_lock = new omp_lock_t[1U<<20];
		for(int i=0;i<(1U<<20);++i)
			omp_init_lock(&hash_lock[i]);
		if(INTEST)
		{
			cout << "start to count kmer" << endl;
			system("date");
			system("ps ux");
		}
#pragma omp parallel for schedule(dynamic)
		for(int i=0;i<GenomeNum;++i)
			count_Occ_In_Str(hash_lock, ctgstr[i],i);
		cerr << "Finished counting occurences in strings. " << endl;
		if(INTEST)
		{
			cout << "end of counting kmer" << endl;
			system("date");
			system("ps ux");
			cout << endl;
		}
		////////////////////////////////////////////////////////
	/*	if(INTEST)
		{
			map<long long,long long>M;
			for(int i=0;i<DBHASHSIZE;++i)
			{
				unsigned ct = 0;
				for(DBKmerNode* p = KmerMap[i];p!=NULL;p=p->next)
					++ct;
				unsigned cbase = 1;
				while(cbase*10 < ct)cbase*=10;
				ct/=cbase;ct*=cbase;
				++M[ct];
			}
			cerr << "Hash table list-length statistics: " <<endl;
			for(map<int,int>::const_iterator itr=M.begin();itr!=M.end();++itr)
				cerr << itr->first << '\t' << itr->second << endl;

			M.clear();
			for(int i=0;i<DBHASHSIZE;++i)
			{
				for(DBKmerNode* p = KmerMap[i];p!=NULL;p=p->next)
				{
				unsigned ct = p->VSize;
				unsigned cbase = 1;
				while(cbase*10 < ct)cbase*=10;
				ct/=cbase;ct*=cbase;
				++M[ct];
				}
			}
			cerr << "Hash table array-length statistics: " <<endl;
			for(map<int,int>::const_iterator itr=M.begin();itr!=M.end();++itr)
				cerr << itr->first << '\t' << itr->second << endl;
			exit(-1);
		}*/
		////////////////////////////////////////////////////////
	/*	for(int i=0;i<NodePool.getNodeNum();++i)
			NodePool.mallocVector(i);
		*/
		for(int i=0;i<DBHASHSIZE;++i)
			for(DBKmerNode* p = KmerMap[i];p!=NULL;p=p->next)
				p->mallocVector();
		cerr << "Finished mallocing vectors in nodes. " << endl;
		//////////////////////////////////////////////////////////
		system("ps ux");
#pragma omp parallel for schedule(dynamic)
		for(int i=0;i<GenomeNum;++i)
			mapto(hash_lock, ctgstr[i], i);
		cerr << "Finished mapping k-mers in strings. " << endl;
#pragma omp parallel for schedule(dynamic)
		for(int i=0;i<NodePool.getNodeNum();++i)
			NodePool.sortVector(i);
		cerr << "Finished sorting vectors in nodes. " << endl;
/////////////////////////////////////////////////////////////////////////
		if(INTEST)
		{
		for(int i=0;i<NodePool.getNodeNum();++i)
			if(NodePool.getCapacity(i) != NodePool.getVSize(i))
				assert(false);
		}
//////////////	/////////////////////////////////////////////////////////////////////

#pragma omp parallel for schedule(dynamic)
		for(int i=0;i<NodePool.getNodeNum();++i)
			NodePool.Capacity2UniqueIDNum(i);
		cerr << "Finished turning capacity to number of unique id num. " << endl;
		for(int i=0;i<(1U<<20);++i)
			omp_destroy_lock(&hash_lock[i]);
		delete[] hash_lock;
		delete[]ctgstr;
		cerr << "Genome DB Initialization is finished."<<endl;
	}
	void getStatis()
	{
		cerr << "Start to get statistics. " << endl;
		cout << "hash vector size:\t" << endl;
		map<int,int>vecSize;
		for(int i=NodePool.getNodeNum()-1;i>=0;--i)
		{
			int ct = NodePool.getVSize(i);
			if(ct > 10000)ct = (ct/10000)*10000;
			else if(ct > 1000)ct = (ct/1000)*1000;
			else if(ct > 100)ct = (ct/100)*100;
			else if(ct > 10)ct = (ct/10)*10;
			++vecSize[ct];
			if(ct > 1000)
				cout <<hex << NodePool.getKmer(i).U64 << '\t' << dec << ct << endl;
		}
		for(map<int,int>::const_iterator itr = vecSize.begin();itr!=vecSize.end();++itr)
			cout << (itr->first) << '\t' << itr->second << endl;

		cout << "hash clision size:\t" << endl;
		vecSize.clear();
		for(int i=0;i<DBHASHSIZE;++i)
		{
			int ct = 0;
			for(DBKmerNode* p = KmerMap[i];p!=NULL;p=p->next)
				++ct;
			++vecSize[ct];
		}
		for(map<int,int>::const_iterator itr = vecSize.begin();itr!=vecSize.end();++itr)
			cout << (itr->first) << '\t' << itr->second << endl;
	}
	//////////////////////////////////////////////////////////////////////
//	void calTaxoForCtg(const string& str, unsigned& taxid, unsigned& taxLevel,double& score)
	void calTaxoForCtg(const string& str, map<int,double>& sp_score)
	{
//		map<unsigned, double>genoScore;
		sp_score.clear();
		DBKmerType kmer0(str,false,true);
		DBKmerType kmer1(str,true,true);
		DBKmerType kmerR;
		int tlength = str.length();
		int TotalScore = 0;
		//////////////////////////////////////////////////////////////////////////
		map<DBKmerType, unsigned> kmerInCtg;
		for(int i=PARA_KMER-1;i<tlength;++i)
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
				++kmerInCtg[kmer0];
			else
				++kmerInCtg[kmer1];
		}
		for(map<DBKmerType,unsigned>::const_iterator itr=kmerInCtg.begin();itr!=kmerInCtg.end();++itr)
		{
			DBKmerNode* p = findKmer(itr->first);
			if(p==NULL || p->VSize <= 0)continue;
			unsigned* pend = (p->GeID) + (p->VSize);
			set<unsigned>candsp;
			for(unsigned* itr2 = (p->GeID); itr2<pend;)
			{
				unsigned* itst = itr2++;
				while(itr2 < pend && (*itst) == (*itr2))
					++itr2;
				candsp.insert(Taxo[*itst][6]);
			}
			double unit = 1.0/candsp.size();
			for(set<unsigned>::const_iterator itr2=candsp.begin();itr2!=candsp.end();++itr2)
				sp_score[*itr2] += unit;

/*			double unit = 1.0/(p->VSize);
			for(unsigned* itr2 = (p->GeID); itr2<pend;)
			{
				unsigned* itst = itr2++;unsigned num = 1;
				while(itr2 < pend && (*itst) == (*itr2))
					++itr2,++num;
				assert(num==(itr2-itst));
				sp_score[Taxo[*itst][6]] += unit*min(num,itr->second);
				/*to be tested:  only consider occurence*/
/*			}*/
		}
		//////////////////////////////////////////////////////////////////////////
		/*
		if(TotalScore==0)
		{
			maxScore = -1;
			return;
		}
		entropy.clear();
		entropy.resize(7);
		/////////////////////////
		for(map<unsigned,double>::iterator itr = genoScore.begin();itr != genoScore.end();++itr)
			(itr->second) *= (3e6/GenomeLength[itr->first]);
		//	(itr->second) *= (9e12/GenomeLength[itr->first]/GenomeLength[itr->first]);
		simulateId = 0;simulateScore = 0;
		maxId = 0;maxScore = 0;
		for(map<unsigned,double>::const_iterator itr = genoScore.begin();itr != genoScore.end();++itr)
		{
			if(itr->second > maxScore){
				maxScore = itr->second;
				maxId = itr->first;
			}
			if(Taxo[itr->first][4]==119060 && itr->second > simulateScore)
			{
				simulateScore = itr->second;
				simulateId = itr->first;
			}
		}
		//////////////////////////////////
		vector<unsigned>maxTaxo(7);
		for(int i=0;i<7;++i)
			maxTaxo[i] = Taxo[maxId][i];
		////////////////////////////////////////////////////////
		map<unsigned,pair<unsigned, double> > Lev[2];
		for(map<unsigned,double>::const_iterator itr = genoScore.begin();itr != genoScore.end();++itr)
			if(itr->second > (Lev[0][itr->first].second))
				(Lev[0][itr->first].second) = itr->second;
		for(map<unsigned,pair<unsigned, double> >::iterator itr=Lev[0].begin();itr!=Lev[0].end();++itr)
			(itr->second).first = itr->first;

		for(int i=6;i>=0;--i)
		{
			int idx1 = i&0x1, idx2=1-idx1;
			Lev[idx2].clear();
			for(map<unsigned,pair<unsigned,double> >::const_iterator itr = Lev[idx1].begin();itr != Lev[idx1].end();++itr)
			{
				unsigned pa = Taxo[(itr->second).first][i];
				if((Lev[idx2][pa].second) < (itr->second.second))
					Lev[idx2][pa] = itr->second;
			}
			vector<double>tscores;
			for(map<unsigned,pair<unsigned,double> >::const_iterator itr = Lev[idx1].begin();itr != Lev[idx1].end();++itr)
				if(Taxo[itr->second.first][i]==maxTaxo[i])
					tscores.push_back(itr->second.second);
			double tsum = 0;
			for(vector<double>::const_iterator itr=tscores.begin();itr!=tscores.end();++itr)
				tsum += *itr;
			double tentropy = 0;
			for(vector<double>::const_iterator itr=tscores.begin();itr!=tscores.end();++itr)
				tentropy += (*itr/tsum)*log(*itr/tsum);
			entropy[i] = tentropy/log(2);
		}*/
	}
	int getGenomeNum()
	{
		return GenomeNum;
	}
	void clearDB()
	{
		if(KmerMap!=NULL)
		{
			for(int i=0;i<DBHASHSIZE;++i)
				if(KmerMap[i]!=NULL)
					delete[]KmerMap[i];
			delete[]KmerMap;
			KmerMap = NULL;
		}
		NodePool.clear();
		if(GenomeLength != NULL)
		{
			delete[] GenomeLength;
			GenomeLength = 0;
		}
	}
	void clear()
	{
		clearDB();
		if(Taxo!=NULL)
		{
			for(int i=0;i<GenomeNum;++i)
				delete[]Taxo[i];
			delete[]Taxo;
			Taxo = NULL;
		}
	}

private:
	void count_Occ_Hash(const DBKmerType& kmer_, const unsigned hashidx)
	{
		unsigned idx = hashidx;
		if(KmerMap[idx] == NULL)
			KmerMap[idx] = NodePool.getNew(kmer_);
		else
		{
			DBKmerNode* p = KmerMap[idx];
			for(;p->next!=NULL;p=p->next)
				if(p->kmer == kmer_)
				{
					p->incVSize();
					return;
				}
			if(p->kmer == kmer_)
				p->incVSize();
			else
				p->next = NodePool.getNew(kmer_);
		}
	}
	void insertHash(const DBKmerType& kmer_, const unsigned hashidx, int id)
	{
		assert(KmerMap[hashidx]!=NULL);
		DBKmerNode* p = KmerMap[hashidx];
		for(;p != NULL;p = p->next)
			if(p->kmer == kmer_)
			{
				p->push_back(id);
				return;
			}
		assert(false);
	}
	DBKmerNode* findKmer(const DBKmerType & kmer_)
	{
		DBKmerNode* ans = KmerMap[kmer_.hash()];
		for(;ans!=NULL;ans = ans->next)
			if(ans->kmer == kmer_)
				return ans;
		return NULL;
	}
	void count_Occ_In_Str(omp_lock_t* hash_lock,const string& str, int idx)
	{
		assert(KmerMap !=NULL);
		DBKmerType kmer0(str,false,true);
		DBKmerType kmer1(str,true,true);
		int tlength = str.length()-1;
		for(int i=PARA_KMER-1;i<str.length();++i)
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
				count_Occ_Hash(kmer0, hashidx);
				omp_unset_lock(&hash_lock[hashidx&0xfffff]);
			}
			else
			{
				unsigned hashidx = kmer1.hash();
				omp_set_lock(&hash_lock[hashidx&0xfffff]);
				count_Occ_Hash(kmer1, hashidx);
				omp_unset_lock(&hash_lock[hashidx&0xfffff]);
			}
		}
	}
	void mapto(omp_lock_t* hash_lock,const string& str, int idx)
	{
		assert(KmerMap !=NULL);
		DBKmerType kmer0(str,false,true);
		DBKmerType kmer1(str,true,true);
		int tlength = str.length()-1;
		for(int i=PARA_KMER-1;i<str.length();++i)
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
				insertHash(kmer0,hashidx, idx);
				omp_unset_lock(&hash_lock[hashidx&0xfffff]);
			}
			else
			{
				unsigned hashidx = kmer1.hash();
				omp_set_lock(&hash_lock[hashidx&0xfffff]);
				insertHash(kmer1,hashidx, idx);
				omp_unset_lock(&hash_lock[hashidx&0xfffff]);
			}
		}
	}
};
#endif
