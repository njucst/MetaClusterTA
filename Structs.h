/*
 * last fixed: 2012.11.10.
 * by Wang Yi.
 * */
#ifndef MCH_HYBRID_STRUCTS_H_

#define MCH_HYBRID_STRUCTS_H_

#include <cmath>
#include <cassert>
#include <vector>
#include <iostream>
#include <algorithm>
#include <omp.h>
#include "KMER.h"

//typedef unsigned long long INodeRef;
//const unsigned long long INULL = ~0ULL;
typedef KMER KmerType;

//const int TRANLEN = 25;
//const unsigned long long TRANMASK = 0x3ffffffffffff;
//const unsigned HASHMASK = 0x3fffffff;
//const unsigned HASHSIZE = 1U<<30;

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////functions

////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////structure & function
struct TaxoInfo
{
	int gid;
	int max,sum;
};
struct SimpleTaxoInfoClass
{
	int taxid;
	double score;
	SimpleTaxoInfoClass():taxid(0),score(0){}
	void set(int id_,double score_)
	{
		taxid = id_;
		score = score_;
	}
};
struct ClustTaxoInfoClass
{
	vector<SimpleTaxoInfoClass> taxo;
	vector<double> alignscore;
	ClustTaxoInfoClass()
	{
		taxo.resize(40);
		alignscore.resize(40);
	}
};
////////////////////////////////////////////////////////////////////////
//#pragma pack(1)
struct NodeIDPosi
{
	unsigned id;
	unsigned posi;//bit 0: isRev;bit 1: isCtg;
	void set(unsigned id_, unsigned posi_)
	{
		id = id_;
		posi = posi_;
	}
	bool isRev()
	{
		return posi & 0x1;
	}
	bool isCtg()
	{
		return posi & 0x2;
	}
	bool operator<(const NodeIDPosi& t)const
	{
		return posi<t.posi;
	}
};

#pragma pack()
class KmerNode
{
public:
	friend class KmerNodeAloc;
	unsigned VSize,Capacity,CtgNum;
	KmerNode* next;
	KmerType kmer;
	NodeIDPosi* myvector;
	~KmerNode()
	{
		if(myvector)
			delete[] myvector;
	}
	void push_back(unsigned ele,unsigned posi)
	{
		assert(Capacity>0);
		if(VSize==Capacity)
		{
			Capacity <<= 1;
			NodeIDPosi* ori = myvector;
			myvector = new NodeIDPosi[Capacity];
			copy(ori,ori+VSize,myvector);
			delete[]ori;
		}
		myvector[VSize].id=ele;
		myvector[VSize].posi=posi;
		++VSize;
	}
	void naiveset(KmerType kmer_,unsigned id, unsigned posi)
	{
		kmer = kmer_;
		VSize = 1;
		Capacity = 2;
		CtgNum=0;
		myvector = new NodeIDPosi[2];
		myvector[0].set(id,posi);
		next = NULL;
	}
	void mallocVector()
	{
		myvector = new NodeIDPosi[VSize];
	}

	void clear()
	{
		VSize = 0;
		if(myvector!=NULL)
		{
			delete[]myvector;
			myvector = NULL;
		}
	}
	void moveFrom(KmerNode& node)
	{
		myvector = node.myvector;
		next = node.next;
		kmer = node.kmer;
		VSize = node.VSize;
		Capacity = node.Capacity;
		CtgNum = node.CtgNum;

		node.VSize = 0;
		node.myvector = NULL;
	}
private:
/*	KmerNode(KmerType kmer_,unsigned VSize_,KmerNode* next_)
	{
		kmer = kmer_;
		VSize = VSize_;
		next = next_;
		myvector = NULL;
	}*/
	KmerNode()
	{
		VSize = 0;
		Capacity = 0;
		CtgNum=0;
		myvector = NULL;
		next = NULL;
	}
	void fixCtgNum()
	{
		CtgNum = VSize;
	}
	KmerNode(const KmerNode &node){}
	const KmerNode &operator=(const KmerNode &node){return *this;}
};

//have to process memory myself.
class KmerNodeAloc
{
public:
	const static int RowSizeBit = 20;
	const static int RowSize = 1U<<(RowSizeBit);
	const static int RowNum = 1U<<19;
	const static unsigned MASK = 0xfffff;
	KmerNode** AllNodes;
	omp_lock_t getnew_lock;

	KmerNodeAloc()
	{
		NodeNum = 0;
		AllNodes = new KmerNode*[RowNum];
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
		}
		delete[] AllNodes;
		NodeNum = 0;
		AllNodes = NULL;
	}
/*	void reIni()
	{
		for(unsigned i=0;i<NodeNum;++i)
			AllNodes[i>>RowSizeBit][i&MASK].clear();
		for(unsigned i=0;i<RowNum;++i)
		{
			if(AllNodes[i] != NULL)
				delete[] AllNodes[i];
			AllNodes[i] = NULL;
		}
		NodeNum = 0;
	}*/
	virtual ~KmerNodeAloc()
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

	KmerNode* getNew(KmerType kmer_,unsigned id_,unsigned posi_)
	{
		omp_set_lock(&getnew_lock);
		unsigned rowNo = NodeNum>>RowSizeBit;
		if(AllNodes[rowNo]==NULL)
			AllNodes[rowNo] = new KmerNode[RowSize];
		unsigned curNum = NodeNum++;
		omp_unset_lock(&getnew_lock);
		AllNodes[rowNo][curNum&MASK].naiveset(kmer_,id_,posi_);
		return &AllNodes[rowNo][curNum&MASK];
	}
	void fixCtgNum()
	{
		for(int i=0;i<NodeNum;++i)
			AllNodes[i>>RowSizeBit][i&MASK].fixCtgNum();
	}

	unsigned long long getNodeNum() const
	{
		return NodeNum;
	}

	KmerNode* getRef(unsigned id) const
	{
		return &AllNodes[id>>RowSizeBit][id&MASK];
	}

	void shrinkSize(KmerNode** KmerMap, const int Thresh=2)//shrinkSize1(INodeRef* KmerMap)
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
				KmerNode* result = updateNext(KmerMap, getRef(srcid), getKmer(srcid).hash(), getRef(tgtid));
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
	NodeIDPosi* getVector(unsigned id)const
	{
		return AllNodes[id>>RowSizeBit][id&MASK].myvector;
	}
	KmerType getKmer(unsigned id) const
	{
		return AllNodes[id>>RowSizeBit][id&MASK].kmer;
	}
private:
	inline KmerNode* updateNext(KmerNode** KmerMap,const KmerNode* next_,unsigned hashid, KmerNode* nextnew)
	{
		if(KmerMap[hashid] == next_)
		{
			KmerMap[hashid] = nextnew;
			return KmerMap[hashid];
		}
		for(KmerNode* curr = KmerMap[hashid];curr != NULL;curr = curr -> next)
		{
			if(curr->next == next_)
			{
				curr->next = nextnew;
				return curr;
			}
		}
		return NULL;
	}
	KmerNode* getNext(unsigned id)
	{
		return &AllNodes[id>>RowSizeBit][id&MASK];
	}

	void setVSize(unsigned id,unsigned VSize_)
	{
		AllNodes[id>>RowSizeBit][id&MASK].VSize = VSize_;
	}
	void incVSize(unsigned id)
	{
		++AllNodes[id>>RowSizeBit][id&MASK].VSize;
	}
	unsigned setNext(unsigned id,KmerNode* next_)
	{
		AllNodes[id>>RowSizeBit][id&MASK].next = next_;
	}
	unsigned setKmer(unsigned id,KmerType kmer_)
	{
		AllNodes[id>>RowSizeBit][id&MASK].kmer = kmer_;
	}
	void mallocVector(unsigned id)
	{
		AllNodes[id>>RowSizeBit][id&MASK].mallocVector();
	}
	void push_back(unsigned id,unsigned ele,unsigned posi)
	{
		AllNodes[id>>RowSizeBit][id&MASK].push_back(ele,posi);
	}
	void clearVect(unsigned id)
	{
//		delete[] AllNodes[id>>RowSizeBit][id&MASK].myvector;
		AllNodes[id>>RowSizeBit][id&MASK].clear();
	}

	unsigned long long NodeNum;
	KmerNodeAloc(const KmerNodeAloc &uset){}
	const KmerNodeAloc &operator=(const KmerNodeAloc &uset){return *this;}
};//INodePool;
////////////////////////////////////////////////////
#endif
