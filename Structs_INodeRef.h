/*
 * last fixed: 2012.11.10.
 * by Wang Yi.
 * */
#ifndef MCH_HYBRID_STRUCTS_H_

#define MCH_HYBRID_STRUCTS_H_

#include <cmath>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <omp.h>
#include "KMER.h"

typedef unsigned long long INodeRef;
const unsigned long long INULL = ~0ULL;
typedef KMER KmerType;

//const int TRANLEN = 25;
//const unsigned long long TRANMASK = 0x3ffffffffffff;
//const unsigned HASHMASK = 0x3fffffff;
//const unsigned HASHSIZE = 1U<<30;

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////functions

////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////structure & function
//#pragma pack(1)
struct NodeIDPosi
{
	unsigned id;
	unsigned posi;
	bool operator<(const NodeIDPosi& t)const
	{
		return posi<t.posi;
	}
};

#pragma pack()
struct KmerNode
{
	unsigned VSize;
	INodeRef next;
	KmerType kmer;
	NodeIDPosi* myvector;
	KmerNode(KmerType kmer_,unsigned VSize_,INodeRef next_)
	{
		kmer = kmer_;
		VSize = VSize_;
		next = next_;
		myvector = NULL;
	}
	KmerNode()
	{
		VSize = 0;
		myvector = NULL;
		next = INULL;
	}
	void push_back(unsigned ele,unsigned posi)
	{
		myvector[VSize].id=ele;
		myvector[VSize].posi=posi;
		++VSize;
	}
	void naiveset(KmerType kmer_,unsigned VSize_)
	{
		kmer = kmer_;
		VSize = VSize_;
		myvector = NULL;
		next = INULL;
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

		node.VSize = 0;
		node.myvector = NULL;
	}
private:
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
			AllNodes[i] = NULL;
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

	INodeRef getNew(KmerType kmer_,unsigned VSize_)
	{
		omp_set_lock(&getnew_lock);
		unsigned rowNo = NodeNum>>RowSizeBit;
		if(AllNodes[rowNo]==NULL)
			AllNodes[rowNo] = new KmerNode[RowSize];
		unsigned curNum = NodeNum++;
		omp_unset_lock(&getnew_lock);
		AllNodes[rowNo][curNum&MASK].naiveset(kmer_,VSize_);
		return curNum;
	}

	unsigned long long getNodeNum() const
	{
		return NodeNum;
	}
	unsigned getVSize(unsigned id)
	{
		return AllNodes[id>>RowSizeBit][id&MASK].VSize;
	}
	INodeRef getNext(unsigned id)
	{
		return AllNodes[id>>RowSizeBit][id&MASK].next;
	}
	KmerType getKmer(unsigned id)
	{
		return AllNodes[id>>RowSizeBit][id&MASK].kmer;
	}
	NodeIDPosi* getVector(unsigned id)
	{
		return AllNodes[id>>RowSizeBit][id&MASK].myvector;
	}

	void setVSize(unsigned id,unsigned VSize_)
	{
		AllNodes[id>>RowSizeBit][id&MASK].VSize = VSize_;
	}
	void incVSize(unsigned id)
	{
		++AllNodes[id>>RowSizeBit][id&MASK].VSize;
	}
	unsigned setNext(unsigned id,INodeRef next_)
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
	KmerNode* getRef(unsigned id)
	{
		return &AllNodes[id>>RowSizeBit][id&MASK];
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

	inline INodeRef updateNext(INodeRef* KmerMap,const INodeRef next_,unsigned hashid, INodeRef nextnew)
	{
		if(KmerMap[hashid] == next_)
		{
			KmerMap[hashid] = nextnew;
			return hashid;
		}
		for(INodeRef curr = KmerMap[hashid];curr != INULL;curr = getNext(curr))
		{
			if(getNext(curr)==next_)
			{
				AllNodes[curr>>RowSizeBit][curr&MASK].next = nextnew;
				return curr;
			}
		}
		return INULL;
	}
	void shrinkSize(INodeRef* KmerMap, const int Thresh)//shrinkSize1(INodeRef* KmerMap)
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
				unsigned result = updateNext(KmerMap, srcid, getKmer(srcid).hash(), tgtid);
				assert(result==INULL);
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
		std::cout<<"deleted rows:\t"<<oldnum-((NodeNum-1)>>RowSizeBit)<<std::endl;
	}
private:
	unsigned long long NodeNum;
	KmerNodeAloc(const KmerNodeAloc &uset){}
	const KmerNodeAloc &operator=(const KmerNodeAloc &uset){return *this;}
};//INodePool;
////////////////////////////////////////////////////
////////for same NodeAloc
struct SameTNode
{
	unsigned count,type;
	KmerType kmer;
	SameTNode(){}
	void set(KmerType kmer_,unsigned count_,unsigned type_)
	{
		kmer = kmer_;
		count = count_;
		type = type_;
	}
	void getData(KmerType &kmer_, unsigned &count_, unsigned &type_)
	{
		kmer_ = kmer;
		count_ = count;
		type_ = type;
	}
};
class SameTNAloc
{
public:
	const static int RowSize = 1U<<20;
	const static int RowNum = 1U<<12;
	const static int RowSizeBit = 20;
	const static unsigned MASK = 0xfffff;
	SameTNode** NodeArray;
	omp_lock_t getnew_lock;

	SameTNAloc()
	{
		NodeNum = 0;
		NodeArray = new SameTNode*[RowNum];
		for(unsigned i=0;i<RowNum;++i)
		{
			NodeArray[i] = NULL;
		}
		omp_init_lock(&getnew_lock);
	}
	/*
	virtual ~SameTNAloc()
	{
		clear();
		omp_destroy_lock(&getnew_lock);
	}
	*/

	unsigned getNew(KmerType kmer_, unsigned count_, unsigned type_)
	{
		omp_set_lock(&getnew_lock);
		unsigned rowNo = NodeNum>>RowSizeBit;
		if(NodeArray[rowNo]==NULL)
			NodeArray[rowNo] = new SameTNode[RowSize];
		unsigned curNum = NodeNum++;
		omp_unset_lock(&getnew_lock);
		NodeArray[rowNo][curNum&MASK].set(kmer_, count_, type_);
		return curNum;
	}
	unsigned getNodeNum()
	{
		return NodeNum;
	}
	/*
	unsigned getTest()
	{
		return test;
	}
	*/

	unsigned getType(unsigned id)
	{
		return NodeArray[id>>RowSizeBit][id&MASK].type;
	}
	unsigned getCount(unsigned id)
	{
		return NodeArray[id>>RowSizeBit][id&MASK].count;
	}
	KmerType getKmer(unsigned id)
	{
		return NodeArray[id>>RowSizeBit][id&MASK].kmer;
	}
	void getData(unsigned id,KmerType &kmer_, unsigned &count_, unsigned &type_)
	{
		NodeArray[id>>RowSizeBit][id&MASK].getData(kmer_,count_,type_);
	}
	void clear()
	{
		for(unsigned i=0;i<RowNum;++i)
		{
			if(NodeArray[i]!=NULL)
				delete[]NodeArray[i];
			NodeArray[i] = NULL;
		}
		delete[]NodeArray;
	}
private:
	unsigned NodeNum;
//	unsigned test;
	SameTNAloc(const SameTNAloc &uset){}
	const SameTNAloc &operator=(const SameTNAloc &uset){return *this;}
};
#endif
