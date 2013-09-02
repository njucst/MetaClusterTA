/*
 * Reason: small arrays cannot be recycled efficiently by glibc. we need to manamge them manually.
 * Note: better not to use it when L is large. If there's no choice, SmallArrayRowSize should be smaller for large L.
 * by ywang@hku
 * 2013.05.05
 */
#ifndef MCH_HYBRID_SMALLARRAY_H_

#define MCH_HYBRID_SMALLARRAY_H_
#include <omp.h>
#include <iostream>
using namespace std;
template <class T, int L>
class SmallArray
{
public:
	SmallArray()
	{
		AllNodes = new T*[RowNum];
		for(unsigned i=0;i<RowNum;++i)
			AllNodes[i] = NULL;
		omp_init_lock(&getnew_lock);
		NodeNum = 0;
	}
	T* getNew()
	{
		omp_set_lock(&getnew_lock);
		unsigned rowNo = NodeNum>>RowSizeBit;
		if(AllNodes[rowNo]==NULL)
			AllNodes[rowNo] = new T[RowSize*L];
		unsigned curNum = NodeNum++;
		omp_unset_lock(&getnew_lock);
		return &AllNodes[rowNo][(curNum&MASK)*L];
	}
	void test()
	{
		unsigned rowNo = NodeNum>>RowSizeBit;
		for(int i=0;i<rowNo;++i)
			for(int j=0;j<(RowSize*L);++j)
				AllNodes[i][j] = i*(RowSize*L)+j;
	}
	bool isEmpty()
	{
		return NULL==AllNodes;
	}
	void clear()
	{
		if(NULL==AllNodes)return;
		for(unsigned i=0;i<RowNum;++i)
			if(AllNodes[i] != NULL)
				delete[] AllNodes[i];
		delete[]AllNodes;
		AllNodes = NULL;
	}
	virtual ~SmallArray()
	{
		omp_destroy_lock(&getnew_lock);
		if(NULL==AllNodes)return;
		for(unsigned i=0;i<RowNum;++i)
			if(AllNodes[i] != NULL)
				delete[] AllNodes[i];
		delete[]AllNodes;
		AllNodes = NULL;
	}
private:
	const static int RowSizeBit = 16;
	const static int RowSize = 1U<<(RowSizeBit);
	const static int RowNum = 1U<<20;
	const static unsigned MASK = 0xffff;
	T **AllNodes;
	unsigned long long NodeNum;
	omp_lock_t getnew_lock;

	SmallArray(const SmallArray &uset){}
	const SmallArray &operator=(const SmallArray &uset){return *this;}
};
#endif
