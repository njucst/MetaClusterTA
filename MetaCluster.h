#ifndef __MetaCluster_H_

#define __MetaCluster_H_
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include "TomAlgorithm.h"

using namespace std;
class MetaCluster
{
public:
//	explicit MetaCluster(int KmerLen_,int size_,int**spear_,int**KmerDistri_,int GenoNum_,int** Component_,int kmeansize,int MaxSpecies_,int MinSpecies_)
	explicit MetaCluster(int KmerLen_,int size_,int**ReverKmerDistri_,int kmeansize,int MaxSpecies_=0,int MinSpecies_=2,int GenoNum_=0,int**Component_=NULL)
	{
//////////////////////////////////////////////////////
		KmerLen = KmerLen_;
		ReverPosition = getReverPosition(KmerLen);
		ReverSize = tomReverSize(KmerLen);
//////////////////////////////////////////////////////
		Size = size_;
		KmerDistri = ReverKmerDistri_;
		CtgLen = new int[Size];

		MaxSpecies = MaxSpecies_;
		MinSpecies = MinSpecies_;

		cerr << "Size:\t" << Size << endl;
		cerr << "ReverSize:\t" << ReverSize << endl;

			int sumlen = 0;
			for(int i=0;i<Size;++i)
			{
				int tlen = 0;
				for(int j=0;j<ReverSize;++j)
					tlen += KmerDistri[i][j];
				CtgLen[i] = tlen;
				sumlen += tlen;
			}
			cerr << "sumlen:\t" << sumlen << endl;

		if(MaxSpecies==0)
		{
			MaxSpecies = max(sumlen/500000,2);
			MinSpecies = max(sumlen/3000000,2);
	//		MaxSpecies = max(MinSpecies,MaxSpecies);
		}
		cerr << "MaxSpecies:\t" << MaxSpecies << endl;
		cerr << "MinSpecies:\t" << MinSpecies << endl;

		GenoNum = GenoNum_;
		if(GenoNum > 0)
		{
			Component = new int*[Size];
			for(int i=0;i<Size;++i)
			{
				Component[i] = new int[GenoNum];
				for(int j=0;j<GenoNum;++j)
					Component[i][j] = Component_[i][j];
			}
		}
		else Component = NULL;

//		int sixall=0;for(int i=0;i<Size;++i)sixall+=SixTMerSize[i];
		//classes = Size/100;
		if(Size < MaxSpecies)
		{
			MaxSpecies = Size/2;
			cerr << "MaxSpecies is too large. We will start with half of the group number." << endl;
		}
			
		if(kmeansize==0)
			classes = MaxSpecies;//sixall/600000+1;
		else classes = kmeansize;

//		cerr<<"Size:"<<Size<<",\t"<<sixall<<endl;
//		cerr<<"classes:"<<classes<<",\t"<<sixall/1200000<<endl;
		cerr<<"Size:"<<Size<<",\tclasses:"<<classes<<endl;
		cerr<<"KmerDistri last:"<<KmerDistri[0][0]<<"\t"<<KmerDistri[Size-1][2]<<endl;//<<",\t"<<SixTMerSize[Size-1]<<endl;
		distSum = new int[classes];
		distNum = new int[classes];
		distSumSquare = new int[classes];

		distMean = new int[classes];
		distSD = new int[classes];
//////////////////////////////////////////////////////
		isOutlier = new bool[Size];
		type = new int[Size];
		best = new int[Size];
		aux = new int[Size];
//////////////////////////////////////////////////////
		SpearRank = new int*[Size];
		for(int i=0;i<Size;++i)
			SpearRank[i] = toSpear(KmerDistri[i],ReverSize);
//			SpearRank[i] = toSpear(tomNormalize_rever(KmerDistri[i],ReverPosition,KmerLen,false),ReverSize);
//////////////////////////////////////////////////////
		centerDistri = new int*[classes];
		for(int i=0;i<classes;++i)
		{
			centerDistri[i] = new int[ReverSize];
		}
		centerRank = new int*[classes];
		for(int i=0;i<classes;++i)
		{
			centerRank[i] = new int[ReverSize];
		}
	}
	virtual ~MetaCluster(){}
////////////////////////////////////////////////////////
	double muiltkmeans(const int ROUND,int classes_);
	int MergeClusters(double shre_t);
	int iterMeta(const int ROUND,double shre_t);
///////////////////////////////////////////////////////
	int** getComp()const{return Component;}
	int getGenoNum()const{return GenoNum;}
	
	int Size;
	int classes;
	int **SpearRank;
	int **KmerDistri;
//	int *SixTMerSize;
///////////////////////////////////////////////////////
	int *type;
	int *best;
	bool *isOutlier;
private:
	int *CtgLen;
	//static const int KmerLen = 4;
	//static const int DisSize = 1<<(KmerLen<<1);
	int KmerLen;
//	int DisSize;
	int ReverSize;

	int MaxSpecies;
	int MinSpecies;

	int *distSum;
	int *distNum;
	int *distSumSquare;

	int *distMean;
	int *distSD;

	int **centerRank;
	int **centerDistri;
//////////////////////////////////////////////////////
	int* ReverPosition;
//////////////////////////////////////////////////////
	int *aux;
//////////////////////////////////////////////////////
	int GenoNum;
	int** Component;
//////////////////////////////////////////////////////
	static const int Times=200;
//////////////////////////////////////////////////////
	void ComputeRank(int id);
	double Distance(const int* rank1,const int* rank2);
////////////////////////////////////////////////////////
	double DistributePoints();
	bool CalCenter();
	double kmeanscluster(int classes_);
	double DistanceAverage(vector<int> &kv1,vector<int> &kv2);
	double IntraDistance(vector<int> &kv1);

//////////////////////////////////////////////////////
	MetaCluster(){}
	MetaCluster(const MetaCluster &kmeans){}
	const MetaCluster &operator=(const MetaCluster &kmeans){return *this;}
};
#endif
