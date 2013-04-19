/* tester for simulated data;
 * by Wang Yi
 * wangyi.tom@gmail.com
 * April 19,2013
 */
#ifndef MCH_HYBRID_ACCTESTER_H_

#define MCH_HYBRID_ACCTESTER_H_

#include <cstdio>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include "GLOBAL.h"
#include "USet.h"

using namespace std;
class AccTester
{
public:
	vector<vector<int> >CtgComp;
	int GenomeNum;
	int CtgNum;
	AccTester()
	{
		minId = 1;
		CtgNum = 0;
		GenomeNum = 0;
		TCounter = 0;;
	}
	void init_total(int TotalNum_)
	{
		All_Reads_GenomeId.resize(TotalNum_);
		All_Reads_Posi.resize(TotalNum_);
		Reads_CtgId.resize(TotalNum_);
	}
	void getGId(int idx, char*Buf)
	{
		int gid=-1,tmp,posi=-1;
		sscanf(Buf,">read%d_%d/%d\t%d",&gid,&tmp,&tmp,&posi);
		All_Reads_GenomeId[idx] = gid;
		All_Reads_Posi[idx] = posi;
		if(gid < minId)	minId = gid;
		if(gid > GenomeNum)GenomeNum = gid;
	}
	void init(const int TotalNum_, int*MatchId_, const int ReadNum_,int*NewIdToOldId_)
	{
		for(int i=0;i<TotalNum_;++i)
		{
			Reads_CtgId[i] = MatchId_[i];
			if(MatchId_[i] > CtgNum)
				CtgNum = MatchId_[i];
		}
		++CtgNum;
		++GenomeNum;
		if(minId == 1)
		{
			for(int i=0;i<TotalNum_;++i)
				--All_Reads_GenomeId[i];
			--GenomeNum;
		}
		CtgComp.resize(CtgNum);
		for(int i=0;i<CtgNum;++i)
			CtgComp[i].resize(GenomeNum);
		for(int i=0;i<TotalNum_;++i)
			if(Reads_CtgId[i]>=0)
				++CtgComp[Reads_CtgId[i]][All_Reads_GenomeId[i]];

		{
		}

		UnmapedNum = ReadNum_;
		Unmaped_Reads_GenomeId.resize(UnmapedNum);
		for(int i=0;i<UnmapedNum;++i)
			Unmaped_Reads_GenomeId[i] = All_Reads_GenomeId[NewIdToOldId_[i]];
	}

	//pre-condition: the first CtgNum elements in uset are contigs and the rest elements are reads.
	void calAcc(USet& uset, const int Thresh = 100,bool isOutputMtx=false)
	{
		assert(uset.size()==CtgNum+UnmapedNum);
		cout << "calAcc: " << isOutputMtx << endl;

		map<int,int> toNewClustId;
		int idx = 0;
		for(int i=0;i<uset.size();++i)
		{
			int cur = uset.find(i);
			if(toNewClustId.find(cur)==toNewClustId.end())
				toNewClustId[cur] = idx++;
		}
		/////////////////////////////
		////initialize comp
		vector<vector<int> >Comp;
		Comp.resize(idx);
		for(int i=0;i<idx;++i)
			Comp[i].resize(GenomeNum);
		
		for(int i=0;i<CtgNum;++i)
			for(int j=0;j<GenomeNum;++j)
				Comp[toNewClustId[uset.find(i)]][j] += CtgComp[i][j];
		for(int i=CtgNum;i<uset.size();++i)
			++Comp[toNewClustId[uset.find(i)]][Unmaped_Reads_GenomeId[i-CtgNum]];

		int Asum = 0, Amax=0;;
		for(int i=0;i<idx;++i)
		{
			int sum = 0;
			for(int j=0;j<GenomeNum;++j)
				sum += Comp[i][j];
			Asum += sum;

			int max = 0;
			for(int j=0;j<GenomeNum;++j)
				if(Comp[i][j] > max)
					max = Comp[i][j];
			Amax += max;

			if(sum < Thresh)continue;
			if(isOutputMtx)
			{
				for(int j=0;j<GenomeNum;++j)
					cout << Comp[i][j] << '\t';
				cout << max/(double)sum << '\t' << sum << endl;
			}
		}
		cout << TCounter << " precision:\t" << Amax/(double)Asum << '\t' << Amax << '\t' << Asum << endl;
		cerr << TCounter << " precision:\t" << Amax/(double)Asum << '\t' << Amax << '\t' << Asum << endl;

		Amax=0;
		for(int i=0;i<GenomeNum;++i)
		{
			int max=0;
			for(int j=0;j<idx;++j)
				if(max < Comp[j][i])
					max = Comp[j][i];
			Amax += max;
		}
		cout << TCounter << " sensitivity:\t" << Amax/(double)Asum << '\t' << Amax << '\t' << Asum << endl;
		cerr << TCounter << " sensitivity:\t" << Amax/(double)Asum << '\t' << Amax << '\t' << Asum << endl;
		++TCounter;
		/////////////////////////////////////////////////////
		map<unsigned,unsigned> Tsize;
		for(int i=0;i<uset.size();++i)
		{
			if(uset.find(i)!=i)continue;
			if(uset.getCtgLen(i)/1000 < 20)
				++Tsize[uset.getCtgLen(i)/1000];
			else
				++Tsize[20];
		}
		cout << " ctg lengths after merging: " << endl;
		for(map<unsigned,unsigned>::const_iterator itr=Tsize.begin();itr!=Tsize.end();++itr)
			cout << itr->first << '\t' << itr->second << endl;
	}
	void getPreSen4Other(int**Component_,int Size,int*best, const int Thresh=100)
	{
	/*	int idx = Size;
		int Asum = 0, Amax=0;;
		for(int i=0;i<idx;++i)
		{
			int sum = 0;
			for(int j=0;j<GenomeNum;++j)
				sum += Comp[i][j];
			Asum += sum;

			int max = 0;
			for(int j=0;j<GenomeNum;++j)
				if(Comp[i][j] > max)
					max = Comp[i][j];
			Amax += max;

			if(sum < Thresh)continue;
			for(int j=0;j<GenomeNum;++j)
				cerr << Comp[i][j] << '\t';
			cerr << max/(double)sum << '\t' << sum << endl;
		}
		cerr << "precision:\t" << Amax/(double)Asum << '\t' << Amax << '\t' << Asum << endl;

		Amax=0;
		for(int i=0;i<GenomeNum;++i)
		{
			int max=0;
			for(int j=0;j<idx;++j)
				if(max < Comp[j][i])
					max = Comp[j][i];
			Amax += max;
		}
		cerr << "sensitivity:\t";
		cerr << Amax/(double)Asum << '\t' << Amax << '\t' << Asum << endl;
		*/
		map<int,int> toNewClustId;
		int idx = 0;
		{
			set<int>allids;
			for(int i=0;i<Size;++i)
				if(best[i]>=0)
					allids.insert(best[i]);
			for(set<int>::const_iterator itr=allids.begin();itr!=allids.end();++itr)
				toNewClustId[*itr] = idx++;
		}
/*
		for(int i=0;i<Size;++i)
		{
			int cur = best[i];
			if(cur>=0 && toNewClustId.find(cur)==toNewClustId.end())
				toNewClustId[cur] = idx++;
		}*/
		/////////////////////////////
		////initialize comp
		vector<vector<int> >Comp;
		Comp.resize(idx);
		for(int i=0;i<idx;++i)
			Comp[i].resize(GenomeNum);
		
		int outlierctg = 0;
		int outlierread= 0;
		for(int i=0;i<Size;++i)
		{
			if(best[i]>=0)
			{
				for(int j=0;j<GenomeNum;++j)
					Comp[toNewClustId[best[i]]][j] += Component_[i][j];
			}
			else
			{
			   	++outlierctg;
				for(int j=0;j<GenomeNum;++j)
					outlierread += Component_[i][j];
			}
		}
		cerr << "outlierctg: \t" << outlierctg << endl;
		cerr << "outlierread: \t"<< outlierread<< endl;

		cout << ">clust matrix." << endl;
		int Asum = 0, Amax=0;;
		for(int i=0;i<idx;++i)
		{
			int sum = 0;
			for(int j=0;j<GenomeNum;++j)
				sum += Comp[i][j];
			Asum += sum;

			int max = 0,maxid=0;
			for(int j=0;j<GenomeNum;++j)
				if(Comp[i][j] > max)
				{
					max = Comp[i][j];
					maxid = j;
				}
			Amax += max;

//			if(sum < Thresh)continue;
			for(int j=0;j<GenomeNum;++j)
				cout << Comp[i][j] << '\t';
			cout << max/(double)sum << '\t' << maxid << '\t' << sum << endl;
		}
		cerr << "precision:\t" << Amax/(double)Asum << '\t' << Amax << '\t' << Asum << endl;
		cout << "precision:\t" << Amax/(double)Asum << '\t' << Amax << '\t' << Asum << endl;

		Amax=0;
		for(int i=0;i<GenomeNum;++i)
		{
			int max=0;
			for(int j=0;j<idx;++j)
				if(max < Comp[j][i])
					max = Comp[j][i];
			Amax += max;
		}
		cerr << "sensitivity:\t" << Amax/(double)Asum << '\t' << Amax << '\t' << Asum << endl;
		cout << "sensitivity:\t" << Amax/(double)Asum << '\t' << Amax << '\t' << Asum << endl;
	}
private:
	int minId;
	int UnmapedNum;
	int TCounter;
	vector<int> All_Reads_GenomeId;
	vector<int> All_Reads_Posi;

	vector<int> Unmaped_Reads_GenomeId;
	vector<int> Reads_CtgId;
};

#endif
