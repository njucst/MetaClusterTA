#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <map>
#include <vector>
#include <string>
#include "Methods.h"
using namespace std;

bool cmpReads(const KmerNode* p1, const KmerNode* p2)
{
	return (p1->VSize - p1->CtgNum) < (p2->VSize - p2->CtgNum);
}

void calTaxoForCtg(const BWTs& bwts, const string& str, map<int,double>& sp_score)
{
	set<string>kmers;
	{
		string rstr = str;
		for(int i=0;i<str.length();++i)
		{
			char ch = str[str.length()-1-i];
			switch(ch)
			{
				case 'A':rstr[i]='T';break;
				case 'C':rstr[i]='G';break;
				case 'G':rstr[i]='C';break;
				case 'T':rstr[i]='A';break;
				default:rstr[i]='N';
			}
		}
		int strend = (str.length()-(PARA_KMER-1));
		for(int i=0;i<strend;++i)
		{
			kmers.insert( str.substr(i,PARA_KMER));
			kmers.insert(rstr.substr(i,PARA_KMER));
		}
	}
	for(set<string>::const_iterator itr = kmers.begin();itr!=kmers.end();++itr)
	{
		set<int> cand_taxid;
	   	bwts.search(*itr, cand_taxid);
		if(cand_taxid.size()==0)continue;
		double unit = 1.0/cand_taxid.size();
		for(set<int>::const_iterator itr2 = cand_taxid.begin();itr2!=cand_taxid.end();++itr2)
			sp_score[*itr2] += unit;
	}
}

void compa_read(int ReadLen,const int &position1,const int &position2,int &match,int &mismatch,const ULLN &read1,const ULLN &read2)
{
	if(position1 > position2)
	{
		int k=position1-position2;
		ULLN differ= (read2>>(2*k));
//		differ ^= (read1 & AndMask[ReadLen-k]);
		differ ^= read1; differ.keeptopk(ReadLen-k);
		mismatch = differ.non0base();
		match = ReadLen-k-mismatch;
		//Note: can change == to ^, so that can return a value;
	}
	else if(position1 < position2)
	{
		int k=position2-position1;
		ULLN differ= (read1>>(2*k));
//		differ ^= (read2 & AndMask[ReadLen-k]);
		differ ^= read2; differ.keeptopk(ReadLen-k);
		mismatch = differ.non0base();
		match = ReadLen-k-mismatch;
		//Note: can change == to ^, so that can return a value;
	}
	else
	{
		ULLN differ= read1;
		differ ^= read2;
		mismatch = differ.non0base();
		match = ReadLen-mismatch;
	}
}

void MergeReads(const KmerNodeAloc& NodePool, const ReadsClass& Reads, USet& uset, int CtgNum)
{
	////////////////////const variables
//	const int StrLen = 50;
	const int LARGELENGTH = 8000;
	const int MINSCORE = 50;
	const int FragLen = 1000;
	const double Cdf_Cl = 1.65;

	////////////////////const variables
	int confidence[PARA_READ+1];
	{
		double r = 0.01;
		double p = 1+r*r-5*r/3;
		double var = p*(1-p);
		for(int i=0;i<=PARA_READ;++i)
			confidence[i] = i*p-sqrt(var*i)*Cdf_Cl;
	}
	if(INTEST)
	{
		cerr << "confidence levels: " << endl;
		for(int i=0;i<=PARA_READ;++i)
			cerr << i << ":\t" << confidence[i] << endl;
	}

	////////////////////sort KmerNodes
	int NodeNum = NodePool.getNodeNum();
	KmerNode** SortedNodes = new KmerNode*[NodeNum];
	for(int i=0; i<NodeNum; ++i)
		SortedNodes[i] = NodePool.getRef(i);
	sort(SortedNodes, SortedNodes+NodeNum, cmpReads);
	////////////////////merge nodes only with reads
	omp_lock_t type_lock;
	omp_init_lock(&type_lock);
	int ReadLen = Reads.ReadLen;
#pragma omp parallel for schedule(dynamic)
	for(int i=0; i<NodeNum; ++i)
	{
		NodeIDPosi* vect = SortedNodes[i]->myvector;
		int ctgnum = SortedNodes[i]->CtgNum;
		int sub_size = SortedNodes[i]->VSize - ctgnum;
		int* indexes = new int[sub_size];
		int* positions = new int[sub_size];
		bool* isRev = new bool[sub_size];
		for(int j=0;j<sub_size;++j)
		{
			indexes[j] = vect[j+ctgnum].id;
			positions[j] = vect[j+ctgnum].posi >> 2;
			isRev[j] = vect[j+ctgnum].isRev();
		}
	
		ULLN read1,read2;
		for(int i2=0;i2<sub_size;++i2)
		{
			int indexes1 = indexes[i2];
			int uidx1 = indexes1 + CtgNum;

			if(isRev[i2])	read1 = (*Reads.R2[indexes1]);
			else	read1 = (*Reads.R1[indexes1]);

	/*		int j=i+1;
			while(j<sub_size && positions[j]-positions[i2]<=ReadLen-StrLen)
				++j;*/
			for(int j=i2+1;j<sub_size;++j)
			{
				int indexes2 = indexes[j];
				int uidx2 = indexes2 + CtgNum;
				if(uset.find(uidx1)==uset.find(uidx2))
					continue;
			/*	if(positions[j]-positions[i2]>ReadLen-StrLenLowCover)
					break;*/
				/////////////////////////////////////////////////////////////////////////////////////
	//			if((uset.getReadNum(uidx1)<LARGESIZE && uset.getReadNum(uidx2)<LARGESIZE )||(uset.getReadNum(uidx1)<Frag) || uset.getReadNum(uidx2)<Frag)
				if((uset.getCtgLen(indexes1)<LARGELENGTH && uset.getCtgLen(indexes2)<LARGELENGTH )||(uset.getCtgLen(indexes1)<FragLen) || uset.getCtgLen(indexes2)<FragLen)
				{
					if(isRev[j])	read2 = (*Reads.R2[indexes2]);
					else read2 = (*Reads.R1[indexes2]);
					int match,mismatch;
					compa_read(ReadLen, positions[i2],positions[j],match,mismatch,read1,read2);
					assert(match >= 32);
				/////////////////bin small ones
			//		if(match-StrLenLowCover>=0 && (match-StrLenLowCover >= confidence[match+mismatch-StrLenLowCover]))
					if(match>=MINSCORE && (match-MINSCORE >= confidence[match+mismatch-MINSCORE]))
					{
						omp_set_lock(&type_lock);
						uset.Union(uidx1,uidx2);
						omp_unset_lock(&type_lock);
					}
				}
			}
		}
		delete[]indexes;delete[]positions;delete[]isRev;
	}
	omp_destroy_lock(&type_lock);
	delete[]SortedNodes;
	///////////////////////////////////////////////////////////////////////
}

void MergeCtgReads(const KmerNodeAloc& NodePool, const ReadsClass& Reads, const ContigsClass& Ctgs, USet& uset)
{
	////////////////////const variables
//	const int StrLenLowCover = 32;
//	const int StrLen = 50;
	const int LARGELENGTH = 8000;
	const int FragLen = 1000;
	const int MINSCORE = 50;
	const double Cdf_Cl = 1.65;

	////////////////////const variables
	int confidence[PARA_READ+1];
	{
		double r = 0.01;
		double p = 1+r*r-5*r/3;
		double var = p*(1-p);
		for(int i=0;i<=PARA_READ;++i)
			confidence[i] = i*p-sqrt(var*i)*Cdf_Cl;
	}
	/////////////////////////////////////////
	int CtgNum = Ctgs.CtgNum;

	////////////////////sort KmerNodes
	int NodeNum = NodePool.getNodeNum();
	KmerNode** SortedNodes = new KmerNode*[NodeNum];
	for(int i=0; i<NodeNum; ++i)
		SortedNodes[i] = NodePool.getRef(i);
	sort(SortedNodes, SortedNodes+NodeNum, cmpReads);
	////////////////////merge nodes only with reads
	omp_lock_t type_lock;
	omp_init_lock(&type_lock);
	int ReadLen = Reads.ReadLen;
#pragma omp parallel for schedule(dynamic)
	for(int i=0; i<NodeNum; ++i)
	{
		NodeIDPosi* vect = SortedNodes[i]->myvector;
		int ctgnum = SortedNodes[i]->CtgNum;
		int sub_size = SortedNodes[i]->VSize;
		int* indexes = new int[sub_size];
		int* positions = new int[sub_size];
		bool* isRev = new bool[sub_size];
		for(int j=0;j<sub_size;++j)
		{
			indexes[j] = vect[j].id;
			positions[j] = vect[j].posi >> 2;
			isRev[j] = vect[j].isRev();
		}
		for(int i2=0;i2<ctgnum;++i2)
		{
			int indexes1 = indexes[i2];
			assert(indexes1<CtgNum);
			BaseStr* ctg1 = Ctgs.contigs[indexes1];

	//		ULLN read2;
			for(int j=ctgnum;j<sub_size;++j)
			{
				int indexes2 = indexes[j];
				assert(indexes2<Reads.ReadNum);
				int uidx2 = indexes2 + CtgNum;
				if(uset.find(indexes1)==uset.find(uidx2) || Reads.Score[indexes2]<MINSCORE)
					continue;
	/*			if(positions[j]-positions[i2]>ReadLen-StrLenLowCover)
					break;*/
				/////////////////////////////////////////////////////////////////////////////////////
				if((uset.getCtgLen(indexes1)<LARGELENGTH && uset.getCtgLen(uidx2)<LARGELENGTH )||(uset.getCtgLen(indexes1)<FragLen) || uset.getCtgLen(uidx2)<FragLen)
				{
	/*				if(isRev[j])	read2 = (*Reads.R2[indexes2]);
					else read2 = (*Reads.R1[indexes2]);*/
					int match,mismatch;
					if(isRev[i2])
					{
						if(isRev[j])
							ctg1->checkRead(positions[i2],(*Reads.R1[indexes2]),ReadLen-1-positions[j],ReadLen,match,mismatch);
						else
							ctg1->checkRead(positions[i2],(*Reads.R2[indexes2]),ReadLen-1-positions[j],ReadLen,match,mismatch);
					}
					else
					{
						if(isRev[j])
							ctg1->checkRead(positions[i2],(*Reads.R2[indexes2]),PARA_KMER+positions[j]-1,ReadLen,match,mismatch);
						else
							ctg1->checkRead(positions[i2],(*Reads.R1[indexes2]),positions[j]+PARA_KMER-1,ReadLen,match,mismatch);
					}
/*					if(INTEST && match < 32)
					{
						omp_set_lock(&type_lock);
						cerr << dec << sub_size << '\t' << ctgnum << '\t' << ((vect[i2].posi)&0x3) << '\t' << ((vect[j].posi)&0x3) << endl;
						cerr << SortedNodes[i]->kmer << endl;
						cerr << dec << "match1 < 32:\t" << match << '\t' <<  indexes1 <<'\t' << positions[i2]<<'\t' << indexes2<<'\t'<<positions[j]+PARA_KMER-1 << endl;
						cerr << ctg1->subStr(positions[i2]-positions[j]+1-PARA_KMER, ReadLen) << endl;
						cerr << ctg1->subStr(positions[i2]+1-PARA_KMER, ReadLen) << endl;
						omp_unset_lock(&type_lock);
					}*/
					assert(match >= 32);
					if(match>=MINSCORE && (match-MINSCORE >= confidence[match+mismatch-MINSCORE]))
					{
						omp_set_lock(&type_lock);
						uset.Union(indexes1,uidx2);
						omp_unset_lock(&type_lock);
					}
				}
			}
		}
		delete[]indexes;delete[]positions;delete[]isRev;
	}
	omp_destroy_lock(&type_lock);
	delete[]SortedNodes;
	///////////////////////////////////////////////////////////////////////
}

void MergePE(const ReadsClass& Reads, USet& uset, AccTester& acc_tester)
{
	/////////////////////////////////////////get pair-end graph (PEGraph);
	map<unsigned, map<unsigned, unsigned> >PEGraph;
	SparseMatrix SPMTX;
	int SPMTXSIZE = SPMTX.Size;
	omp_lock_t* hash_lock = new omp_lock_t[SPMTXSIZE];
	for(int i=0;i<SPMTXSIZE;++i)
		omp_init_lock(&hash_lock[i]);

#pragma omp parallel for
	for(int i=0;i<Reads.TotalNum;i+=2)
	{
		int idx1 = Reads.MatchId[i];
		int idx2 = Reads.MatchId[i+1];
		if(idx1!=idx2 && (idx1>=0 && idx2 >= 0))
		{
			unsigned id = SPMTX.hash(idx1,idx2);
			omp_set_lock(&hash_lock[id]);
			SPMTX.insert(idx1,idx2);
			SPMTX.insert(idx2,idx1);
			omp_unset_lock(&hash_lock[id]);
		}
	}
	for(int i=0;i<SPMTXSIZE;++i)
		omp_destroy_lock(&hash_lock[i]);
	delete[] hash_lock;
	SPMTX.toNeighbor(PEGraph);
	////////////////////////////////////////////

	//////////////////////////////////////merge pair-end reads
	map<unsigned,unsigned> weightsum;
	map<unsigned,unsigned> friends;
	for(map<unsigned, map<unsigned, unsigned> >::const_iterator itr1 = PEGraph.begin(); itr1 != PEGraph.end(); ++itr1)
	{
		unsigned tsum = 0;
		for(map<unsigned, unsigned>::const_iterator itr2 = itr1->second.begin(); itr2 != itr1->second.end(); ++itr2)
			tsum += itr2->second;
		weightsum[itr1->first] = tsum;
		friends[itr1->first] = itr1->second.size();
	}

	for(map<unsigned, map<unsigned, unsigned> >::const_iterator itr1 = PEGraph.begin(); itr1 != PEGraph.end(); ++itr1)
		for(map<unsigned, unsigned>::const_iterator itr2 = itr1->second.begin(); itr2 != itr1->second.end(); ++itr2)
			if((itr2->second/(double)(max(weightsum[itr1->first],weightsum[itr2->first])))>=0.25)
				uset.Union(itr1->first, itr2->first);

	if(INTEST)acc_tester.calAcc(uset);
	for(map<unsigned, map<unsigned, unsigned> >::const_iterator itr1 = PEGraph.begin(); itr1 != PEGraph.end(); ++itr1)
	{
		for(map<unsigned, unsigned>::const_iterator itr2 = itr1->second.begin(); itr2 != itr1->second.end(); ++itr2)
		{
			if(uset.getCtgLen(itr1->first)>10000 && uset.getCtgLen(itr2->first)>10000)
				continue;
			if((itr2->second/(double)(max(weightsum[itr1->first],weightsum[itr2->first])))>=0.20)
				uset.Union(itr1->first, itr2->first);
		}
	}
	if(INTEST)acc_tester.calAcc(uset);
	for(map<unsigned, map<unsigned, unsigned> >::const_iterator itr1 = PEGraph.begin(); itr1 != PEGraph.end(); ++itr1)
	{
		for(map<unsigned, unsigned>::const_iterator itr2 = itr1->second.begin(); itr2 != itr1->second.end(); ++itr2)
		{
			if(uset.getCtgLen(itr1->first)>10000 && uset.getCtgLen(itr2->first)>10000)
				continue;
			if((itr2->second/(double)(max(weightsum[itr1->first],weightsum[itr2->first])))>=0.15)
				uset.Union(itr1->first, itr2->first);
		}
	}
	//////////////////////////////////////////////////////////////////////////////
	////////////////
	/*
	cerr << "start to get purity of each ctg." << endl;
	vector<double> purity(acc_tester.CtgNum,0);
	vector<int>majorid(acc_tester.CtgNum,0);
	for(int i=0;i<acc_tester.CtgNum;++i)
	{
		int tmajor = 0, tmax = 0, sum = 0;
		for(int j=0;j<acc_tester.GenomeNum;++j)
		{
			int t = acc_tester.CtgComp[i][j];
			sum += t;
			if(t > tmax)
			{
				tmajor = j;
				tmax = t;
			}
		}
		purity[i] = tmax/(double)sum;
		majorid[i] = tmajor;
	}
	////////////////
	cerr << "pe-graph\t" << endl;
	for(map<unsigned, map<unsigned, unsigned> >::const_iterator itr1 = PEGraph.begin(); itr1 != PEGraph.end(); ++itr1)
	{
		for(map<unsigned, unsigned>::const_iterator itr2 = itr1->second.begin(); itr2 != itr1->second.end(); ++itr2)
		{
			if((itr2->second/(double)(max(weightsum[itr1->first],weightsum[itr2->first])))<0.1)
				continue;

			if(majorid[itr1->first]==majorid[itr2->first])
				cerr << "true\t";
			else
				cerr << "false\t";
			cerr << itr2->second/(double)weightsum[itr1->first] << '\t';
			cerr << itr2->second/(double)weightsum[itr2->first] << '\t';
			cerr << friends[itr1->first] <<'\t'<<friends[itr2->first]<<'\t';
			cerr << purity[itr1->first] << '\t' << purity[itr2->first] << '\t';
			cerr << itr2->second << '\t';
			cerr << weightsum[itr1->first] << '\t';
			cerr << weightsum[itr2->first] << '\t';
			cerr << uset.getCtgLen(itr1->first) << '\t';
			cerr << uset.getCtgLen(itr2->first) << '\t';
			cerr << itr2->second/(double)(max(weightsum[itr1->first],weightsum[itr2->first]))<<'\t';
			cerr << endl;
		}
	}
	cerr << "end of pe-graph\t" << endl;
	*/
}

void MergeAsStep1(const ContigsClass& Ctgs,const ReadsClass& Reads, KmerNodeAloc& NodePool, USet& uset, AccTester& acc_tester, int pair_end_merge_threshold)
{
	///////////////////////testing
	if(INTEST)
	{
		map<int,int>halfmapscore;
		int unmatch = 0;
		int bothunmapped = 0;
		for(int i=0;i<Reads.TotalNum;i+=2)
		{
			int idx1 = Reads.MatchId[i];
			int idx2 = Reads.MatchId[i+1];
			if(idx1!=idx2 && (idx1>=0 && idx2 >= 0))
				++unmatch;
			if(idx1<0 && idx2<0)
				++bothunmapped;
			if(idx1 >=0 && idx2<0)
				++halfmapscore[Reads.Score[i+1]];
			if(idx2 >=0 && idx1<0)
				++halfmapscore[Reads.Score[i]];
		}
		cerr << "Mixed mapped pair-end reads: \t" << unmatch/(double)Reads.TotalNum << '\t' << unmatch << '\t' << Reads.TotalNum << endl;
		cerr << "buth unmaped pair-end reads: \t" << bothunmapped/(double)Reads.TotalNum << '\t' << bothunmapped << '\t' << Reads.TotalNum << endl;
		cerr << "half mapped scores:\t"<< endl;
		for(map<int,int>::const_iterator itr=halfmapscore.begin();itr!=halfmapscore.end();++itr)
			cerr << itr->first << '\t' << itr->second << endl;;
	}
	/////////////////////merge pair-end reads
/*	for(int i=0;i<Reads.TotalNum;i+=2)
	{
		if(Reads.MatchId[i]<0 && Reads.MatchId[i+1]<0)
			uset.Union(Ctgs.CtgNum+Reads.toNewId(i+1),Ctgs.CtgNum+Reads.toNewId(i));
	}*/
	if(INTEST)acc_tester.calAcc(uset);
	///////////////////////////////////////////////////////////////////////
	MergePE(Reads, uset, acc_tester);
	if(INTEST)acc_tester.calAcc(uset);
/*	MergeCtgReads(NodePool, Reads, Ctgs, uset);
	if(INTEST)acc_tester.calAcc(uset);
	MergeReads(NodePool, Reads, uset, Ctgs.CtgNum);
	if(INTEST)acc_tester.calAcc(uset);
	*/
	///////////////////////////////////////////////////////////////////////
}

void getScoreV(int** GenoDBTaxo, const vector<int>&CtgId, const int taxoLevel, int& AnsLevel, double& AnsScore)
{
	map<int,int> TaxoStatis;
	for(vector<int>::const_iterator itr=CtgId.begin();itr!=CtgId.end();++itr)
		++TaxoStatis[GenoDBTaxo[*itr][taxoLevel]];
	double ans = 0;
	for(map<int,int>::const_iterator itr=TaxoStatis.begin();itr!=TaxoStatis.end();++itr)
	{
		double p = itr->second/(double)CtgId.size();
		ans += p*log(p);
	}
	int maxId = 0,maxC = 0;
	for(map<int,int>::const_iterator itr=TaxoStatis.begin();itr!=TaxoStatis.end();++itr)
		if(itr->second > maxC)
		{
			maxId = itr->first;
			maxC = itr->second;
		}
	AnsLevel = maxId;
	AnsScore = -ans;
}

void getTaxo4Clust(int** GenoDBTaxo, const vector<int>&CtgId, double Thresh, int& AnsLevel, double& AnsScore)
{
	for(int i=6;i>=0;--i)
	{
		getScoreV(GenoDBTaxo, CtgId, i, AnsLevel, AnsScore);
		if(AnsScore >= Thresh)
			return;
	}
}

double getScoreEntropy(const map<int,long long> &TaxoComp)
{
	long long sum = 0;
	for(map<int,long long>::const_iterator itr = TaxoComp.begin(); itr!=TaxoComp.end(); ++itr)
		sum += itr->second;
	if(sum==0) return 1e9;
	double ans = 0;
	for(map<int,long long>::const_iterator itr = TaxoComp.begin(); itr!=TaxoComp.end(); ++itr)
	{
		double tp = (itr->second)/(double)sum;
		ans += tp*log(tp);
	}
	return ans/(log(2));
}
double getScoreMax(const map<int,long long> &TaxoComp)
{
	long long sum = 0,max=0;
	for(map<int,long long>::const_iterator itr = TaxoComp.begin(); itr!=TaxoComp.end(); ++itr)
	{
		sum += itr->second;
		if(itr->second > max)
			max = itr->second;
	}
	if(sum==0)return 1e9;
	return max/(double)sum;
}
double getScoreMax2(const map<int,long long> &TaxoComp)
{
	long long max1=0,max2=0,sum=0;
	for(map<int,long long>::const_iterator itr = TaxoComp.begin(); itr!=TaxoComp.end(); ++itr)
	{
		sum += itr->second;
		if(itr->second > max1)
		{
			max2 = max1;
			max1 = itr->second;
		}
		else if(itr->second > max2)
			max2 = itr->second;
	}
	if(sum==0)return 1e9;
	return (max1-max2)/(double)sum;
}

void outputResult(int* MatchId, long long ReadN, USet& uset,int* metabest, MCPara& mcpara,string infile,vector<ClustTaxoInfoClass>& taxoofclust)
{
	map<int,vector<unsigned> > R;
	vector<unsigned> Gid(ReadN,0);
	for(long long i=0;i<ReadN;++i)
	{
		int ctgid = MatchId[i];
		if(ctgid < -1)
			ctgid = -ctgid-2;
		if(ctgid==-1)
		{
			R[-1].push_back(i);
			continue;
		}
		ctgid = uset.find(ctgid);
		long long metaid;
		if((metaid=mcpara.getNewId(ctgid))<0)
		{
			R[-1].push_back(i);
			continue;
		}
		R[metabest[metaid]].push_back(i);
	}
	ifstream ifs(infile.c_str());
	vector<string> readtitle(ReadN,"");
	if(!ifs.fail())
	{
		const int BMax = 10000;
		char Buf[BMax];
		for(long long i=0;i<ReadN;++i)
		{
			if(ifs.eof())break;
			ifs.getline(Buf,BMax);
			readtitle[i] = Buf;
			ifs.getline(Buf,BMax);
		}
	}
	
	ofstream ofs((infile+".clust").c_str());
	for(map<int,vector<unsigned> >::const_iterator itr1 = R.begin();itr1!=R.end();++itr1)
	{
		ofs << "cluster " << itr1->first << '\t' << itr1->second.size() << endl;
		ofs << "Taxo:\t";
		if(itr1->first >=0 && itr1->first < taxoofclust.size())
			for(int i=0;i<28;++i)
				ofs << taxoofclust[itr1->first].taxo[i].taxid << ':' << taxoofclust[itr1->first].taxo[i].score << '\t';
		ofs << endl;

		for(vector<unsigned>::const_iterator itr2 = itr1->second.begin();itr2!=itr1->second.end();++itr2)
		{
			if(readtitle[*itr2].length()==0)
				ofs << "< read " << (*itr2) << endl;
			else 
				ofs << readtitle[*itr2] << endl;
		}
	}
}
