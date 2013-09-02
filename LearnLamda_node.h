/*
 * last fixed: 2012.01.05.
 * by Wang Yi.
 * */
#ifndef __LEARNLAMDA_H_

#define __LEARNLAMDA_H_
#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include "Ncbi_nodes_dmp.h"

using namespace std;
/*bool isSameFold(const string&p1,const string&p2){
	int idx1 = p1.length()-1, idx2 = p2.length()-1;
	while(idx1>0 && p1[idx1]!='/')--idx1;
	while(idx2>0 && p2[idx2]!='/')--idx2;
	if(idx1!=idx2)return false;
	for(;idx1>=0;--idx1)
		if(p1[idx1]!=p2[idx1])
			return false;
	return true;
}*/

bool pairCmp(const pair<double,bool>& t1,const pair<double,bool>& t2){
	if(t1.first==t2.first && t1.second!=t2.second)return !t1.second;
	return t1.first < t2.first;
}
double calLamda(const vector<double>& intra,const vector<double>& inter,int&OTP,int&OTN){
	vector<pair<double,bool> > V(inter.size()+intra.size());
	for(int i=0;i<inter.size();++i)
		V[i] = make_pair(inter[i],true);
	for(int i=0;i<intra.size();++i)
		V[i+inter.size()] = make_pair(intra[i],false);
	sort(V.begin(),V.end(),pairCmp);
	int max = intra.size(), maxid = -1;
	int TP = 0, TN = intra.size();
	OTP = TP; OTN = TN;
	for(int i=0;i<V.size();++i){
		if(V[i].second){
			++TP;
			if(TP+TN>max){
				max = TP+TN;
				maxid = i;
				OTP = TP;
				OTN = TN;
			}
		}
		else  --TN;
	}
	if(maxid==-1)return -1;
	else if(maxid==V.size()-1)return 1e99;
	else return (V[maxid].first + V[maxid+1].first)/2;
}

void getLamda(NCBI_nodes_dmp &Dmp, const string bac_info, const string mtx_path,map<int,vector<double> >&TaxonLamda,map<int,map<int,double> >&TaxonPrec){
//	map<int,map<int,double> >TaxonLamda;
//	map<int,map<int,double> >TaxonPrec;
	const int MAX=1000000;
	char* Buf;
	Buf = new char[MAX];

	vector<int>Tid;
	vector<int>Len;
	vector<string>Path;
	char tBuf[1000];
	ifstream ifs(bac_info.c_str());
	assert(!ifs.fail());
	while(!ifs.eof()){
		ifs.getline(Buf,MAX);
		int gid,tid=-1,len=0;
		sscanf(Buf,"%d\t%d\t%d\t%s",&gid,&tid,&len,tBuf);
		if(tid==-1)break;
		Tid.push_back(tid);
		Len.push_back(len);
		Path.push_back(tBuf);
	}
	ifs.close();
	delete[]Buf;
	int N = Tid.size();
	cerr << "# seq: " << N << endl;

	vector<vector<double> >Mtx(N);
	for(int i=0;i<N;++i)
		Mtx[i].resize(N,0);
	ifs.open(mtx_path.c_str());
	assert(!ifs.fail());
	for(int i=0;i<N;++i)
		for(int j=0;j<N;++j)
			ifs >> Mtx[i][j];
	ifs.close();
//	cerr << Mtx[0][0] << '\t' << Mtx[N-1][N-2] << '\t' << Mtx[N-1][N-1] << endl;

	vector<vector<int> >Taxon(N);
	for(int i=0;i<N;++i)
		Taxon[i] = Dmp.getTaxo(Tid[i]);

//	int L7[] = {26,22,19,15,10,6,4};
	int L7[] = {4,6,10,15,19,22,26};
	const int NL = sizeof(L7)/sizeof(int);
	cerr << "NL :" << NL << endl;
	vector<int>TP(NL,0);
	vector<int>TN(NL,0);
	vector<int>PN(NL,0);
	for(int i=0;i<N;++i){
		vector<double> Intra[NL+1];
		for(int j=0;j<N;++j){
			int k = 0;
			for(;k<NL && Taxon[i][L7[k]]!=Taxon[j][L7[k]];++k);
			Intra[k].push_back(Mtx[i][j]);
		}
		vector<double> lamdas(NL);
		for(int j=0;j<NL;++j){
			int OTP=0,OTN=0;
			lamdas[j]=calLamda(Intra[j], Intra[j+1],OTP,OTN);
			TP[j] += OTP;
			TN[j] += OTN;
			PN[j] += Intra[j].size() + Intra[j+1].size();
			double tall = Intra[j].size()+Intra[j+1].size();
			cout<<i<<'\t'<<j<<'\t'<<(OTP+OTN)/tall<<'\t' << OTP<<'\t'<<OTN<<'\t'<<tall<<endl;
		}
		TaxonLamda[Taxon[i][L7[0]]] = lamdas;
	}
	for(int i=0;i<NL;++i){
		cerr << i << '\t' << (TP[i]+TN[i])/(double)PN[i] << '\t' << TP[i] << '\t' << TN[i] << '\t' << PN[i] << endl;
	}

	/*
	for(int l=0;l<6;++l){
		int corC = 0, incC = 0, sigC = 0;
		int tTP=0,tTN=0, tAll = 0;
		int L1=L7[6-l],L2 = L7[5-l];
		map<int,map<int,vector<int> > > Comp;
		for(int i=0;i<N;++i)
			if(Taxon[i][L1]>0 && Taxon[i][L2]>0)
				Comp[Taxon[i][L1]][Taxon[i][L2]].push_back(i);
		for(map<int,map<int,vector<int> > >::const_iterator itr=Comp.begin();itr!=Comp.end();++itr){
			if(itr->second.size()<2)continue;
			int size = 0;
			for(map<int,vector<int> >::const_iterator itr2=itr->second.begin();itr2!=itr->second.end();++itr2)
				size += itr2->second.size();

			for(map<int,vector<int> >::const_iterator itr2=itr->second.begin();itr2!=itr->second.end();++itr2){
				vector<double>intra;
				vector<double>inter;
				for(vector<int>::const_iterator itr3=itr2->second.begin();itr3!=itr2->second.end();++itr3){
					double maxintra = 0;int maxid1=-1;
					for(vector<int>::const_iterator itr4=itr2->second.begin();itr4!=itr2->second.end();++itr4)
						if(itr4 != itr3 && Len[*itr4]>=0.5*Len[*itr3] 
								&& Taxon[*itr3][Taxon[*itr3].size()-1]!=10239
								&& Taxon[*itr4][Taxon[*itr4].size()-1]!=10239
		//						&& !isSameFold(Path[*itr3],Path[*itr4])
								){
							double tmax = max(Mtx[*itr3][*itr4],Mtx[*itr4][*itr3]);
							if(tmax > maxintra){
								maxintra = tmax;
								maxid1 = *itr4;
							}
						}
		//			intra.push_back(maxintra);
					intra.push_back(maxintra/Len[*itr3]);

					double maxinter = 0;int maxid2 = -1;
					for(map<int,vector<int> >::const_iterator itr4=itr->second.begin();itr4!=itr->second.end();++itr4){
						if(itr4->first == itr2->first)
							continue;
						for(vector<int>::const_iterator itr5=itr4->second.begin();itr5!=itr4->second.end();++itr5){
							double tmax = max(Mtx[*itr3][*itr5],Mtx[*itr5][*itr3]);
							if(tmax > maxinter){
								maxinter = tmax;
								maxid2 = *itr5;
							}
						}
					}
		//			inter.push_back(maxinter);
					inter.push_back(maxinter/Len[*itr3]);
					if(maxid1 == -1)
						++sigC;
					else if(maxintra < maxinter)
						++incC;
					else ++corC;
				}
				int TP,TN;
				TaxonLamda[5-l][itr2->first] = calLamda(intra,inter,TP,TN);
				TaxonPrec[5-l][itr2->first] = (TP+TN)/(double)(intra.size()+inter.size());
				tTP += TP; tTN += TN; tAll += intra.size()+inter.size();
				if(intra.size()>1)
					cout << 5-l << '\t' << (itr2->first) << '\t' << TP << '\t' << TN << '\t' << intra.size() << '\t' << inter.size() <<'\t' << TaxonLamda[5-l][itr2->first]<< '\t' << TaxonPrec[5-l][itr2->first] << endl;
			}
		}
		cerr << l << '\t'<<(tTP+tTN) <<'\t' << (tAll-tTP-tTN) <<'\t' << corC << '\t' << sigC << '\t' << incC << endl;
	}*/
}
#endif
