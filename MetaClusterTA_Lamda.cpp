/*
 * (1) use scores of contigs to predict taxon of contigs
 * and (2) use taxon of contigs to predict taxon of clusters
 */
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <omp.h>

#include "Methods.h"
#include "LearnLamda.h"
using namespace std;
/////////////////////////////////////////////////////////////////////
//parameters
//static const unsigned UNASSIGNED  = ~(0U);
//const double ScoreThresh = -1.0;
int ReadLen = 80;
int MetaKmerLen = 5;
int ClusterSize = 0;
int MaxSpecies = 200;
int MinSpecies = 2;
int StrLen = 50;
int StrLenLowCover = 25;
double MC3_Thresh = 0.94;
int OrphanLen = 400;

int Num_Thread = 0;
int CtgLenThresh = 500;
int AlignThresh = 76;
/////////////////////////////////////////////////////////////////////
//global variables
KmerNode** KmerMap;
ContigsClass Ctgs;
ReadsClass Reads;
KmerNodeAloc NodePool;
USet uset;
NCBI_nodes_dmp NodesDmp;
vector<ClustTaxoInfoClass> TaxoOfClust;

AccTester acc_tester;
////////////////////////////////////////////////////////////////////////////
/////////for test
int Forbid100[] = {1573,1398,1997,75,1490,1666,585,659,37,1995,1121,1862,863,873,683,693,1615,582,1352,2056,1332,2012,1364,22,1504,773,232,193,1404,781,1466,753,962,2037,1652,2002,851,266,201,1672,280,1158,311,1494,1284,823,642,884,347,1731,838,1078,1520,1958,946,1468,665,1688,1689,1382,1945,12,1475,631,620,964,1101,1899,342,128,1091,1326,1748,24,932,1090,854,79,879,1474,1125,1142,531,858,1679,1669,2023,1274,26,1698,206,968,1623,1627,412,56,411,1874,1339,926};
struct SpInfoOfCtgs
{
	int length;
	int NReads;

	int matchedLength;
	int spTaxId;
	int refLevel;
	int spFId;
	SpInfoOfCtgs():length(0),matchedLength(0),spTaxId(0),refLevel(0),spFId(0),NReads(0){}
};
vector<SpInfoOfCtgs> spInfoOfCtgs;
////////////////////////////////////////////////////////////////////////////
void usage()
{
	cerr << "usage:\t MetaClusterTA contig.fa read.fa bwt.idx nodes.dmp mtx.csv [options]" << endl;
	cerr << endl;
	exit(-1);
}

void printtime(string str = "")
{
	static time_t rawtime;
	static tm* timeinfo;
	cerr << str << endl;
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	cerr<<asctime(timeinfo);
}
//////////////////////////////
void init()
{
	{
		unsigned long long u64Exp = ((((1ULL << 50)%HASHSIZE)<<14))%HASHSIZE;
		unsigned long long tExp = 1;
		for(int i=0;i<10;++i)
		{
			U64HASH[i] = tExp;
			tExp *= u64Exp;
			tExp %= HASHSIZE;
		}
	}
	KmerMap = new KmerNode*[HASHSIZE];
	for(unsigned i=0;i<HASHSIZE;++i)
		KmerMap[i] = NULL;
}

void input(int argc, char* argv[])
{
	if(argc<4)
		usage();
	//Important! readlength should belong to (50,128]. Default value = 75.
	AddParameter("ReadLen", &ReadLen, INTEGER);
	//We use q-mer distribution in Phase 3. Here you can set qmer = 4 or 5. Default value = 4.
	AddParameter("qmer", &MetaKmerLen, INTEGER);
	//You can fix the #. of species. Default value = 0. 
	//In this case, the program will predict the #. of species.
	AddParameter("Species", &ClusterSize, INTEGER);
	//We use binary search to predict #. of species, this is max value for the first iteration.
	//Set this number properly will improve the performance and speed up the program. Default = 600.
	AddParameter("MaxSpecies", &MaxSpecies, INTEGER);
	//We use binary search to predict #. of species, this is min value for the first iteration. 
	//Set this number properly will improve the performance and speed up the program. Default = 2.
	AddParameter("MinSpecies", &MinSpecies, INTEGER);
	//wmer in phase 1, default = 50;
	AddParameter("wmer", &StrLen, INTEGER);
	//wmer for low coverage read , default = 25; if a smaller value is needed, mapinto() should be modifed to mappint more 16-mers.
	AddParameter("lowwmer", &StrLenLowCover, INTEGER);
	//Threshold for further-merge in MetaCluster 3.0. Default valule = 0.9. If you set this value larger, program will conceptually predict more species. 
	AddParameter("MetaC3Threshold", &MC3_Thresh, FLOAT);
	//groups of size smaller than Orphan16merThresh will be considered as orphans
	AddParameter("Orphan", &OrphanLen, INTEGER);
	//
	AddParameter("TN_WEIGHT", &TN_WEIGHT, FLOAT);
	//////////////////////////////////////////////////////////////////////////
	//set minimum length of contigs.
	AddParameter("CtgLenThresh", &CtgLenThresh, INTEGER);
	AlignThresh = ReadLen*0.95;
	AddParameter("AlignThresh", &AlignThresh, INTEGER);

	ProcessParameters(argc,argv);
	/////////////////////////////////////////////////////////////////
	cerr << dec << "ReadLen:\t " << ReadLen << endl;
	cerr << "CtgLenThresh:\t " << CtgLenThresh << endl;
	cerr << "AlignThresh:\t " << AlignThresh << endl;
	cerr << "MC3_Thresh:\t" << MC3_Thresh << endl;
	cerr << "TN_WEIGHT:\t" << TN_WEIGHT << endl;
	/////////////////////////////////////////////////////////////////
	if(INTEST){
		printtime("before loading contigs. ");
		system("ps ux");
	}
	Ctgs.init(argv[1], CtgLenThresh,KmerMap, NodePool);
	if(INTEST){
		cerr << "# of Contigs: " << Ctgs.CtgNum << "\tNodeNum:" << NodePool.getNodeNum()<< endl;
		printtime("Contigs loaded. ");
	}
	NodePool.fixCtgNum();
	Reads.init(argv[2], ReadLen, AlignThresh, KmerMap, NodePool, Ctgs, acc_tester);
	if(INTEST){
		cerr << "NodeNum:" << NodePool.getNodeNum()<< endl;
		printtime("initializing uset: ");
		system("ps ux");
	}
	{
		int CtgNum = Ctgs.CtgNum;
		unsigned* ctglen = new unsigned[CtgNum];
		for(int i=0;i<CtgNum;++i)
			ctglen[i] = Ctgs.contigs[i]->length;
		uset.init(Ctgs.CtgNum + Reads.ReadNum, CtgNum, ctglen);
		delete[]ctglen;
	}
	delete[] KmerMap;
	KmerMap = NULL;
	///////////////////////////////////////////////////////////
	if(INTEST){
		cerr << dec << "Total Reads:\t"<< Reads.TotalNum << endl;
		cerr << "Unmaped Reads:\t" << Reads.ReadNum << endl;
		cerr << "AlignThresh:\t" << Reads.AlignThresh << endl;
		printtime("initializing acc_tester:\t");
		system("ps ux");
		acc_tester.init(Reads.TotalNum, Reads.MatchId, Reads.ReadNum,Reads.NewIdToOldId);
	}
//	NodePool.shrinkSize(NULL,2);
}

bool tCtgInfoSort(const vector<int>& t1,const vector<int>&t2)
{
	if(t1[5]==t2[5])
	{
		if(t1[6]==t2[6])
			return t1[1]>t2[1];
		else return t1[6]<t2[6];
	}
	return t1[5]<t2[5];
}

//const int L8[] = {26,22,19,15,10,6,4};
const int L8[] = {4,6,10,15,19,22,26};
const int NL = 7;
void anaCluster(MCPara& mcpara, MetaCluster& metacluster,vector<unsigned>& ReadNumInCtg,char* argv[])
{
	if(INTEST){
		///////////////////////////////////////////////
		//initial Forbid
		spInfoOfCtgs.resize(Ctgs.CtgNum);
		for(int i=0;i<Reads.TotalNum;++i)
			if(Reads.MatchId[i]>=0)
				++spInfoOfCtgs[Reads.MatchId[i]].NReads;
		ifstream ifs("/home/ywang/DB/list/No_Taxid.list");
		assert(!ifs.fail());
		const int MAXFile = 5000;
		int Fid2Taxid[MAXFile];
		for(int i=0;i<MAXFile;++i)
			Fid2Taxid[i] = 0;
		while(!ifs.eof()){
			int tn=-1,ttax;
			ifs >> tn >> ttax;
			if(tn>=0)Fid2Taxid[tn]=ttax;
		}
		ifs.close();
		vector<int> Forbid(Forbid100,Forbid100 + 100);
		if(acc_tester.GenomeNum == 2){
			int Forbid2[] = {716,753};
			Forbid.assign(Forbid2,Forbid2+2);
		}
		static int ForbidNum = Forbid.size();
		cerr << ForbidNum << '\t' << acc_tester.GenomeNum << endl;
		assert(ForbidNum == acc_tester.GenomeNum);
		for(int i=0;i<ForbidNum;++i)
			ForbidTaxid.insert(Fid2Taxid[Forbid[i]]);
		cerr << "# Forbidden Taxids:" << ForbidTaxid.size() << endl;
		for(int i=0;i<Ctgs.CtgNum;++i){
			char tbuf[100];
			sscanf(Ctgs.info[i].c_str(),"%s\t%d/%d\t%d", tbuf,&spInfoOfCtgs[i].matchedLength,&spInfoOfCtgs[i].length,&spInfoOfCtgs[i].spFId);
			spInfoOfCtgs[i].spTaxId = Fid2Taxid[spInfoOfCtgs[i].spFId];
		}
		map<int,int>Fid2Lid;
		for(int i=0;i<100;++i)
			Fid2Lid[Forbid100[i]] = i/25;
		for(int i=0;i<Ctgs.CtgNum;++i)
			spInfoOfCtgs[i].refLevel = Fid2Lid[spInfoOfCtgs[i].spFId];
		cerr << "Taxoid of forbid: "<<Fid2Taxid[Forbid[0]] << '\t'<<Fid2Taxid[Forbid[1]]<<endl;
	}

	int VCtgNum = metacluster.Size;
	int ClustNum = 0;
	for(int i=0;i<VCtgNum;++i)
		ClustNum = max(metacluster.best[i],ClustNum);
	++ClustNum;
	vector<int>toClustId(Ctgs.CtgNum,0);
	for(int i=0;i<Ctgs.CtgNum;++i)
		toClustId[i] = metacluster.best[mcpara.getNewId(uset.find(i))];
	vector<vector<int> >CtgIdInCluster(ClustNum);
	for(int i=0;i<Ctgs.CtgNum;++i)
		CtgIdInCluster[toClustId[i]].push_back(i);
	TaxoOfClust.resize(ClustNum);
	vector<int>ClustLength(ClustNum);
	for(int i=0;i<CtgIdInCluster.size();++i)
		for(vector<int>::const_iterator itr2 = CtgIdInCluster[i].begin();itr2!=CtgIdInCluster[i].end();++itr2)
			ClustLength[i] += Ctgs.contigs[*itr2]->str.length();
	///////////////////////////////////////////////////////////////////////
	if(INTEST){printtime("Before initiating bwts. ");system("ps ux");}
	BWTs bwts(argv[3]);
	NodesDmp.init(argv[4]);
	if(INTEST){printtime("Before Annotating contigs:");system("ps ux");}
	vector<map<pair<int,int>,double> >taxid_score(Ctgs.CtgNum);
#pragma omp parallel for schedule(dynamic)
	for(int i=0;i<Ctgs.CtgNum;++i)
		calTaxoForCtg(bwts, Ctgs.contigs[i]->str, taxid_score[i]);
//		calTaxoForCtg(bwts, Ctgs.contigs[i]->str, TaxoInfo[i], NodesDmp,taxid_score[i]);
	bwts.clear();
	///////////////////////////////////////////////////
	////cal lamda
	map<int,double> TaxonLamda;
	map<int,double> TaxonPrec;
	getLamda(NodesDmp,(string(argv[3])+".info"),argv[5],TaxonLamda,TaxonPrec);
	///////////////////////////////////////////////////
	//anotate clusters
	//////////////annotate clusters from score
for(double P_CLUST = 0;P_CLUST <=0.9;P_CLUST+=0.1){
	for(int i=0;i<CtgIdInCluster.size();++i){
		map<pair<int,int>,double> ClustTaxidScore;
		for(vector<int>::const_iterator itr1 = CtgIdInCluster[i].begin(); itr1!=CtgIdInCluster[i].end();++itr1)
			for(map<pair<int,int>,double>::const_iterator itr = taxid_score[*itr1].begin();itr!=taxid_score[*itr1].end();++itr)
				ClustTaxidScore[itr->first] += itr->second;
		double maxScore = 0;
		for(map<pair<int,int>,double>::const_iterator itr = ClustTaxidScore.begin();itr!=ClustTaxidScore.end();++itr)
			maxScore = max(maxScore,itr->second);
		set<int>cand_taxo;
		for(map<pair<int,int>,double>::const_iterator itr = ClustTaxidScore.begin();itr!=ClustTaxidScore.end();++itr)
			if(itr->second >= maxScore*(1.0-P_CLUST))
				cand_taxo.insert(itr->first.second);
		vector<vector<int> >V;
		for(set<int>::const_iterator itr=cand_taxo.begin();itr!=cand_taxo.end();++itr)
			V.push_back(NodesDmp.getTaxo(*itr));
		if(V.size()==0){
			cout << "cluster info:\t";
			cout << i << "\t" << 0 << '\t' << 9 << endl;
			continue;
		}
		int taxon = 0,level = 99;
		for(int k=0;k<NL;++k){
			bool allSame = true;
			for(int j=1;j<V.size() && allSame;++j){
				if(V[j][L8[k]] != V[0][L8[k]])
					allSame = false;
			}
			if(allSame){
				taxon = V[0][L8[k]];
				level = k;
				break;
			}
		}
		//////////////////////////////////////
		//use lamda to decide taxon level.
		double candScore = maxScore/ClustLength[i];
		for(int k=level;k<NL-1;++k){
			if(TaxonLamda.find(taxon)!=TaxonLamda.end()){
		//		if(candScore > TaxonLamda[taxon])
				if(candScore >= 0)
					break;
				else{
					taxon = V[0][L8[k+1]];
					level = k+1;
				}
			}
		}
		/////////////////////////////////////
		TaxoOfClust[i].taxon = taxon;
		TaxoOfClust[i].level = level;
//		cout << i << "\tcluster taxon:\t" << taxon << '\t' << level << endl;
	}
	{
		cout << "P_CLUST = " << P_CLUST << endl;
		cerr << "cal taxonomy for clusters" << endl;
		int Evaluate[10][10][10];
		for(int i=0;i<10;++i)
			for(int j=0;j<10;++j)
				for(int k=0;k<10;++k)
					Evaluate[i][j][k] = 0;
		for(int i=0;i<CtgIdInCluster.size();++i){
			int taxon =  TaxoOfClust[i].taxon;
			int level =  TaxoOfClust[i].level;
			vector<int> predTaxon = NodesDmp.getTaxo(taxon);
			int predLevel = 0;
			for(;predLevel<NL;++predLevel)
				if(predTaxon[L8[predLevel]]>0)
					break;
	
			for(int i1=0;i1<CtgIdInCluster[i].size();++i1){
				vector<int> correctTaxon = NodesDmp.getTaxo(spInfoOfCtgs[CtgIdInCluster[i][i1]].spTaxId);
				int correctLevel = 0;
				for(;correctLevel<NL;++correctLevel)
					if(predTaxon[L8[correctLevel]]==correctTaxon[L8[correctLevel]] && predTaxon[L8[correctLevel]])
						break;
				Evaluate[spInfoOfCtgs[CtgIdInCluster[i][i1]].refLevel][predLevel][correctLevel]+= spInfoOfCtgs[CtgIdInCluster[i][i1]].NReads;
			}
		}
		cout << "evaluate: " << endl;
		cerr << "eva cluster: " << endl;
		for(int l=0;l<10;++l){
			long long sum = 0;
			for(int i=0;i<=NL;++i)
				for(int j=0;j<=NL;++j)
					sum += Evaluate[l][i][j];
			if(sum==0)continue;
			cout << "reference level:" << l << endl;
			for(int i=0;i<=NL;++i){
				for(int j=0;j<=NL;++j)
					cout << Evaluate[l][i][j] << '\t';
				cout << endl;
			}
			cout << endl<< endl;

			int lower = 0,high=0,corr=Evaluate[l][l][l];
			for(int i=0;i<l;++i)
				for(int j=0;j<=l;++j)
					lower += Evaluate[l][i][j];
			for(int i=l+1;i<NL;++i)
				high += Evaluate[l][i][i];
			cerr << lower << '\t' << high << '\t' << corr << '\t' <<sum << endl;
		}
	}
}
	//////////////annotate clusters from LCA contigs
/*	for(int i=0;i<CtgIdInCluster.size();++i){
		set<int>cand_taxo;
		for(int i1=0;i1<CtgIdInCluster[i].size();++i1){
			int i2 = CtgIdInCluster[i][i1];
			double maxScore = 0;int maxtaxon = 0;
			for(map<pair<int,int>,double>::const_iterator itr = taxid_score[i2].begin();itr!=taxid_score[i2].end();++itr){
				maxScore = max(maxScore,itr->second);
				maxtaxon = itr->first.second;
			}
			cand_taxo.insert(maxtaxon);
		}
		if(cand_taxo.size()==0)
			continue;
		vector<vector<int> >V;
		for(set<int>::const_iterator itr=cand_taxo.begin();itr!=cand_taxo.end();++itr)
			V.push_back(NodesDmp.getTaxo(*itr));
		int taxon = 0,level = 99;
		for(int k=0;k<NL;++k){
			int non0taxon = 0;
			for(int j=0;j<V.size();++j)
				non0taxon = max(0,V[j][L8[k]]);
			bool allSame = true;
			for(int j=0;j<V.size() && allSame;++j){
				if(V[j][L8[k]] && V[j][L8[k]] != non0taxon)
					allSame = false;
			}
			if(allSame){
				taxon = non0taxon;
				level = k;
				break;
			}
		}
	//	TaxoOfClust[i].taxo[itr2->first].set(maxTaxoid,maxScore/ClustLength[cidx]);
	//	TaxoOfClust[i].alignscore[itr2->first]=maxScore;
		TaxoOfClust[i].taxon = taxon;
		TaxoOfClust[i].level = level;
		cout << i << "\tcluster taxon:\t" << taxon << '\t' << level << endl;
	}*/
	if(INTEST){
		////////////////////////////////////////////////////////////////////
		///////////////////evaluate contigs
for(double P1 = 0;P1 <=0.9;P1+=0.1){
	cout << "P1 = " << P1 << endl;
			int Evaluate[10][10][10];
			for(int i=0;i<10;++i)
				for(int j=0;j<10;++j)
					for(int k=0;k<10;++k)
						Evaluate[i][j][k] = 0;
			cerr << "cal taxonomy for contigs" << endl;
			for(int i=0;i<Ctgs.CtgNum;++i){
				double maxScore = 0;
				for(map<pair<int,int>,double>::const_iterator itr = taxid_score[i].begin();itr!=taxid_score[i].end();++itr)
					maxScore = max(maxScore,itr->second);
				set<int>cand_taxo;
				for(map<pair<int,int>,double>::const_iterator itr = taxid_score[i].begin();itr!=taxid_score[i].end();++itr){
					if(itr->second >= maxScore*(1.0-P1))
						cand_taxo.insert(itr->first.second);
				}
				vector<vector<int> >V;
				for(set<int>::const_iterator itr=cand_taxo.begin();itr!=cand_taxo.end();++itr)
					V.push_back(NodesDmp.getTaxo(*itr));
				if(V.size()==0){
					cout << "ctg info:\t" << spInfoOfCtgs[i].matchedLength << '\t' << spInfoOfCtgs[i].length << '\t' << spInfoOfCtgs[i].spFId << '\t' << spInfoOfCtgs[i].spTaxId <<'\t';
					cout << i << "\t" << 0 << '\t' << 9 << endl;
					continue;
				}
				int taxon = -1,level = 99;
				for(int k=0;k<NL;++k){
					bool allSame = true;
					for(int j=1;j<V.size() && allSame;++j){
						if(V[j][L8[k]] != V[0][L8[k]])
							allSame = false;
					}
					if(allSame){
						taxon = V[0][L8[k]];
						level = k;
						break;
					}
				}
				cout << "ctg info:\t" << spInfoOfCtgs[i].matchedLength << '\t' << spInfoOfCtgs[i].length << '\t' << spInfoOfCtgs[i].spFId << '\t' << spInfoOfCtgs[i].spTaxId <<'\t';
				cout << i << "\t" << taxon << '\t' << level << endl;
				///////////////////////////////
				vector<int> predTaxon = NodesDmp.getTaxo(taxon);
				int predLevel = 0;
				for(;predLevel<NL;++predLevel)
					if(predTaxon[L8[predLevel]]>0)
						break;
				vector<int> correctTaxon = NodesDmp.getTaxo(spInfoOfCtgs[i].spTaxId);
				int correctLevel = 0;
				for(;correctLevel<NL;++correctLevel)
					if(predTaxon[L8[correctLevel]]==correctTaxon[L8[correctLevel]] && predTaxon[L8[correctLevel]])
						break;
				Evaluate[spInfoOfCtgs[i].refLevel][predLevel][correctLevel]+= spInfoOfCtgs[i].NReads;
			}
			cout << "evaluate contigs: " << endl;
			for(int l=0;l<10;++l){
				long long sum = 0;
				for(int i=0;i<=NL;++i)
					for(int j=0;j<=NL;++j)
						sum += Evaluate[l][i][j];
				if(sum==0)continue;
				cout << "contig reference level by count:" << l << endl;
				for(int i=0;i<=NL;++i){
					for(int j=0;j<=NL;++j)
						cout << Evaluate[l][i][j] << '\t';
					cout << endl;
				}
				cout << endl<< endl;
			}cout << endl << endl << endl;
			for(int l=0;l<10;++l){
				long long sum = 0;
				for(int i=0;i<=NL;++i)
					for(int j=0;j<=NL;++j)
						sum += Evaluate[l][i][j];
				if(sum==0)continue;
				cout << "contig reference level by perc:" << l << endl;
				for(int i=0;i<=NL;++i){
					double tsum = 0;
					for(int j=0;j<=NL;++j)
						tsum += Evaluate[l][i][j];
					cout << tsum << "\t:\t";
					for(int j=0;j<=NL;++j)
						cout << Evaluate[l][i][j]/tsum << '\t';
					cout << endl;
				}
				cout << endl<< endl;
			}
		}
		////////////////////////////////////////////////////////////////////
		///////////////////evaluate clusters
	}
}

void printVirtualContigs(map<int,vector<int> >& cluster_ctg,string filename)
{
	const int insertionN = 100;
	char buf[1000];
	for(int i=0;i<insertionN;++i)
		buf[i]='N';
	buf[insertionN] = 0;
	string insertion(buf);

	ofstream ofs(filename.c_str());
	assert(!ofs.fail());
	for(map<int,vector<int> >::const_iterator itr = cluster_ctg.begin();itr!=cluster_ctg.end();++itr)
	{
		ofs << Ctgs.info[*(itr->second.begin())] <<'\t';
		for(vector<int>::const_iterator itr2 = itr->second.begin();itr2!=itr->second.end();++itr2)
			ofs << *itr2 << '_';
		ofs << endl;
		for(vector<int>::const_iterator itr2 = itr->second.begin();itr2!=itr->second.end();++itr2)
		{
			if(itr2 != itr->second.begin())
			   	ofs << insertion;
			ofs << (Ctgs.contigs[*itr2]->str);
		}
		ofs << endl;
	}
	ofs.close();
}

void printClusterIdOfCtgs(string filename)
{
	if(INTEST){
		map<int,vector<int> > cluster_ctg;
		for(int i=0;i<Ctgs.CtgNum;++i)
			cluster_ctg[uset.find(i)].push_back(i);
		/*
		cout << "cluster id after step1: " << Ctgs.CtgNum << '\t' << cluster_ctg.size()<<endl;
		for(map<int,vector<int> >::const_iterator itr = cluster_ctg.begin();itr!=cluster_ctg.end();++itr)
		{
			cout << (itr->first) << '\t' << (itr->second.size()) << ":\t";
			for(vector<int>::const_iterator itr2 = itr->second.begin();itr2!=itr->second.end();++itr2)
				cout << *itr2 << '\t';
			cout << endl;
		}
		cout << "<end> cluster id after step1: " << Ctgs.CtgNum<<endl;
		*/
		printVirtualContigs(cluster_ctg,filename);
	}
}

int main(int argc, char* argv[])
{
	if(argc < 6)
		usage();
	init();
	input(argc,argv);
	if(INTEST)system("ps ux");

	printtime("Main: before MergeAsStep1. ");
	MergeAsStep1(Ctgs,Reads,NodePool,uset,acc_tester,AlignThresh);
	NodePool.clear();
	if(INTEST)acc_tester.calAcc(uset);

	printtime("Main: before MetaCluster. ");
	if(INTEST)system("ps ux");
	////////////////////////////////////////
	if(INTEST)
	{
		printClusterIdOfCtgs((string(argv[1])+".step1.contig.fa").c_str());
/*		MCPara mcpara(MetaKmerLen, Ctgs, uset, acc_tester);
		MetaCluster metacluster(mcpara.KmerLen, mcpara.Size, mcpara.ReverKmerDistri, ClusterSize, MaxSpecies, MinSpecies, mcpara.GenoNum , mcpara.Component);
		ClusterSize? metacluster.muiltkmeans(10,ClusterSize) : metacluster.muiltkmeans(10,500);
		metacluster.clear();//best[]&Component[][] are not cleared.
		*/
	}
	////////////////////////////////////////
	MCPara mcpara(MetaKmerLen, Ctgs, uset, acc_tester);
	MetaCluster metacluster(mcpara.KmerLen, mcpara.Size, mcpara.ReverKmerDistri, ClusterSize, MaxSpecies, MinSpecies, mcpara.GenoNum , mcpara.Component);
	if(INTEST)system("ps ux");
	ClusterSize ? (metacluster.muiltkmeans(10,ClusterSize)) : (metacluster.iterMeta(10,MC3_Thresh));
	metacluster.clear();//best[]&Component[][] are not cleared.
	if(INTEST)system("ps ux");
	printtime("Main: before anaCluster. ");
	anaCluster(mcpara, metacluster, Reads.ReadNumInCtg,argv);
	printtime("Main: before output. ");
	outputResult(Reads.MatchId, Reads.TotalNum, uset,metacluster.best, mcpara,argv[2], TaxoOfClust);
	if(INTEST)
	{
		system("ps ux");
		printtime("Clustering finished. ");
		acc_tester.getPreSen4Other(metacluster.getComp(),metacluster.Size,metacluster.best);
/*		for(int i=0;i<TaxoOfClust.size();++i)
		{
			cout << i << '\t';
			for(int level=0;level<7;++level)
			{
				cout << TaxoOfClust[idx].taxo[level].taxid<<'\t';
				cout << TaxoOfClust[idx].taxo[level].score<<'\t';
				cout << TaxoOfClust[idx].alignscore[level]<<'\t';
			}
		}*/
		map<int,vector<int> > cluster_ctg;
		for(int i=0;i<Ctgs.CtgNum;++i)
			cluster_ctg[metacluster.best[mcpara.getNewId(uset.find(i))]].push_back(i);
		cout << "cluster id after step3: " << Ctgs.CtgNum << '\t' << cluster_ctg.size()<<endl;
		for(map<int,vector<int> >::const_iterator itr = cluster_ctg.begin();itr!=cluster_ctg.end();++itr)
		{
			cout << (itr->first) << '\t' << (itr->second.size()) << ":\t";
			for(vector<int>::const_iterator itr2 = itr->second.begin();itr2!=itr->second.end();++itr2)
				cout << *itr2 << '\t';
			cout << endl;
		}
		cout << "<end> cluster id after step3: " << Ctgs.CtgNum<<endl;
		printVirtualContigs(cluster_ctg,(string(argv[1])+".step3.contig.fa"));
	}
	cerr << "Finished." << endl;
	return 0;
}
