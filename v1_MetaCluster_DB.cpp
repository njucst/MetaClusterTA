/*
 * use scores of contigs to predict scores of clusters directly
 */
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <omp.h>

/*
#include "GLOBAL.h"
#include "Utils.h"
#include "CtgsReads.h"
#include "BaseStr.h"
#include "AccTester.h"
#include "MetaCluster.h"
#include "USet.h"
*/

//#include "GenomeDB_KmerDistri.h"
#include "Methods.h"
#include "MCPara.h"
#include "Ncbi_nodes_dmp.h"
//#include "DBEntropy.h"
using namespace std;
/////////////////////////////////////////////////////////////////////
static const unsigned UNASSIGNED  = ~(0U);
const double ScoreThresh = -1.0;
//parameters
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
time_t rawtime;
struct tm* timeinfo;

KmerNode** KmerMap;
ContigsClass Ctgs;
ReadsClass Reads;
KmerNodeAloc NodePool;
USet uset;
NCBI_nodes_dmp NodesDmp;
//GenomesClass GenoDB;
//DBEntropy GenoDB;
vector<map<int,double> >TaxoInfo;
vector<ClustTaxoInfoClass> TaxoOfClust;

AccTester acc_tester;

////////////////////////////////////////////////////////////////////////////
void usage()
{
	cerr << "usage:\t MetaCluster_HB contig.fa read.fa bwt.idx nodes.dmp [options]" << endl;
	cerr << endl;
	exit(-1);
}

void printtime(string str = "")
{
	cerr << str << endl;
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	cerr<<asctime(timeinfo);
}

void init()
{
	///////////////////////////////////////////////////////
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
	///////////////////////////////////////////////////////
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
	/////////////////////////////////////////////////////////////////
	if(INTEST)printtime("before loading contigs. ");
	Ctgs.init(argv[1], CtgLenThresh,KmerMap, NodePool);
	if(INTEST)
	{
		cerr << "# of Contigs: \t" << Ctgs.CtgNum << endl;
		/*
		map<int,int>MId;
		map<int,int>MOc;
		omp_lock_t tlock;
		omp_init_lock(&tlock);
#pragma omp parallel for
		for(unsigned i=0;i<HASHSIZE;++i)
		{
			for(KmerNode* p = KmerMap[i];p!=NULL;p=p->next)
			{
				set<int>ids;
				for(int j=0;j<p->VSize;++j)
					ids.insert(p->myvector[j].id);
				omp_set_lock(&tlock);
				++MId[ids.size()];
				++MOc[p->VSize];
				omp_unset_lock(&tlock);
			}
		}
		cout << "Id Occurence: " << endl;
		for(map<int,int>::const_iterator itr=MId.begin();itr!=MId.end();++itr)
			cout << itr->first << '\t' << itr->second << endl;
		cout << "Kmer Occurence: " << endl;
		for(map<int,int>::const_iterator itr=MOc.begin();itr!=MOc.end();++itr)
			cout << itr->first << '\t' << itr->second << endl;
		*/
		printtime("Contigs loaded. ");
		///////////////////////////////////////////////////////////
	}
	NodePool.fixCtgNum();
	if(INTEST)cerr<<"Ctg Num fixed." << endl;
	////////////////////////////////
//	GenoDB.init(argv[3], KmerMap, NodePool);
/*	GenoDB.init(argv[3]);
	system("ps ux");
	if(INTEST)cerr<<"Genome DB loaded: "<< Ctgs.CtgNum << endl;
	*/
/*	TaxoInfo.resize(Ctgs.CtgNum);
	if(INTEST)cerr<<"Before Annotating contigs:" << endl;
#pragma omp parallel for
	for(int i=0;i<Ctgs.CtgNum;++i)
	{
		unsigned taxid=UNASSIGNED;double score;
		unsigned simulateId;double simulateScore;
		vector<double>entropy;
		GenoDB.calTaxoForCtg(Ctgs.contigs[i]->str, taxid, score,entropy,simulateId,simulateScore);
		if(score < 0)taxid = UNASSIGNED;
		TaxoInfo[i].set(taxid,score);
		if(score>=0)
		{
			cout << i << '\t' << taxid << '\t' << score <<'\t'<<simulateScore << '\t' << GenoDB.GenomeLength[taxid] << '\t';
			cout << GenoDB.GenomeLength[simulateId] << '\t';
			for(int j=0;j<7;++j)cout << GenoDB.Taxo[taxid][j] << '\t';
			for(int j=0;j<7;++j)cout << entropy[j] << '\t';
			cout << Ctgs.contigs[i]->length <<'\t'<< Ctgs.info[i] << endl;
		}
	}*/
	/////////////////////////////////////////////////////////////////
	KmerDistriPara Para5mer(5);
	int** ctgSpear = Ctgs.getSpear(Para5mer);
	for(int i=0;i<Ctgs.CtgNum;++i)
	{
		ctgSpear[i] = NULL;
	}
	////////////////////////////////
	printtime("initializing reads. ");
	Reads.init(argv[2], ReadLen, AlignThresh, KmerMap, NodePool, Ctgs, acc_tester);
	{
		printtime("initializing uset: ");
		int CtgNum = Ctgs.CtgNum;
		unsigned* ctglen = new unsigned[CtgNum];
		for(int i=0;i<CtgNum;++i)
			ctglen[i] = Ctgs.contigs[i]->length;
		uset.init(Ctgs.CtgNum + Reads.ReadNum, CtgNum, ctglen);
		delete[]ctglen;
	}

	if(INTEST)
	{
		printtime("Reads loaded. ");
		cerr << dec << "Total Reads:\t"<< Reads.TotalNum << endl;
		cerr << "Unmaped Reads:\t" << Reads.ReadNum << endl;
		cerr << "AlignThresh:\t" << Reads.AlignThresh << endl;
		///////////////////////////////////////////////////////////////////
		printtime("initializing acc_tester:\t");
		acc_tester.init(Reads.TotalNum, Reads.MatchId, Reads.ReadNum,Reads.NewIdToOldId);
	/*	acc_tester.calAcc(uset);*/
	}
	delete[] KmerMap;
	KmerMap = NULL;
//	NodePool.shrinkSize(NULL,2);
	//////////////////////////////////////////////////////////
}

void anaCluster(MCPara& mcpara, MetaCluster& metacluster,vector<unsigned>& ReadNumInCtg,char* argv[])
{
	int VCtgNum = metacluster.Size;
	int ClustNum = 0;
	for(int i=0;i<VCtgNum;++i)
		ClustNum = max(metacluster.best[i],ClustNum);
	++ClustNum;
	vector<int>toClustId(Ctgs.CtgNum,0);
	for(int i=0;i<Ctgs.CtgNum;++i)
		toClustId[i] = metacluster.best[mcpara.getNewId(uset.find(i))];
	cout << "Cluster Id of Ctgs:\t"<< endl; 
	for(int i=0;i<Ctgs.CtgNum;++i)
		cout << toClustId[i] << endl;
	vector<vector<int> >CtgIdInCluster(ClustNum);
	for(int i=0;i<Ctgs.CtgNum;++i)
		CtgIdInCluster[toClustId[i]].push_back(i);
	TaxoOfClust.resize(ClustNum);
	vector<int>ClustLength(ClustNum);
	for(int i=0;i<CtgIdInCluster.size();++i)
	{
		int len = 0;
		for(vector<int>::const_iterator itr2 = CtgIdInCluster[i].begin();itr2!=CtgIdInCluster[i].end();++itr2)
			len += Ctgs.contigs[*itr2]->str.length();
		ClustLength[i] = len;
	}
	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////
	TaxoInfo.resize(Ctgs.CtgNum);
	if(INTEST)printtime("Before initiating contigs. ");
	BWTs bwts(argv[3]);
	NodesDmp.init(argv[4]);
	if(INTEST)printtime("Before Annotating contigs:");
#pragma omp parallel for
	for(int i=0;i<Ctgs.CtgNum;++i)
		calTaxoForCtg(bwts, Ctgs.contigs[i]->str, TaxoInfo[i]);

	/////////////////////////////////////////////////////////////////
	if(INTEST)
		for(int i=0;i<Ctgs.CtgNum;++i)
		{
			cout << " Contig score: " << i << '\t' << Ctgs.contigs[i]->str.length() << '\t';
			for(map<int,double>::const_iterator itr=TaxoInfo[i].begin();itr!=TaxoInfo[i].end();++itr)
				cout << itr->first << ':' << itr->second << '\t';
			cout << endl;
		}
	//calculate scores
	{
		vector<double>ClusterTotalScore(Ctgs.CtgNum,0);
		for(int i=0;i<CtgIdInCluster.size();++i)
			for(int j=0;j<CtgIdInCluster[i].size();++j)
				for(map<int,double>::const_iterator itr=TaxoInfo[CtgIdInCluster[i][j]].begin();itr!=TaxoInfo[CtgIdInCluster[i][j]].end();++itr)
					ClusterTotalScore[i] += itr->second;
		//////////////////////////////////////////////////////////////////
		if(INTEST)printtime("Before calculating scores:");
		//////////////////////////////////////////////////////////////////
		int cidx = 0;
		for(vector<vector<int> >::const_iterator itr1= CtgIdInCluster.begin();itr1!=CtgIdInCluster.end();++itr1,++cidx)
		{
			map<int,map<int,double> >Node_Score;//TaxonLevel,TaxId,Score
			for(vector<int>::const_iterator itr2 = itr1->begin();itr2!=itr1->end();++itr2)
				for(map<int,double>::const_iterator itr3 = TaxoInfo[*itr2].begin();itr3!=TaxoInfo[*itr2].end();++itr3)
					NodesDmp.addScore2Path(itr3->first, itr3->second, Node_Score);
	
			for(map<int,map<int,double> >::const_iterator itr2= Node_Score.begin();itr2!=Node_Score.end();++itr2)
			{
				int maxTaxoid = 0;double maxScore = 0;
				double TotalScore = ClusterTotalScore[cidx];
				for(map<int,double>::const_iterator itr3 = itr2->second.begin();itr3!=itr2->second.end();++itr3)
				{
					if(itr3->second > maxScore && itr3->first > 0)
					{
						maxScore = itr3->second;
						maxTaxoid = itr3->first;
					}
		//			if(itr3->first > 0)TotalScore += itr3->second;
				}
				TaxoOfClust[cidx].taxo[itr2->first].set(maxTaxoid,maxScore/TotalScore);
				TaxoOfClust[cidx].alignscore[itr2->first]=maxScore;

				cout << "cluster temp score: " << cidx <<'\t';
				cout <<  maxTaxoid << '\t' << maxScore << '\t' << TotalScore << '\t' << ClustLength[cidx] <<'\t';
				cout << NodesDmp.Id2LevelName[itr2->first] << endl;
			}
		}
		cidx = 0;
		for(vector<ClustTaxoInfoClass>::const_iterator itr=TaxoOfClust.begin();itr!=TaxoOfClust.end();++itr,++cidx)
		{
			cout <<"cluster score: " << cidx << '\t';
			for(int i=0;i<40;++i)
				if(itr->taxo[i].score > 1e-9)
					cout <<NodesDmp.Id2LevelName[i]<<":"<<(itr->taxo[i].score)<<'\t';
			cout << endl;
		}
		if(INTEST)printtime("End of calculating scores:");
	}
}

int main(int argc, char* argv[])
{
	if(argc < 4)
		usage();
	init();
	input(argc,argv);

	printtime("Main: before MergeAsStep1. ");
	MergeAsStep1(Ctgs,Reads,NodePool,uset,acc_tester,AlignThresh);

	if(INTEST)acc_tester.calAcc(uset);

	printtime("Main: before MetaCluster. ");
	MCPara mcpara(MetaKmerLen, Ctgs, uset, acc_tester);
	MetaCluster metacluster(mcpara.KmerLen, mcpara.Size, mcpara.ReverKmerDistri, ClusterSize, MaxSpecies, MinSpecies, mcpara.GenoNum , mcpara.Component);
	ClusterSize ? (metacluster.muiltkmeans(10,ClusterSize)) : (metacluster.iterMeta(10,MC3_Thresh));
	printtime("Main: before anaCluster. ");
	anaCluster(mcpara, metacluster, Reads.ReadNumInCtg,argv);
	/* ToDo: clear memory
	 * because outputResult needs O(n)memory, where n is # reads.
	 */
	printtime("Main: before output. ");
	outputResult(Reads.MatchId, Reads.TotalNum, uset,metacluster.best, mcpara,argv[2], TaxoOfClust);

	if(INTEST)
	{
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
	}
	return 0;
}
