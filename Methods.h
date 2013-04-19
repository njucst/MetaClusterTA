/*short methods for MetaCuster
 * by Wang Yi
 * wangyi.tom@gmail.com
 * Nov 9,2012
 */
#ifndef MCH_HYBRID_METHODS_H_

#define MCH_HYBRID_METHODS_H_

#include <iostream>
#include "CtgsReads.h"
#include "MetaCluster.h"
#include "AccTester.h"
#include "SparseMatrix.h"
//#include "DBEntropy.h"
#include "MCPara.h"
#include "BWTs.h"
using namespace std;

bool cmp(const KmerNode* p1, const KmerNode* p2);

//void calTaxoForCtg(const BWTs &bwts, const string& str, map<int,double>& sp_score);
void calTaxoForCtg(const BWTs& bwts, const string& str, map<int,double>& sp_score);
void compa_read(int ReadLen,const int &position1,const int &position2,int &match,int &mismatch,const ULLN &read1,const ULLN &read2);
void MergeReads(const KmerNodeAloc& NodePool, const ReadsClass& Reads, USet& uset, int CtgNum);
void MergeCtgReads(const KmerNodeAloc& NodePool, const ReadsClass& Reads, const ContigsClass& Ctgs, USet& uset);
void MergeAsStep1(const ContigsClass& Ctgs,const ReadsClass& Reads, KmerNodeAloc& NodePool, USet& uset, AccTester& acc_tester, int pair_end_merge_threshold);

void getScoreV(int** GenoDBTaxo, const vector<int>&CtgId, const int taxoLevel, int& AnsLevel, double& AnsScore);
void getTaxo4Clust(int** GenoDBTaxo, const vector<int>&CtgId, double Thresh, int& AnsLevel, double& AnsScore);
double getScoreEntropy(const map<int,long long> &TaxoComp);
double getScoreMax(const map<int,long long> &TaxoComp);
double getScoreMax2(const map<int,long long> &TaxoComp);
void outputResult(int* MatchId, long long ReadN, USet& uset,int* metabest, MCPara& mcpara,string infile,vector<ClustTaxoInfoClass>& taxoofclust);
#endif
