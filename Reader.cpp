#include "Reader.h"
//#include "global.h"
#include "TomAlgorithm.h"

#include <iostream>
#include <string>
#include <cstring>
#include <cmath>
#include <time.h>
//#include <algorithm>
#include <cstdio>

using namespace std;

GenomeReader::GenomeReader()
{
	Total5merSpear = NULL;
	digits = NULL;
	fp= NULL;
	allDistri=NULL;
	distriNum=0;
	glength = 0;
	line = NULL;
}
GenomeReader::GenomeReader(const string &filename)
{
	////////////////////////////////
	Total5merSpear = NULL;
	///////////////////////////////
	digits = NULL;
	fp = fopen(filename.c_str(),"rt");

	allDistri=NULL;
	distriNum=0;
	glength = 0;
	line = NULL;

	if(fp == NULL)
	{
		cerr << "can not open " << filename << endl;
		throw exception();
	}
	line = new char[MAXLENGTH];
	char* buf = new char[MAXLINE];
	fgets(buf,MAXLINE,fp);
	if(buf[0]!='>')
	{
		cerr << "First line is in wrong format:\t"<<buf<<endl;
		throw exception();
	}
	description = string(buf);
	while(fgets(buf, MAXLINE, fp) != NULL)
	{
		if(buf[0] != '>')
		{
			if(buf[strlen(buf)-1]=='\n')
				buf[strlen(buf)-1]=0;
			strcpy(line+glength,buf);
			glength += strlen(buf);
		}
		else
		{
			cerr << "Wrong format:\t"<<buf<<endl;
			throw exception();
		}
	}
	delete[] buf;
	fclose(fp);
	fp=NULL;
}

GenomeReader::~GenomeReader()
{
	if(fp != NULL)
		fclose(fp);
	if(line!=NULL)
		delete[]line;
	if(digits!=NULL)
		delete[]digits;
	if(distriNum>0)
	{
	  for(int i=0;i<distriNum;++i)
	  {
  			delete[]allDistri[i];
	  }delete[]allDistri;
	  distriNum=0;
	}
}

void GenomeReader::clearString()
{
	if(line!=NULL)
	{
		delete[]line;
		line = NULL;
	}
}
void GenomeReader::clear()
{
	if(fp!=NULL)
		fclose(fp);
	if(line!=NULL){
		delete[]line;
		line=NULL;}
	if(digits!=NULL){
		delete[]digits;
		digits=NULL;}
	glength = 0;
	description = "";
	if(distriNum>0)
	{
	  for(int i=0;i<distriNum;++i)
	  {
  		delete[]allDistri[i];
	  }delete[]allDistri;
	  allDistri = NULL;
	  distriNum=0;
	}
}

bool GenomeReader::reRead(const string &filename)
{
	clear();
	fp = fopen(filename.c_str(),"rt");
	if(fp == NULL)
	{
		cerr << "can not open " << filename << endl;
		throw exception();
	}
	line = new char[MAXLENGTH];
	char* buf = new char[MAXLINE];
	glength = 0;
	fgets(buf,MAXLINE,fp);
	if(buf[0]!='>')
	{
		cerr << "First line is in wrong format:\t"<<buf<<endl;
		throw exception();
	}
	description = string(buf);
	while(fgets(buf, MAXLINE, fp) != NULL)
	{
		if(buf[0] != '>')
		{
			if(buf[strlen(buf)-1]=='\n')
					buf[strlen(buf)-1]=0;
			strcpy(line+glength,buf);
			glength += strlen(buf);
		}
		else
   		{
			cerr << "Wrong format:\t"<<buf<<endl;
			throw exception();
		}
	}
	delete[] buf;
	fclose(fp);
	fp=NULL;
	return true;
}

bool GenomeReader::createDigits()
{
	if(glength==0)
		return false;
	if(digits!=NULL)return false;
	digits = new int[glength];
	for(int i=0;i<glength;++i)
	{
		switch(line[i])
		{
			case 'A':digits[i]=0;break;
			case 'C':digits[i]=1;break;
			case 'G':digits[i]=2;break;
			case 'T':digits[i]=3;break;
			default :digits[i]=0-line[i];break;
		}
	}
	return true;
}
int** GenomeReader::getAllDistri(const string &seed,const int fragLen)
{
	//////////////////////////
	if(distriNum==0)
	{
		for(int i=0;i<distriNum;++i)
			delete[]allDistri[i];
		delete[]allDistri;
	}
	distriNum = glength/fragLen;
	allDistri = new int*[distriNum];
	for(int i=0;i<distriNum;++i)
	{
		allDistri[i]=getSeedKDistr(seed,i*fragLen,(i+1)*fragLen,true);
	}
	return allDistri;
	/////////remove duplicate caused by reverse complement
  /*  int *positions = new int[distriSize];
	spearSize = 0;
	for(int i=0;i<distriSize;++i)
	{
		if(tomReverComple(i,seedWeight)>=i)
		{
			positions[i] = spearSize++;
		}
		else
			positions[i] = -1;
	}
	for(int i=0;i<fragNum;++i)
	{
		for(int j=0;j<distriSize;++j)
		{
			if(positions[j]!=-1)
			{
				allDistri[i][positions[j]]=allDistri[i][j];
			}
		}
	}
	delete[] positions;*/
}

double GenomeReader::selfDistance(const string &seed,const int fragLen,double &average,double &stdev,bool newspear=false)
{
	////////////////get time
  //  time_t ltime;char tmpbuf[128];
 //   time(&ltime);cerr<<ctime(&ltime);
	//////////////////////////
	getAllDistri(seed,fragLen);
	int fragNum = distriNum;
	int seedWeight = 0;
	for(int i=0;i<seed.size();++i)
	{
		if(seed[i]=='1')++seedWeight;
	}
//	int distriSize = pow(4,seedWeight);
	int distriSize = 1<<(seedWeight<<1);
	int spearSize = distriSize;
	/////Use spearman here
	int **allSpear = new int *[fragNum];
	if(newspear)
	{
	for(int i=0;i<fragNum;++i)
	{
		allSpear[i]=toNewSpear(allDistri[i],spearSize);
		delete[]allDistri[i];
	}delete[]allDistri;
	}
	else{
	for(int i=0;i<fragNum;++i)
	{
		allSpear[i]=toSpear(allDistri[i],spearSize);
		delete[]allDistri[i];
	}delete[]allDistri;
	}
 //   time(&ltime);cerr<<ctime(&ltime);

	double disSum=0;
	double disSqareSum=0;
 //   cout<<fragNum<<'\t'<<distriSize;
	for(int i=0;i<fragNum-1;++i)
	{
		for(int j=i+1;j<fragNum;++j)
		{
			double distance = 0;
			for(int k=0;k<spearSize;++k)
			{
				distance+=abs(allSpear[i][k]-allSpear[j][k]);
			}
			disSum+=distance;
			disSqareSum+=distance*distance;
		}
	}
 //   time(&ltime);cerr<<ctime(&ltime);
	for(int i=0;i<fragNum;++i)
	{
		delete[]allSpear[i];
	}delete[]allSpear;
	
	double n = fragNum*(fragNum-1)/2;
	average = disSum/n;//(fragNum*(fragNum-1));
	stdev = sqrt((disSqareSum*n-disSum*disSum)/(n*(n-1)));
	return average;
}
///////////kmerLen is the weight of seeda
//////////get k seed distribution of fragment ragion[startpositon,endposition)
int* GenomeReader::getSeedKDistr(const string &seed,int startposition,int endposition,bool revercom)
{
	if(digits==NULL)createDigits();
	if(endposition-startposition<seed.size())return NULL;
	const int seedLen=seed.size();
	int kmerLen = 0;
	for(int i=0;i<seedLen;++i)
		if(seed[i]=='1')++kmerLen;
	/////////////////get index for the seed
	int* seedindex = new int[kmerLen];
	int index = 0;
	for(int i=0;i<seedLen;++i)
	{
		if(seed[i]=='1')
			seedindex[index++]=i;
	}
	/////////////////
//	int distrisize = pow(4,kmerLen);
	int distrisize = 1<<(kmerLen<<1);
	int* distribution = new int[distrisize];
	for(int i=0;i<distrisize;++i)
	{
		distribution[i]=0;
	}
	for(int i=startposition;i<=endposition-seedLen;++i)
	{
		int kmer = 0;//digits[i];
		for(int j=0;j<kmerLen;++j)
		{
			int curDig = digits[i+seedindex[j]];
			if(curDig < 0)
			{
				kmer = -1;
				break;
			}
			kmer <<= 2;
			kmer |= curDig;
		}
		if(kmer<0)continue;
		distribution[kmer]++;
	}
	///////reverse compliment
	if(revercom)
	{
//		const int all_3 = pow(4,kmerLen)-1;
		const int all_3 = (1<<(kmerLen<<1))-1;
		for(int i=endposition-1;i>=startposition+seedLen;--i)
		{
			int kmer = 0;//digits[i];
			for(int j=0;j<kmerLen;++j)
			{
				int curDig = digits[i-seedindex[j]];
				if(curDig < 0)
				{
					kmer = -1;
					break;
				}
				kmer <<= 2;
				kmer |= curDig;
			}
			if(kmer<0)continue;
			distribution[all_3-kmer]++;
		}
	}
	////////////////////////
	delete[]seedindex;
	return distribution;
}

//reverse complement is considered
int* GenomeReader::getKmerSpear(const int startposition, const int fragLen, const KmerDistriPara& Para)
{
	createDigits();
	int KmerLen = Para.KmerLen;
	int ReverSize = Para.ReverSize;
	int NonReverSize = Para.NonReverSize;
	unsigned KmerMask = Para.KmerMask;
	///////////////////////////////////////
	vector<int>nonRever(NonReverSize,0);
	int* ReverDis = new int[ReverSize];
	for(int i=0;i<ReverSize;++i)
		ReverDis[i] = 0;
	unsigned prek = 0;
	int* line_ = digits + startposition;
	for(int j=0;j<KmerLen-1;++line_)
	{
		prek <<= 2;
		if((*line_)>=0)
		{
			prek |= *line_;
			++j;
		}
	}
	int* lend =  digits+startposition+fragLen;
	for(;line_<lend;++line_)
	{
		if(*line >= 0)
		{
			prek <<= 2;
			prek |= *line_; 
			++nonRever[prek & KmerMask];
		}
	}
	for(int i=0;i<NonReverSize;++i)
	{
		if(tomReverComple(i,KmerLen)==i)
			ReverDis[Para.RCIdx[i]] = nonRever[i]*2;
		else
			ReverDis[Para.RCIdx[i]] += nonRever[i];
	}
	return ReverDis;
}

//reverse complement is considered
void GenomeReader::getAllKmerSpear(const KmerDistriPara& Para, const int FragLen, const int Overlap, int** &OutDis,int &FragNum )
{
	int KmerLen = Para.KmerLen;
	int ReverSize = Para.ReverSize;
	int NonReverSize = Para.NonReverSize;
	unsigned KmerMask = Para.KmerMask;
	//////////////////////////////////////////////////////////
	FragNum = (glength-Overlap)/(FragLen-Overlap);
	int PreNum = FragNum;
	if(FragNum*(FragLen-Overlap) < glength-Overlap)
		++FragNum;
	int** ans = new int*[FragNum];
	for(int i=0;i<PreNum;++i)
		ans[i] = getKmerSpear((FragLen-Overlap)*i, FragLen, Para);
	if(FragNum > PreNum)
		ans[PreNum] = getKmerSpear(glength-FragLen, FragLen, Para);
	OutDis = ans;
}

void GenomeReader::getDiffAllKmerSpear(const KmerDistriPara& Para)
{
	for(int i=0;i<GN;++i)
		getAllKmerSpear(Para, 10000*(i+1), 5000*(i+1),All5merSpear[i], SpearSize[i]);
	if(Total5merSpear==NULL)
		Total5merSpear = getKmerSpear(0,glength,Para);
}

double GenomeReader::intraDis(const KmerDistriPara& Para, const int FragLen)
{
	int **OutDis = NULL;
	int FragNum = 0;
	double ans = 0;
	getAllKmerSpear(Para, FragLen, 0, OutDis, FragNum);
	for(int i=0;i<FragNum;++i)
		for(int j=i+1;j<FragNum;++j)
		{
			double tmp = 0;
			for(int k=0;k<Para.ReverSize;++k)
				tmp += abs(OutDis[i][k]-OutDis[j][k]);
			ans += tmp;
		}
	for(int i=0;i<FragNum;++i)
		delete[]OutDis[i];
	delete[]OutDis;
	return ans*2/((FragNum-1)*FragNum);
}
