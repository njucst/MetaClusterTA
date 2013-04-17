#include <iostream>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <map>
#include <string>
using namespace std;
const int MAXL = 10000;
char Buf[MAXL];
const int MAXN = 200;
const int MAXG = 200;
const double taxoThresh = 0.7;
int toOriId[] = {25,26,27,0,1,75,6,7,8,9,10,11,12,13,4,5,14,15,16,18,19,28,29,30,31,32,33,34,35,17,36,37,38,39,76,20,77,78,79,21,22,23,40,41,24,42,43,44,50,51,52,53,54,55,56,57,58,59,60,61,62,2,3,45,46,47,64,63,49,48,65,66,67,99,68,98,69,97,70,71,72,73,74,95,96,87,88,89,90,91,92,93,94,85,86,84,83,82,80,81};
/*struct Clust
{
	int clustId,spId,numMaxRead,numReads;
	int taxoId[7];
	double taxoProb[7];
}cluster[MAXN];*/
int genomeTaxo[MAXG][7];
void usage()
{
	cerr << "./anaMegan treeTaxo.txt genomeTaxo.txt inputClust.txt inputmegan.txt" << endl;
	exit(-1);
}
int main(int argc, char* argv[])
{
	if(argc<4)
		usage();
	ifstream ifs(argv[1]);
	assert(!ifs.fail());
	int idx = 0;
	map<int,int>taxo2level;
	while(!ifs.eof())
	{
		int taxoid=0;
		for(int i=0;i<7;++i)//from kingdom to species.
		{
			ifs >> taxoid;
			taxo2level[taxoid] = i;
		}
		ifs.getline(Buf,MAXL);
		++idx;
	}
	ifs.close();
	cerr << idx << " genomes have been loaded." << endl;
	ifs.open(argv[2]);
	assert(!ifs.fail());
	int tmp;
	int gidx = 0;
	while(!ifs.eof())
	{
		ifs >> tmp;
		for(int i=0;i<7;++i)
			ifs >> genomeTaxo[gidx][i];
		ifs.getline(Buf,MAXL);
		++gidx;
	}
	ifs.close();
	cerr << gidx << " gnomes have been loaded. " << endl;
	//////////////////////////////////////////////////////
	const int RMAX = 100000000;
	int* RGid_clust = new int[RMAX];
	for(int i=0;i<RMAX;++i)
		RGid_clust[i] = -1;
	{
	int C2[4][8];
	int A2[4][8];
	for(int i=0;i<4;++i)
		for(int j=0;j<8;++j)
			A2[i][j] = C2[i][j] = 0;
		///////////////////////////////////////////
	ifs.open(argv[3]);
	assert(!ifs.fail());
	int tgid[7];
	double tscore[7];

	int cid = 0,csize=10000000;
	int curgid = -1,curlevel=7;
	while(!ifs.eof())
	{
		ifs.getline(Buf,MAXL);
		if(Buf[0]=='>')
		{
			if(cid >=0 && csize<100000)continue;
			int spid,rid,peid;
			sscanf(Buf,">read%d_%d/%d",&spid,&rid,&peid);
			--spid;
			rid += peid-1;
			RGid_clust[rid] = curgid;
			///////////////////////////////////////////////////
			int x = toOriId[spid]/25;
			++C2[x][curlevel];
			if(curgid == genomeTaxo[spid][curlevel])
				++A2[x][curlevel];
			///////////////////////////////////////////////////
		}
		else if(Buf[0]=='c')
		{
			sscanf(Buf,"cluster %d\t%d",&cid,&csize);
		}
		else if(Buf[0]=='T')
		{
			sscanf(Buf,"Taxo:\t%d:%lf\t%d:%lf\t%d:%lf\t%d:%lf\t%d:%lf\t%d:%lf\t%d:%lf",&tgid[0],&tscore[0],&tgid[1],&tscore[1],&tgid[2],&tscore[2],&tgid[3],&tscore[3],&tgid[4],&tscore[4],&tgid[5],&tscore[5],&tgid[6],&tscore[6]);
			////////////////////////////////////////////////////////
		/*	cerr << cid << '\t';
			for(int i=0;i<7;++i)
				cerr << tgid[i] << '\t';
			for(int i=0;i<7;++i)
				cerr << tscore[i] << '\t';
			cerr << endl;*/
			///////////////////////////////////////////////////////
			for(int i=6;i>=0;--i)
			{
				if(tscore[i] >= taxoThresh)
				{
					curgid = tgid[i];
					curlevel = i;
					break;
				}
				else  if(i==0)
				{
					curgid = -1;
					curlevel = 7;
				}
			}
		}
	}
	/////////////////////////////////////////////////
	for(int i=0;i<4;++i)
	{
		for(int j=0;j<8;++j)
			cout << A2[i][j] << '\t';
		cout << endl;
	}cout << endl << endl;
	for(int i=0;i<4;++i)
	{
		for(int j=0;j<8;++j)
			cout << C2[i][j] << '\t';
		cout << endl;
	}cout << endl << endl;
	ifs.close();
	}
	//////////////////////////////////////////////////////
	ifs.open(argv[4]);
	assert(!ifs.fail());
	int C[4][7];
	int A[4][7];
	int CL[4][7][8];
	int AL[4][7][8];
	int AL2[4][7][8];
	for(int i=0;i<4;++i)
		for(int j=0;j<7;++j)
			A[i][j] = C[i][j] = 0;
	for(int i=0;i<4;++i)
		for(int j=0;j<7;++j)
			for(int k=0;k<8;++k)
				AL[i][j][k] = CL[i][j][k] = 0;
	for(int i=0;i<4;++i)
		for(int j=0;j<7;++j)
			for(int k=0;k<8;++k)
				AL2[i][j][k] = 0;
	int error = 0;
	int fidx = 0;
	while(!ifs.eof())
	{
		ifs.getline(Buf,MAXL);
		int spid,rid,peid,ttaxo;
		sscanf(Buf,"read%d_%d/%d\t%d",&spid,&rid,&peid,&ttaxo);
		--spid;
		rid += peid-1;

	/*	char* str = Buf+4;
		int spid = atoi(str)-1;
		while((*str)!='\t')++str;
		++str;
		int ttaxo = atoi(str);*/
		int x = toOriId[spid]/25;
		map<int,int>::const_iterator itr = taxo2level.find(ttaxo);
		int y = 0;
		if(itr != taxo2level.end())
		{
			y = itr->second;

			int cttaxo = RGid_clust[rid];
			map<int,int>::const_iterator itr2 = taxo2level.find(cttaxo);
			if(itr2!=taxo2level.end())
				++CL[x][y][itr2->second];
			else
				++CL[x][y][7];

			++C[x][y];
			if(ttaxo == genomeTaxo[spid][y])
			{
				++A[x][y];
				if(cttaxo==genomeTaxo[spid][itr2->second])
				{
					if(itr2!=taxo2level.end())
						++AL[x][y][itr2->second];
					else
						++AL[x][y][7];
				}
			}
			if(cttaxo==genomeTaxo[spid][itr2->second])
			{
				if(itr2!=taxo2level.end())
					++AL2[x][y][itr2->second];
				else
					++AL2[x][y][7];
			}
		}
		else 
			++error;

		++fidx;
	}
	for(int i=0;i<4;++i)
	{
		for(int j=0;j<7;++j)
			cout << A[i][j] << '\t';
		cout << endl;
	}cout << endl << endl;
	for(int i=0;i<4;++i)
	{
		for(int j=0;j<7;++j)
			cout << C[i][j] << '\t';
		cout << endl;
	}
	cout << "error:\t" << error << endl;

	for(int i=0;i<4;++i)
	{
		for(int j=0;j<7;++j)
		{
			for(int k=0;k<8;++k)
				cout << CL[i][j][k] <<'\t';
			cout << endl;
		}
		cout << endl << endl;
	}
	cout << endl << endl;
	for(int i=0;i<4;++i)
	{
		for(int j=0;j<7;++j)
		{
			int tsum = 0;
			for(int k=0;k<8;++k)
				tsum += AL[i][j][k];
			for(int k=0;k<8;++k)
				cout << AL[i][j][k] <<'\t';
			cout << tsum << endl;
		}
		cout << endl << endl;
	}
	cout << endl << endl;
	for(int i=0;i<4;++i)
	{
		for(int j=0;j<7;++j)
		{
			int tsum = 0;
			for(int k=0;k<8;++k)
				tsum += AL2[i][j][k];
			for(int k=0;k<8;++k)
				cout << AL2[i][j][k] <<'\t';
			cout << tsum << endl;
		}
		cout << endl << endl;
	}

	ifs.close();
	return 0;
}
