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
int toOriId[] = {25,26,27,0,1,75,6,7,8,9,10,11,12,13,4,5,14,15,16,18,19,28,29,30,31,32,33,34,35,17,36,37,38,39,76,20,77,78,79,21,22,23,40,41,24,42,43,44,50,51,52,53,54,55,56,57,58,59,60,61,62,2,3,45,46,47,64,63,49,48,65,66,67,99,68,98,69,97,70,71,72,73,74,95,96,87,88,89,90,91,92,93,94,85,86,84,83,82,80,81};
struct Clust
{
	int clustId,spId,numMaxRead,numReads;
	int taxoId[7];
	double taxoProb[7];
}cluster[MAXN];
int genomeTaxo[MAXG][7];
void usage()
{
	cerr << "./anaTaxoThresh clustertaxo.txt genomeTaxo.txt " << endl;
	exit(-1);
}
int main(int argc, char* argv[])
{
	if(argc<3)
		usage();
	ifstream ifs(argv[1]);
	assert(!ifs.fail());
	int idx = 0;
	while(!ifs.eof())
	{
		ifs >> cluster[idx].clustId >> cluster[idx].spId >> cluster[idx].numMaxRead >> cluster[idx].numReads;
		for(int i=0;i<7;++i)
			ifs >> cluster[idx].taxoId[i];
		for(int i=0;i<7;++i)
			ifs >> cluster[idx].taxoProb[i];
		ifs.getline(Buf,MAXL);
		++idx;
	}
	ifs.close();
	cerr << idx << " clusters have been loaded." << endl;
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
	for(int i=0;i<idx;++i)
	{
		for(int j=0;j<7;++j)
		{
			if(cluster[i].taxoId[j]==genomeTaxo[cluster[i].spId][6-j])
			{
				cout << toOriId[cluster[i].spId] << '\t';
				cout << j << '\t' << cluster[i].taxoId[j] << '\t' << cluster[i].taxoProb[j] << '\t';
				if(j==0)cout << '0' << endl;
				else cout << cluster[i].taxoProb[j] << endl;
				break;
			}
			else if(j==6)
			{
				cout << toOriId[cluster[i].spId] << '\t';
				cout << 7 << '\t' << cluster[i].taxoId[j] << '\t' << cluster[i].taxoProb[j] << '\t';
				cout << cluster[i].taxoProb[5] << endl;
			}
		}
	}
	return 0;
}

