#include <iostream>
#include <fstream>
#include <cstdio>
#include <cassert>
#include <cstdlib>
#include <set>
#include <map>
using namespace std;
void usage()
{
	cerr << "./anaH_Anotation clusterComp.txt genomeTaxo.txt Y_X_taxo.txt" << endl;
	exit(-1);
}
const int MAXG = 200;
int genomeTaxo[MAXG][7];
const int MAXL = 20000;
char Buf[MAXL];
int Comp[MAXG][MAXG];

int main(int argc, char* argv[])
{
	if(argc<4)usage();
	ifstream ifs(argv[1]);
	assert(!ifs.fail());
	int cidx = 0;
	while(!ifs.eof())
	{
		for(int i=0;i<100;++i)
			ifs >> Comp[cidx][i];
		++cidx;
	}
	cerr<< "cidx:\t" << cidx << endl;
	ifs.close();

	ifs.open(argv[2]);
	assert(!ifs.fail());
	int tmp,gidx=0;
	while(!ifs.eof())
	{
		ifs >> tmp;
		for(int i=0;i<7;++i)
			ifs >> genomeTaxo[gidx][6-i];
		ifs.getline(Buf,MAXL);
		++gidx;
	}
	cerr<< "gidx:\t" << gidx << endl;
	ifs.close();

	ifs.open(argv[3]);
	assert(!ifs.fail());
	int C[4][8];
	int A[4][8];
	for(int i=0;i<4;++i)
		for(int j=0;j<7;++j)
			C[i][j] = A[i][j] = 0;
	int idx = 0;
	while(!ifs.eof())
	{
		int y,x,taxo;
		ifs >> y >> x >> taxo;
		if(y>6)y=6;
		for(int i=0;i<100;++i)
		{
			if(genomeTaxo[i][y]==taxo)
				C[x][y] += Comp[idx][i];
			A[x][y] += Comp[idx][i];
		}
		++idx;
	}
	cerr << idx << " lines are loaded." << endl;
	
	for(int i=0;i<4;++i)
	{
		for(int j=0;j<7;++j)
			cout << C[i][j] << '\t';
		cout << endl;
	}
	cout << endl << endl;
	for(int i=0;i<4;++i)
	{
		for(int j=0;j<7;++j)
			cout << A[i][j] << '\t';
		cout << endl;
	}
	for(int i=0;i<4;i+=3)
	{
		int sum1 = 0 ,sum2=0;
		for(int j=0;j<7;++j)
		{
			sum1 += C[i][j];
			sum2 += A[i][j];
		}
		cerr << sum1 << '\t' << sum2 << '\t' << sum1/(double)sum2 << endl << endl;
	}
	for(int i=0;i<4;i+=3)
	{
		for(int j=0;j<7;++j)
			cerr << C[i][j]<<'/'<<A[i][j]<<'\t';
		cerr << endl;
	}
	return 0;
}
