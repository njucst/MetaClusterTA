#include <iostream>
#include <cassert>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <string>
#include <dirent.h>
#include <cstdlib>
#include "SmallArray.h"
using namespace std;
struct Node{
	int id,posi;
};
int main()
{
	cerr << log(1.0e-9) << endl;;
	const int N=250000;

	const int L = 2; 
	system("ps ux");
	SmallArray <int,L> SA2;
	int ** V=new int*[N];
	for(int i=0;i<N;++i)
		V[i] = SA2.getNew();

	SA2.test();
	for(int i=0;i<10;++i)
		for(int j=0;j<L;++j)
			cout << V[i][j] << endl;

	system("ps ux");
//	delete[] V[0];
	SA2.clear();
	delete[]V;
	system("ps ux");

/*	for(int i=0;i<N;++i)
		V[i] = (Node*)malloc(2*sizeof(Node));
	system("ps ux");
	for(int i=0;i<N;++i)
		free(V[i]);
	system("ps ux");
	delete[]V;
	system("ps ux");*/

	return 0;
}
