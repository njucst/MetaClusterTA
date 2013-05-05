#include <iostream>
#include <cassert>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <string>
#include <dirent.h>
#include <cstdlib>
using namespace std;
struct Node{
	int id,posi;
};
int main()
{
	const int N=250000000;
	Node ** V=new Node*[N];
/*	for(int i=0;i<N;++i)
		V[i] = (Node*)malloc(2*sizeof(Node));
	system("ps ux");
	for(int i=0;i<N;++i)
		free(V[i]);
	system("ps ux");
	delete[]V;
	system("ps ux");*/

	Node p[1000][2];
	for(int i=0;i<1000;++i)
		V[i] = p[i];
	cerr << sizeof(p) << '\t' << sizeof(V)<<endl;
	return 0;
}
