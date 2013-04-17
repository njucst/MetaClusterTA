#include<iostream>
#include<cstdio>
#include<string>
#include<cstring>
#include<cstdlib>
#include "BWTs.h"
using namespace std;
int main(int argc, char* argv[])
{
	string S[] = {
		"ATTTTTAAGACTCTT",
		"AACATACCGGATCATATTTGAGTCTGTTATTTAGG",
		"GATCATATTTGAGT",
		"TGCATATATCGAGTCAAGATATAA"
	};
	cerr << "# tests:" << sizeof(S)/sizeof(string) << endl;
	BWTs bwts(argv[1]);
	system("ps ux");
	for(int i=0;i<(sizeof(S)/sizeof(string));++i)
	{
		set<int> ans;
		bwts.search(S[i], ans);
		cerr << ans.size() << endl;
		int idx=0;
		for(set<int>::const_iterator itr = ans.begin();itr!=ans.end()&&idx<50;++itr,++idx)
			cerr << *itr << '\t';
		cerr << endl << endl << endl;
	}
	return 0;
}
