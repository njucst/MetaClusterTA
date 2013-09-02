#include<iostream>
#include<cstdio>
#include<string>
#include<cstring>
#include<cstdlib>
#include "BWTs.h"
using namespace std;
int main(int argc, char* argv[])
{
	BWTs bwts(argv[1]);
	set<pair<int,int> > ans;
   	bwts.search("AA",ans);
	cerr << "#. ans: " << ans.size() << endl;
	for(set<pair<int,int> >::const_iterator itr = ans.begin();itr!=ans.end();++itr)
		cerr << (itr->first) << '\t' << (itr->second) << endl;
	return 0;
}
