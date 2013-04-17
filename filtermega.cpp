#include <iostream>
#include <cstdio>
#include <cassert>
#include <cstring>
#include <fstream>
#include <cstdlib>
using namespace std;
const int MAX = 1000;
char Buf[MAX];
int main(int argc,char* argv[])
{
	assert(argc >= 2);
	ifstream ifs(argv[1]);
	assert(!ifs.fail());
	while(!ifs.eof())
	{
		ifs.getline(Buf,MAX);
		int a,b,c,d;
		sscanf(Buf,"read%d_%d/%d\t%d",&a,&b,&c,&d);
		if(d>=0)
			cout << Buf << endl;
	}
}
