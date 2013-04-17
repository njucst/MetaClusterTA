#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>

using namespace std;
void usage()
{
	cerr << "getSubCtgs Ctgs_ClusterId.txt contig.fa targetId" << endl << endl;
	exit(-1);
}

int main(int argc, char* argv[])
{
	if(argc<4)
		usage();
	ifstream ifs(argv[1]);
	if(ifs.fail())
	{
		cerr << "File open failed: " << argv[1] << endl;
		exit(-1);
	}
	vector<int> ClusterId;
	int tmp;
	while(!ifs.eof())
	{
		tmp = -1;
		ifs >> tmp;
		ClusterId.push_back(tmp);
	}
	ifs.close();
	ifs.open(argv[2]);
	if(ifs.fail())
	{
		cerr << "File open failed: " << argv[1] << endl;
		exit(-1);
	}
	const int MAXL = 3000000;
	char* Buf = new char[MAXL];
	string* Ctg = new string[ClusterId.size()];
	string* Discrib = new string[ClusterId.size()];
	int idx = 0;
	while(!ifs.eof() && idx < ClusterId.size())
	{
		ifs.getline(Buf,MAXL);
		Discrib[idx]=Buf;
		ifs.getline(Buf,MAXL);
		Ctg[idx++] = Buf;
	}
	ifs.close();
	int tId = atoi(argv[3]);
	for(int i=0;i<ClusterId.size();++i)
	{
		if(ClusterId[i] == tId)
		{
			cout << Discrib[i] << endl;
			cout << Ctg[i] << endl;
		}
	}

	map<string,int>M;
	for(int i=0;i<ClusterId.size();++i)
		if(ClusterId[i] == tId)
		{
			const char* itr = Discrib[i].c_str();
			while(*itr != '/')++itr;
			++itr;int len = atoi(itr);
			while(*itr != '\t')++itr;
			++itr;string sid = itr;
			M[sid] += len;
		}
	long long sum = 0;
	for(map<string,int>::const_iterator itr=M.begin();itr!=M.end();++itr)
		sum += itr->second;
	for(map<string,int>::const_iterator itr=M.begin();itr!=M.end();++itr)
		cerr<< itr->first << '\t' << (itr->second) << '\t' << (itr->second)/(double)sum << endl;
	return 0;
}
