#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <cstring>
#include <cassert>
#include <cstdio>
#include <cstdlib>

using namespace std;
const int MAXL = 10000000;
char* Buf;
int* gid2treeid;

const int MAXGID = 500000000;

void usage()
{
	cerr<<"getDBStatictic "<<endl;
	exit(-1);
}
void ini_gid2treeid()
{
	gid2treeid = new int[MAXGID];
	for(int i=0;i<MAXGID;++i)
		gid2treeid[i] = -1;

	cerr<<"start to load gi_taxid." << endl;
	char path1[] = "/home/ywang/DB/gi_taxid/gi_taxid_nucl.dmp";
	char path2[] = "/home/ywang/DB/gi_taxid/gi_taxid_prot.dmp";
	ifstream ifs(path1); 
	assert(!ifs.fail());
	int idx = 0;
	while(!ifs.eof())
	{
		ifs.getline(Buf,MAXL);
		int tgid,ttreeid;
		sscanf(Buf,"%d\t%d",&tgid,&ttreeid);
		if(!(tgid>0 && tgid < MAXGID && ttreeid>=0 && tgid<MAXGID))
		{
			cerr << tgid << '\t' << idx << endl;;
		}
		gid2treeid[tgid] = ttreeid;
	//	assert(tgid>=0 && tgid < MAXGID && ttreeid>0 && tgid<MAXGID);
	}
	ifs.close();
	cerr << "gi_taxid_nucl.dmp is loaded." << endl;
	ifs.open(path2);
	assert(!ifs.fail());
	while(!ifs.eof())
	{
		ifs.getline(Buf,MAXL);
		int tgid,ttreeid;
		sscanf(Buf,"%d\t%d",&tgid,&ttreeid);
		assert(tgid>0 && tgid < MAXGID && ttreeid>=0 && tgid<MAXGID);
		gid2treeid[tgid] = ttreeid;
	}
	ifs.close();
	cerr << "gi_taxid_prot.dmp is loaded." << endl;
}

void ini_taxoTree()
{
	ifstream ifs("/home/ywang/DB/gi_taxid/nodes.dmp");
	assert(!ifs.fail());
	while(!ifs.eof())
	{
		ifs.getline(Buf,MAXL);
		int taxoid,parent;
		sscanf(Buf,"%d\t|\t%d",&taxoid,&parent);
		char* str = Buf;int tabcnt=0;
		for(;*str && tabcnt<4;++str)
			if(*str == '\t')++tabcnt;
	}
}

long long round10(long long len)
{
	long long step = 1;
	long long inp = len;
	while(inp > 10)
	{
		inp  /= 10;
		step *= 10;
	}
	return len/step*step;
}

void getStatis_nt(const string path, map<long long,long long>ans)
{
	ans.clear();
	ifstream ifs(path.c_str());
	assert(!ifs.fail());

//	map<long long, long long>SpLen;
	long long readlen=-1;
	map<long long,pair<long long,long long> >gidnum;

	long long idx = 0;
	long long lastspnum = 0;
	while(!ifs.eof())
	{
		ifs.getline(Buf,MAXL);
		if(Buf[0]=='>')
		{
			if(readlen > 0)
				gidnum[lastspnum].second += readlen;
			readlen = 0;
			char* strend = Buf+strlen(Buf);
			long long tgid = 0;
			set<long long>ttreeids;
			for(char* str=Buf+3;*str;++str)
			{
				if(*str == '|' && *(str-1)=='i' && *(str-2)=='g')
				{
					tgid = atoi(++str);
					ttreeids.insert(gid2treeid[tgid]);
				}
			}
			++gidnum[ttreeids.size()].first;
			lastspnum = ttreeids.size();

			++idx;
		}
		else readlen += strlen(Buf);
	}
	cerr << idx << " reads are processed." << endl;
	for(map<long long,pair<long long,long long> >::const_iterator itr=gidnum.begin();itr!=gidnum.end();++itr)
		cerr << itr->first << " -> " << itr->second.first << " -> " << itr->second.second << endl;
	ifs.close();
}

int main(int argc, char* argv[])
{
	if(argc < 1)
		usage();
	Buf = new char[MAXL];

	ini_gid2treeid();
	string str1 = "/home/ywang/DataBase/nt/nt";
	string str2 = "/home/ywang/DataBase/nr/nr";
	map<long long, long long>ans;
	getStatis_nt(str1.c_str(), ans);
	getStatis_nt(str2.c_str(), ans);
	return 0;
}
