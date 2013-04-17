#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
using namespace std;

const int MAX = 10000000;
char Buf[MAX];

void Usage()
{
	cerr << "Turn a contig file from assembly to a simulate-test contig file. " << endl;
	cerr << "Usage:\t ToTestCtgFa input_contig_fa psl_file out_name " << endl;
	cerr << "e.g. ToTestCtgFa contig-100.fa out_100.psl out_100_test.fa " << endl;
	cerr << endl;
}


int main(int argc, char* argv[])
{
	if(argc < 4)
	{
		Usage();
		exit(-1);
	}
	ifstream ifs_ctg(argv[1]);
	if(ifs_ctg.fail())
	{
		cerr << "Open failed: " << argv[1] << endl;
		exit(-1);
	}
	ifstream ifs_psl(argv[2]);
	if(ifs_psl.fail())
	{
		cerr << "Open failed: " << argv[2] << endl;
		exit(-1);
	}
//////////////////////////////////////////////////////
	map<string,pair<int,string> > M;
	int pslcount = 0;
	char tbuf[1000];
	while(!ifs_psl.eof())
	{
		ifs_psl.getline(Buf,MAX);
		if(pslcount++ < 5)
			continue;

		string tstr;int t;
		int match,mismatch,ctgLen,refLen;
		char ctgname[1000];
		char refname[1000];
		sscanf(Buf,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%s",
				&match,&mismatch,&t,&t,&t,&t,&t,&t,tbuf,ctgname,&ctgLen,&t,&t,refname);
		string ctgName = ctgname;
		if(M.find(ctgName)==M.end() || M[ctgName].first < match)
		{
			M[ctgName] = make_pair(match,string(refname));
		}
	}
	ifs_psl.close();

	ofstream ofs(argv[3]);
	if(ofs.fail())
	{
		cerr << "Open failed: " << argv[3] << endl;
		exit(-1);
	}
	while(!ifs_ctg.eof())
	{
		ifs_ctg.getline(Buf,MAX);
		if(Buf[0]!='>')
		{
			cerr << "hit non- '>'. Over." << endl;
			break;
		}
		ofs << Buf;
		map<string,pair<int,string> >::iterator itr = M.find(string(Buf+1));

		ifs_ctg.getline(Buf,MAX);
		if(itr != M.end())
			sprintf(tbuf,"\t%d/%d\t%s\0",(itr->second).first,strlen(Buf),(itr->second).second.c_str());
		else
			sprintf(tbuf,"\t0/%d\tunknown\0",strlen(Buf));
		ofs << tbuf << endl;
		ofs << Buf << endl;
	}
	ifs_ctg.close();
	ofs.close();
}
