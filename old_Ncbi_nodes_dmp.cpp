#ifndef MCH_HYBRID_NCBINODESDMP_H_

#define MCH_HYBRID_NCBINODESDMP_H_
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cassert>
struct NcbiNodesDmpStruct
{
	int father;
	int levelId;
	NcbiNodesDmpStruct()
	{
		father = 0;
		levelId = 0;
	}
	NcbiNodesDmpStruct(int father_, int levelId_)
	{
		father = father_;
		levelId = levelId_;
	}
	void set(int father_, int levelId_)
	{
		father = father_;
		levelId = levelId_;
	}
};

using namespace std;
class NCBI_nodes_dmp
{
public:
	const static int MAXTAXID = 2000000;
	vector<NcbiNodesDmpStruct>TreeNodes;
	map<string,int> LevelName2Id;
	map<int,string> Id2LevelName;
	vector<string>LevelName;

	explicit NCBI_nodes_dmp(){}
	void init(string path)
	{
		ifstream ifs(path.c_str());
		assert(!ifs.fail());
		TreeNodes.resize(MAXTAXID);
		const int MAX = 1000;
		char buf[MAX];

		int Tidx=1;
	/*	string levels[]={"species","genus","family","order","class","phylum","kingdom","no"};
		for(int i=0;i<(sizeof(levels)/sizeof(string));++i)
		{
			Id2LevelName[Tidx] = levels[i];
			LevelName2Id[levels[i]] = Tidx++;
		}*/

		while(!ifs.eof())
		{
			int tson=0,tfather=0;
			ifs >> tson >> buf >> tfather >> buf >> buf;
			if(tson && LevelName2Id.find(buf) == LevelName2Id.end())
			{
				Id2LevelName[Tidx]=buf;
				LevelName2Id[buf] = Tidx++;
			}
			if(tson)
				TreeNodes[tson].set(tfather,LevelName2Id[buf]);
			ifs.getline(buf,MAX);
		}
		ifs.close();

		int Id_Of_NoRank = LevelName2Id[string("no")];
		for(vector<NcbiNodesDmpStruct>::iterator itr=TreeNodes.begin();itr!=TreeNodes.end();++itr)
		{
			int father = itr->father;
			int fatherlevel = TreeNodes[father].levelId;
			while(father && father!=1 && (fatherlevel==Id_Of_NoRank))
			{
			   	assert(father!=TreeNodes[father].father);
				father = TreeNodes[father].father;
				fatherlevel = TreeNodes[father].levelId;
			}
			itr->father = father;
		}
		//////////////////////////////
	/*	for(map<string,int>::const_iterator itr=LevelName2Id.begin();itr!=LevelName2Id.end();++itr)
			cerr << "Node LevelName: " << (itr->first) << '\t' << (itr->second) << endl;
		map<int,int>lcnt;
		for(vector<NcbiNodesDmpStruct>::const_iterator itr=TreeNodes.begin();itr!=TreeNodes.end();++itr)
			++lcnt[TreeNodes[itr->father].levelId];
		for(map<int,int>::const_iterator itr=lcnt.begin();itr!=lcnt.end();++itr)
			cerr << "level count: " << itr->first << '\t' << itr->second << endl;
		map<int,map<int,int> >flcnt;int t1=0;
		for(vector<NcbiNodesDmpStruct>::const_iterator itr=TreeNodes.begin();itr!=TreeNodes.end();++itr,++t1)
		{
			++flcnt[TreeNodes[itr->father].levelId][itr->levelId];
			if(TreeNodes[itr->father].levelId == Id_Of_NoRank && itr->father!=1)
				cerr << "exception " << t1 << '\t' << itr->father << endl;
		}
		for(map<int,map<int,int> >::const_iterator itr=flcnt.begin();itr!=flcnt.end();++itr)
		{
			cerr << "father is " << Id2LevelName[itr->first] << endl;
			for(map<int,int>::const_iterator itr2 = itr->second.begin();itr2!=itr->second.end();++itr2)
				cerr<< Id2LevelName[itr2->first] << '\t' << itr2->second << endl;
		}*/
		//////////////////////////////
	}
	virtual ~NCBI_nodes_dmp(){}
private:
	NCBI_nodes_dmp(const NCBI_nodes_dmp& ncbi_nodes_dmp){}
	const NCBI_nodes_dmp& operator=(const NCBI_nodes_dmp& ncbi_nodes_dmp){return*this;}
};
#endif

int main(int argc, char* argv[])
{
	NCBI_nodes_dmp tdmp;
	tdmp.init(argv[1]);
	return 0;
}
