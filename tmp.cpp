#ifndef MCH_HYBRID_NCBINODESDMP_H_

#define MCH_HYBRID_NCBINODESDMP_H_
//This line is to test git.
//This line is to test git2.
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
		////////////////////////////////////////////////////
		for(vector<NcbiNodesDmpStruct>::const_iterator itr=TreeNodes.begin();itr!=TreeNodes.end();++itr)
			assert(!(TreeNodes[itr->father].levelId == Id_Of_NoRank && itr->father!=1));
		////////////////////////////////////////////////////
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
