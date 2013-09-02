#ifndef MCH_HYBRID_NCBINODESDMP_H_

#define MCH_HYBRID_NCBINODESDMP_H_
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <cassert>
#include <cstdlib>
using namespace std;

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

class NCBI_nodes_dmp
{
public:
	const static int MAXTAXID = 2000000;
	vector<NcbiNodesDmpStruct>TreeNodes;
	map<string,int> LevelName2Id;
	map<int,string> Id2LevelName;
	int MaxTaxId;

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
		MaxTaxId = Tidx;
		////////////////////////////////////////////////////
		//remove "no rank"
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
		//topological sort
		map<int,int>Topology;
		set<pair<int,int> >Edge;
		for(vector<NcbiNodesDmpStruct>::const_iterator itr=TreeNodes.begin();itr!=TreeNodes.end();++itr)
			if(itr->levelId != Id_Of_NoRank && itr->levelId != TreeNodes[itr->father].levelId)
				Edge.insert(make_pair(itr->levelId,TreeNodes[itr->father].levelId));
		topologicalSort(Edge, Topology);
		for(vector<NcbiNodesDmpStruct>::iterator itr=TreeNodes.begin();itr!=TreeNodes.end();++itr)
			itr->levelId = Topology[itr->levelId];
		Id2LevelName.clear();
		for(map<string,int>::iterator itr=LevelName2Id.begin();itr!=LevelName2Id.end();++itr)
		{
			itr->second = Topology[itr->second];
			Id2LevelName[itr->second] = itr->first;
		}
		Id_Of_NoRank = LevelName2Id[string("no")];
		for(vector<NcbiNodesDmpStruct>::const_iterator itr=TreeNodes.begin();itr!=TreeNodes.end();++itr)
			assert(!(TreeNodes[itr->father].levelId == Id_Of_NoRank && itr->father!=1));
		////////////////////////////////////////////////////
		{
			for(map<int,string>::const_iterator itr=Id2LevelName.begin();itr!=Id2LevelName.end();++itr)
				cerr << "Taxon name: " << itr->first << '\t' << itr->second << endl;
	/*	map<int,int>rTopology;
		for(map<int,int>::const_iterator itr=Topology.begin();itr!=Topology.end();++itr)
			rTopology[itr->second] = itr->first;
		for(map<int,int>::const_iterator itr=rTopology.begin();itr!=rTopology.end();++itr)
			cerr << "Topology: " << Id2LevelName[itr->second] << '\t' << itr->first << endl;
		*/
		}
		////////////////////////////////////////////////////
	}

	void topologicalSort(const set<pair<int,int> >& Edge, map<int,int>& result)
	{
		result.clear();
		map<int,int>degree;
		for(set<pair<int,int> >::const_iterator itr=Edge.begin();itr!=Edge.end();++itr)
		{
			degree[itr->first];
			++degree[itr->second];
		}
		int ridx=1;
		while(!degree.empty())
		{
			int id0=-1;
			for(map<int,int>::const_iterator itr=degree.begin();itr!=degree.end();++itr)
			{
				if(itr->second == 0)
				{
					id0 = itr->first;
					break;
				}
			}
			assert(id0 >= 0);
			result[id0]=ridx++;
			degree.erase(id0);
			for(set<pair<int,int> >::const_iterator itr=Edge.begin();itr!=Edge.end();++itr)
				if(itr->first==id0)
					--degree[itr->second];
		}
	}

	template<class T>
	void addScore2Path(int leaf,T score,map<int,map<int,T> >&Node_Score)const
	{
//		assert(leaf>=0 && leaf < MAXTAXID);
		if(leaf<0 || leaf>=MAXTAXID)
			return;
		while(leaf != TreeNodes[leaf].father)
		{
			Node_Score[TreeNodes[leaf].levelId][leaf] += score;
			leaf = TreeNodes[leaf].father;
		}
		Node_Score[TreeNodes[leaf].levelId][leaf] += score;
	}

	vector<int> getTaxo(int leaf)const
	{
		int vsize = Id2LevelName.size();
		vector<int> ans(vsize,0);
		if(leaf<0 || leaf>=MAXTAXID)
			return ans;
		while(leaf != TreeNodes[leaf].father)
		{
			ans[TreeNodes[leaf].levelId] = leaf;
			leaf = TreeNodes[leaf].father;
		}
		return ans;
	}

	int getMaxTaxId()const
	{ 
		return MaxTaxId;
	}
	int getLevelId(int taxid)const
	{
		return TreeNodes[taxid].levelId;
	}
	virtual ~NCBI_nodes_dmp(){}
private:
	NCBI_nodes_dmp(const NCBI_nodes_dmp& ncbi_nodes_dmp){}
	const NCBI_nodes_dmp& operator=(const NCBI_nodes_dmp& ncbi_nodes_dmp){return*this;}
};
#endif
