/*
 * last fixed: 2013.04.17.
 * by Wang Yi.
 * add fid for testing.
 * */
#ifndef MCH_HYBRID_BWTS_H_

#define MCH_HYBRID_BWTS_H_
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <set>
#include "BWTDtStr.h"
using namespace std;
class BWTs
{
public:
	void search(const string& str,set<pair<int,int> >& ans)const;

	explicit BWTs(string filepath);
	void clear();
	virtual ~BWTs();
private:
	vector<BWTDtStr> bwts;
	BWTs(){}
};
#endif
