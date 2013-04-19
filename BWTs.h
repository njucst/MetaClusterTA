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
	void search(const string& str,set<int>& ans)const;

	explicit BWTs(string filepath);
	virtual ~BWTs(){}
	/////////////////////////////////////////////
	const static int ForbidNum = 100;
	static int Forbid[ForbidNum];
	vector<short*>fids;
	/////////////////////////////////////////////
private:
	vector<BWTDtStr> bwts;
	BWTs(){}
};
#endif
