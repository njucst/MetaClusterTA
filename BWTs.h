/*
 * last fixed: 2013.01.21.
 * by Wang Yi.
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
private:
	vector<BWTDtStr> bwts;
	BWTs(){}
};
#endif
