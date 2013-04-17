/*
 * last fixed: 2013.04.16.
 * by Wang Yi.
 * */
#ifndef MCH_HYBRID_BWTFIDGENERATOR_H_

#define MCH_HYBRID_BWTFIDGENERATOR_H_
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <string>
#include <vector>
#include "BWTPARA.h"
/////////////////////////////////////////////
////////////////////////////////////////////
//DC3
//http://www.cppblog.com/superKiki/archive/2010/05/15/115421.html
class BWT_FidGenerator
{
public:
	long long calBWT(char* input);
	void calBWTfromCompleGenome(std::string FoldPath,std::string gi_taxid_nucl_path, std::string outfile);

	unsigned* getSA() const {return SA;}
	char* getBWT() const {return bwt;}

	BWT_FidGenerator()
	{
		str = new unsigned[MAXN*3];
		SA = new unsigned[MAXN*3];
		wa = new unsigned[MAXN];
		wb = new unsigned[MAXN];
		ws = new unsigned[MAXN];
		wv = new unsigned[MAXN];

		bwt = new char[MAXN+2];
		inputstr = new char[MAXN+2];
	}
	~BWT_FidGenerator()
	{
		delete[]wv;
		delete[]ws;
		delete[]wb;
		delete[]wa;
		delete[]SA;
		delete[]str;
		
		delete[] bwt;
		delete[] inputstr;
	}
private:
	const static long long MAXM = 256; 

	const static char chA = 4;
	const static char chC = 5;
	const static char chG = 6;
	const static char chT = 7;
	const static char chES = 2; //extra end sigh attached to string.(e.g. '$');
	const static char chUDef = 3; //undefined char (e.g. 'N'); OR, char inserted between strings.

	unsigned* wa,*wb,*ws;
	unsigned* str;
	unsigned* wv;
	char* inputstr;
	unsigned* SA;

	char* bwt; 
	std::vector<std::pair<int,int> >BinarySearchTable;
	std::vector<std::pair<int,int> >FBinarySearchTable;

	bool c0(unsigned *r,long long a,long long b);
	bool c12(long long k,unsigned *r,long long a,long long b);
	void sort(unsigned *r,unsigned *a,unsigned *b,long long n,long long maxm);
	void dc3(unsigned *r,unsigned *sa,long long n,long long m);

	long long inputToStr(char* input);
	long long calSA(char* input);

	int binary_search_taxid(int sa_idx);
	int Fbinary_search_taxid(int sa_idx);
	void outputSA_BWT_Taxid(std::string outfile, long long bwtstridx);

};
#endif
