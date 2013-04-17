#ifndef __READER_H_

#define __READER_H_
#include <cstdio>
#include <iostream>
#include <string>
#include "KmerDistriPara.h"

class GenomeReader
{
public:
	GenomeReader();
	explicit GenomeReader(const std::string &filename);
	virtual ~GenomeReader();

	bool reRead(const std::string &filename);
	void clearString();
	void clear();

	std::size_t size(){return glength;}
	std::size_t length(){return glength;}
	std::string info(){return description;}
	int *getdigits(){if(digits==NULL)createDigits();return digits;}

	void reLocate(){std::fseek(fp, 0, SEEK_SET);}

	char* getLine(){return line;}
	bool createDigits();

	int getDistriNum(){return distriNum;}//number of fragments : glength/fragLen

	int** getAllDistri(const std::string &seed,const int fragLen);//get all fragments of a certain Length
	int* getSeedKDistr(const std::string &seed,int startposition,int endposition,bool revercom=true);
	double selfDistance(const std::string &seed,const int fragLen,double &average, double &stdev,bool newspear);//default value for newspear is false
	//////////////////////////////////////////////////////////////////////
	//added at 2012.11.30
	int* getKmerSpear(const int startposition, const int fragLen, const KmerDistriPara& Para);
	void getDiffAllKmerSpear(const KmerDistriPara& Para);
	double intraDis(const KmerDistriPara& Para, const int FragLen); 
	//////////////////////////////////////////////////////////////////////
	////////////////////////////////////
	//added at 2012.11.30
	const static int GN = 11;
	int** All5merSpear[GN];
	int SpearSize[GN];
	int* Total5merSpear;

protected:
	std::FILE* fp;
	char* line;
	std::size_t glength;
	std::string description;
	int* digits;
	int**allDistri;//get all fragments of a certain Length
	int distriNum;

private:
	//////////////////////////////////////////////////////////////////////
	//added at 2012.11.30
//	int* getKmerSpear(const int startposition, const int fragLen, const KmerDistriPara& Para);
	void getAllKmerSpear(const KmerDistriPara& Para, const int FragLen, const int Overlap, int** &OutDis,int &FragNum );
	//////////////////////////////////////////////////////////////////////
	static const int MAXLENGTH = 20000000;
	static const int MAXLINE = 20000000;
	GenomeReader(const GenomeReader &reader){}
	const GenomeReader &operator =(const GenomeReader &reader){return *this;}
};
#endif

