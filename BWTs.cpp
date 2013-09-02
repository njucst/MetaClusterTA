#include "BWTs.h"
#include <iostream>
#include <fstream>
using namespace std;
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
void to_format(char* source,unsigned length, char*target)
{
	for(unsigned i=0;i<length;++i)
	{
		switch(source[i])
		{
			case 'A': source[i] = 4;break;
			case 'C': source[i] = 5;break;
			case 'G': source[i] = 6;break;
			case 'T': source[i] = 7;break;
			case 'N': source[i] = 3;break;
			default: source[i] = 2;
		}
	}
	source[length]=0;
	for(unsigned i=0;i<length;i+=2)
		target[i/2] = (source[i]<<4)|(source[i+1]);
	target[(length+1)/2]=0;
}
BWTs::BWTs(string filepath)
{
	ifstream ifs(filepath.c_str());
	if(ifs.fail())
	{
		cerr<< "idx file open failed: "<<filepath<<endl;
		return;
	}
	static const long long MAXLEN = MAXN+3;
	char* bwt_in_char = new char[MAXLEN];
	int* NTaxon = new int[1];
	char* Id2Taxon= new char[MAXLEN/10];
	int* Taxid_of_SA = new int[MAXLEN];
	while(!ifs.eof())
	{
		ifs.getline(bwt_in_char,MAXLEN);
		if(ifs.eof())
			break;
		long long length = strlen(bwt_in_char);

		long long nuc_length = length*2;
		if((bwt_in_char[length-1]&0xf)==0)
			--nuc_length;
		//////////////////////////////////
		ifs.read((char*)Taxid_of_SA,nuc_length*4);
		ifs.read((char*)NTaxon,4);
		ifs.read(Id2Taxon,NTaxon[0]*4);
		BWTDtStr bwtdtstr;
		bwtdtstr.preprocess(bwt_in_char, Taxid_of_SA, nuc_length, (int*)Id2Taxon, NTaxon[0]);

		bwts.push_back(bwtdtstr);
	}
	ifs.close();
	///////////////////////////////////////
/*	for(unsigned i=0;i<bwts.size();++i)
	{
		cout << "bwt id: " << i <<'\t' << bwts[i].getlength()<<'\t' << bwts[i].Taxid_of_SA[0] << '\t' << bwts[i].Taxid_of_SA[1] << '\t'
			<< bwts[i].Taxid_of_SA[bwts[i].getlength()-2] << '\t'
			<< bwts[i].Taxid_of_SA[bwts[i].getlength()-1] << endl;
	}*/
	////////////////////////////////////////
	delete[] Taxid_of_SA;
	delete[] bwt_in_char;
	delete[] NTaxon;
	delete[] Id2Taxon;
	return;
}

void BWTs::search(const string& str,set<pair<int,int> >& ans)const
{
	int NSeq = 0;
	for(vector<BWTDtStr>::const_iterator itr = bwts.begin();itr!=bwts.end();++itr)
	{
		pair<long long,long long> interval = itr->search(str);
		if(interval.first >= 0)
			for(long long i=interval.first ; i<= interval.second;++i){
				pair<int,int> t = itr->getSidTaxidOfPosition(i);
				t.first += NSeq;
				ans.insert(t);
			}
		NSeq += itr->getNSeq();
	}
}
void BWTs::clear()
{
	for(vector<BWTDtStr>::iterator itr=bwts.begin();itr!=bwts.end();++itr)
		itr->clear();
	bwts.clear();
}
BWTs::~BWTs()
{
	for(vector<BWTDtStr>::iterator itr=bwts.begin();itr!=bwts.end();++itr)
		itr->clear();
}
