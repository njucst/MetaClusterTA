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
	static const long long MAXLEN = MAXN+3;
	char* bwt_in_char = new char[MAXLEN];
	ifstream ifs(filepath.c_str());
	if(ifs.fail())
	{
		cerr<< "idx file open failed: "<<filepath<<endl;
		return;
	}
	while(!ifs.eof())
	{
		int* Taxid_of_SA = new int[MAXLEN];
		ifs.getline(bwt_in_char,MAXLEN);
		if(ifs.eof())
			break;
		long long length = strlen(bwt_in_char);

		long long nuc_length = length*2;
		if((bwt_in_char[length-1]&0xf)==0)
			--nuc_length;
		ifs.read((char*)Taxid_of_SA,nuc_length*4);
		BWTDtStr bwtdtstr;
		bwtdtstr.preprocess(bwt_in_char, Taxid_of_SA);
/*		char* bwt_in_DESIGN_format = new char[length+2];
		strcpy(bwt_in_DESIGN_format,bwt_in_char);
		int nuc_length = length*2;
		if((bwt_in_DESIGN_format[length-1]&0xf)==0)
			--nuc_length;
		ifs.read((char*)Taxid_of_SA,nuc_length*4);
		BWTDtStr bwtdtstr;
		bwtdtstr.preprocess(bwt_in_DESIGN_format, Taxid_of_SA);
		*/
		bwts.push_back(bwtdtstr);

/*		cout << length << '\t' << nuc_length << endl;
		for(int i=0;i<length;++i)
			cout << (int)bwt_in_char[i]<<' ';
		cout << endl;
		for(int i=0;i<length;++i)
			cout<<Taxid_of_SA[i]<<' ';
		cout << endl;*/

//		ifs.getline(bwt_in_char,MAXLEN);
	}
	ifs.close();
	///////////////////////////////////////
	for(unsigned i=0;i<bwts.size();++i)
	{
		cout << "test: " << i <<'\t' << bwts[i].getlength()<<'\t' << bwts[i].Taxid_of_SA[0] << '\t' << bwts[i].Taxid_of_SA[1] << '\t'
			<< bwts[i].Taxid_of_SA[bwts[i].getlength()-2] << '\t'
			<< bwts[i].Taxid_of_SA[bwts[i].getlength()-1] << endl;
	}
	////////////////////////////////////////
	delete[] bwt_in_char;
	return;
}

void BWTs::search(const string& str,set<int>& ans)const
{
	for(vector<BWTDtStr>::const_iterator itr = bwts.begin();itr!=bwts.end();++itr)
	{
		pair<long long,long long> interval = itr->search(str);
		if(interval.first >= 0)
			for(long long i=interval.first ; i<= interval.second;++i)
				ans.insert(itr->Taxid_of_SA[i]);
			//	ans.push_back(itr->Taxid_of_SA[i]);
	}
//	set<int>sans(ans.begin(),ans.end());
//	vector<int>vans(sans.begin(),sans.end());
//	return ans;
}

////////////////////////////////////////////////////////////////////////////
