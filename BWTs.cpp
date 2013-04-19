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
	ifstream ifs_fid((filepath+".fid").c_str());
	if(ifs_fid.fail())cerr<<"No fid file found."<<endl;
	else cerr << "fid file is found." << endl;
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

		/////////////////////////////////
		if(!ifs_fid.fail())
		{
			ifs_fid.read((char*)Taxid_of_SA,nuc_length*4);
			short* Fid_of_SA = new short[nuc_length];
			for(int i=0;i<nuc_length;++i)
				Fid_of_SA[i] = Taxid_of_SA[i];
			fids.push_back(Fid_of_SA);
		}
		//////////////////////////////////
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
	ifs_fid.close();
	///////////////////////////////////////
	for(unsigned i=0;i<bwts.size();++i)
	{
		cout << "bwt id: " << i <<'\t' << bwts[i].getlength()<<'\t' << bwts[i].Taxid_of_SA[0] << '\t' << bwts[i].Taxid_of_SA[1] << '\t'
			<< bwts[i].Taxid_of_SA[bwts[i].getlength()-2] << '\t'
			<< bwts[i].Taxid_of_SA[bwts[i].getlength()-1] << endl;
	}
	////////////////////////////////////////
	delete[] bwt_in_char;
	return;
}

void BWTs::search(const string& str,set<int>& ans)const
{
	int itridx = 0;
	for(vector<BWTDtStr>::const_iterator itr = bwts.begin();itr!=bwts.end();++itr,++itridx)
	{
		pair<long long,long long> interval = itr->search(str);
		if(interval.first >= 0)
			for(long long i=interval.first ; i<= interval.second;++i)
			{
				bool notbanned = true;
				if(fids.size()>0)
					for(int j=0;j<ForbidNum;++j)
						if(fids[itridx][i]==Forbid[j])
						{
							notbanned = false;
							break;
						}
				if(notbanned)
					ans.insert(itr->Taxid_of_SA[i]);
			}
	}
}
int BWTs::Forbid[ForbidNum] = {1573,1398,1997,75,1490,1666,585,659,37,1995,1121,1862,863,873,683,693,1615,582,1352,2056,1332,2012,1364,22,1504,773,232,193,1404,781,1466,753,962,2037,1652,2002,851,266,201,1672,280,1158,311,1494,1284,823,642,884,347,1731,838,1078,1520,1958,946,1468,665,1688,1689,1382,1945,12,1475,631,620,964,1101,1899,342,128,1091,1326,1748,24,932,1090,854,79,879,1474,1125,1142,531,858,1679,1669,2023,1274,26,1698,206,968,1623,1627,412,56,411,1874,1339,926};

////////////////////////////////////////////////////////////////////////////
