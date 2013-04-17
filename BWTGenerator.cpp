#include  "BWTGenerator.h"
#include <dirent.h>
#include <fstream>
#include <iostream>

#define F(x) ((x)/3+((x)%3==1?0:tb))
#define G(x) ((x)<tb?(x)*3+1:((x)-tb)*3+2)

using namespace std;
bool BWTGenerator::c0(unsigned *r,long long a,long long b)
{
	return r[a]==r[b]&&r[a+1]==r[b+1]&&r[a+2]==r[b+2];
}
bool BWTGenerator::c12(long long k,unsigned *r,long long a,long long b)
{
	if(k==2) return r[a]<r[b]||r[a]==r[b]&&c12(1,r,a+1,b+1);
	else return r[a]<r[b]||r[a]==r[b]&&wv[a+1]<wv[b+1];
}
void BWTGenerator::sort(unsigned *r,unsigned *a,unsigned *b,long long n,long long m)
{
    long long i;
    for(i=0;i<n;i++) wv[i]=r[a[i]];
    for(i=0;i<m;i++) ws[i]=0;
    for(i=0;i<n;i++) ws[wv[i]]++;
    for(i=1;i<m;i++) ws[i]+=ws[i-1];
    for(i=n-1;i>=0;i--) b[--ws[wv[i]]]=a[i];
    return;
}
void BWTGenerator::dc3(unsigned *r,unsigned *sa,long long n,long long m)
{
    long long i,j;unsigned *rn=r+n;unsigned *san=sa+n;long long ta=0,tb=(n+1)/3,tbc=0,p;
    r[n]=r[n+1]=0;
    for(i=0;i<n;i++) if(i%3!=0) wa[tbc++]=i;
    sort(r+2,wa,wb,tbc,m);
    sort(r+1,wb,wa,tbc,m);
    sort(r,wa,wb,tbc,m);
    for(p=1,rn[F(wb[0])]=0,i=1;i<tbc;i++)
        rn[F(wb[i])]=c0(r,wb[i-1],wb[i])?p-1:p++;
    if(p<tbc) 
		dc3(rn,san,tbc,p);
    else 
		for(i=0;i<tbc;i++) 
			san[rn[i]]=i;
    for(i=0;i<tbc;i++) 
		if(san[i]<tb) 
			wb[ta++]=san[i]*3;
    if(n%3==1) 
		wb[ta++]=n-1;
    sort(r,wb,wa,ta,m);
    for(i=0;i<tbc;i++) 
		wv[wb[i]=G(san[i])]=i;
    for(i=0,j=0,p=0;i<ta && j<tbc;p++)
        sa[p]=c12(wb[j]%3,r,wa[i],wb[j])?wa[i++]:wb[j++];
    for(;i<ta;p++) 
		sa[p]=wa[i++];
    for(;j<tbc;p++) 
		sa[p]=wb[j++];
    return;
}
/*
bool BWTGenerator::c0(unsigned *r,long long a,long long b)
{
	return r[a]==r[b]&&r[a+1]==r[b+1]&&r[a+2]==r[b+2];
}
bool BWTGenerator::c12(long long k,unsigned *r,long long a,long long b)
{
	if(k==2) return r[a]<r[b]||r[a]==r[b]&&c12(1,r,a+1,b+1);
	else return r[a]<r[b]||r[a]==r[b]&&wv[a+1]<wv[b+1];
}
void BWTGenerator::sort(unsigned *r,unsigned *a,unsigned *b,long long n,long long maxm)
{
    long long i;
    for(i=0;i<n;i++) wv[i]=r[a[i]];
    for(i=0;i<maxm;i++) ws[i]=0;
    for(i=0;i<n;i++) ws[wv[i]]++;
    for(i=1;i<maxm;i++) ws[i]+=ws[i-1];
    for(i=n-1;i>=0;i--) b[--ws[wv[i]]]=a[i];
    return;
}
void BWTGenerator::dc3(unsigned *r,unsigned *sa,long long n,long long maxm)
{
    long long i,j,ta=0,tb=(n+1)/3,tbc=0,p;
	unsigned* san = sa+n;
	unsigned *rn = r+n;
    r[n]=r[n+1]=0;
    for(i=0;i<n;i++) if(i%3!=0) wa[tbc++]=i;
    sort(r+2,wa,wb,tbc,maxm);
    sort(r+1,wb,wa,tbc,maxm);
    sort(r,wa,wb,tbc,maxm);
    for(p=1,rn[F(wb[0])]=0,i=1;i<tbc;i++)
        rn[F(wb[i])]=c0(r,wb[i-1],wb[i])?p-1:p++;
    if(p<tbc) 
		dc3(rn,san,tbc,p);
    else 
		for(i=0;i<tbc;i++) 
			san[rn[i]]=i;
    for(i=0;i<tbc;i++) 
		if(san[i]<tb) 
			wb[ta++]=san[i]*3;
    if(n%3==1) 
		wb[ta++]=n-1;
    sort(r,wb,wa,ta,maxm);
    for(i=0;i<tbc;i++) 
		wv[wb[i]=G(san[i])]=i;
    for(i=0,j=0,p=0;i<ta && j<tbc;p++)
        sa[p]=c12(wb[j]%3,r,wa[i],wb[j])?wa[i++]:wb[j++];
    for(;i<ta;p++) 
		sa[p]=wa[i++];
    for(;j<tbc;p++) 
		sa[p]=wb[j++];
    return;
}*/
long long BWTGenerator::inputToStr(char* input)
{
	long long Len = strlen(input);
	assert(Len < MAXN);//not allowed to be equal, because we will add an extra chES at end;
	for(long long i=0;i<Len;++i)
	{
		switch(input[i])
		{
			case 'A':str[i]=chA;break;
			case 'C':str[i]=chC;break;
			case 'G':str[i]=chG;break;
			case 'T':str[i]=chT;break;
			case 'a':str[i]=chA;break;
			case 'c':str[i]=chC;break;
			case 'g':str[i]=chG;break;
			case 't':str[i]=chT;break;
			default: str[i]=chUDef;
		}
	}
	str[Len]=chES;
	str[Len+1] = 0;
	return Len+1;
}
long long BWTGenerator::calSA(char* input)
{
	long long Len = inputToStr(input);
	dc3(str,SA,Len,MAXM);
	return Len;
}
long long BWTGenerator::calBWT(char* input)
{
	long long Len = calSA(input);
	for(long long i=0;i<Len;++i)
		if(SA[i])
			bwt[i] = input[SA[i]-1];
		else
			bwt[i] = '$';
	bwt[Len] = 0;
	return Len;
}

int BWTGenerator::binary_search_taxid(int sa_idx)
{
	long long min_ = 0, max_ = BinarySearchTable.size()-1;
	while(min_ < max_)
	{
		long long mid_ = (min_+max_)/2;
		if(BinarySearchTable[mid_].first <= sa_idx)
			min_ = mid_+1;
		else
			max_ = mid_;
	}
	if(BinarySearchTable[min_].first >= sa_idx)
		return BinarySearchTable[min_].second;
	return -2;
}

inline char to3bits(char c1)
{
	switch(c1)
	{
		case 'A': return 4;
		case 'C': return 5;
		case 'G': return 6;
		case 'T': return 7;
		case 'N': return 3;
		case '\0': return 0;
		default: return 2;
	}
	return 2;
}
inline char combine_chars(char c1,char c2)
{
	return ((to3bits(c1)<<4)|(to3bits(c2)));
}
void BWTGenerator::outputSA_BWT_Taxid(string outfile, long long bwtstridx )
{
	ofstream ofs(outfile.c_str(), ios_base::out|ios_base::app|ios_base::binary);
	assert(!ofs.fail());

	for(long long i=0;i<bwtstridx;i+=2)
		ofs << combine_chars(bwt[i],bwt[i+1]);
	ofs << endl;
	for(long long i=0;i<bwtstridx;++i)
		SA[i]=binary_search_taxid(SA[i]);
	ofs.write((char*)SA,bwtstridx*4); 
	ofs.close();
//	cerr << "length: " << bwtstridx << '\t' << (int)bwt[((bwtstridx-1)>>1)<<1] << '\t' << (int)combine_chars(bwt[0],bwt[1])<<'\t'<<(int)combine_chars(bwt[2],bwt[3])<< '\t'<<(int)combine_chars(bwt[((bwtstridx-3)>>1)<<1],bwt[(((bwtstridx-3)>>1)<<1)+1])<<'\t' <<(int) combine_chars(bwt[((bwtstridx-1)>>1)<<1],bwt[(((bwtstridx-1)>>1)<<1)+1])<< endl;
}
//	ofs << bwt << endl;
	//	ofs << binary_search_taxid(SA[i])<<' ';
	/*
	for(long long i=0;i<bwtstridx;++i)
		ofs << bwt[i];
	ofs << endl;
	for(long long i=0;i<bwtstridx;++i)
		ofs << binary_search_taxid(SA[i]) << ' ';
	ofs << endl;
	ofs << inputstr << endl;

	for(vector<pair<int,int> >::const_iterator itr=BinarySearchTable.begin();itr!=BinarySearchTable.end();++itr)
		ofs << itr->first <<','<<itr->second<<';';
	ofs << endl;
	for(long long i=0;i<bwtstridx;++i)
		ofs << SA[i] << ' ';
	ofs << endl;
	ofs.close();
	*/
/////////////////////////////////////////////////////////////////////
//input from complete genomes
void getFilePaths(std::string FoldPath, std::vector<std::string>& filepaths)
{
	if(FoldPath[FoldPath.length()-1]!='/')
		FoldPath += "/";
	struct dirent *ent = NULL;
	DIR *pDir;
	if((pDir = opendir(FoldPath.c_str())) != NULL)
	{
		while(NULL != (ent = readdir(pDir)))
		{
			if(ent->d_type == 8)					// d_type:8-file,4-fold;
			{
				std::string filepath = ent->d_name;
				if(filepath.substr(filepath.length()-4)==".fna")
					filepaths.push_back(FoldPath + ent->d_name);
				//get fna file
			}
			else if(ent->d_name[0] != '.')
				getFilePaths(FoldPath + ent->d_name, filepaths);
		}
		closedir(pDir);
	}
}
void BWTGenerator::calBWTfromCompleGenome(std::string FoldPath,std::string gi_taxid_nucl_path,string outfile)
{
	////get gi-taxid array
	const long long MAX_GI = 1<<30;
	ifstream gifs(gi_taxid_nucl_path.c_str());
	if(gifs.fail())
	{
		cerr << "gi_taxid_nucl file open failed: "<<gi_taxid_nucl_path << endl;
		exit(-1);
	}
	const int BufMax=500;
	char Buf[BufMax];
	int* Gi_Taxid = new int[MAX_GI];
	for(long long i = 0;i<MAX_GI;++i)
		Gi_Taxid[i] = -1;
	cerr << " start to load gi_taxid_nucl.dmp." << endl;
	{
		long long tgi,ttaxid,txtidx = 0;
		while(!gifs.eof())
		{
			tgi=-1,ttaxid = -1;
			gifs.getline(Buf,BufMax);
//			if(txtidx<10)cerr <<Buf<<endl;
			sscanf(Buf,"%lld\t%lld",&tgi,&ttaxid);
			if(tgi>=MAX_GI)cerr<<"Array too small. current gid: "<<tgi<<"\t"<<txtidx<<"\t"<<Buf<<endl;
			else if(tgi>=0)
				Gi_Taxid[tgi] = ttaxid;
//			if(txtidx%10000000==0)cerr<<"line tag: "<<txtidx<<"\t"<<Buf<<endl;
			++txtidx;
		}
		gifs.close();
		cerr << " gi_taxid_nucl.dmp is loaded. # of lines:"<<txtidx<<"; max gi:"<<tgi<<"; max taxid:"<<ttaxid << endl;
	}
	//////////////////////////////////////////////////////////////
	////get string
	outfile += ".idx";
	system((string("rm ")+outfile).c_str());
	std::vector<std::string> filepaths;
	getFilePaths(FoldPath,filepaths);
	if(INTEST)
	{
		cerr << "There are " << filepaths.size() <<" files:"<< endl;
		for(vector<string>::const_iterator itr=filepaths.begin();itr!=filepaths.end();++itr)
			cerr << *itr << endl;
		cerr << endl;
	}
	BinarySearchTable.clear();
	char* genome = new char[40000000];
	long long bwtstridx = 0;
	int fileidx = 0;
	for(vector<string>::const_iterator itr = filepaths.begin();itr!=filepaths.end();++itr)
	{
		cerr << "Processing file " << fileidx++ << " : " << (*itr) << endl;
		ifstream ifs(itr->c_str());
		if(ifs.fail())
		{
			std::cerr << "File open failed: " << (*itr) << endl;
			continue;
		}
		ifs.getline(Buf,BufMax);
		int current_gi;
		sscanf(Buf,">gi|%d",&current_gi);
		if(Gi_Taxid[current_gi]<0)
		{
			cerr<<"gi<0: "<<fileidx << " " << current_gi << '\t'<<Gi_Taxid[current_gi]<<endl;
			ifs.close();
			continue;
		}
		unsigned charidx = 0;
		while(!ifs.eof())
		{
			ifs.getline(Buf,BufMax);
			strcpy(genome+charidx,Buf);
			charidx += strlen(Buf);
		}
		if(bwtstridx + charidx >= MAXN-1)
		{
			long long bwtlen=calBWT(inputstr);
			outputSA_BWT_Taxid(outfile, bwtlen);
			///////
			bwtstridx = 0;
			BinarySearchTable.clear();
		}
		else
		{
			genome[charidx++]='M';
			genome[charidx]=0;
		}
		strcpy(inputstr+bwtstridx,genome);
		bwtstridx += charidx;
		BinarySearchTable.push_back(make_pair(bwtstridx,Gi_Taxid[current_gi]));
		ifs.close();
	}
	if(bwtstridx > 0)
	{
		long long bwtlen = calBWT(inputstr);
		outputSA_BWT_Taxid(outfile, bwtlen);
	}
	delete[] genome;
	delete[] Gi_Taxid;
}
//
//////////////////////////////////////////////////////////////////////////
