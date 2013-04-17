//DC3
//http://www.cppblog.com/superKiki/archive/2010/05/15/115421.html
#include <iostream>
#define F(x) ((x)/3+((x)%3==1?0:tb))
#define G(x) ((x)<tb?(x)*3+1:((x)-tb)*3+2)
using namespace std;
//typedef unsigned unsigned;
//typedef unsigned char unsigned;
const int maxn = 100000;
unsigned wa[maxn],wb[maxn],wv[maxn],ws2[maxn];
unsigned sa[maxn],rank[maxn];
/*int c0(int *r,int a,int b)
{
	return r[a]==r[b]&&r[a+1]==r[b+1]&&r[a+2]==r[b+2];
}
int c12(int k,int *r,int a,int b)
{
	if(k==2) return r[a]<r[b]||r[a]==r[b]&&c12(1,r,a+1,b+1);
	else return r[a]<r[b]||r[a]==r[b]&&wv[a+1]<wv[b+1];
}
void sort(int *r,int *a,int *b,int n,int m)
{
    int i;
    for(i=0;i<n;i++) wv[i]=r[a[i]];
    for(i=0;i<m;i++) ws2[i]=0;
    for(i=0;i<n;i++) ws2[wv[i]]++;
    for(i=1;i<m;i++) ws2[i]+=ws2[i-1];
    for(i=n-1;i>=0;i--) b[--ws2[wv[i]]]=a[i];
    return;
}
void dc3(int *r,int *sa,int n,int m)
{
    int i,j,*rn=r+n,*san=sa+n,ta=0,tb=(n+1)/3,tbc=0,p;
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
*/
/*
bool c0(char *r,long long a,long long b)
{
	return r[a]==r[b]&&r[a+1]==r[b+1]&&r[a+2]==r[b+2];
}
bool c12(long long k,char *r,long long a,long long b)
{
	if(k==2) return r[a]<r[b]||r[a]==r[b]&&c12(1,r,a+1,b+1);
	else return r[a]<r[b]||r[a]==r[b]&&wv[a+1]<wv[b+1];
}
void sort(char *r,unsigned *a,unsigned *b,long long n,long long maxm)
{
    long long i;
    for(i=0;i<n;i++) wv[i]=r[a[i]];
    for(i=0;i<maxm;i++) ws2[i]=0;
    for(i=0;i<n;i++) ws2[wv[i]]++;
    for(i=1;i<maxm;i++) ws2[i]+=ws2[i-1];
    for(i=n-1;i>=0;i--) b[--ws2[wv[i]]]=a[i];
    return;
}
void dc3(char *r,unsigned *sa,long long n,long long maxm)
{
    long long i,j,ta=0,tb=(n+1)/3,tbc=0,p;
	unsigned* san = sa+n;
	char *rn = r+n;
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
bool c0(unsigned *r,int a,int b)
{
	return r[a]==r[b]&&r[a+1]==r[b+1]&&r[a+2]==r[b+2];
}
bool c12(int k,unsigned *r,int a,int b)
{
	if(k==2) return r[a]<r[b]||r[a]==r[b]&&c12(1,r,a+1,b+1);
	else return r[a]<r[b]||r[a]==r[b]&&wv[a+1]<wv[b+1];
}
void sort(unsigned *r,unsigned *a,unsigned *b,int n,int m)
{
    int i;
    for(i=0;i<n;i++) wv[i]=r[a[i]];
    for(i=0;i<m;i++) ws2[i]=0;
    for(i=0;i<n;i++) ws2[wv[i]]++;
    for(i=1;i<m;i++) ws2[i]+=ws2[i-1];
    for(i=n-1;i>=0;i--) b[--ws2[wv[i]]]=a[i];
    return;
}
void dc3(unsigned *r,unsigned *sa,int n,int m)
{
    int i,j;unsigned *rn=r+n;unsigned *san=sa+n;int ta=0,tb=(n+1)/3,tbc=0,p;
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

int main()
{
	char str[] = "ATGGTTCCGATGCCATTGCCGGAGTCGTCAATATCATTTTGAAATCTGACAATCATGGTGGCAGTGCCAGAAGTCAGATCGGACAGACCTATGCCGGAGATGGATTGGTCGGACAGGCCGGTTTTAATAAGGGATTCAAAATTGGCCATTCCGGTTTCTTTGATGTCGCCTTCGATTTTCGTCATCAAAATCATACCAGCCGTGATGGTATCGATAGCCGAACCCAGCGCCATAGTTTGAAAGTTGTCGGTGATCCGATGGCGACCCGCTATAATCTAGCGATCAATGCCGGTTATGATTTTGGAAATGGCATCGAAATTTATACGACTGATAGCTATTCGCATCGCAATMATAATATAACATTATTTAATAATTTATATCCAATAAATAAATACTTTATACCCATAAATTATATTTAAATTTTATAATTAAATAAAATAATTATNNNAAAAATTATCGAGGAGAAGCCAAATTATCTAAAAATAAAGAGGGGGCGTATGAAAAATTTTATTAAAAANNNNNNNTTTCTTTTTTCCTTTTCATGGTTTTTCTTGCAATATCGATCTTGATTTGATTGCACCTGAATTGATTGACCATAMATCAGAGCTTTATCAGGATGAAGCCTTTGTCCAGCATCAATAAACGGCAATTTCTTTCGCTTTCTGCTCTGACCGCTTTTTTACCAATGGCAACCCAAGTCTTGGGTAAAAATTCCGGTAAAGAGGCCTCGAAGAATTTGGMACGCACAAGGGGTGACCGTAGCTTCGACGCTTGAGGCTCAAAATATATCGGCATGGTTCTTTACCAATGGTGCCAGCACCCGCACCCGCGGTTTGGATTTTACM";
//	int str2[] = {1,2,1,2,3,4,0};
	unsigned* str2 = new unsigned[sizeof(str)];
	for(int i=0;i<sizeof(str);++i)
		str2[i] = str[i];
	cerr << sizeof(str) << '\t' << sizeof(str2) << endl;
	dc3(str2,sa,sizeof(str),'Z');
//	dc3(str,sa,sizeof(str),300);

	for(int i=0;i<sizeof(str);++i)
		cout << sa[i] << ' ';
	cout << endl;
	return 0;
}
