#include <iostream>
#include "USet.h"
unsigned USet::getReadNum(const unsigned x)
{
	return readnum[find(x)];
}

unsigned USet::getCtgLen(const unsigned x)
{
	return ctglen[find(x)];
}

unsigned USet::find(const unsigned x)
{
	//std::cerr<<"Finding "<<x<<std::endl;
	unsigned y,temp,root;
	root=x;
	while(parent[root]!=root)
		root=parent[root];
	y=x;
	while(parent[y]!=y)//compress path
	{
		temp=parent[y];
		parent[y]=root;
		y=temp;
	}
	return root;
}

void USet::Union(unsigned x,unsigned y)
{
	unsigned u,v;
	u=find(x);v=find(y);
	if(u==v)return;
	if(rank[u]<=rank[v])//Union operation
	{
		parent[u]=v;
		readnum[v] += readnum[u];
		ctglen[v] += ctglen[u];
		if(rank[u]==rank[v])
			rank[v]++;
	}
	else
	{
		parent[v]=u;
		readnum[u] += readnum[v];
		ctglen[u] += ctglen[v];
	}
}

/*void USet::incctglen( unsigned i)
{
	++ctglen[find(i)];
}
void USet::make_set(const unsigned &x)
{
	parent[x]=x;
	readnum[x]=1;
	rank[x]=0;
}*/
/*void USet::ReInitial()
{
		for(unsigned i=0;i<Size;++i)
			parent[i]=i;
		for(unsigned i=0;i<Size;++i)
			rank[i]=0;
		for(unsigned i=0;i<Size;++i)
			readnum[i]=1;
		for(unsigned i=0;i<Size;++i)
			ctglen[i]=0;
}*/
