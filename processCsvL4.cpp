#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>
using namespace std;
struct clustClass
{
	int gid[7];//species->kingdom
//	double taxoScore[7];
	double maxScore[7];
	double allScore[7];
	int ReadN,maxRN;
	int WrongMajorSp;
	int RealMajorSp;

	int correct_taxo_level;
};
struct genomeClass
{
	int id;
	int taxoId[7];//species->kingdom
};
void usage()
{
	cerr << "processCsv tobeProcessed.txt " << endl;
	exit(-1);
}
const int GenoN = 100;
const int CluN = 37;
const int MAX = 100000;
const double GetTaxoThresh = 0.7;
char Buf[MAX];
/////////////////////////////
clustClass cluster[CluN];
genomeClass genome[GenoN];
int cluMtx[CluN][GenoN];

int main(int argc, char* argv[])
{
	if(argc < 2)
		usage();
	ifstream ifs(argv[1]);
	assert(!ifs.fail());
	////////////////////////////////////
	/*{
		ifs.getline(Buf,MAX);
		int id1,id2;
		for(int i=0;i<GenoN;++i)
		{
			ifs.getline(Buf,MAX);
			sscanf(Buf,"%d\t%d",&id1,&id2);
			i = id2;
		}
	}*/
	////////////////////////////////////
	{
		ifs.getline(Buf,MAX);
		for(int i=0;i<GenoN;++i)
		{
			ifs.getline(Buf,MAX);
			sscanf(Buf,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",&genome[i].id,&genome[i].taxoId[6],&genome[i].taxoId[5],&genome[i].taxoId[4],&genome[i].taxoId[3],&genome[i].taxoId[2],&genome[i].taxoId[1],&genome[i].taxoId[0]);
		}
	}
	////////////////////////////////////
	{
		for(int level=6;level>=0;--level)
		{
			ifs.getline(Buf,MAX);
			int tmp;double dtmp;
			for(int i=0;i<CluN;++i)
			{
				ifs.getline(Buf,MAX);
				sscanf(Buf,"%d\t%d\t%lf\t%lf",&tmp,&cluster[i].gid[level],&cluster[i].maxScore[level],&cluster[i].allScore[level]);
			}
		}
	}
	////////////////////////////////////
	{
		ifs.getline(Buf,MAX);
		for(int i=0;i<CluN;++i)
		{
			for(int j=0;j<GenoN;++j)
				ifs >> cluMtx[i][j];
			ifs.getline(Buf,MAX);
		}
	}
	////////////////////////////////////
	int refLevel_assignLevel_All[4][8];
	int refLevel_assignLevel_Correct[4][8];
	for(int i=0;i<4;++i)
		for(int j=0;j<8;++j)
			refLevel_assignLevel_All[i][j] = refLevel_assignLevel_Correct[i][j] = 0;
	for(int i=0;i<CluN;++i)
	{
		int maxR = 0,sumR = 0, majorId = 0;
		for(int j=0;j<GenoN;++j)
		{
			if(cluMtx[i][j] > maxR)
			{
				maxR = cluMtx[i][j];
				majorId = j;
			}
			sumR += cluMtx[i][j];
		}
		int refLevel = majorId/25;
		cout << i << '\t' << majorId << '\t' << majorId << '\t' << refLevel << '\t' << maxR << '\t' << sumR << '\t';
		for(int j=0;j<7;++j)
			cout << cluster[i].gid[j] << '\t';
		for(int j=0;j<7;++j)
			cout << cluster[i].maxScore[j] << '\t';
		for(int j=0;j<7;++j)
		//	cout << cluster[i].allScore[j] << '\t';
			cout << genome[majorId].taxoId[j] << '\t';
		int correctT = 0;
		for(;correctT<6;++correctT)
			if(cluster[i].gid[correctT] == genome[majorId].taxoId[correctT])
				break;
		cout << correctT << '\t';
		int AssignT = 0;
		for(;AssignT<7;++AssignT)
			if(cluster[i].maxScore[AssignT]/cluster[i].allScore[AssignT]>=GetTaxoThresh)
				break;
		cout << AssignT << '\t';
		cout << cluster[i].allScore[AssignT] << '\t';
		cout << endl;

		for(int j=0;j<GenoN;++j)
		{
			int trefLevel = j/25;
			refLevel_assignLevel_All[trefLevel][AssignT] += cluMtx[i][j];
			if(cluster[i].gid[AssignT] == genome[j].taxoId[AssignT])
				refLevel_assignLevel_Correct[trefLevel][AssignT] += cluMtx[i][j];
		}
	}

	{	
		cout << endl << endl;
		int asumR = 0,amaxR=0;
		for(int i=0;i<CluN;++i)
		{
			int maxR = 0,sumR = 0;
			for(int j=0;j<GenoN;++j)
			{
				if(cluMtx[i][j] > maxR)
					maxR = cluMtx[i][j];
				sumR += cluMtx[i][j];
			}
			amaxR += maxR;
			asumR += sumR;
		}
		cout << "precision:\t" << amaxR/(double)asumR << '\t' << amaxR << '\t' << asumR << endl;
		amaxR = 0;
		for(int j=0;j<GenoN;++j)
		{
			int maxR = 0;
			for(int i=0;i<CluN;++i)
				if(cluMtx[i][j] > maxR)
					maxR = cluMtx[i][j];
			amaxR += maxR;
		}
		cout << "sensitivity:\t" << amaxR/(double)asumR << '\t' << amaxR << '\t' << asumR << endl;
	}

	{
		cout << endl << endl;
		for(int i=0;i<4;++i)
		{
			for(int j=0;j<8;++j)
				cout << refLevel_assignLevel_Correct[i][j] << '\t';
			cout << endl;
		}
		cout << endl ;
		for(int i=0;i<4;++i)
		{
			for(int j=0;j<8;++j)
				cout << refLevel_assignLevel_All[i][j] << '\t';
			cout << endl;
		}
		cout << endl ;
		for(int i=0;i<4;++i)
		{
			for(int j=0;j<8;++j)
				cout << refLevel_assignLevel_Correct[i][j]/(double)refLevel_assignLevel_All[i][j] << '\t';
			cout << endl;
		}

		cout << endl ;
		for(int i=0;i<4;++i)
		{
			int tsum = 0;
			for(int j=0;j<8;++j)
				tsum += refLevel_assignLevel_All[i][j];
			for(int j=0;j<8;++j)
				cout << refLevel_assignLevel_All[i][j]/(double)tsum << '\t';
			cout << endl;
		}
	}
	////////////////////////////////////
	return 0;
}
