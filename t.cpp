#include "LearnLamda.h"
using namespace std;
int main(int argc, char*argv[]){
	string bac_info = "/home/ywang/SA/Bac_Sa.idx.info";
	string node_dmp = "/home/ywang/DB/gi_taxid/nodes.dmp";
	string mtx_path = argv[1];

	NCBI_nodes_dmp Dmp;
	Dmp.init(node_dmp.c_str());

//	map<int,map<int,double> >Lamda;
	map<int,vector<double> >Lamda;
	map<int,vector<double> >Prec;
	getLamda(Dmp, bac_info, mtx_path,Lamda,Prec);
	return 0;
}

