#include <iostream>
#include  "BWT_FidGenerator.h"
using namespace std;
void usage()
{
	cerr << "testBWT_FidGenerator CompleteGenomeFolderPath gi_taxid_nucl.dmp outfilename" << endl;
	exit(-1);
}
int main(int argc, char* argv[])
{
	if(argc < 4)
		usage();
	BWT_FidGenerator bwtgenerator;
	cerr << "initiation of genorator is finished." << endl;
	bwtgenerator.calBWTfromCompleGenome(argv[1],argv[2],argv[3]);
	return 0;
}
