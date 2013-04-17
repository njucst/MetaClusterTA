#include <cstdlib>
#include <iostream>
#include <fstream>

#include "CtgsReads.h"

using namespace std;

void usage()
{
	cerr << "TestDB DB_Path_File.txt " << endl;
	exit(-1);
}
int main(int argc, char* argv[])
{
	if(argc<2)
		usage();

	return 0;
}

