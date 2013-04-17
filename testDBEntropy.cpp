#include <iostream>
#include "DBEntropy.h"

using namespace std;

int main()
{
	string path = "/home/ywang/DB/list/BacTaxo.csv";
	DBEntropy DB(path);
	DB.getStatis();
	return 0;
}
