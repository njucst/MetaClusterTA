#include <iostream>
#include "BWTGenerator.h"
#include <cstdio>
#include <cstring>
#include <string>
#include <dirent.h>
#include <cstdlib>
using namespace std;
void List(const char *path, int level)
{
	struct dirent *ent = NULL;
	DIR *pDir;
	if((pDir = opendir(path)) != NULL)
	{
		while(NULL != (ent = readdir(pDir)))
		{
			if(ent->d_type == 8)					// d_type:8-??,4-??
			{
				for(int i = 0; i < level; i++)
					printf("--");
				printf("%s/%s\n", path,ent->d_name);
			}
			else if(ent->d_name[0] != '.')
			{
				for(int i = 0; i < level; i++)
					printf("+ ");
				char *p = new char[strlen(path) + strlen(ent->d_name) + 2];
				strcpy(p, path);
				strcat(p, "/");
				strcat(p, ent->d_name);
				printf("%s/%s\n", path,ent->d_name);
				List(p, level+1);					// ???????
				free(p);
			}
		}
		closedir(pDir);
	}
}
int main()
{
//	cout << sizeof(unsigned short) << endl;
//	List("/home/ywang/Hybrid",1); 
//	char str[] = "ACACGT";
	char str[] = "ATGGTTCCGATGCCATTGCCGGAGTCGTCAATATCATTTTGAAATCTGACAATCATGGTGGCAGTGCCAGAAGTCAGATCGGACAGACCTATGCCGGAGATGGATTGGTCGGACAGGCCGGTTTTAATAAGGGATTCAAAATTGGCCATTCCGGTTTCTTTGATGTCGCCTTCGATTTTCGTCATCAAAATCATACCAGCCGTGATGGTATCGATAGCCGAACCCAGCGCCATAGTTTGAAAGTTGTCGGTGATCCGATGGCGACCCGCTATAATCTAGCGATCAATGCCGGTTATGATTTTGGAAATGGCATCGAAATTTATACGACTGATAGCTATTCGCATCGCAATMATAATATAACATTATTTAATAATTTATATCCAATAAATAAATACTTTATACCCATAAATTATATTTAAATTTTATAATTAAATAAAATAATTATNNNAAAAATTATCGAGGAGAAGCCAAATTATCTAAAAATAAAGAGGGGGCGTATGAAAAATTTTATTAAAAANNNNNNNTTTCTTTTTTCCTTTTCATGGTTTTTCTTGCAATATCGATCTTGATTTGATTGCACCTGAATTGATTGACCATAMATCAGAGCTTTATCAGGATGAAGCCTTTGTCCAGCATCAATAAACGGCAATTTCTTTCGCTTTCTGCTCTGACCGCTTTTTTACCAATGGCAACCCAAGTCTTGGGTAAAAATTCCGGTAAAGAGGCCTCGAAGAATTTGGMACGCACAAGGGGTGACCGTAGCTTCGACGCTTGAGGCTCAAAATATATCGGCATGGTTCTTTACCAATGGTGCCAGCACCCGCACCCGCGGTTTGGATTTTACM";
	BWTGenerator sfa;

/*	unsigned str2[10000];
	for(int i=0;i<sizeof(str);++i)
		str2[i] = str[i];
	unsigned SA[10000];
	sfa.dc3(str2,SA,sizeof(str),'Z');*/


	sfa.calBWT(str);//result);
	char* result = sfa.getBWT();
	cout << result << endl;
	unsigned *SA = sfa.getSA();

	for(int i=0;i<=strlen(str);++i)
		cout <<SA[i] << '_';
	cout << endl;
	return 0;
}
