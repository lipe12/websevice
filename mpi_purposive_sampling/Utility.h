#pragma once


#include <string>


using namespace std;


#define VERYSMALL 0.00005
#define	VERYBIG 999999999999999


#define NOINTVALUE -9999
#define NOFLOATVALUE -9999
#define NODOUBLEVALUE -9999


#define MAX_ITEM_NUM 10000



void connectPatch(int irow, int icol, float *resultBuf, int *patchBuf, int patchID, int nrow, int ncol);

void makeDir(string name);

char *extract_filename(char *filespec);

char *strncpynm(char *s1, char *s2, int n, int m);


