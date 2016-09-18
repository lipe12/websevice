#include "Utility.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef __linux__
#include <sys/stat.h>
#endif

#ifdef _WIN32 || _WIN64
#include <direct.h>
#endif


/*
void connectPatch(int irow, int icol, float *resultBuf, int *patchBuf, int patchID, int nrow, int ncol)
{
	int deltaRow[8] = {1, -1, 1, -1, 1, 0, -1, 0};

	int deltaCol[8] = {1, -1,-1, 1, 0, 1, 0, -1};

	for(int i=0; i<8; i++)
	{
		if(irow + deltaRow[i] >=0 && irow + deltaRow[i]< nrow && icol+deltaCol[i] >=0 && icol+deltaCol[i]<ncol)
		{
			if(resultBuf[(irow + deltaRow[i])*ncol+icol + deltaCol[i]] > NOFLOATVALUE && patchBuf[(irow + deltaRow[i])*ncol+ icol + deltaCol[i]] == 0)  //8 neighborhood
			{
				patchBuf[(irow + deltaRow[i])*ncol + icol+deltaCol[i]] = patchID;
				connectPatch(irow + deltaRow[i],icol+deltaCol[i], resultBuf, patchBuf, patchID, nrow, ncol);

			}
		}
	}

}*/



void makeDir(string name)
{
    #ifdef __linux__
	  mkdir(name.c_str() , S_IRWXU | S_IRWXG | S_IRWXO);
	#endif
	#ifdef _WIN32 || _WIN64
	  _mkdir(name.c_str());
	#endif
}



char *extract_filename(char *filespec)
{
 int len=0, j=0;
 char *tmpstr;

 len=strlen(filespec);
 if ((tmpstr=(char *)malloc(sizeof(char)*len))==NULL)
    {
    printf("\n\tOut of memory in extract_filename, program aborted");
    exit(-1);
    }
 filespec[len]='\0';
 strcpy(tmpstr, filespec);
 for (j=len-1; j>=0; j--)
    {
    if (filespec[j]=='/' || filespec[j]=='\\' || filespec[j]==':')
       {
       j++;
       break;
       }
    }
 if (j<0)
    j=0;
 strncpynm(tmpstr, filespec, j, len-1);
 return(tmpstr);

}



char* strncpynm(char *s1, char *s2, int n, int m){
	int i, j=0;
	for (i=n;i<=m && s2[i]!='\0';i++){
		s1[j]=s2[i];
		j++;
	}
	s1[j]='\0';
	return(s1);
}

