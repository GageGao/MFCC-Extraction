#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <memory.h>
#include "mfcc.h"

int usage()
{
	printf(".exe 提取特征的列表 保存特征文件名 特征阶数\n");
}

int main(int argc,char *argv[])
{
	FILE *fp,*lp,*mp;
	int nSize,i,j;
	char *pBuf;
	short *pAudio;
	float *pMfcc;
	int nFrames,order,dim;
	char fn[256];

	order = atoi(argv[3]); 
	dim = MFCCNUM*(order+1);

	lp = fopen(argv[1],"rb");
	mp = fopen(argv[2],"wb");
	if ( NULL == lp )
	{
		printf("can't open file %s!",argv[1]);
		getchar();
		return 0;
	}

	InitParam();

	while (!feof(lp))
	{
		memset(fn,0,256);
		fscanf(lp,"%s",fn);
		fp = fopen(fn,"rb");
		if ( NULL == fp )
		{
			printf("can't open file!\n");
			getchar();
			return 0;
		}
		fseek(fp,0,SEEK_END);
		nSize = ftell(fp);
		fseek(fp,0,SEEK_SET);
		pBuf = (char *)malloc(nSize);
		fread(pBuf,1,nSize,fp);
		fclose(fp);

		nFrames = (int)((nSize/2-(FRAME_SIZE-FRAME_STEP))/FRAME_STEP);
		pMfcc = (float *)malloc(nFrames*dim*sizeof(float));
		if ( NULL == pMfcc )
			return 0;

		pAudio = (short *)pBuf;
		mfcc(pAudio,nFrames,pMfcc,order);
		for ( i = 0; i < nFrames; i++ )
		{
			for ( j = 0; j < dim; j++ )
			{
				fprintf(mp,"%f\t",pMfcc[dim*i+j]);
			}
			fprintf(mp,"\n");
		}

		free(pMfcc);
	}

	UInitParam();

	fclose(mp);
	fclose(lp);

	return 0;
}