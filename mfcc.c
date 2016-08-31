#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <memory.h>
#include "mfcc.h"

#define PI   3.14159265358979
#define TPI  6.28318530717959     /* PI*2 */

float* g_pmfpoint;	//BANKNUM+2
float* g_pfpoint;	//BANKNUM+2
float* g_pftomf;	//FFT_NUM/2
float* g_pcos;		//BANKNUM

void Realft (float *s,int l);
void GenParam();

//初始化相关参数
int InitParam()
{
	g_pmfpoint = (float *)malloc((BANKNUM+2)*sizeof(float));
	g_pfpoint = (float *)malloc((BANKNUM+2)*sizeof(float));
	g_pftomf = (float *)malloc((FFT_NUM/2)*sizeof(float));
	g_pcos = (float *)malloc(BANKNUM*MFCCNUM*sizeof(float));

	if ( NULL == g_pmfpoint || NULL == g_pfpoint || NULL == g_pftomf || NULL == g_pcos )
		return 0;

	GenParam();

	return 1;
}

int UInitParam()
{
	free(g_pmfpoint);
	free(g_pfpoint);
	free(g_pftomf);
	free(g_pcos);

	return 1;
}

//
int ExtractStaticMfcc(short *pAudio,int l,float *pMfcc)
{
	float Data[FFT_NUM+1] = {0};
	float ef[FFT_NUM/2] = {0};
	int i;

	if ( l > FRAME_SIZE )
		return 0;

	for ( i = 0; i < l; i++ )
		Data[i+1] = (float)pAudio[i];

	Realft(Data,FFT_NUM);
	for ( i = 0; i < FFT_NUM/2; i++ )
		ef[i] = sqrt(Data[2*i+1]*Data[2*i+1] + Data[2*(i+1)]*Data[2*(i+1)]);
	GetMfcc(ef,pMfcc);

	return 1;

}

//delta
int mfcc(short *pAudio,int nFrames,float *pMfcc,int order)
{
	int dim,i;
	short *pTmp;

	dim = MFCCNUM*(order+1);
	pTmp = pAudio;
	for ( i = 0; i < nFrames; i++ )
	{
		ExtractStaticMfcc(pTmp,FRAME_SIZE,pMfcc+i*dim);
		pTmp += FRAME_STEP;
	}

	for ( i = 1; i < order+1; i++ )
	{
		CalcDelta(pMfcc,nFrames,dim,i);
	}

	return 1;
}


int CalcDelta(float *pMfcc,int nFrames,int dim,int order)
{
	int i,j;
	float Tmp;

	for ( i = 0; i < nFrames; i++ )
	{
		if ( i < DELTA_K )
			for ( j = 0; j < MFCCNUM; j++ )
				pMfcc[i*dim+MFCCNUM*order+j] = pMfcc[(i+1)*dim+MFCCNUM*(order-1)+j]-pMfcc[i*dim+MFCCNUM*(order-1)+j];
		else if ( i >= nFrames-DELTA_K)
			for ( j = 0; j < MFCCNUM; j++ )
				pMfcc[i*dim+MFCCNUM*order+j] = pMfcc[i*dim+MFCCNUM*(order-1)+j]-pMfcc[(i-1)*dim+MFCCNUM*(order-1)+j];
		else
			for ( j = 0; j < MFCCNUM; j++ )
			{
				Tmp = pMfcc[(i+1)*dim+MFCCNUM*(order-1)+j]-pMfcc[(i-1)*dim+MFCCNUM*(order-1)+j]+2*(pMfcc[(i+2)*dim+MFCCNUM*(order-1)+j]-pMfcc[(i-2)*dim+MFCCNUM*(order-1)+j]);
				pMfcc[i*dim+MFCCNUM*order+j] = Tmp/10;
			}
	}

	return 1;
}

int GetMfcc(float *pfft, float *pmfcc)
{
	float pbf[BANKNUM];
	int i,k,j;

	for( i=0;i<BANKNUM;i++)
	{
		float startf=g_pfpoint[i];
		float midf=g_pfpoint[i+1];
		float endf=g_pfpoint[i+2];
		
		pbf[i]=0.0;
		for(k=0;k<(FFT_NUM/2);k++)
		{
			float tf=2.0*k*FWIDTH/FFT_NUM;
			if((tf>=startf)&&(tf<midf))
				pbf[i]=pbf[i]+pfft[k]*(tf-g_pfpoint[i])/(g_pfpoint[i+1]-g_pfpoint[i]);

			if((tf>=midf)&&(tf<endf))
				pbf[i]=pbf[i]+pfft[k]*(g_pfpoint[i+2]-tf)/(g_pfpoint[i+2]-g_pfpoint[i+1]);
		}
	}

	for(i=0;i<MFCCNUM;i++)
	{
		pmfcc[i]=0.0;
		for(j=0;j<BANKNUM;j++)
		{
			if(pbf[i]<0.0000001)
				pbf[i]=0.0000001;
		    pmfcc[i]=pmfcc[i]+log10(pbf[j])*g_pcos[i*BANKNUM+j];
		}
	}

	return 1;
}

void GetParamForMel()
{
	float mffrom=1125.0*log(FFROM/700.0+1.0);
	float mfto=1125.0*log(FTO/700.0+1.0);
	int i,k;

	float mintval=(mfto-mffrom)/(BANKNUM+1);
	for(i=0;i<(BANKNUM+2);i++)
	{
		g_pmfpoint[i]=mffrom+i*mintval;
		g_pfpoint[i]=700.0*(exp(g_pmfpoint[i]/1125.0)-1.0);
	}
	for(k=0;k<(FFT_NUM/2);k++)
	{
		float f=2.0*k*FWIDTH/FFT_NUM;
		g_pftomf[k]=1125.0*log(f/700.0+1.0);
	}
}

void GenParam()
{
	int i,j;

	GetParamForMel();

	for(i=0;i<MFCCNUM;i++)
	{
		for(j=0;j<BANKNUM;j++)
			g_pcos[i*BANKNUM+j]=cos(PI*i*(2.0*j+1.0)/(2.0*BANKNUM));
	}
}

/* EXPORT-> FFT: apply fft/invfft to complex s */
void FFT(float *s, int l,int invert)
{
	int ii,jj,n,nn,limit,m,j,inc,i;
	float wx,wr,wpr,wpi,wi,theta;
	float xre,xri,x;

	n=l;
	nn=n / 2; j = 1;
	for (ii=1;ii<=nn;ii++) {
		i = 2 * ii - 1;
		if (j>i) {
			xre = s[j]; xri = s[j + 1];
			s[j] = s[i];  s[j + 1] = s[i + 1];
			s[i] = xre; s[i + 1] = xri;
		}
		m = n / 2;
		while (m >= 2  && j > m) {
			j -= m; m /= 2;
		}
		j += m;
	};
	limit = 2;
	while (limit < n) {
		inc = 2 * limit; theta = TPI / limit;
		if (invert) theta = -theta;
		x = sin(0.5 * theta);
		wpr = -2.0 * x * x; wpi = sin(theta); 
		wr = 1.0; wi = 0.0;
		for (ii=1; ii<=limit/2; ii++) {
			m = 2 * ii - 1;
			for (jj = 0; jj<=(n - m) / inc;jj++) {
				i = m + jj * inc;
				j = i + limit;
				xre = wr * s[j] - wi * s[j + 1];
				xri = wr * s[j + 1] + wi * s[j];
				s[j] = s[i] - xre; s[j + 1] = s[i + 1] - xri;
				s[i] = s[i] + xre; s[i + 1] = s[i + 1] + xri;
			}
			wx = wr;
			wr = wr * wpr - wi * wpi + wr;
			wi = wi * wpr + wx * wpi + wi;
		}
		limit = inc;
	}
	if (invert)
		for (i = 1;i<=n;i++) 
			s[i] = s[i] / nn;

}

/* EXPORT-> Realft: apply fft to real s */
void Realft (float *s,int l)
{
	int n, n2, i, i1, i2, i3, i4;
	float xr1, xi1, xr2, xi2, wrs, wis;
	float yr, yi, yr2, yi2, yr0, theta, x;

	n=l / 2; n2 = n/2;
	theta = PI / n;
	FFT(s,l,0);
	x = sin(0.5 * theta);
	yr2 = -2.0 * x * x;
	yi2 = sin(theta); yr = 1.0 + yr2; yi = yi2;
	for (i=2; i<=n2; i++) {
		i1 = i + i - 1;      i2 = i1 + 1;
		i3 = n + n + 3 - i2; i4 = i3 + 1;
		wrs = yr; wis = yi;
		xr1 = (s[i1] + s[i3])/2.0; xi1 = (s[i2] - s[i4])/2.0;
		xr2 = (s[i2] + s[i4])/2.0; xi2 = (s[i3] - s[i1])/2.0;
		s[i1] = xr1 + wrs * xr2 - wis * xi2;
		s[i2] = xi1 + wrs * xi2 + wis * xr2;
		s[i3] = xr1 - wrs * xr2 + wis * xi2;
		s[i4] = -xi1 + wrs * xi2 + wis * xr2;
		yr0 = yr;
		yr = yr * yr2 - yi  * yi2 + yr;
		yi = yi * yr2 + yr0 * yi2 + yi;
	}
	xr1 = s[1];
	s[1] = xr1 + s[2];
	s[2] = 0.0;
}
