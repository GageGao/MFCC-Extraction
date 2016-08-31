#ifndef _MFCC_H
#define _MFCC_H

#define FFROM	(300)
#define FTO		(7600)
#define BANKNUM	(24)
#define FFT_NUM	(512)
#define MFCCNUM	(13)
#define FWIDTH	(8000)

#define DELTA_K	(2)

#define FRAME_SIZE	(400)
#define FRAME_STEP	(160)

int InitParam();
int UInitParam();

//pAudio [in]: one frame data
//l [in]: the length of data 
//pMfcc [out]: save static mfcc
int ExtractStaticMfcc(short *pAudio,int l,float *pMfcc);

//pAudio [in]: raw input
//nFrames [in]: frames of raw input
//pMfcc [out]: save mfcc
//e.g. order=0 represents static, order=1 represents delta
int mfcc(short *pAudio,int nFrames,float *pMfcc,int order);

#endif //_MFCC_H