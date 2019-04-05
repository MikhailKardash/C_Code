#include <stdlib.h>
#include <stdio.h>
#include <float.h>
//#include "wave.h"
#include <sndfile.h>
#include <math.h>

#define Fs 8000
#define PI 3.14159265

void SampledSinusoid(float f, float L, float* x);
void iir(float* x, float* y);


int main(int argc, char *argv[])
{

	//Require 2 arguments: input file and output file
	if(argc < 2)
	{
		printf("Not enough arguments \n");
		return -1;
	}

	SF_INFO sndInfoOut;
	sndInfoOut.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
	sndInfoOut.channels = 1;
	sndInfoOut.samplerate = Fs;
	SNDFILE *sndFileOut = sf_open(argv[1], SFM_WRITE, &sndInfoOut);

	//Insert your code here
	float F = 8000;
	float L = 0.5;
	int samplenum = 8000/2;
	float* x;
	x = calloc(samplenum*32 + 32, sizeof(float));
    SampledSinusoid(F, L, x);
    float* y = calloc(samplenum*32 + 32, sizeof(float));
    iir(x,y);

    for(int i = 0; i < samplenum*32 + 32; i++) {
        //Use y+i or x+i for parts 8/7 respectively
        sf_writef_float(sndFileOut, y+i, 1);
    }

	sf_write_sync(sndFileOut);
	sf_close(sndFileOut);

	for (int k = 0; k < 10 ;k++) {
        printf("%f\n",*(y+k));
	}
	printf("\n");
    for (int k = 0; k <10; k++) {
        printf("%f\n",*(x+k));
    }
	free(x);
	free(y);

	return 1;
}

void SampledSinusoid(float f, float L, float* x)
{
    float notes[] = {0, 440 ,495 ,550 ,587 ,660 ,733, 825};
	float notes2[] = {0, 206 ,220, 248 ,275 ,293 ,330 ,367};
	int nu=0;
	int do1=1;
	int re=2;
	int mi=3;
	int fa=4;
	int so=5;
	int la=6;
	int ti=7;
	int dob=1;
	int reb=2;
	int mib=3;
	int fab=4;
	int sob=5;
	int lab=6;
	int tib=7;
	int melody[]={do1,nu,do1,nu,so,nu,so,nu,la,nu,la,nu,so,nu,nu,nu,fa,nu,fa,nu,mi,nu,mi,nu,re,nu,re,nu,do1,nu,nu,nu};
	int chorus[]={dob,sob,mib,sob,dob,sob,mib,sob,dob,lab,fab,lab,dob,sob,mib,sob,tib,sob,fab,sob,dob,sob,mib,sob,tib,sob,fab,sob,dob,sob,mib,sob};
    float q = f*L;
    int max_iter = (int) q;
    float temp;
    for (int i = 0; i < 32; i++) {
        for (int j = 0; j < max_iter; j++) {
            double f1 = (double) notes[melody[i]];
            double f2 = (double) notes2[chorus[i]];
            int totindex = i*max_iter + j;
            double time = ((double) totindex)/((double) f);
            float x1 = (float) 0.6*sin(2*PI*f1*time);
            float x2 = (float) 0.4*sin(2*PI*f2*time);
            *(x + totindex) = x1 + x2;
        }
    }
}

void iir(float* x, float* y)
{
    float b[] = {1, -3.1820023, 3.9741082, -2.293354, 0.52460587}; //y[n] coeffs
    float a[] = {0.62477732, -2.444978, 3.64114, -2.444978, 0.62477732}; //x[n] coeffs
    //initial conditions.
    float temp[5] = {*(x),*(x+1),*(x+2),*(x+3),*(x+4)};
    float temp2[4];
    *(y) = a[0]/b[0]*temp[0];
    temp2[0] = *(y);
    *(y+1) = (a[1]*temp[0] + a[0]*temp[1] - b[1]*temp2[0])/b[0];
    temp2[1] = *(y+1);
    *(y+2) = (a[2]*temp[0] + a[1]*temp[1]+ a[0]*temp[2] - b[2]*temp2[0]- b[1]*temp2[1])/b[0];
    temp2[2] = *(y+2);
    *(y+3) = (a[3]*temp[0] + a[2]*temp[1]+ a[1]*temp[2] + a[0]*temp[3] - b[3]*temp2[0]- b[2]*temp2[1] - b[1]*temp2[2])/b[0];
    for (int i = 4; i < 4000*32 + 32; i++) {
        temp[0] = *(x+i-4);
        temp[1] = *(x+i-3);
        temp[2] = *(x+i-2);
        temp[3] = *(x+i-1);
        temp[4] = *(x+i);
        temp2[0] = *(y+i-4);
        temp2[1] = *(y+i-3);
        temp2[2] = *(y+i-2);
        temp2[3] = *(y+i-1);
        *(y+i) = (a[4]*temp[0] + a[3]*temp[1]+ a[2]*temp[2] + a[1]*temp[3] + a[0]*temp[4] - b[4]*temp2[0]- b[3]*temp2[1] - b[2]*temp2[2]- b[1]*temp2[3])/b[0];
    }//and now we have y[i];
}
