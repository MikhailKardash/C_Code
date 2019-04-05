#include <stdlib.h>
#include <stdio.h>
#include <float.h>
//#include "wave.h"
#include <sndfile.h>
#include <math.h>

#define Fs 48000
#define pi 3.14159265
#define r 0.8

int main(int argc, char *argv[])
{

	//Require 2 arguments: input file and output file
	if(argc < 3)
	{
		printf("Not enough arguments \n");
		return -1;
	}

	SF_INFO sndInfoOut;
    SF_INFO sndInfo;
	SNDFILE *sndFile = sf_open(argv[1], SFM_READ, &sndInfo);
    SNDFILE *sndFilecorrupt = sf_open(argv[2], SFM_READ, &sndInfo);
    
	sndInfoOut.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
	sndInfoOut.channels = 1;
	sndInfoOut.samplerate = Fs;
	SNDFILE *sndFileOut = sf_open(argv[3], SFM_WRITE, &sndInfoOut);
    SNDFILE *sndFileOut2 = sf_open(argv[4], SFM_WRITE, &sndInfoOut);
    SNDFILE *sndFileOut3 = sf_open(argv[5], SFM_WRITE, &sndInfoOut);
    
	//Insert your code here
	int F = 48000;
	float L = 0.5;
	int samplenum = 48000/2;
	float* x;
    float* noise;
    float* noisex;
    x = calloc(960000, sizeof(float));
    noise = calloc(960000, sizeof(float));
    noisex = calloc(960000, sizeof(float));
    int ii;
    for(ii=0; ii < 960000; ii++)
    {
        sf_readf_float(sndFile, x+ii, 1);
        sf_readf_float(sndFilecorrupt, noise+ii, 1);
        *(noisex+ii)= *(x+ii) + *(noise+ii);
    }
    
    float* y = calloc(960000 + F/2, sizeof(float));
    
    for (int k=0; k<960000; k++)
    {
            *(y+k) = *(y+k) + 0.7*(*(x+k));
            *(y+k+F/2) = *(y+k+F/2) + 0.3*(*(x+k));       
    }
    
    int n[] = {-7,-6,-5,-4,-3,-2,-1,10,1,2,3,4,5,6,7};
    float hn[15]; 
    for(int l=0; l<15;l++)
    {
            hn[l] = 2*sin(0.1*pi*n[l])/(pi*n[l]);
    }
    hn[8]= -0.8;
    float* z = calloc(960000 + 14, sizeof(float));
    for (int k = 0; k < 960000; k++)
    {
        *(z+k+7) = *(noisex+k);
    }
   
    float* outfir = calloc(960000, sizeof(float));
    float temp;
    
    for (int k = 0; k < 960000; k++)
    {
        temp = 0;
        for (int j=0;j <= 14 ; j++) 
        {
            temp = temp + hn[j]*(*(noisex+k+j));
        }
        *(outfir+k) = temp;
    }
    
    float* outiir = calloc(960000, sizeof(float));
    float k1=-2*cos(pi/10);
    float k2= -2*r*cos(pi/10);
    float k3=r*r;
    *(outiir) = *(z+7) + k1*(*(z+6)) + *(z+5);
    *(outiir+1) = -k2*(*(outiir)) + *(z+8) + k1*(*(z+7)) + *(z+6);
    for(int ll=2; ll < 960000 ; ll++) 
    {
    *(outiir+ll)=-k2*(*(outiir+ll-1))+*(z+7+ll)+k1*(*(z+6+ll))+*(z+5+ll)-k3*(*(outiir+ll-2));
    }
   
    for(int i = 0; i < 960000 + F/2 ; i++) {
       
        sf_writef_float(sndFileOut, y+i, 1);
    }
    sf_write_sync(sndFileOut);
	sf_close(sndFileOut);
    
    for(int i = 0; i < 960000; i++) {
        sf_writef_float(sndFileOut2, outfir+i, 1);
    }
   
    sf_write_sync(sndFileOut2);
	sf_close(sndFileOut2);
    
  for(int i = 0; i < 960000; i++) {
        sf_writef_float(sndFileOut3, outiir+i, 1);
    }
    sf_write_sync(sndFileOut3);
	sf_close(sndFileOut3);
    
	for (int k = 0; k < 10 ;k++) {
        printf("%f\n",*(y+k));
	}
	printf("\n");
    for (int k = 0; k <10; k++) {
        printf("%f\n",*(x+k));
    }
	free(x);
	free(y);
    free(z);
    free(outfir);
    free(noisex);
    free(noise);
    free(outiir);
    
	return 1;
}