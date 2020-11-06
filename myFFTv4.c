// mducng, SoC, D2D, g2touch
// FFT version 4
// Usage:
// to complie: gcc -lm myFFTv4.c -o fft4
// to run:     ./fft test32.txt 32 2
// FFT output in outFFT_real, outFFT_image
// Radix 2
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
////////////////////////////////////////////////////////////
#define PI 3.14159265359
#define MAXPOW 24
struct complex
{
    double r;
    double i;
};
int pow_2[MAXPOW];
////////////////////////////////////////////////////////////
void twiddle(struct complex *W, int N, double stuff)
{
    W->r=cos(stuff*2.0*PI/(double)N);
    W->i=-sin(stuff*2.0*PI/(double)N);
}
////////////////////////////////////////////////////////////
void bit_reverse_reorder(struct complex *W, int N)
{
    int 	bits, i, j, k;
    double 	tempr, tempi;
    
    for (i=0; i<MAXPOW; i++)
	{
		if (pow_2[i]==N) 
			bits=i;
	}

    for (i=0; i<N; i++)
    {
		j=0;
		for (k=0; k<bits; k++)
		{
			if (i&pow_2[k]) 
			{
				j+=pow_2[bits-k-1];
			}
		}
		if (j>i)  /** Only make "up" swaps */
		{
			tempr	=	W[i].r;
			tempi	=	W[i].i;
			W[i].r	=	W[j].r;
			W[i].i	=	W[j].i;
			W[j].r	=	tempr;
			W[j].i	=	tempi;
		}
    }
}
////////////////////////////////////////////////////////////
/** RADIX-2 FFT ALGORITHM */
void radix2(struct complex *data, int N)
{
    int    n2, k1, N1, N2;
    struct complex W, bfly[2];
    N1	=	2;
    N2	=	N/2;
    /** Do 2 Point DFT */
    for (n2=0; n2<N2; n2++)
    {
		/** Don't hurt the butterfly */
		twiddle(&W, N, (double)n2);
		bfly[0].r = (data[n2].r + data[N2 + n2].r);
		bfly[0].i = (data[n2].i + data[N2 + n2].i);
		bfly[1].r = (data[n2].r - data[N2 + n2].r) * W.r 
					- ((data[n2].i - data[N2 + n2].i) * W.i); 
		bfly[1].i = (data[n2].i - data[N2 + n2].i) * W.r 
					+ ((data[n2].r - data[N2 + n2].r) * W.i);
		/** In-place results */
		for (k1=0; k1<N1; k1++)
		{
			data[n2 + N2*k1].r = bfly[k1].r;
			data[n2 + N2*k1].i = bfly[k1].i;
		}
    }
    /** Don't recurse if we're down to one butterfly */
    if (N2!=1)
	{
		for (k1=0; k1<N1; k1++)
			radix2(&data[N2*k1], N2);
	}
}
////////////////////////////////////////////////////////////
void main(int argc, char *argv[])
{
    FILE   *infile;
	FILE   *outfileR,*outfileI;
    int    N, radix, numsamp;
    int    i;
    struct complex *data;
	// get data from file
#ifndef GEN
    if (argc<3)
    {
		printf("usage:\n");
		printf("    fft [input file] [sequence length]\n");
		printf("        output: DFT\n");
		exit(1);
    }
    else if ((infile=fopen(argv[1], "r"))==NULL)
    {
		printf("Error reading input sequence file: %s\n", argv[1]);
		exit(1);
    }
    sscanf(argv[2], "%d", &N);
#endif
    /** Set up power of two arrays */
    pow_2[0]=1;
    for (i=1; i<MAXPOW; i++)
	{
		pow_2[i]=pow_2[i-1]*2;
	}
    if ((data=malloc(sizeof(struct complex)*(size_t)N))==NULL)
    {
		fprintf(stderr, "Out of memory!\n");
		exit(1);
    }
#ifndef GEN
    for (i=0; i<N; i++)
    {
		fscanf(infile, "%lf", &data[i].r);
		data[i].i=0.0;
    }
#endif    
    // radix-2 algorithm 
	radix2(data, N);
	// bit ordering
	bit_reverse_reorder(data, N);
	// to output FFT data
	//// real part
	outfileR=fopen("outFFT_real.txt", "w");
	fprintf(outfileR,"%%%%%%%%%%%%%%%%\n");
	//// image part
	outfileI=fopen("outFFT_image.txt", "w");
	fprintf(outfileI,"%%%%%%%%%%%%%%%%\n");
	for (i=0; i<N; i++)
	{
	    fprintf(outfileR,"%f\n", data[i].r);
	    fprintf(outfileI,"%f\n", data[i].i);		
	}
#ifndef GEN
    fclose(infile);
	fclose(outfileR);
	fclose(outfileI);
#endif
}
