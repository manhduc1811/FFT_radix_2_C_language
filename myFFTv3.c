// mducng, SoC, D2D, g2touch
// FFT version 3
// Usage:
// to complie: gcc -lm myFFTv3.c -o fft3
// to run:     ./fft test32.txt 32 2
// FFT output in outFFT_real, outFFT_image

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#define PI 3.14159265359
#define MAXPOW 24

struct complex
{
    double r;
    double i;
};

int pow_2[MAXPOW];
int pow_4[MAXPOW];

void twiddle(struct complex *W, int N, double stuff)
{
    W->r=cos(stuff*2.0*PI/(double)N);
    W->i=-sin(stuff*2.0*PI/(double)N);
}

void bit_reverse_reorder(struct complex *W, int N)
{
    int bits, i, j, k;
    double tempr, tempi;
    
    for (i=0; i<MAXPOW; i++)
	if (pow_2[i]==N) bits=i;

    for (i=0; i<N; i++)
    {
		j=0;
		for (k=0; k<bits; k++)
			if (i&pow_2[k]) j+=pow_2[bits-k-1];

		if (j>i)  /** Only make "up" swaps */
		{
			tempr=W[i].r;
			tempi=W[i].i;
			W[i].r=W[j].r;
			W[i].i=W[j].i;
			W[j].r=tempr;
			W[j].i=tempi;
		}
    }
}
void bit_r4_reorder(struct complex *W, int N)
{
    int bits, i, j, k;
    double tempr, tempi;
    
    for (i=0; i<MAXPOW; i++)
	if (pow_2[i]==N) bits=i;

    for (i=0; i<N; i++)
    {
	j=0;
	for (k=0; k<bits; k+=2)
	{
	    if (i&pow_2[k]) j+=pow_2[bits-k-2];
	    if (i&pow_2[k+1]) j+=pow_2[bits-k-1];
	}

	if (j>i)  /** Only make "up" swaps */
	{
	    tempr=W[i].r;
	    tempi=W[i].i;
	    W[i].r=W[j].r;
	    W[i].i=W[j].i;
	    W[j].r=tempr;
	    W[j].i=tempi;
	}
    }
}

/** RADIX-4 FFT ALGORITHM */
void radix4(struct complex *x, int N)
{ 
    int    n2, k1, N1, N2;
    struct complex W, bfly[4];

    N1=4;
    N2=N/4;
    
    /** Do 4 Point DFT */ 
    for (n2=0; n2<N2; n2++)
    {
	/** Don't hurt the butterfly */
	bfly[0].r = (x[n2].r + x[N2 + n2].r + x[2*N2+n2].r + x[3*N2+n2].r);
	bfly[0].i = (x[n2].i + x[N2 + n2].i + x[2*N2+n2].i + x[3*N2+n2].i);

	bfly[1].r = (x[n2].r + x[N2 + n2].i - x[2*N2+n2].r - x[3*N2+n2].i);
	bfly[1].i = (x[n2].i - x[N2 + n2].r - x[2*N2+n2].i + x[3*N2+n2].r);

	bfly[2].r = (x[n2].r - x[N2 + n2].r + x[2*N2+n2].r - x[3*N2+n2].r);
	bfly[2].i = (x[n2].i - x[N2 + n2].i + x[2*N2+n2].i - x[3*N2+n2].i);

	bfly[3].r = (x[n2].r - x[N2 + n2].i - x[2*N2+n2].r + x[3*N2+n2].i);
	bfly[3].i = (x[n2].i + x[N2 + n2].r - x[2*N2+n2].i - x[3*N2+n2].r);


	/** In-place results */
	for (k1=0; k1<N1; k1++)
	{
	    twiddle(&W, N, (double)k1*(double)n2);
	    x[n2 + N2*k1].r = bfly[k1].r*W.r - bfly[k1].i*W.i;
	    x[n2 + N2*k1].i = bfly[k1].i*W.r + bfly[k1].r*W.i;
	}
    }
    
    /** Don't recurse if we're down to one butterfly */
    if (N2!=1)
	for (k1=0; k1<N1; k1++)
	{
	    radix4(&x[N2*k1], N2);
	}
}

/** RADIX-2 FFT ALGORITHM */
void radix2(struct complex *data, int N)
{
    int    n2, k1, N1, N2;
    struct complex W, bfly[2];

    N1=2;
    N2=N/2;
    
    /** Do 2 Point DFT */
    for (n2=0; n2<N2; n2++)
    {
	/** Don't hurt the butterfly */
	twiddle(&W, N, (double)n2);
	bfly[0].r = (data[n2].r + data[N2 + n2].r);
	bfly[0].i = (data[n2].i + data[N2 + n2].i);
	bfly[1].r = (data[n2].r - data[N2 + n2].r) * W.r - 
	    ((data[n2].i - data[N2 + n2].i) * W.i); 
	bfly[1].i = (data[n2].i - data[N2 + n2].i) * W.r +
	    ((data[n2].r - data[N2 + n2].r) * W.i);

	/** In-place results */
	for (k1=0; k1<N1; k1++)
	{
	    data[n2 + N2*k1].r = bfly[k1].r;
	    data[n2 + N2*k1].i = bfly[k1].i;
	}
    }
    
    /** Don't recurse if we're down to one butterfly */
    if (N2!=1)
	for (k1=0; k1<N1; k1++)
	    radix2(&data[N2*k1], N2);
}

void main(int argc, char *argv[])
{
    FILE   *infile;
	FILE   *outfileR,*outfileI;
    int    N, radix, numsamp;
    int    i;
    struct complex *data;
    
#ifndef GEN
    if (argc<4)
    {
	printf("usage:\n");
	printf("    fft [input file] [sequence length] [radix]\n");
	printf("        output: DFT\n");
	exit(1);
    }
    else if ((infile=fopen(argv[1], "r"))==NULL)
    {
	printf("Error reading input sequence file: %s\n", argv[1]);
	exit(1);
    }
    
    sscanf(argv[2], "%d", &N);
    sscanf(argv[3], "%d", &radix);
#endif

    /** Set up power of two arrays */
    pow_2[0]=1;
    for (i=1; i<MAXPOW; i++)
	pow_2[i]=pow_2[i-1]*2;
    pow_4[0]=1;
    for (i=1; i<MAXPOW; i++)
	pow_4[i]=pow_4[i-1]*4;
    
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

    if (radix==2) radix2(data, N);
    if (radix==4) radix4(data, N);

    if (radix==2) bit_reverse_reorder(data, N);
    if (radix==4) bit_r4_reorder(data, N);

	outfileR=fopen("outFFT_real.txt", "w");
	outfileI=fopen("outFFT_image.txt", "w");
	
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
