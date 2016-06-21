#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// rtlsdr-to-gqrx Copyright 2014 Paul Brewer KI6CQ
// License: CC BY-SA 3.0 or GNU GPL 3.0
// IQ file converter 
// from rtl_sdr recording format -- interleaved unsigned char
// to gqrx/gnuradio .cfile playback format -- complex64 

/*	TODO
 *	- Templatize for various input file formats (complex, int8 etc)
 *	- command line args
 *	- put on website
 *	- figure out exact normalization desired
 */
int main(int argc, char *argv[])
{
  int byte1, byte2;  			// int -- not unsigned char -- see fgetc man page
  signed char fc;			// TODO difference between char and int?
  const size_t fc_size = sizeof(fc);
  FILE *infile,*outfile;
  const char *infilename = argv[1];
  const char *outfilename = argv[2];
  if (argc<3){
    printf("usage:  rtlsdr-to-gqrx infile outfile\n");
    exit(1);
  }
  infile=fopen(infilename,"rb");
  outfile=fopen(outfilename,"wb");
  if ((infile==NULL) || (outfile==NULL)){
    printf("Error opening files\n");
    exit(1);
  }
  while ((byte1=fgetc(infile)) != EOF){
    if ((byte2=fgetc(infile)) == EOF){
      exit(0);
    }

/*	c  < 0.04  -> no satellites			USING PARK DATA
 *      c == 0.04  -> 1 satellite
 *	c == 0.06  -> 2 satellites           2 preamble candidates
 *	c == 0.075 -> 4 satellite segfault
 *      c == 0.09  -> 4 satellite segfault
 *      c == 0.099 -> 4 satellite segfault
 *	c == 0.1   -> 4 satellite           3 preamble candidates
 *      c == 0.11  -> 4 satellite segfault
 *	c == 0.125 -> 4 satellite segfault
 *	c == 0.250 -> 4 satellite segfault
 *	c == 0.375 -> 4 satellite no craash 2 preamble candidates
 *	c  > 0.5 -> crash
 */
    float coeff = 0.1;				// 0.5 crashes, 0.1 works < 0.03 no acq
    float i = (byte1-127);			// -1 to 1
    float q = (byte2-127);			// -1 to 1
    float A = sqrt(i*i + q*q);			// amplitude -> -1.414 to 1.414
    float p = atan2( q, i);
    fc = A * cos(p) / (1.414 * coeff);		// 1.414 normalizes sqrt. coeff to scale.

    fwrite(&fc,fc_size,1,outfile);
  }
  return 0;
}
