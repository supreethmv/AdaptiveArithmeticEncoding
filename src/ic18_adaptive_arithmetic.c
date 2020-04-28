/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                  PROGRAMMING EXERCISE FOR THE LECTURE                    */
/*                            IMAGE COMPRESSION                             */
/*                   ADAPTIVE INTEGER ARITHMETIC CODING                     */
/*                  (Copyright by Pascal Peter, 5/2018)                     */
/*                                                                          */
/*--------------------------------------------------------------------------*/

/* global includes */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <getopt.h>
#include <omp.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

/* local includes */
#include "alloc.h"              /* memory allocation */
#include "image_io.h"           /* reading and writing pgm and ppm images */
#include "bfio.h"               /* writing and reading of bitfiles */

/* defines */
/* version */
#define VERSION 1.1-05-18
/* supported input formats */
#define FORMAT_PGM 0
#define FORMAT_PPM 1
/* auxiliary functions */
#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))
/* maximum number of channels */
#define MAXCHANNELS 3
/* define maximum grey value */
#define MAXGREYVALUE 255
/* console formatting */
#define ONE_UP "\033[5D\033[1A"
#define CLRLINE "\033[K"

/* definition of compressed image datatype and struct */
typedef struct ImageData ImageData;
struct ImageData {
  long*** orig_rgb;    /* original image (RGB) */
  long size_orig;      /* size of raw input image */
  long nx, ny, nc;     /* image dimensions and channels */
  long encoder;        /* encoder to apply */
  long** orig_hist;    /* channel histograms of original image */
  double** orig_prob;  /* channel probabilities of original image */ 
};

long log2(long x) {
  int logx = 0;
  while (x >>= 1) ++logx;
  return logx;
}


/*--------------------------------------------------------------------------*/
void init_image (ImageData* img) {
  /* set image variables to default values */
  img->nx=img->ny=img->nc=0;
  img->encoder=0;
  img->size_orig=0;
  img->orig_rgb=0;
  img->orig_hist=0;
  img->orig_prob=0;
}

/*--------------------------------------------------------------------------*/
void alloc_image (ImageData* img,long nx,long ny) {
  /* allocate all image data arrays */
  if (img->orig_rgb == 0)
    alloc_long_cubix(&img->orig_rgb,MAXCHANNELS,nx+2,ny+2);
  if (img->orig_hist == 0)
    alloc_long_matrix(&img->orig_hist,MAXCHANNELS,MAXGREYVALUE);
  if (img->orig_prob == 0)
    alloc_matrix(&img->orig_prob,MAXCHANNELS,MAXGREYVALUE);
}

/*--------------------------------------------------------------------------*/
void destroy_image (ImageData* img) {
  /* disalloc all image data arrays */
  long nx, ny;
  nx = img->nx;
  ny = img->ny;
  /* disalloc images */
  if (img->orig_rgb != 0)
    disalloc_long_cubix(img->orig_rgb,MAXCHANNELS,nx+2,ny+2);
  if (img->orig_hist != 0)
    disalloc_long_matrix(img->orig_hist,MAXCHANNELS,MAXGREYVALUE);
  if (img->orig_prob != 0)
    disalloc_matrix(img->orig_prob,MAXCHANNELS,MAXGREYVALUE);
}

/*--------------------------------------------------------------------------*/
void write_comment_string(ImageData* img, char* additional, char* comments)
/* writes all important information of an R-EED image into a comment string
 * parameter additional allows to add more custom text, set to 0 if not needed.
 * make sure that the comments array is preallocated with sufficient space. */
{
  /* write parameter values in comment string */
  comments[0] = '\0';
  comment_line (comments, (char*)"# IMAGE COMPRESSION - PROGRAMMING EXERCISE\n");
}

/*--------------------------------------------------------------------------*/
/* for user interaction */
/* Returns true if it is possible to read from cin */
int canread()
{
  struct timeval tv;
  tv.tv_sec = 0;
  tv.tv_usec = 0;
  fd_set read_set;
  FD_ZERO(&read_set);
  FD_SET(0,&read_set);
  select(1,&read_set,NULL,NULL, &tv);
  return (FD_ISSET(0,&read_set));
}

/*--------------------------------------------------------------------------*/
void copy_cubix(double ***source, double ***target,
		long sc, long nc, long nx, long ny) {
  /* copy cubix from channel sc to channel nc */
  long i,j,m;
  for (m=0;m<nc;m++) {
    for (i=1;i<=nx;i++) {
      for (j=1;j<=ny;j++) {
        target[m][i][j]=source[m][i][j];
      }
    }
  }
}

/*--------------------------------------------------------------------------*/
void copy_cubix_long(long ***source, long ***target,
                    long sc, long nc, long nx, long ny) {
  /* copy cubix (only nc channels, starting with channel sc) */
  long i,j,m;
  for (m=sc;m<nc;m++) {
    for (i=1;i<=nx;i++) {
      for (j=1;j<=ny;j++) {
        target[m][i][j]=source[m][i][j];
      }
    }
  }
}

/*--------------------------------------------------------------------------*/
void copy_matrix_long(long **source, long **target,long nx, long ny) {
  /* copy input matrix to target matrix */
  long i,j;
  for (i=0;i<nx;i++) {
    for (j=0;j<ny;j++) {
      target[i][j]=source[i][j];
    }
  }
}

/*--------------------------------------------------------------------------*/
void copy_vector_long(long *source, long *target, long nx) {
  /* copy input vector to target vector */
  long i;
  for (i=0;i<nx;i++) {
    target[i]=source[i];
  }
}

/*--------------------------------------------------------------------------*/
long get_size_of_file(char* filename) {
  /* compute and return size of file with arbitrary content */
  FILE   *infile;             /* input file */
  long i;
  /* open file */
  infile = fopen(filename,"r");
  if (NULL == infile)
    {
      printf("Could not open file '%s' for reading, aborting (10)\n",filename);
      exit(-1);
    }
  /* Find length of file */
  fseek(infile,0,SEEK_END);
  i = ftell(infile);
  fclose(infile);
  return i;
}

/*--------------------------------------------------------------------------*/
double get_compression_ratio(char* original_file,
                             char* compressed_file) {
  /* computes compression ratio between two files on the hard disk (X:1) */
  long original_size;
  long compressed_size;
  
  original_size = get_size_of_file(original_file);
  compressed_size = get_size_of_file(compressed_file);
  
  return (double)original_size/(double)compressed_size;
}

/*--------------------------------------------------------------------------*/
void convert_image_to_vector (long*** image_matrix,
                              long bx, long by, long nx, long ny, long nc,
                              long* image_vector) {
  long i,j,c;

  /* parse image channel by channel, row by row */
  for (c = 0; c < nc; c++) 
    for (j = by; j < ny+by;j++)
      for (i = bx; i < nx+bx;i++) {
        image_vector[c*nx*ny+(j-by)*nx+(i-bx)]=image_matrix[c][i][j];
      }
}

/*--------------------------------------------------------------------------*/
void convert_vector_to_image (long* image_vector,
                              long bx, long by, long nx, long ny, long nc,
                              long*** image_matrix) {
  long i,j,c;

  /* parse image channel by channel, row by row */
  for (c = 0; c < nc; c++) 
    for (j = by; j < ny+by;j++)
      for (i = bx; i < nx+bx;i++) {
        image_matrix[c][i][j]=image_vector[c*nx*ny+(j-by)*nx+(i-bx)];
        /*  printf("%ld %ld %ld (%ld): %ld\n",
           c,i,j,c*nx*ny+(j-by)*nx+(i-bx),image_matrix[c][i][j]); */ 
      }
}

/*--------------------------------------------------------------------------*/
void print_usage_message() {
  printf("./compress -i input_file -o output_prefix -S target_ratio [optional parameters]\n");
  printf("list of mandatory paramters:\n");
  printf("-i input_file   (string): uncompressed image, e.g. \"image.ppm\"\n");
  printf("-o out_prefix   (string): prefix for output files, e.g. \"my_image\"\n");
  printf("list of optional paramters:\n");
  printf("-d discretisation parameter M (int)   : yields M:=2^(d+2), default: 12\n");
  printf("-r readjustment parameter   r (float) : default: 0.3\n");
  printf("-D debug file                 (string): filename for debug information\n");  
}

/*--------------------------------------------------------------------------*/

/* apply WNC algorithm for adaptive arithmetic integer encoding */
void encode_adaptive_wnc(
  long* sourceword,   /* array containing n numbers from {0,...,s}
                         where s is the end of file symbol */
  long n,             /* length of sourceword */
  long s,             /* size of source alphabet */
  double r,           /* rescaling parameter */
  long  M,            /* WNC discretisation parameter M */
  FILE* debug_file,   /* 0 - no output, 1 - debug output to file */
  BFILE* compressed)  {/* binary file for compressed bitstring */

  long i,j;      /* loop variables */
  long L,H;      /* low and high interval endpoints of current interval */
  long oldL;     /* temporary variable to preserve L for computing new interval*/
  long C;        /* sum of all counters */
  long k;        /* underflow counter */
  long* counter; /* array of counters for adaptive probabilities */
  long M12=M/2, M14=M/4, M34=3*M/4;  /* time savers */
  long symbol;   /* index of current symbol in counter array */
  long csum;     /* sum of counters 0,...,symbol-1 */
                        
  /* allocate memory */
  alloc_long_vector(&counter,s);
  
  /* initialise counters and C */
  for (i=0;i<s;i++) {
    counter[i]=1;
  }
  C = s;

  if (debug_file != 0) {
    fprintf(debug_file,"n: %ld, s: %ld, r: %f, M: %ld\n",n,s,r,M);
  }
  
  if ((double)C>(double)M/4.0+2.0) {
    if (debug_file != 0) {
      fprintf(debug_file,
              "M=%ld is too small (C=%ld), setting M to %ld\n",M,C,8*C);
    }
    M=8*C;
    M12=M/2; M14=M/4; M34=3*M/4;
  }

  /* initialise interval endpoints and underflow counter */

  /* insert your code here */

  /* encode sourceword */
  for (i=0;i<=n;i++) {
    if (debug_file != 0) {
      fprintf(debug_file,"sourceword[%ld]=%ld\n",i,sourceword[i]);
    }
    
    /* underflow expansions/rescaling */
    while (1) {
      /* insert your code here */

      if (debug_file != 0) {
        fprintf(debug_file,"k: %ld, [%ld, %ld)\n",k,L,H);
      }
      break;
    }
    
    /* readjustment */

    /* insert your code here */
    
    /* encode symbol */
    if (i<=n-1) {
      /* insert your code here */
      
      if (debug_file != 0) {
        fprintf(debug_file,"new [L,H) = [%ld,%ld)\n",L,H);
      }
      counter[symbol]++; C++;
    } else {
      /* last step */
      if (debug_file != 0) {
        fprintf(debug_file,"last interval - written bits:");
      }

      /* insert your code here */
    }
  }
}


/*--------------------------------------------------------------------------*/

/* apply WNC algorithm for adaptive arithmetic integer encoding */
void decode_adaptive_wnc(
    BFILE* compressed,  /* binary file with compressed bitstring */
    long n,             /* length of sourceword */
    long s,             /* size of source alphabet */
    double r,           /* rescaling parameter */
    long  M,            /* WNC discretisation parameter M */
    FILE* debug_file,   /* 0 - no output, 1 - debug output to file */
    long* sourceword)  {/* array containing n numbers from {0,...,s}
                           where s is the end of file symbol */
  long i,j;      /* loop variables */
  long L,H;      /* low and high interval endpoints of current interval */
  long oldL;     /* temporary variable to preserve L for computing new interval*/
  long C;        /* sum of all counters */
  long* counter; /* array of counters for adaptive probabilities */
  long M12=M/2, M14=M/4, M34=3*M/4;  /* time savers */
  long symbol;   /* index of current symbol in counter array */
  long csum;     /* sum of counters 0,...,symbol-1 */
  long w;        /* variable for finding correct decoding inverval */
  long v;        /* partial dyadic fraction */
  long b;        /* auxiliary variable for reading individual bits */
  long N;        /* auxiliary variable for determining number of initial bits
                    for v */
 
  /* allocate memory */
  alloc_long_vector(&counter,s);
  
  /* initialise counters and C */
  for (i=0;i<s;i++) {
    counter[i]=1;
  }
  C = s;

  if ((double)C>(double)M/4.0+2.0) {
    if (debug_file != 0) {
      fprintf(debug_file,
              "M=%ld is too small (C=%ld), setting M to %ld\n",C,M,8*C);
    }
    M=8*C;
    M12=M/2; M14=M/4; M34=3*M/4;
  }

  /* initialise interval endpoints*/

  /* insert your code here */
 
  /* read first bits of codeword to obtain initival v */

  /* insert your code here */

  if (debug_file != 0) {
    fprintf(debug_file,"initial v: %ld (%ld first bits from coded file)\n",
           v,N);
  }

  /* decode sourceword */
  for (i=0;i<n;i++) {
  
    /* underflow expansions/rescaling */
    while (1) {
      
      /* insert your code here */

      if (debug_file != 0) {
        fprintf(debug_file,"v: %ld, [%ld, %ld)\n",v,L,H);
      }
      break;
    }
    
    /* readjustment */

    /* insert your code here */

    /* decode symbol */

    /* insert your code here */

    if (debug_file != 0) {
      fprintf(debug_file,"[c_i,c_i-1) = [%ld %ld) ",csum,csum+counter[symbol]);
      fprintf(debug_file,"w: %ld symbol[%ld]: %ld, new [L,H)=[%ld,%ld)\n",
              w,i,symbol,L,H);
    }
  }
}


/*--------------------------------------------------------------------------*/

int main (int argc, char **args)

{
  /* user interaction */
  char   ch;                   /* for reading single chars */
  char   used[256];            /* flag for all possible chars, used to check for
                                  multiple occurences of input parameters */
  long len;
  
  /* parameters */

  /* file IO */
  char *output_file = 0;       /* file name of output image */
  char *input_file = 0;        /* file name of uncompressed image */
  char tmp_file[1000];         /* string for intermediate filenames */
  char total_file[1000];       /* file name of total compressed output */
  char comments[10000];        /* string for comments */
  char *program_call;          /* call of the compression program */
  char extension[5];           /* extension of input file */
  long format;                 /* format of input file */

  /* image information */
  long nx, ny, nc;                /* image dimensions */

  /* loop variables */
  long i;

  /* image struct that contains all information
   * about original image, compressed image, and parameters */
  ImageData image;

  /* compression/decompression */
  long flag_compress=0;
  double r = 0.3;             /* readjustment parameter */
  long   M = (long)pow(2,12);  /* discretisation parameter */
  long   s = 256;             /* length of source alphabet */
  BFILE* binary_file=0;       /* binary file for compressed bitstring */
  char*  debug_file=0;        /* filename for writing debug information */
  FILE*  dfile=0;             /* file for writing debug information */
  long*  sourceword;          /* sourceword for compression 
                                 (image in vector format) */
  
  printf ("\n");
  printf ("PROGRAMMING EXERCISE FOR IMAGE COMPRESSION\n\n");
  printf ("**************************************************\n\n");
  printf ("    Copyright 2018 by Pascal Peter                \n");
  printf ("    Dept. of Mathematics and Computer Science     \n");
  printf ("    Saarland University, Saarbruecken, Germany    \n\n");
  printf ("    All rights reserved. Unauthorised usage,      \n");
  printf ("    copying, hiring, and selling prohibited.      \n\n");
  printf ("    Send bug reports to                           \n");
  printf ("    peter@mia.uni-saarland.de                     \n\n");
  printf ("**************************************************\n\n");

  /* ARGUMENT PROCESSING ******************************************************/
  init_image(&image); /* initialise image with standard parameters */

  for (i = 0; i <= 255; i++) {
    used[i] = 0;
  }

  while ((ch = getopt(argc,args,"i:o:d:r:D:")) != -1) {
    used[(long)ch]++;
    if (used[(long)ch] > 1) {
      printf("Duplicate parameter: %c\n",ch);
      printf("Please check your input again.\n");
    }
    
    switch(ch) {
    case 'i': input_file = optarg;break;
    case 'o': output_file = optarg;break;
    case 'd': M = pow(2,atoi(optarg));break;
    case 'D': debug_file = optarg;break;
    case 'r': r = atof(optarg);break;
    default:
      printf("Unknown argument.\n");
      print_usage_message();
      return 0;
    }
  }

  if (output_file == 0 || input_file == 0) {
    printf("ERROR: Missing mandatory parameter, aborting.\n");
    print_usage_message();
    return 0;
  }


  /* prepare file names */
  sprintf(total_file,"%s.coded",output_file);

  /* create reboot string */
  len = 0;
  for (i=0;i<argc;i++) {
    len += strlen(args[i]) + 1;
  }

  program_call = (char*)malloc( sizeof(char) * ( len + 3 ) );
  sprintf(program_call," ");
  for (i=0;i<argc;i++) {
    sprintf(program_call, "%s%s ",program_call,args[i]);
  }

  /* DETERMINE COMPRESSION/DECOMPRESSION MODE**********************************/

  /* try to identify the file format via the extension */
  strcpy(extension, input_file+(strlen(input_file)-4));
  if      (!strcmp(extension, ".pgm")) format = FORMAT_PGM;
  else if (!strcmp(extension, ".ppm")) format = FORMAT_PPM;
  else {
    printf("ERROR: Extension %s not supported for input file, aborting.\n",
           extension);
    print_usage_message();
    return 0;
  }

  /* determine if we are in compression or decompression mode */
  if (format==FORMAT_PPM || format==FORMAT_PGM) {
    flag_compress = 1;
  } else {
    flag_compress = 0;
  }
  
  if (flag_compress == 1) {
    /* COMPRESS ***************************************************************/

    /* read input image */
    if (format==FORMAT_PPM) {
      read_ppm_and_allocate_memory(input_file,&image.nx,&image.ny,
                                   &image.orig_rgb);
      image.nc = 3;
      printf("Image %s loaded (PPM).\n",input_file);
    } else if (format==FORMAT_PGM) {
      read_pgm_header(input_file,&image.nx,&image.ny);
      alloc_long_cubix(&image.orig_rgb,MAXCHANNELS,image.nx+2,image.ny+2);
      read_pgm_and_allocate_memory(input_file,&image.nx,&image.ny,
                                   &image.orig_rgb[0]);
      image.nc = 1;
      printf("Image %s loaded (PGM).\n",input_file);
    }
    nx = image.nx; ny = image.ny; nc = image.nc;
    image.size_orig=get_size_of_file(input_file);

    /* allocate memory */
    alloc_image(&image,nx,ny);

    alloc_long_vector(&sourceword,nc*nx*ny);
    convert_image_to_vector(image.orig_rgb,1,1,nx,ny,nc,sourceword);

    /* compress with WNC */
    printf("Compressing with WNC (M=%ld, r=%f)\n",M,r);
    sprintf(tmp_file,"%s.wnc",output_file);
    binary_file = bfopen(tmp_file,"w");
    if (debug_file != 0) {
      dfile = fopen(debug_file,"w");
    }
    encode_adaptive_wnc(sourceword,nx*ny*nc,s,r,M,dfile,binary_file);
    bfclose(binary_file);
    printf("Resulting compression ratio: %f:1\n",
           get_compression_ratio(input_file,tmp_file));

    binary_file = bfopen(tmp_file,"r");
    printf("Decompressing with WNC (M=%ld, r=%f)\n",M,r);
    decode_adaptive_wnc(binary_file,nx*ny*nc,s,r,M,dfile,sourceword);
    convert_vector_to_image(sourceword,1,1,nx,ny,nc,image.orig_rgb);
    bfclose(binary_file);
    if (debug_file !=0) {
      fclose(dfile);
    }
     
    /* write image data */
    write_comment_string(&image,0,comments);
    if (format==FORMAT_PPM) {
      sprintf(tmp_file,"%s.ppm",output_file);
      write_ppm(image.orig_rgb, nx, ny, tmp_file, comments);
    } else {
      sprintf(tmp_file,"%s.pgm",output_file);
      write_pgm(image.orig_rgb[0], nx, ny, tmp_file, comments);
    }
  } else {
    /* DECOMPRESS *************************************************************/

  }
  
  /* ---- free memory  ---- */
  destroy_image(&image);
  free(program_call);

  return(0);
}
