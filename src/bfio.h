#ifndef MIAPDE_BFIO_H
#define MIAPDE_BFIO_H


/* ATTENTION: This is a custom lib modification by Pascal Peter (2015) */
/* changes: new setter and getter functions for bytes */

/*--------------------------------------------------------------------------*/
/*                                                                          */
/*  File:       bfio.h                                                      */
/*                                                                          */
/*  Purpose:    Definitions and functions for bitwise reading and writing   */
/*              of files. Implementation e.g. includes a new datatype       */
/*              wrapping the FILE datatype and io-functions wrapping        */
/*              stdio-functions trying to provide a similar thus familiar   */
/*              interface.                                                  */
/*                                                                          */
/*  Author:     Leif Bergerhoff, bergerhoff@mia.uni-saarland.de             */
/*  Date:       May 19, 2015                                                */
/*                                                                          */
/*--------------------------------------------------------------------------*/

#include <stdio.h>
#include "defines.h"

/* table of BITFILE-modes and corresponding FILE-modes
 * BITFILE      underlying FILE
 * r            r
 * r+           r+
 * w            w+
 * w+           w+
 * a            r+/w+
 * a+           r+/w+ */

/* definition of the datatype BITFILE for bitfiles */
typedef struct _BITFILE {
    FILE *file;                     /* the underlying FILE pointer */
    long file_size;                 /* size of FILE in bytes, number of bytes */
    char *path;                     /* location of (bit)FILE */
    MFILE_MODE mode;                /* the file access mode */
    unsigned char *buffer_r;        /* reading buffer */
    unsigned char *buffer_w;        /* wrinting buffer */
    long *pos_byte_r;               /* current byteposition of reading buffer
                                       in FILE */
    long *pos_byte_w;               /* current byteposition of writing buffer
                                       in FILE */
    int *pos_r;                     /* current position in reading buffer */
    int *pos_w;                     /* current position in writing buffer */
    int valid;                      /* number of valid bit in last byte */
    unsigned char header;           /* first byte of FILE,
                                       first bit contains the bit used for
                                       filling up the last byte of FILE,
                                       bit 2 - 8 contain the beginning
                                       bitstream */
    unsigned char fill;             /* the bit used for filling up the last
                                       byte '0' or '1' */
} BITFILE;

/* definition of the datatype BFILE */
typedef BITFILE BFILE;

/* setup reading/writing buffers AND positions according to the bitfile's
   mode */
int bfsetbuf( BFILE *fp );

/* Open a bitstream */
BFILE *bfopen( const char *path, const char *mode );

/* close a bitstream */
int bfclose( BFILE *fp );

/* flush a bitstream */
int bfflush( BFILE *stream );

/* print out information about bitstream */
int bfinfo( BFILE *fp );

/* set the bitfile pointer to a specific bitposition in bitstream */
int bfseek( BFILE *stream, long offset, int whence );

/* tell bitposition in bitstream (not including header bit) */
long bftell( BFILE *stream );

/* write one bit to the bitstream */
int bfputb( const int c, BFILE *stream );

/* write multiple bits to the bitstream */
size_t bfwrite( const void *ptr, size_t size, size_t nmemb, BFILE *stream );

/* write multiple bits to the bitstream using a c-string as input */
size_t bfstrtobf( const char *ptr, BFILE *stream );

/* read one bit from the bitstream */
int bfgetb( BFILE *stream );

/* read multiple bits from the bitstream */
size_t bfread( void *ptr, size_t size, size_t nmemb, BFILE *stream );

/* added by Pascal Peter: custom lib! */
void set_byte(BFILE *binary_file, long b);
long get_byte(BFILE *binary_file);
void set_bit(BFILE *binary_file, long b);
long get_bit(BFILE *binary_file);
void set_bits(BFILE *binary_file, long b, long bitdepth);
long get_bits(BFILE *binary_file, long bitdepth);
#endif
