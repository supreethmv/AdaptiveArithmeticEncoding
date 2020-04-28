/*--------------------------------------------------------------------------*/
/*                                                                          */
/*  File:       bfio.c                                                      */
/*                                                                          */
/*  Purpose:    Functions for bitwise reading and writing of files.         */
/*              Implementation includes io-functions wrapping               */
/*              stdio-functions trying to provide a similar thus familiar   */
/*              interface.                                                  */
/*                                                                          */
/* Author:      Leif Bergerhoff, bergerhoff@mia.uni-saarland.de             */
/* Date:        May 19, 2015                                                */
/*                                                                          */
/*--------------------------------------------------------------------------*/

#include "bfio.h"
#include "math.h"
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <stdio.h>

//#define BFIO_WRITE_CHECKS

/*--------------------------------------------------------------------------*/

static int initBitfile

    ( BFILE *fp      ,          /* pointer to bitfile structure */
      const char *path,         /* filename of bitfile */
      const MFILE_MODE mode  )  /* mode of bitfile */

    /* setup the bitfile structure
       return 0 on success, 1 else. */

{
    size_t len = 0;             /* length of path string */

    /* verify pointer */
    if( fp == NULL )
        return 1;

    /* initialize structure */
    memset( fp, 0, sizeof( BFILE ) );
    fp->file = NULL;
    fp->path = NULL;
    fp->buffer_r = NULL;
    fp->buffer_w = NULL;
    fp->pos_byte_r = NULL;
    fp->pos_byte_w = NULL;
    fp->pos_r = NULL;
    fp->pos_w = NULL;

    /* open the underlying FILE in correct mode AND
       estimate the FILEsize AND set fp->file_size accordingly */
    switch( mode ) {

        case MF_R:
            fp->file = fopen( path, "r" );
            if( fp->file == NULL )
                goto FreeAndFail;

            if( fseek( fp->file, 0, SEEK_END ) == -1 )
                goto FreeAndFail;

            if( ( fp->file_size = ftell( fp->file ) ) == -1 )
                goto FreeAndFail;

            if( fseek( fp->file, 0, SEEK_SET ) == -1 )
                goto FreeAndFail;

            break;

        case MF_RP:
            fp->file = fopen( path, "r+" );
            if( fp->file == NULL )
                goto FreeAndFail;

            if( fseek( fp->file, 0, SEEK_END ) == -1 )
                goto FreeAndFail;

            if( ( fp->file_size = ftell( fp->file ) ) == -1 )
                goto FreeAndFail;

            if( fseek( fp->file, 0, SEEK_SET ) == -1 )
                goto FreeAndFail;

            break;

        case MF_W:
            fp->file = fopen( path, "w+" );
            if( fp->file == NULL )
                goto FreeAndFail;

            fp->file_size = 0;

            break;

        case MF_WP:
            fp->file = fopen( path, "w+" );
            if( fp->file == NULL )
                goto FreeAndFail;

            fp->file_size = 0;

            break;

        case MF_A:
            /* check for existence of FILE */
            fp->file = fopen( path, "r+" );
            if( fp->file == NULL ) {
                /* FILE does NOT exist yet -> create a new FILE */
                fp->file = fopen( path, "w+" );
                if( fp->file == NULL )
                    goto FreeAndFail;

                fp->file_size = 0;
            }
            else {
                /* FILE exists */
                if( fseek( fp->file, 0, SEEK_END ) == -1 )
                    goto FreeAndFail;

                if( ( fp->file_size = ftell( fp->file ) ) == -1 )
                    goto FreeAndFail;

                if( fseek( fp->file, 0, SEEK_SET ) == -1 )
                    goto FreeAndFail;
            }

            break;

        case MF_AP:
            /* check for existence of FILE */
            fp->file = fopen( path, "r+" );
            if( fp->file == NULL ) {
                /* FILE does NOT exist yet -> create a new FILE */
                fp->file = fopen( path, "w+" );
                if( fp->file == NULL )
                    goto FreeAndFail;

                fp->file_size = 0;
            }
            else {
                /* FILE exists */
                if( fseek( fp->file, 0, SEEK_END ) == -1 )
                    goto FreeAndFail;

                if( ( fp->file_size = ftell( fp->file ) ) == -1 )
                    goto FreeAndFail;

                if( fseek( fp->file, 0, SEEK_SET ) == -1 )
                    goto FreeAndFail;
            }

            break;

        default:
            return 1;

    }

    /* keep the FILE's path */
    len = strlen( path ) + 1;
    fp->path = (char*)malloc( len );
    if( fp->path == NULL )
        goto FreeAndFail;

    strncpy( fp->path, path, len );

    /* keep the bitfiles mode */
    fp->mode = mode;

    /* setup buffer and position pointer */
    if( bfsetbuf( fp ) != 0 )
        goto FreeAndFail;

    /* for readable files: setup header AND setup filling bit AND setup valid
       bits AND set FILE pointer to correct position */
    if( ( fp->mode == MF_R || fp->mode == MF_RP || fp->mode == MF_A ||
          fp->mode == MF_AP ) && fp->file_size > 0 ) {

        /* read header byte */
        if( fseek( fp->file, 0, SEEK_SET ) == -1 )
            goto FreeAndFail;

        if( fread( &fp->header, 1, 1, fp->file ) != 1 )
            goto FreeAndFail;

        /* estimate the filling bit */
        fp->fill = fp->header >> 7;

        /* read last byte */
        if( fseek( fp->file, -1, SEEK_END ) == -1 )
            goto FreeAndFail;

        if( fread( fp->buffer_r, 1, 1, fp->file ) != 1 )
            goto FreeAndFail;

        /* estimate the number of valid bit in last byte */
        while(
            ( ( ( *fp->buffer_r >> fp->valid ) & 1 ) &&
            fp->fill == 1 ) ||
            ( !( ( *fp->buffer_r >> fp->valid ) & 1 ) &&
            fp->fill == 0 )
        ) {
            fp->valid++;

            if( fp->file_size == 1 && fp->valid == 7 )
                break;
        }

        fp->valid = 8 - fp->valid;

        /* setup reading buffer */
        *fp->buffer_r = fp->header;

        if( fp->mode == MF_A || fp->mode == MF_AP ) {
            /* set buffer to last (potentially incomplete) byte */

            if( fp->valid == 8 ) {
                /* complete */
                *fp->pos_byte_w = fp->file_size;
                *fp->buffer_w = 0;
                *fp->pos_w = 0;
            }
            else {
                /* incomplete */
                if( fseek( fp->file, -1, SEEK_END ) == -1 )
                    goto FreeAndFail;

                if( fread( fp->buffer_w, 1, 1, fp->file ) != 1 )
                    goto FreeAndFail;

                *fp->pos_w = fp->valid;

                if( fseek( fp->file, -1, SEEK_END ) == -1 )
                    goto FreeAndFail;

                *fp->pos_byte_w = fp->file_size - 1;
            }

        }
        else {
            /* set FILE pointer to first byte again */
            if( fseek( fp->file, 0, SEEK_SET ) == -1 )
                goto FreeAndFail;
        }

    }

    return 0;

FreeAndFail:
    if( fp->file != NULL )
        fclose( fp->file );

    free( fp->path );

    memset( fp, 0, sizeof( BFILE ) );
    fp->file = NULL;
    fp->path = NULL;

    return 1;
}

/*--------------------------------------------------------------------------*/

int bfsetbuf

    ( BFILE *fp )

    /* setup reading/writing buffers AND positions according to the bitfile's
       mode.
       returns 0 on success, 1 else. */

{

    /* check pointer validity */
    if( fp == NULL )
        return 1;

    if( fp->buffer_r != NULL || fp->buffer_w != NULL ||
        fp->pos_byte_r != NULL || fp->pos_byte_w != NULL ||
        fp->pos_r != NULL || fp->pos_w != NULL )
        return 1;

    /* setup reading/writing buffers AND positions according to the bitfile's
       mode */
    switch( fp->mode ) {

        case MF_R:
            /* initialize reading buffer */
            fp->buffer_r = (unsigned char*)malloc( 1 );
            if( fp->buffer_r == NULL )
                goto FreeAndFail;

            *fp->buffer_r = 0;

            /* initialize buffer position */
            fp->pos_byte_r = (long*)malloc( sizeof( long ) );
            if( fp->pos_byte_r == NULL )
                goto FreeAndFail;

            *fp->pos_byte_r = 0;

            /* initialize reading position */
            fp->pos_r = (int*)malloc( sizeof( int ) );
            if( fp->pos_r == NULL )
                goto FreeAndFail;

            *fp->pos_r = 1;

            /* initialize writing buffer */
            fp->buffer_w = fp->buffer_r;

            /* initialize buffer position */
            fp->pos_byte_w = fp->pos_byte_r;

            /* initialize writing position */
            fp->pos_w = fp->pos_r;

            break;

        case MF_RP:
            /* initialize reading buffer */
            fp->buffer_r = (unsigned char*)malloc( 1 );
            if( fp->buffer_r == NULL )
                goto FreeAndFail;

            *fp->buffer_r = 0;

            /* initialize buffer position */
            fp->pos_byte_r = (long*)malloc( sizeof( long ) );
            if( fp->pos_byte_r == NULL )
                goto FreeAndFail;

            *fp->pos_byte_r = 0;

            /* initialize reading position */
            fp->pos_r = (int*)malloc( sizeof( int ) );
            if( fp->pos_r == NULL )
                goto FreeAndFail;

            *fp->pos_r = 1;

            /* initialize writing buffer */
            fp->buffer_w = fp->buffer_r;

            /* initialize buffer position */
            fp->pos_byte_w = fp->pos_byte_r;

            /* initialize writing position */
            fp->pos_w = fp->pos_r;

            break;

        case MF_W:
            /* initialize writing buffer */
            fp->buffer_w = (unsigned char*)malloc( 1 );
            if( fp->buffer_w == NULL )
                goto FreeAndFail;

            *fp->buffer_w = 0;

            /* initialize buffer position */
            fp->pos_byte_w = (long*)malloc( sizeof( long ) );
            if( fp->pos_byte_w == NULL )
                goto FreeAndFail;

            *fp->pos_byte_w = 0;

            /* initialize writing position */
            fp->pos_w = (int*)malloc( sizeof( int ) );
            if( fp->pos_w == NULL )
                goto FreeAndFail;

            *fp->pos_w = 1;

            /* initialize reading buffer */
            fp->buffer_r = fp->buffer_w;

            /* initialize buffer position */
            fp->pos_byte_r = fp->pos_byte_w;

            /* initialize reading position */
            fp->pos_r = fp->pos_w;

            /* set number of valid bits */
            fp->valid = 1;

            break;

        case MF_WP:
            /* initialize writing buffer */
            fp->buffer_w = (unsigned char*)malloc( 1 );
            if( fp->buffer_w == NULL )
                goto FreeAndFail;

            *fp->buffer_w = 0;

            /* initialize buffer position */
            fp->pos_byte_w = (long*)malloc( sizeof( long ) );
            if( fp->pos_byte_w == NULL )
                goto FreeAndFail;

            *fp->pos_byte_w = 0;

            /* initialize writing position */
            fp->pos_w = (int*)malloc( sizeof( int ) );
            if( fp->pos_w == NULL )
                goto FreeAndFail;

            *fp->pos_w = 1;

            /* initialize reading buffer */
            fp->buffer_r = fp->buffer_w;

            /* initialize buffer position */
            fp->pos_byte_r = fp->pos_byte_w;

            /* initialize reading position */
            fp->pos_r = fp->pos_w;

            /* set number of valid bits */
            fp->valid = 1;

            break;

        case MF_A:
            /* initialize writing buffer */
            fp->buffer_w = (unsigned char*)malloc( 1 );
            if( fp->buffer_w == NULL )
                goto FreeAndFail;

            *fp->buffer_w = 0;

            /* initialize buffer position */
            fp->pos_byte_w = (long*)malloc( sizeof( long ) );
            if( fp->pos_byte_w == NULL )
                goto FreeAndFail;

            *fp->pos_byte_w = 0;

            /* initialize writing position */
            fp->pos_w = (int*)malloc( sizeof( int ) );
            if( fp->pos_w == NULL )
                goto FreeAndFail;

            *fp->pos_w = 1;

            /* initialize reading buffer */
            fp->buffer_r = fp->buffer_w;

            /* initialize buffer position */
            fp->pos_byte_r = fp->pos_byte_w;

            /* initialize reading position */
            fp->pos_r = fp->pos_w;

            break;

        case MF_AP:
            /* initialize writing buffer */
            fp->buffer_w = (unsigned char*)malloc( 1 );
            if( fp->buffer_w == NULL )
                goto FreeAndFail;

            *fp->buffer_w = 0;

            /* initialize buffer position */
            fp->pos_byte_w = (long*)malloc( sizeof( long ) );
            if( fp->pos_byte_w == NULL )
                goto FreeAndFail;

            *fp->pos_byte_w = 0;

            /* initialize writing position */
            fp->pos_w = (int*)malloc( sizeof( int ) );
            if( fp->pos_w == NULL )
                goto FreeAndFail;

            *fp->pos_w = 1;

            /* initialize reading buffer */
            fp->buffer_r = (unsigned char*)malloc( 1 );
            if( fp->buffer_r == NULL )
                goto FreeAndFail;

            *fp->buffer_r = 0;

            /* initialize buffer position */
            fp->pos_byte_r = (long*)malloc( sizeof( long ) );
            if( fp->pos_byte_r == NULL )
                goto FreeAndFail;

            *fp->pos_byte_r = 0;

            /* initialize reading position */
            fp->pos_r = (int*)malloc( sizeof( int ) );
            if( fp->pos_r == NULL )
                goto FreeAndFail;

            *fp->pos_r = 1;

            break;


        default:
            goto FreeAndFail;

    }

    return 0;

FreeAndFail:
    if( fp->buffer_r == fp->buffer_w )
        free( fp->buffer_r );
    else {
        free( fp->buffer_r );
        free( fp->buffer_w );
    }

    if( fp->pos_byte_r == fp->pos_byte_w )
        free( fp->pos_byte_r );
    else {
        free( fp->pos_byte_r );
        free( fp->pos_byte_w );
    }

    if( fp->pos_r == fp->pos_w )
        free( fp->pos_r );
    else {
        free( fp->pos_r );
        free( fp->pos_w );
    }

    fp->buffer_r = NULL;
    fp->buffer_w = NULL;
    fp->pos_byte_r = NULL;
    fp->pos_byte_w = NULL;
    fp->pos_r = NULL;
    fp->pos_w = NULL;

    return 1;
}


/*--------------------------------------------------------------------------*/

BFILE *bfopen

    ( const char *path,     /* 'path'-argument of fopen() */
      const char *mode )    /* 'mode'-argument of fopen() */

    /* Open a bitstream.
       for possible choices for mode refer to 'man fopen'.
       returns 0 on success, 1 else. */

{
    BFILE *new_bf = NULL;   /* new bitfile */

    /* allocate memory for bitfile structure */
    new_bf = (BFILE*)malloc( sizeof( BFILE ) );
    if( new_bf == NULL )
        return NULL;

    /* init bitfile in appropriate file mode */
    if( strncmp( mode, "r+", 2 ) == 0 ) {
        if( initBitfile( new_bf, path, MF_RP ) != 0 )
            goto free;
    }
    else if( strncmp( mode, "w+", 2 ) == 0 ) {
        if( initBitfile( new_bf, path, MF_WP ) != 0 )
            goto free;
    }
    else if( strncmp( mode, "a+", 2 ) == 0 ) {
        if( initBitfile( new_bf, path, MF_AP ) != 0 )
            goto free;
    }
    else if( *mode == 'r' ) {
        if( initBitfile( new_bf, path, MF_R ) != 0 )
            goto free;
    }
    else if( *mode == 'w' ) {
        if( initBitfile( new_bf, path, MF_W ) != 0 )
            goto free;
    }
    else if( *mode == 'a' ) {
        if( initBitfile( new_bf, path, MF_A ) != 0 )
            goto free;
    }

    return new_bf;

free:
    free( new_bf );

    return NULL;
}

/*--------------------------------------------------------------------------*/

int bfclose

    ( BFILE *fp )           /* pointer to bitstream */

    /* close a bitstream */

{

    /* validate argument */
    if( fp == NULL ) {
        errno = EBADF;
        return EOF;
    }

    /* flush buffer */
    if( fp->mode != MF_R && fp->mode != MF_CLOSED ) {
        if( bfflush( fp ) == EOF )
            return EOF;
    }

    if( fp->file != NULL )
        fclose( fp->file );

    free( fp->path );

    if( fp->buffer_r == fp->buffer_w )
        free( fp->buffer_r );
    else {
        free( fp->buffer_r );
        free( fp->buffer_w );
    }

    if( fp->pos_byte_r == fp->pos_byte_w )
        free( fp->pos_byte_r );
    else {
        free( fp->pos_byte_r );
        free( fp->pos_byte_w );
    }

    if( fp->pos_r == fp->pos_w )
        free( fp->pos_r );
    else {
        free( fp->pos_r );
        free( fp->pos_w );
    }

    memset( fp, 0, sizeof( BFILE ) );
    fp->file = NULL;
    fp->path = NULL;
    fp->buffer_r = NULL;
    fp->buffer_w = NULL;
    fp->pos_byte_r = NULL;
    fp->pos_byte_w = NULL;
    fp->pos_r = NULL;
    fp->pos_w = NULL;

    free( fp );

    return 0;
}

/*--------------------------------------------------------------------------*/

int bfflush

    ( BFILE *stream )

    /* flush a bitstream.
       Upon successful completion 0 is returned.  Otherwise, EOF  is  returned
       and errno is set to indicate the error. */

{
    unsigned char header = 0;       /* new header byte */

    /* verify argument */
    if( stream == NULL ) {
        errno = EBADF;
        return EOF;
    }

    /* verify stream is writable */
    if( stream->mode == MF_CLOSED ||
        stream->mode == MF_R ) {
        errno = EMEDIUMTYPE;
        return EOF;
    }

    /* write current buffer to file, if current bitposition is not zero */
    if( *stream->pos_w != 0 ) {
        /* ensure correct position of FILE pointer */
        if( fseek( stream->file, *stream->pos_byte_w, SEEK_SET ) == -1 )
            return EOF;

        if( fwrite( stream->buffer_w, 1, 1, stream->file ) != 1 ) {
            errno = EIO;
            return EOF;
        }

        /* increase file size, if neccessary */
        if( *stream->pos_byte_w == stream->file_size ) {
            stream->file_size++;

        }

        /* reset header, if current byte position is 0 */
        if( *stream->pos_byte_w == 0 )
            stream->header = *stream->buffer_w;

        /* if reading buffer is set to same byte, update reading buffer */
        if( *stream->pos_byte_r == *stream->pos_byte_w )
            *stream->buffer_r = *stream->buffer_w;

    }

    /* write a new header bit */
    /* write the filling bit being used in the last byte to the
     b eginning of* the file, therefor use 1 bits and put them
     into the header */

    /* set first bit of header */
    header = stream->fill;

    /* shift filling bit to the most left position in new header byte */
    header = header << 7;

    /* cut old heading bit from old header */
    stream->header = stream->header << 1;
    stream->header = stream->header >> 1;

    /* combine old and new header */
    stream->header |= header;

    if( fseek( stream->file, 0, SEEK_SET ) == -1 )
        return EOF;

    if( fwrite( &stream->header, 1, 1, stream->file ) != 1 ) {
        errno = EIO;
        return EOF;
    }

    /* restore position */
    if( fseek( stream->file, *stream->pos_byte_w + 1, SEEK_SET ) == -1 )
        return EOF;

    return 0;
}

/*--------------------------------------------------------------------------*/

int bfinfo

    ( BFILE *fp )   /* pointer to bitstream */

    /* print out information about bitstream */

{

    if( fp == NULL )
        return 1;

    fprintf( stdout, "Bitfile file ptr:\t\t" );
    if( fp->file == NULL )
        fprintf( stdout, "NULL\n" );
    else
        fprintf( stdout, "NOT NULL\n" );

    fprintf( stdout ,"Bitfile mode:\t\t\t" );
    switch( fp->mode ) {

        case MF_CLOSED:
                fprintf( stdout, "closed\n" );
            break;

        case MF_R:
                fprintf( stdout, "r\n" );
            break;

        case MF_RP:
                fprintf( stdout, "r+\n" );
            break;

        case MF_W:
                fprintf( stdout, "w\n" );
            break;

        case MF_WP:
                fprintf( stdout, "w+\n" );
            break;

        case MF_A:
                fprintf( stdout, "a\n" );
            break;

        case MF_AP:
                fprintf( stdout, "a+\n" );
            break;

    }

    fprintf( stdout, "Bitfile path:\t\t\t%s\n", fp->path );

    if( fp->buffer_r != NULL )
        fprintf( stdout, "Bitfile buffer (read):\t\t%u\n", *fp->buffer_r );
    if( fp->buffer_w != NULL )
        fprintf( stdout, "Bitfile buffer (write):\t\t%u\n", *fp->buffer_w );

    fprintf( stdout, "Bitfile byte position (read):\t%ld\n", *fp->pos_byte_r );
    fprintf( stdout, "Bitfile byte position (write):\t%ld\n", *fp->pos_byte_w );

    if( fp->pos_r != NULL )
        fprintf( stdout, "Bitfile position (read):\t%d\n", *fp->pos_r );
    if( fp->pos_w != NULL )
        fprintf( stdout, "Bitfile position (write):\t%d\n", *fp->pos_w );

    fprintf( stdout, "Bitfile valid:\t\t\t%u\n", fp->valid );
    fprintf( stdout, "Bitfile header:\t\t\t%u\n", fp->header );
    fprintf( stdout, "Bitfile fill:\t\t\t%u\n", fp->fill );

    return 0;
}

/*--------------------------------------------------------------------------*/

int bfseek

    ( BFILE *stream,    /* pointer to bitstream */
      long offset,      /* bitposition in bitstream (without header byte) */
      int whence )      /* reference position as passed to fseek() */

    /* set the bitfile reading buffer and position to a specific
       bitposition in bitstream */

{
    long offset_cor = 0;    /* corrected offset (corrected for header bit) */
    long byte = 0;          /* byte position in FILE */
    long bit = 0;           /* bit position in byte */

    /* check for pointer validity */
    if( stream == NULL || stream->file == NULL ) {
        errno = EINVAL;
        return -1;
    }

    /* check for incorrect filemode */
    if( stream->mode == MF_A ) {
        errno = EMEDIUMTYPE;
        return -1;
    }

    /* flush bitfile */
    if( stream->mode != MF_R  ) {
        if( bfflush( stream ) == EOF )
            return -1;
    }

    /* behave according to reference position, i.e. estimate bit- and
       byte-position w.r.t. reference position */
    switch( whence ) {

        case SEEK_SET:
            if( offset < 0 ) {
                errno = EINVAL;
                return -1;
            }

            /* correct offset for header bit */
            offset_cor = offset + 1;

            /* estimate position in bytes */
            byte = offset_cor / 8;

            /* estimate bitnumber asked for in estimated byte */
            bit = offset_cor % 8;

            break;

        case SEEK_CUR:

            /* estimate position in bytes, correct for available bits
               in current byte */
            /* available bits in current byte:
               8 - ( stream->pos + 2 ) = 6 - stream->pos */
            byte = ( offset + *stream->pos_r - 7 ) / 8;

            /* estimate bitposition asked for in estimated byte */
            bit = ( offset + *stream->pos_r ) % 8;

            break;

        case SEEK_END:
            if( offset > 0 ) {
                errno = EINVAL;
                return -1;
            }

            /* estimate byte position, starting from last VALID bit,
             t o read incomplete bytes, subtract 1,                       *
             furthermore convert to offset counting from beginning of file */
            byte = stream->file_size + ( offset + stream->valid ) / 8 - 1;

            /* estimate bitnumber asked for in estimated byte */
            bit = ( offset + stream->valid ) % 8;

            /* set to offset from beginning */
            whence = SEEK_SET;

            break;

        default:
            errno = EINVAL;
            return -1;

    }

    /* set FILE pointer right before estimated byte and load byte to buffer */
    if( fseek( stream->file, byte, whence ) == -1 )
        return -1;

    if( fread( stream->buffer_r, 1, 1, stream->file ) != 1 ) {
        errno = ENODATA;
        return -1;
    }

    /* set current byte position */
    *stream->pos_byte_r = byte;
    *stream->pos_byte_w = byte;

    /* setup buffer usage */
    *stream->pos_r = bit;
    *stream->pos_w = bit;

    return 0;
}

/*--------------------------------------------------------------------------*/

long bftell

    ( BFILE *stream )   /* pointer to bitstream */

    /* tell bitposition in bitstream (not including header bit), referring to
       the reading buffer and the bitposition within */

{
    long retVal = 0;        /* return value, i.e. current bitposition */

    if( stream == NULL ) {
        errno = EINVAL;
        return -1;
    }

    if( stream->pos_r == NULL || stream->pos_byte_r == NULL ) {
        errno = EINVAL;
        return -1;
    }

    retVal = *stream->pos_byte_r * 8 + *stream->pos_r;

    return retVal;
}

/*--------------------------------------------------------------------------*/

int bfputb

    ( const int c,          /* bit value to be written */
      BFILE *stream )       /* pointer to bitstream */

    /* write one bit to the bitstream.
       return the bit value being written on success, EOF else. */

{
    int retVal = 0;             /* return value */
    unsigned char buffer = 0;   /* temporary buffer (remaining bits) */
    unsigned char bit = 0;      /* new bit, written to buffer */
    int i = 0;                  /* loop variable */
    int n_fill = 0;             /* number of filling bit */

    /* check for validity of pointer */
    if( stream == NULL )
        return EOF;

    if( stream->buffer_w == NULL || stream->pos_w == NULL ||
        stream->pos_byte_w == NULL )
        return EOF;

    /* check for correct filemode */
    if( stream->mode == MF_CLOSED ||
        stream->mode == MF_R )
        return EOF;

    /* save unused bits in buffer, i.e. the bits following the current,
       but only do this, if not in last byte, if so fill byte using the new
       filling bit */
    if( *stream->pos_byte_w == stream->file_size &&
        *stream->pos_w == stream->valid ) {
        /* we're in last bitposition of the whole file */

        /* estimate the new filling bit */
        if( c != 0 )
            stream->fill = 0;
        else
            stream->fill = 1;

        /* fill temporary buffer using the filling bit */
        n_fill = 7 - *stream->pos_w;

        for( i = 0; i < n_fill; i++ )
            buffer = ( buffer << 1 ) | stream->fill;

    }
    else {
        /* bitposition within the bitfile */

        /* keep the bits in buffer, following the current bit */
        buffer = *stream->buffer_w << ( *stream->pos_w + 1 );
        buffer >>= ( *stream->pos_w + 1 );

    }

    /* prepare writing buffer, i.e. shift out the remaining bit in current
       byte, i.e. all bits following the current position */
    *stream->buffer_w >>= 8 - *stream->pos_w;
    *stream->buffer_w <<= 8 - *stream->pos_w;

    /* prepare bit, written to buffer */
    if( c != 0 ) {
        bit |= 1;
        bit <<= 7 - *stream->pos_w;
        retVal = 1;
    }

    /* get writing buffer by XORing temporary buffer, bit and buffer */
    *stream->buffer_w |= buffer;
    *stream->buffer_w |= bit;

    /* prepare writing buffer */
    (*stream->pos_w)++;

    /* if writing to the last bit,
       keep the number of valid bit in last = current byte */
    if( *stream->pos_byte_w == stream->file_size &&
        *stream->pos_w > stream->valid ) {
        /* we're in last bitposition of the whole file */
        stream->valid = *stream->pos_w;
    }

    /* flush bitBuffer when reaching position first position out of buffer
       byte */
    if( *stream->pos_w == 8 ) {

        /* ensure FILE pointer being correctly positioned */
        if( fseek( stream->file, *stream->pos_byte_w, SEEK_SET ) == -1 )
            return EOF;

        if( fwrite( stream->buffer_w, 1, 1, stream->file ) != 1 )
            retVal = EOF;

        /* update file size, if writing at the end of file */
        if( *stream->pos_byte_w == stream->file_size )
            stream->file_size++;

        /* keep the first byte for header manipulation */
        if( *stream->pos_byte_w == 0 )
            stream->header = *stream->buffer_w;

        /* if reading buffer is set to same byte, update reading buffer */
        if( *stream->pos_byte_r == *stream->pos_byte_w )
            *stream->buffer_r = *stream->buffer_w;

        /* reset writing buffer, increase byte counter */
        *stream->buffer_w = 0;
        *stream->pos_w = 0;

        (*stream->pos_byte_w)++;

        /* set number of valid bits to zero */
        stream->valid = 0;

        /* estimate the new filling bit */
        if( c != 0 )
            stream->fill = 0;
        else
            stream->fill = 1;

    }

    return retVal;
}

/*--------------------------------------------------------------------------*/

size_t bfwrite

    ( const void *ptr,  /* pointer to bitarray */
      size_t size,      /* size of one element in the bitarray */
      size_t nmemb,     /* number of elements in the bitarray */
      BFILE *stream )   /* pointer to bitstream */

    /* write multiple bits to the bitstream */

{
    size_t retVal = 0;      /* return value, i.e. number of items written */
    size_t i = 0;           /* loop variable */
    void *zeros = NULL;     /* 1-bit-full-of-zeros-vector */

    if( ptr == NULL || stream == NULL )
        return 0;

    zeros = (void*)malloc( size );
    if( zeros == NULL )
        return 0;

    memset( zeros, 0, size );

    for( i = 0; i < nmemb; i++ ) {

        if( memcmp( (size_t*)ptr + i * size, zeros, size ) == 0 ) {
            /* add '0' */
            if( bfputb( 0, stream ) != 0 ) {
                free( zeros );
                return retVal;
            }
        }
        else {
            /* add '1' */
            if( bfputb( 1, stream ) != 1 ) {
                free( zeros );
                return retVal;
            }
        }

        retVal++;

    }

    free( zeros );

    return retVal;
}

/*--------------------------------------------------------------------------*/

size_t bfstrtobf

    ( const char *ptr,  /* pointer to c-string */
      BFILE *stream )   /* pointer to bitstream */

    /* write multiple bits to the bitstream using a c-string as input */

{
    size_t retVal = 0;     /* return value */
    size_t i = 0;          /* loop variable */
    size_t len = 0;        /* string's length */


    if( ptr == 0 || stream == NULL )
        return 0;

    len = strlen( ptr );

    for( i = 0; i < len; i++ ) {

        if( ptr[ i ] == '0' ) {
            if( bfputb( 0, stream ) == EOF )
                return retVal;
        }
        else {
            if( bfputb( 1, stream ) == EOF )
                return retVal;
        }

        retVal++;

    }

    return retVal;
}

/*--------------------------------------------------------------------------*/

int bfgetb

    ( BFILE *stream )   /* pointer to bitstream */

    /* read next bit from the bitstream.
       return the bit value on success, EOF else. */

{
    int retVal = 0;         /* return value */

    /* check for validity of pointer */
    if( stream == NULL )
        return EOF;

    if( stream->file == NULL || stream->buffer_r == NULL ||
        stream->pos_r == NULL || stream->pos_byte_r == NULL )
        return EOF;

    /* check for correct filemode */
    if( stream->mode == MF_CLOSED ||
        stream->mode == MF_W ||
        stream->mode == MF_A )
        return EOF;

    /* reached last bit of bitfile */
    if( *stream->pos_byte_r + 1 == stream->file_size &&
        ( *stream->pos_r == 8 || *stream->pos_r == stream->valid ) )
            return EOF;

    /* read next byte from file */
    if( *stream->pos_r == 8 ) {

        /* ensure FILE pointer is positioned correctly */
        if( fseek( stream->file, *stream->pos_byte_r + 1, SEEK_SET ) == -1 )
            return EOF;

        if( ( *stream->buffer_r = fgetc( stream->file ) ) == EOF )
            return EOF;

        *stream->pos_r = 0;

        (*stream->pos_byte_r)++;

    }

    /* shift the leading bit to the most right position in bitBuffer,
       then bitwise AND will give either 0 or 1 */
    retVal = ( *stream->buffer_r >> ( 7 - *stream->pos_r ) ) & 0x01;

    /* set bitposition to next bit */
    (*stream->pos_r)++;

    return retVal;
}

/*--------------------------------------------------------------------------*/

size_t bfread

    ( void *ptr,        /* pointer to bitarray */
      size_t size,      /* size of one element in the bitarray */
      size_t nmemb,     /* number of elements to read from bitstream */
      BFILE *stream )   /* pointer to bitstream */

    /* read multiple bits from the bitstream */

{
    size_t retVal = 0;              /* return value, i.e. number of items
                                       read */
    size_t i = 0;                   /* loop variables */
    int tmp = 0;                    /* temporary bit-value */
    unsigned char* ptr_u = NULL;    /* pointer to bitarray (casted) */

    if( ptr == NULL || stream == NULL )
        return retVal;

    for( i = 0; i < nmemb; i++ ) {

        if( ( tmp = bfgetb( stream ) ) == EOF )
            return retVal;

        memset( (size_t*)ptr + i * size, 0, size );

        if( tmp == 1 ) {
            ptr_u = (unsigned char*)((size_t*)ptr + i * size);
            *ptr_u = 1;
        }

        retVal++;

    }

    return retVal;
}

/*----------------------------------------------------------------------------*/
void set_byte(BFILE *binary_file, long b) {
  long i;

  if (b>255) {
    printf("Warning: writing value %ld as a byte (> 255)\n",b);
  }

  //printf("setting byte %d\n");

  for (i=0; i< 8; i++) {
    //    printf("%ld ",b%2);
    bfputb(b%2,binary_file);
    b=b/2;
  }
}

/*----------------------------------------------------------------------------*/
long get_byte(BFILE *binary_file) {
  long i;
  long p,b;

  b=0;
  p=1;
  long tmp;
  for (i=0; i<8; i++) {
    tmp = bfgetb(binary_file);
    b=b+tmp*p;
    //    printf("b: %ld, bit %ld p: %ld\n",b,tmp,p);
    p=p*2;
  }
  return b;
}

/*----------------------------------------------------------------------------*/
void set_bits(BFILE *binary_file, long b, long bitdepth) {
  long i;
#ifdef BFIO_WRITE_CHECKS
  long max_value;

  max_value = pow(2,bitdepth);

  if (b>max_value) {
    printf("Warning: writing value %ld as a bit depth (> %ld)\n",b,max_value);
  }
#endif

  //printf("setting byte %d\n");

  for (i=0; i<bitdepth; i++) {
    //    printf("%ld ",b%2);
    bfputb(b%2,binary_file);
    b=b/2;
  }
}

/*----------------------------------------------------------------------------*/
long get_bits(BFILE *binary_file, long bitdepth) {
  long i;
  long p,b;

  b=0;
  p=1;
  long tmp;
  for (i=0; i<bitdepth; i++) {
    tmp = bfgetb(binary_file);
    b=b+tmp*p;
     printf("b: %ld, bit %ld p: %ld\n",b,tmp,p);
    p=p*2;
  }
  return b;
}


/*----------------------------------------------------------------------------*/
void set_bit(BFILE *binary_file, long b) {
#ifdef BFIO_WRITE_CHECKS
  if (b>1) {
    printf("Warning: writing value %ld as a bit (> 1)\n",b);
  }
#endif

  bfputb(b,binary_file);
}

/*----------------------------------------------------------------------------*/
long get_bit(BFILE *binary_file) {
  return bfgetb(binary_file);
}
