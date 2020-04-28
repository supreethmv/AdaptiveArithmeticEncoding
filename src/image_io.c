#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "alloc.h"
#include "image_io.h"

/*--------------------------------------------------------------------------*/

void read_pgm_and_allocate_memory
(const char  *file_name,   /* name of pgm file */
 long        *nx,          /* image size in x direction, output */
 long        *ny,          /* image size in y direction, output */
 long        ***u)         /* image, output */

/*
  reads a greyscale image that has been encoded in pgm format P5;
  allocates memory for the image u;
  adds boundary layers of size 1 such that
  - the relevant image pixels in x direction use the indices 1,...,nx
  - the relevant image pixels in y direction use the indices 1,...,ny
*/

{
  FILE   *inimage;    /* input file */
  char   row[80];     /* for reading data */
  long   i, j;        /* loop variables */

  /* open file */
  inimage = fopen (file_name, "rb");
  if (NULL == inimage)
    {
      printf ("could not open file '%s' for reading, aborting.\n", file_name);
      exit (1);
    }

  /* read header */
  fgets (row, 80, inimage);          /* skip format definition */
  fgets (row, 80, inimage);
  while (row[0]=='#')                /* skip comments */
    fgets (row, 80, inimage);
  sscanf (row, "%ld %ld", nx, ny);   /* read image size */
  fgets (row, 80, inimage);          /* read maximum grey value */

  /* allocate memory if necessary */
  if (*u==0)
    alloc_long_matrix (u, (*nx)+2, (*ny)+2);

  /* read image data row by row */
  for (j=1; j<=(*ny); j++)
    for (i=1; i<=(*nx); i++)
      (*u)[i][j] = (long) getc(inimage);

  /* close file */
  fclose(inimage);

  return;

} /* read_pgm_and_allocate_memory */

/*--------------------------------------------------------------------------*/

void read_pgm_header
(const char  *file_name,   /* name of pgm file */
 long        *nx,          /* image size in x direction, output */
 long        *ny)          /* image size in y direction, output */

/*
  reads a greyscale image that has been encoded in pgm format P5;
*/

{
  FILE   *inimage;    /* input file */
  char   row[80];     /* for reading data */

  /* open file */
  inimage = fopen (file_name, "rb");
  if (NULL == inimage)
    {
      printf ("could not open file '%s' for reading, aborting.\n", file_name);
      exit (1);
    }

  /* read header */
  fgets (row, 80, inimage);          /* skip format definition */
  fgets (row, 80, inimage);
  while (row[0]=='#')                /* skip comments */
    fgets (row, 80, inimage);
  sscanf (row, "%ld %ld", nx, ny);   /* read image size */

  /* close file */
  fclose(inimage);

  return;

} /* read_pgm_and_allocate_memory */


/*--------------------------------------------------------------------------*/

void read_ppm_and_allocate_memory

(const char  *file_name,    /* name of pgm file */
 long        *nx,           /* image size in x direction, output */
 long        *ny,           /* image size in y direction, output */
 long       ****u)          /* image, output */

/*
  reads a colour image that has been encoded in pgm format P6;
  allocates memory for the image u;
  adds boundary layers of size 1 such that
  - the relevant image pixels in x direction use the indices 1,...,nx
  - the relevant image pixels in y direction use the indices 1,...,ny
*/

{
  FILE   *inimage;    /* input file */
  char   row[80];     /* for reading data */
  long   i, j, m;     /* loop variables */

  /* open file */
  inimage = fopen (file_name, "rb");
  if (NULL == inimage)
    {
      printf ("could not open file '%s' for reading, aborting.\n", file_name);
      exit (1);
    }

  /* read header */
  fgets (row, 80, inimage);          /* skip format definition */
  fgets (row, 80, inimage);
  while (row[0]=='#')                /* skip comments */
    fgets (row, 80, inimage);
  sscanf (row, "%ld %ld", nx, ny);   /* read image size */
  fgets (row, 80, inimage);          /* read maximum grey value */

  /* allocate memory */
  if (*u==0)
    alloc_long_cubix (u, 3, (*nx)+2, (*ny)+2);

  /* read image data row by row */
  for (j=1; j<=(*ny); j++)
    for (i=1; i<=(*nx); i++)
      for (m=0; m<3; m++)
        (*u)[m][i][j] = (long) getc(inimage);

  /* close file */
  fclose(inimage);

  return;

} /* read_ppm_and_allocate_memory */

/*--------------------------------------------------------------------------*/

void comment_line

(char* comment,       /* comment string (output) */
 char* lineformat,    /* format string for comment line */
 ...)                 /* optional arguments */

/*
  Add a line to the comment string comment. The string line can contain plain
  text and format characters that are compatible with sprintf.
  Example call: print_comment_line(comment,"Text %f %d",double_var,int_var);
  If no line break is supplied at the end of the input string, it is added
  automatically.
*/

{
  char     line[80];
  va_list  arguments;

  /* get list of optional function arguments */
  va_start (arguments, lineformat);

  /* convert format string and arguments to plain text line string */
  vsprintf (line, lineformat, arguments);

  /* add line to total commentary string */
  strncat (comment, line, 80);

  /* add line break if input string does not end with one */
  if (line[strlen(line)-1] != '\n')
    sprintf (comment, "%s\n", comment);

  /* close argument list */
  va_end (arguments);

  return;

} /* comment_line */

/*--------------------------------------------------------------------------*/

void write_pgm

(long  **u,           /* image, unchanged */
 long   nx,           /* image size in x direction */
 long   ny,           /* image size in y direction */
 char   *file_name,   /* name of pgm file */
 char   *comments)    /* comment string (set 0 for no comments) */

/*
  writes a greyscale image into a pgm P5 file;
*/

{
  FILE           *outimage;  /* output file */
  long           i, j;       /* loop variables */
  double          aux;        /* auxiliary variable */
  unsigned char  byte;       /* for data conversion */

  /* open file */
  outimage = fopen (file_name, "wb");
  if (NULL == outimage)
    {
      printf("Could not open file '%s' for writing, aborting\n", file_name);
      exit(1);
    }

  /* write header */
  fprintf (outimage, "P5\n");                  /* format */
  if (comments != 0)
    fprintf (outimage, comments);             /* comments */
  fprintf (outimage, "%ld %ld\n", nx, ny);     /* image size */
  fprintf (outimage, "255\n");                 /* maximal value */

  /* write image data */
  for (j=1; j<=ny; j++)
    for (i=1; i<=nx; i++)
      {
        aux = (double)u[i][j] + 0.499999;    /* for correct rounding */
        if (aux < 0.0)
          byte = (unsigned char)(0.0);
        else if (aux > 255.0)
          byte = (unsigned char)(255.0);
        else
          byte = (unsigned char)(aux);
        fwrite (&byte, sizeof(unsigned char), 1, outimage);
      }

  /* close file */
  fclose (outimage);

  return;

} /* write_pgm */

/*--------------------------------------------------------------------------*/

void write_mask

(long  **u,          /* image, unchanged */
 long   nx,           /* image size in x direction */
 long   ny,           /* image size in y direction */
 char   *file_name,   /* name of pgm file */
 char   *comments)    /* comment string (set 0 for no comments) */

/*
  writes a greyscale image into a pgm P5 file;
*/

{
  FILE           *outimage;  /* output file */
  long           i, j;       /* loop variables */
  unsigned char  byte;       /* for data conversion */

  /* open file */
  outimage = fopen (file_name, "wb");
  if (NULL == outimage)
    {
      printf("Could not open file '%s' for writing, aborting\n", file_name);
      exit(1);
    }

  /* write header */
  fprintf (outimage, "P5\n");                  /* format */
  if (comments != 0)
    fprintf (outimage, comments);             /* comments */
  fprintf (outimage, "%ld %ld\n", nx, ny);     /* image size */
  fprintf (outimage, "255\n");                 /* maximal value */

  /* write image data */
  for (j=1; j<=ny; j++)
    for (i=1; i<=nx; i++)
      {
        if (u[i][j] < 0.5)
          byte = (unsigned char)(255.0);
        else
          byte = (unsigned char)(0.0);
        fwrite (&byte, sizeof(unsigned char), 1, outimage);
      }

  /* close file */
  fclose (outimage);

  return;

} /* write_pgm */

/*--------------------------------------------------------------------------*/

void write_ppm

(long   ***u,         /* image, unchanged */
 long   nx,           /* image size in x direction */
 long   ny,           /* image size in y direction */
 char   *file_name,   /* name of pgm file */
 char   *comments)    /* comment string (set 0 for no comments) */

/*
  writes a greyscale image into a pgm P5 file;
*/

{
  FILE           *outimage;  /* output file */
  long           i, j, m;       /* loop variables */
  double          aux;        /* auxiliary variable */
  unsigned char  byte;       /* for data conversion */

  /* open file */
  outimage = fopen (file_name, "wb");
  if (NULL == outimage)
    {
      printf("Could not open file '%s' for writing, aborting\n", file_name);
      exit(1);
    }

  /* write header */
  fprintf (outimage, "P6\n");                  /* format */
  if (comments != 0)
    fprintf (outimage, comments);             /* comments */
  fprintf (outimage, "%ld %ld\n", nx, ny);     /* image size */
  fprintf (outimage, "255\n");                 /* maximal value */

  /* write image data */
  for (j=1; j<=ny; j++)
    for (i=1; i<=nx; i++)
      for (m=0; m<3; m++)
        {
          aux = (double)u[m][i][j] + 0.499999;    /* for correct rounding */
          if (aux < 0.0)
            byte = (unsigned char)(0.0);
          else if (aux > 255.0)
            byte = (unsigned char)(255.0);
          else
            byte = (unsigned char)(aux);
          fwrite (&byte, sizeof(unsigned char), 1, outimage);
        }

  /* close file */
  fclose (outimage);

  return;

} /* write_ppm */
