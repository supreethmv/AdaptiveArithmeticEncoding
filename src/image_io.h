#ifndef IMAGE_IO_H_
#define IMAGE_IO_H_

/*--------------------------------------------------------------------------*/

void read_pgm_header
(const char  *file_name,   /* name of pgm file */
 long        *nx,          /* image size in x direction, output */
 long        *ny);         /* image size in y direction, output */

/*--------------------------------------------------------------------------*/

void read_pgm_and_allocate_memory
(const char  *file_name,    /* name of pgm file */
 long        *nx,           /* image size in x direction, output */
 long        *ny,           /* image size in y direction, output */
 long        ***u);         /* image, output */

/*
  reads a greyscale image that has been encoded in pgm format P5;
  allocates memory for the image u;
  adds boundary layers of size 1 such that
  - the relevant image pixels in x direction use the indices 1,...,nx
  - the relevant image pixels in y direction use the indices 1,...,ny
*/

/*--------------------------------------------------------------------------*/

void read_ppm_and_allocate_memory

(const char  *file_name,    /* name of pgm file */
 long        *nx,           /* image size in x direction, output */
 long        *ny,           /* image size in y direction, output */
 long       ****u);         /* image, output */

/*
  reads a colour image that has been encoded in pgm format P6;
  allocates memory for the image u;
  adds boundary layers of size 1 such that
  - the relevant image pixels in x direction use the indices 1,...,nx
  - the relevant image pixels in y direction use the indices 1,...,ny
*/

/*--------------------------------------------------------------------------*/

void comment_line

(char* comment,       /* comment string (output) */
 char* lineformat,    /* format string for comment line */
 ...);                /* optional arguments */

/*
  Add a line to the comment string comment. The string line can contain plain
  text and format characters that are compatible with sprintf.
  Example call: print_comment_line(comment,"Text %f %d",double_var,int_var);
  If no line break is supplied at the end of the input string, it is added
  automatically.
*/

/*--------------------------------------------------------------------------*/

void write_pgm

(long  **u,           /* image, unchanged */
 long   nx,           /* image size in x direction */
 long   ny,           /* image size in y direction */
 char   *file_name,   /* name of pgm file */
 char   *comments);   /* comment string (set 0 for no comments) */

/*
  writes a greyscale image into a pgm P5 file;
*/

/*--------------------------------------------------------------------------*/

void write_mask

(long **u,         /* image, unchanged */
 long   nx,           /* image size in x direction */
 long   ny,           /* image size in y direction */
 char   *file_name,   /* name of pgm file */
 char   *comments);   /* comment string (set 0 for no comments) */

/*
  writes a greyscale image into a pgm P5 file;
*/

/*--------------------------------------------------------------------------*/

void write_ppm

(long ***u,           /* image, unchanged */
 long   nx,           /* image size in x direction */
 long   ny,           /* image size in y direction */
 char   *file_name,   /* name of pgm file */
 char   *comments);   /* comment string (set 0 for no comments) */

/*
  writes a greyscale image into a pgm P5 file;
*/


#endif /* IMAGE_IO_H_ */
