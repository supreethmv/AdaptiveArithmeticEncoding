#ifndef ALLOC_H_
#define ALLOC_H_

/*--------------------------------------------------------------------------*/
void alloc_vector

(double **vector,   /* vector */
		long  n1);         /* size */
/* allocates memory for a vector of size n1 */

/*--------------------------------------------------------------------------*/

void alloc_matrix

(double ***matrix,  /* matrix */
		long  n1,         /* size in direction 1 */
		long  n2);         /* size in direction 2 */

/* allocates memory for matrix of size n1 * n2 */

/*--------------------------------------------------------------------------*/

void alloc_cubix

(double ****cubix,  /* cubix */
		long  n1,         /* size in direction 1 */
		long  n2,         /* size in direction 2 */
		long  n3);         /* size in direction 3 */

/* allocates memory for cubix of size n1 * n2 * n3 */

/*--------------------------------------------------------------------------*/

void alloc_quadrix

(double *****quadrix,  /* quadrix */
		long  n1,            /* size in direction 1 */
		long  n2,            /* size in direction 2 */
		long  n3,            /* size in direction 3 */
		long  n4);            /* size in direction 4 */

/* allocates memory for quadrix of size n1 * n2 * n3 * n4 */

/*--------------------------------------------------------------------------*/

void disalloc_vector

(double *vector,    /* vector */
		long  n1);         /* size */

/* disallocates memory for a vector of size n1 */

/*--------------------------------------------------------------------------*/

void disalloc_matrix

(double **matrix,   /* matrix */
		long  n1,         /* size in direction 1 */
		long  n2);         /* size in direction 2 */

/* disallocates memory for matrix of size n1 * n2 */

/*--------------------------------------------------------------------------*/

void disalloc_cubix

(double ***cubix,   /* cubix */
		long  n1,         /* size in direction 1 */
		long  n2,         /* size in direction 2 */
		long  n3);         /* size in direction 3 */

/* disallocates memory for cubix of size n1 * n2 * n3 */


/*--------------------------------------------------------------------------*/

void disalloc_quadrix

(double ****quadrix,   /* quadrix */
		long  n1,            /* size in direction 1 */
		long  n2,            /* size in direction 2 */
		long  n3,            /* size in direction 3 */
		long  n4);            /* size in direction 4 */

/* disallocates memory for cubix of size n1 * n2 * n3 * n4 */


/*--------------------------------------------------------------------------*/

void alloc_string

     (char **vector,  /* vector */
      long  n1);      /* size */

/*
 allocates memory for a string of size n1
*/

/*--------------------------------------------------------------------------*/

void alloc_double_vector

     (double **vector,  /* vector */
      long  n1);        /* size */

/*
 allocates memory for a vector of size n1
*/

/*--------------------------------------------------------------------------*/

void alloc_long_vector

     (long **vector,  /* vector */
      long  n1);      /* size */

/*
 allocates memory for a vector of size n1
*/

/*--------------------------------------------------------------------------*/

void alloc_double_matrix

     (double ***matrix,  /* matrix */
      long  n1,          /* size in direction 1 */
      long  n2);         /* size in direction 2 */

/*
 allocates memory for a matrix of size n1 * n2
*/


/*--------------------------------------------------------------------------*/

void alloc_long_matrix

     (long ***matrix,  /* matrix */
      long  n1,        /* size in direction 1 */
      long  n2);       /* size in direction 2 */

/*
 allocates memory for a matrix of size n1 * n2
*/

/*--------------------------------------------------------------------------*/

void alloc_long_cubix

     (long ****cubix,  /* cubix */
      long  n1,        /* size in direction 1 */
      long  n2,        /* size in direction 2 */
      long  n3);       /* size in direction 3 */

/*
 allocates memory for cubix of size n1 * n2 * n3
*/

/*--------------------------------------------------------------------------*/

void disalloc_string

     (char *vector,  /* vector */
      long  n1);     /* size */

/*
 frees memory for a string of size n1
*/

/*--------------------------------------------------------------------------*/

void disalloc_double_vector

     (double *vector,  /* vector */
      long  n1);       /* size */

/*
 frees memory for a vector of size n1
*/

/*--------------------------------------------------------------------------*/

void disalloc_long_vector

     (long *vector,  /* vector */
      long  n1);     /* size */

/*
 frees memory for a vector of size n1
*/

/*--------------------------------------------------------------------------*/

void disalloc_double_matrix

     (double **matrix,  /* matrix */
      long  n1,         /* size in direction 1 */
      long  n2);        /* size in direction 2 */

/*
 frees memory for a matrix of size n1 * n2
*/

/*--------------------------------------------------------------------------*/

void disalloc_long_matrix

     (long **matrix,  /* matrix */
      long  n1,       /* size in direction 1 */
      long  n2);      /* size in direction 2 */

/*
 frees memory for a matrix of size n1 * n2
*/

/*--------------------------------------------------------------------------*/

void disalloc_long_cubix

     (long ***cubix,  /* cubix */
      long  n1,       /* size in direction 1 */
      long  n2,       /* size in direction 2 */
      long  n3);      /* size in direction 3 */

/*
 frees memory for a cubix of size n1 * n2 * n3
*/

#endif
