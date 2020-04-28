#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "alloc.h"

/*--------------------------------------------------------------------------*/
void alloc_vector

(double **vector,   /* vector */
 long  n1)         /* size */
/* allocates memory for a vector of size n1 */
{
  *vector = (double *) malloc (n1 * sizeof(double));
  if (*vector == NULL)
    {
      printf("alloc_vector: not enough memory available\n");
      exit(1);
    }
  return;
}

/*--------------------------------------------------------------------------*/

void alloc_matrix

(double ***matrix,  /* matrix */
 long  n1,         /* size in direction 1 */
 long  n2)         /* size in direction 2 */

/* allocates memory for matrix of size n1 * n2 */


{
  long i;

  *matrix = (double **) malloc (n1 * sizeof(double *));
  if (*matrix == NULL)
    {
      printf("alloc_matrix: not enough memory available\n");
      exit(1);
    }
  for (i=0; i<n1; i++)
    {
      (*matrix)[i] = (double *) malloc (n2 * sizeof(double));
      if ((*matrix)[i] == NULL)
        {
          printf("alloc_matrix: not enough memory available\n");
          exit(1);
        }
    }
  return;
}

/*--------------------------------------------------------------------------*/

void alloc_cubix

(double ****cubix,  /* cubix */
 long  n1,         /* size in direction 1 */
 long  n2,         /* size in direction 2 */
 long  n3)         /* size in direction 3 */

/* allocates memory for cubix of size n1 * n2 * n3 */


{
  long i, j;

  *cubix = (double ***) malloc (n1 * sizeof(double **));
  if (*cubix == NULL)
    {
      printf("alloc_cubix: not enough memory available\n");
      exit(1);
    }
  for (i=0; i<n1; i++)
    {
      (*cubix)[i] = (double **) malloc (n2 * sizeof(double *));
      if ((*cubix)[i] == NULL)
        {
          printf("alloc_cubix: not enough memory available\n");
          exit(1);
        }
      for (j=0; j<n2; j++)
        {
          (*cubix)[i][j] = (double *) malloc (n3 * sizeof(double));
          if ((*cubix)[i][j] == NULL)
            {
              printf("alloc_cubix: not enough memory available\n");
              exit(1);
            }
        }
    }
  return;
}

/*--------------------------------------------------------------------------*/

void alloc_quadrix

(double *****quadrix,  /* quadrix */
 long  n1,            /* size in direction 1 */
 long  n2,            /* size in direction 2 */
 long  n3,            /* size in direction 3 */
 long  n4)            /* size in direction 4 */

/* allocates memory for quadrix of size n1 * n2 * n3 * n4 */


{
  long i, j, k;

  *quadrix = (double ****) malloc (n1 * sizeof(double ***));
  if (*quadrix == NULL)
    {
      printf("alloc_quadrix: not enough memory available\n");
      exit(1);
    }
  for (i=0; i<n1; i++)
    {
      (*quadrix)[i] = (double ***) malloc (n2 * sizeof(double **));
      if ((*quadrix)[i] == NULL)
        {
          printf("alloc_quadrix: not enough memory available\n");
          exit(1);
        }
      for (j=0; j<n2; j++)
        {
          (*quadrix)[i][j] = (double **) malloc (n3 * sizeof(double *));
          if ((*quadrix)[i][j] == NULL)
            {
              printf("alloc_quadrix: not enough memory available\n");
              exit(1);
            }
          for (k=0; k<n3; k++)
            {
              (*quadrix)[i][j][k] = (double *) malloc (n4 * sizeof(double));
              if ((*quadrix)[i][j][k] == NULL)
                {
                  printf("alloc_quadrix: not enough memory available\n");
                  exit(1);
                }
            }
        }
    }
  return;
}

/*--------------------------------------------------------------------------*/

void disalloc_vector

(double *vector,    /* vector */
 long  n1)         /* size */

/* disallocates memory for a vector of size n1 */

{
  free(vector);
  return;
}

/*--------------------------------------------------------------------------*/

void disalloc_matrix

(double **matrix,   /* matrix */
 long  n1,         /* size in direction 1 */
 long  n2)         /* size in direction 2 */

/* disallocates memory for matrix of size n1 * n2 */

{
  long i;

  for (i=0; i<n1; i++)
    free(matrix[i]);

  free(matrix);

  return;
}

/*--------------------------------------------------------------------------*/

void disalloc_cubix

(double ***cubix,   /* cubix */
 long  n1,         /* size in direction 1 */
 long  n2,         /* size in direction 2 */
 long  n3)         /* size in direction 3 */

/* disallocates memory for cubix of size n1 * n2 * n3 */

{
  long i, j;

  for (i=0; i<n1; i++)
    for (j=0; j<n2; j++)
      free(cubix[i][j]);

  for (i=0; i<n1; i++)
    free(cubix[i]);

  free(cubix);

  return;
}

/*--------------------------------------------------------------------------*/

void disalloc_quadrix

(double ****quadrix,   /* quadrix */
 long  n1,            /* size in direction 1 */
 long  n2,            /* size in direction 2 */
 long  n3,            /* size in direction 3 */
 long  n4)            /* size in direction 4 */

/* disallocates memory for cubix of size n1 * n2 * n3 * n4 */

{
  long i, j, k;

  for (i=0; i<n1; i++)
    for (j=0; j<n2; j++)
      for (k=0; k<n3; k++)
        free(quadrix[i][j][k]);

  for (i=0; i<n1; i++)
    for (j=0; j<n2; j++)
      free(quadrix[i][j]);

  for (i=0; i<n1; i++)
    free(quadrix[i]);

  free(quadrix);

  return;
}


void alloc_string

(char **vector,  /* vector */
 long  n1)       /* size */

/*
  allocates memory for a string of size n1
*/

{
  *vector = (char*) malloc (n1 * sizeof(char));

  if (*vector == NULL)
    {
      printf("alloc_string: not enough memory available\n");
      exit(1);
    }

  return;

}  /* alloc_string */

/*--------------------------------------------------------------------------*/

void alloc_double_vector

(double **vector,  /* vector */
 long  n1)         /* size */

/*
  allocates memory for a vector of size n1
*/

{
  *vector = (double *) malloc (n1 * sizeof(double));

  if (*vector == NULL)
    {
      printf("alloc_double_vector: not enough memory available\n");
      exit(1);
    }

  return;

}  /* alloc_double_vector */

/*--------------------------------------------------------------------------*/

void alloc_long_vector

(long **vector,  /* vector */
 long  n1)       /* size */

/*
  allocates memory for a vector of size n1
*/

{
  *vector = (long *) malloc (n1 * sizeof(long));

  if (*vector == NULL)
    {
      printf("alloc_long_vector: not enough memory available\n");
      exit(1);
    }

  return;

}  /* alloc_long_vector */

/*--------------------------------------------------------------------------*/

void alloc_double_matrix

(double ***matrix,  /* matrix */
 long  n1,          /* size in direction 1 */
 long  n2)          /* size in direction 2 */

/*
  allocates memory for a matrix of size n1 * n2
*/

{
  long i;    /* loop variable */

  *matrix = (double **) malloc (n1 * sizeof(double *));

  if (*matrix == NULL)
    {
      printf("alloc_double_matrix: not enough memory available\n");
      exit(1);
    }

  for (i=0; i<n1; i++)
    {
      (*matrix)[i] = (double *) malloc (n2 * sizeof(double));
      if ((*matrix)[i] == NULL)
        {
          printf("alloc_double_matrix: not enough memory available\n");
          exit(1);
        }
    }

  return;

}  /* alloc_double_matrix */

/*--------------------------------------------------------------------------*/

void alloc_long_matrix

(long ***matrix,  /* matrix */
 long  n1,        /* size in direction 1 */
 long  n2)        /* size in direction 2 */

/*
  allocates memory for a matrix of size n1 * n2
*/

{
  long i;    /* loop variable */

  *matrix = (long **) malloc (n1 * sizeof(long *));

  if (*matrix == NULL)
    {
      printf("alloc_long_matrix: not enough memory available\n");
      exit(1);
    }

  for (i=0; i<n1; i++)
    {
      (*matrix)[i] = (long *) malloc (n2 * sizeof(long));
      if ((*matrix)[i] == NULL)
        {
          printf("alloc_long_matrix: not enough memory available\n");
          exit(1);
        }
    }

  return;

}  /* alloc_long_matrix */

/*--------------------------------------------------------------------------*/

void alloc_long_cubix

(long ****cubix,  /* cubix */
 long  n1,        /* size in direction 1 */
 long  n2,        /* size in direction 2 */
 long  n3)        /* size in direction 3 */

/*
  allocates memory for cubix of size n1 * n2 * n3
*/

{
  long i, j;  /* loop variables */

  *cubix = (long ***) malloc (n1 * sizeof(long **));

  if (*cubix == NULL)
    {
      printf("alloc_long_cubix: not enough memory available\n");
      exit(1);
    }

  for (i=0; i<n1; i++)
    {
      (*cubix)[i] = (long **) malloc (n2 * sizeof(long *));
      if ((*cubix)[i] == NULL)
        {
          printf("alloc_long_cubix: not enough memory available\n");
          exit(1);
        }
      for (j=0; j<n2; j++)
        {
          (*cubix)[i][j] = (long *) malloc (n3 * sizeof(long));
          if ((*cubix)[i][j] == NULL)
            {
              printf("alloc_long_cubix: not enough memory available\n");
              exit(1);
            }
        }
    }

  return;

}  /* alloc_long_cubix */

/*--------------------------------------------------------------------------*/

void disalloc_string

(char *vector,  /* vector */
 long  n1)      /* size */

/*
  frees memory for a string of size n1
*/

{

  free(vector);
  return;

}  /* disalloc_string */

/*--------------------------------------------------------------------------*/

void disalloc_double_vector

(double *vector,  /* vector */
 long  n1)        /* size */

/*
  frees memory for a vector of size n1
*/

{

  free(vector);
  return;

}  /* disalloc_double_vector */

/*--------------------------------------------------------------------------*/

void disalloc_long_vector

(long *vector,     /* vector */
 long  n1)         /* size */

/*
  frees memory for a vector of size n1
*/

{

  free(vector);
  return;

}  /* disalloc_long_vector */

/*--------------------------------------------------------------------------*/

void disalloc_double_matrix

(double **matrix,  /* matrix */
 long  n1,         /* size in direction 1 */
 long  n2)         /* size in direction 2 */

/*
  frees memory for a matrix of size n1 * n2
*/

{
  long i;   /* loop variable */

  for (i=0; i<n1; i++)
    free(matrix[i]);

  free(matrix);

  return;

}  /* disalloc_double_matrix */

/*--------------------------------------------------------------------------*/

void disalloc_long_matrix

(long **matrix,    /* matrix */
 long  n1,         /* size in direction 1 */
 long  n2)         /* size in direction 2 */

/*
  frees memory for a matrix of size n1 * n2
*/

{
  long i;   /* loop variable */

  for (i=0; i<n1; i++)
    free(matrix[i]);

  free(matrix);

  return;

}  /* disalloc_long_matrix */

/*--------------------------------------------------------------------------*/

void disalloc_long_cubix

(long ***cubix,  /* cubix */
 long  n1,       /* size in direction 1 */
 long  n2,       /* size in direction 2 */
 long  n3)       /* size in direction 3 */

/*
  frees memory for a cubix of size n1 * n2 * n3
*/

{
  long i, j;   /* loop variables */

  for (i=0; i<n1; i++)
    for (j=0; j<n2; j++)
      free(cubix[i][j]);

  for (i=0; i<n1; i++)
    free(cubix[i]);

  free(cubix);

  return;

}  /* disalloc_long_cubix */




