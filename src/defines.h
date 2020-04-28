#ifndef MIAPDE_DEFINES_H
#define MIAPDE_DEFINES_H

/*--------------------------------------------------------------------------*/
/*                                                                          */
/*  File:       defines.h                                                   */
/*                                                                          */
/*  Purpose:    Definitions e.g. of datatype being used throughout the      */
/*              whole library.                                              */
/*                                                                          */
/*  Author:     Leif Bergerhoff, bergerhoff@mia.uni-saarland.de             */
/*  Date:       April 13, 2015                                              */
/*                                                                          */
/*--------------------------------------------------------------------------*/

/* definition of access modes, cp. eg. man fopen() */
typedef enum {
    MF_CLOSED = 0,  /* closed */
    MF_R = 1,       /* reading, eqivalent to 'r' */
    MF_RP = 2,      /* eqivalent to 'r+' */
    MF_W = 3,       /* writing, equivalent to 'w' */
    MF_WP = 4,      /* eqivalent to 'w+' */
    MF_A = 5,       /* appending (writing at the end of file), equiv. to 'a' */
    MF_AP = 6,      /* eqivalent to 'a+' */
} MIAPDEFILE_MODE;

/* definition of the datatype MFILE_MODE */
typedef MIAPDEFILE_MODE MFILE_MODE;

#endif
