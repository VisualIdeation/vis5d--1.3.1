/* tokenize.c */

/*
 * Vis5D system for visualizing five dimensional gridded data sets.
 * Copyright (C) 1990 - 2000 Bill Hibbard, Johan Kellum, Brian Paul,
 * Dave Santek, and Andre Battaiola.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * As a special exception to the terms of the GNU General Public
 * License, you are permitted to link Vis5D with (and distribute the
 * resulting source and executables) the LUI library (copyright by
 * Stellar Computer Inc. and licensed for distribution with Vis5D),
 * the McIDAS library, and/or the NetCDF library, where those
 * libraries are governed by the terms of their own licenses.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include "../config.h"


/*
 * Divide an input string into a sequence of tokens similar to how argc and
 * argv pass arguments to main().
 */



#include <stdlib.h>
#include <string.h>
#include "misc_i.h"
/* 31Jan06  Phil McDonald */
#include "memory.h"
/* end PM */



#define MAX_TOKENS  1000
#define MAX_TOKEN_LENGTH 1000



char **tokenize( char *str, int *ntokens )
{
   char temp[MAX_TOKEN_LENGTH];
   char **tokens;
   int i, j, n;


   /* allocate an array [MAX_TOKENS] of pointer to char */
/* 31Jan06  Phil McDonald */
   tokens = (char **) TMP_CALLOC (MAX_TOKENS, sizeof (char *));
/* end PM */
   if (!tokens) {
      return NULL;
   }

   i = j = n = 0;
   for (i=0;;i++) {
      if (str[i]==' ' || str[i]=='\t' || str[i]=='\n' || str[i]==0) {
         if (j>0) {
            /* end of a token */
            temp[j] = 0;
            if (n>=MAX_TOKENS) break;  /*done*/
/* 31Jan06  Phil McDonald */
            tokens[n++] = TMP_STRDUP(temp);
            j = 0;
         }
         if (str[i]==0) break;  /*all done*/
      }
      else {
         if (j<MAX_TOKEN_LENGTH) {
            temp[j] = str[i];
            j++;
         }
      }
   }

   *ntokens = n;

   return tokens;
}



void free_tokens( char **tokens )
{
   int i;

   for (i=0;tokens[i] && i<MAX_TOKENS;i++) {
/* 31Jan06  Phil McDonald */
      TMP_FREE (tokens[i]);
   }

   TMP_FREE (tokens);
/* end PM */
}


