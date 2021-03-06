/* read_grid.c */

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
 * Functions for reading McIDAS GRID files
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "binio.h"
#include "file_i.h"
#include "grid_i.h"
#include "memory.h"
#include "misc_i.h"
#include "proj_i.h"
#include "projlist_i.h"
#include "read_grid_i.h"
#include "v5d.h"


int Debug_i;


/* Meters per degree: */
#define M_PER_DEG  111137.0


/* ADDED BY WLH 5-8-95 */


int uc_i[1000];
int neguc_i[200];
int ttyfd_i[10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

#ifdef HAVE_MCIDAS

#define IGGMAX F77_FUNC(iggmax,IGGMAX)
#define IGGET F77_FUNC(igget,IGGET)
extern int IGGMAX ( int *, int *);
extern int IGGET ( int *, int *, int *, void *, int *, int *, int * );

#ifndef MCIDAS_SIDECAR

#define KLTWIN F77_FUNC(kltwin,KLTWIN)
#define LWINIT F77_FUNC(lwinit,LWINIT)
extern void KLTWIN( void );
extern void LWINIT( void );


static void init_mcidas_kludge( void )
{
  int i;

  for (i=0; i<1000; i++) uc_i[i] = 0;
  for (i=0; i<200; i++) neguc_i[i] = 0;
  KLTWIN();
  LWINIT();
}
#endif /* end of ifndef MCIDAS_SIDECAR */
#endif /* end of ifdef HAVE_MCIDAS */



/*
 * Scan the named file, createing a grid_info struct for each grid.  Store
 * the grid_info structs in the grid data base.
 * Input:  name - name of McIDAS GR3D file.
 *         db - the grid data base
 * Return:  number of grids found in the file
 */
int get_grid_info( char *name, struct grid_db *db )
{
#ifdef HAVE_MCIDAS
   struct grid_info *info;
   int i;
   int grids = 0;
   int len, filenumber;
   int maxgrids, start;
   float minheight, maxheight;
   float args[20];

#ifndef MCIDAS_SIDECAR
   init_mcidas_kludge();
#endif

   /* Extract file number */
   len = strlen(name);
   filenumber = atoi( name + len - 4 );


   /* Call IGGMAX to determine how many grids are in the file */
   maxgrids = IGGMAX( &filenumber, &start );

   if (maxgrids==0) {
      printf("Error in IGGMAX\n");
      return 0;
   }

   /*printf("maxgrids=%d\n", maxgrids);*/

   minheight = 10000.0;
   maxheight = -10000.0;

   grids = 0;

   for (i=1;i<=maxgrids;i++) {
      int k, max=MAXROWS*MAXCOLUMNS, nr, nc, header[64];
      float data[MAXROWS*MAXCOLUMNS];
      int vertlevel, units;
      float height;
      int date, time, validtime, seconds, days;

      /* Call McIDAS function to get 2-D grid info */
      k = IGGET( &filenumber, &i, &max, data, &nr, &nc, header );
      if (k==-1) {
         /* bad grid number (i) */
         continue;
      }
      else if (k==-2) {
         /* bad file number */
         break;
      }
      /* k==0 means OK */

      info = alloc_grid_info();
/* 01Feb06  Phil McDonald */
      info->FileName = TMP_STRDUP (name);
/* end PM */
      info->Format = FILE_GRID;
      info->FileNumber = filenumber;
      info->Position = i;

      info->Nr = nr;
      info->Nc = nc;
      info->Nl = 1;

      /* Date and Time */
      date = header[3];
      time = header[4];
      validtime = header[5];
      /* add validtime to date,time to get actual date,time */
      seconds = v5dHHMMSStoSeconds(time) + validtime * 60 * 60;
      days = v5dYYDDDtoDays(date) + seconds / (24*60*60);
      seconds = seconds % (24*60*60);
      info->DateStamp = v5dDaysToYYDDD( days );
      info->TimeStamp = v5dSecondsToHHMMSS( seconds );

      { /* Get variable name, remove trailing spaces */
         char varname[10];
         int j;
         /* Do this instead of memcpy because it works on Crays */
         varname[0] = (char) ((header[6] >> 24) & 0xff);
         varname[1] = (char) ((header[6] >> 16) & 0xff);
         varname[2] = (char) ((header[6] >> 8) & 0xff);
         varname[3] = (char) ( header[6] & 0xff);
         varname[4] = 0;
         for (j=3;j>=0 && varname[j]==' ';j--) {
            varname[j] = 0;
         }
/* 01Feb06  Phil McDonald */
         info->VarName = TMP_STRDUP (varname);
/* end PM */
      }


      /*** Map Projection ***/
      switch (header[33]) {
         case 1:   /* Pseudo-Mercator */
            args[0] = header[34] / 10000.0;   /* north bound */
            args[1] = header[35] / 10000.0;   /* west bound */
            args[2] = header[38] / 10000.0;   /* row inc */
            args[3] = header[38] / 10000.0;   /* column inc */
            info->Proj = new_projection( db, PROJ_LINEAR, nr, nc, args );
            break;
         case 2:   /* Polar Stereographic or Lambert Conformal */
            if (header[38]==header[39]) {
               /* Polar Stereographic */
               /* TODO: verify these parameters are 100% correct! */
               args[0] = header[38] / 10000.0;   /* center latitude */
               args[1] = header[37] / 10000.0;   /* center longitude */
               args[2] = header[34] / 10000.0;   /* center grid row */
               args[3] = header[35] / 10000.0;   /* center grid column */
               args[4] = header[36] / 10000.0;   /* Col. inc in km */
               info->Proj = new_projection( db, PROJ_STEREO, nr, nc, args );
            }
            else {
               /* Lambert Conformal */
               args[0] = header[38] / 10000.0;   /* Lat1 */
               args[1] = header[39] / 10000.0;   /* Lat2 */
               args[2] = header[34] / 10000.0;   /* pole row */
               args[3] = header[35] / 10000.0;   /* pole column */
               args[4] = header[37] / 10000.0;   /* Cent. Long. */
               args[5] = header[36] / 1000.0;    /* Col. Inc in km */
               info->Proj = new_projection( db, PROJ_LAMBERT, nr, nc, args );
            }
            break;
         case 3:   /* Equidistant */
            args[0] = header[34] / 10000.0;  /* north bound */
            args[1] = header[35] / 10000.0;  /* west bound */
            args[2] = header[38] / M_PER_DEG;  /* row increment */
            args[3] = header[37] / M_PER_DEG;  /* column increment */
            /* header[36] = rotation, unsupported by Vis5D */
            info->Proj = new_projection( db, PROJ_LINEAR, nr, nc, args );
            break;
         case 4:   /* Pseudo-Mercator (more general) */
            /* TODO: is this 100% correct? */
            args[0] = header[34] / 10000.0;   /* north bound */
            args[1] = header[35] / 10000.0;   /* west bound */
            args[2] = header[38] / 10000.0;   /* row inc */
            args[3] = header[39] / 10000.0;   /* column inc */
            info->Proj = new_projection( db, PROJ_LINEAR, nr, nc, args );
            break;
         case 5:   /* No Navigation */
            args[0] = 0.0;   /* north */
            args[1] = 0.0;   /* west */
            args[2] = 1.0;   /* row inc */
            args[3] = 1.0;   /* column inc */
            info->Proj = new_projection( db, PROJ_GENERIC, nr, nc, args );
            break;
         case 6:
            /* Lambert Conformal tangent cone projection */
            args[0] = header[38] / 10000.0;   /* Lat1 */
            args[1] = header[38] / 10000.0;   /* Lat2 */
            args[2] = header[34] / 10000.0;   /* pole row */
            args[3] = header[35] / 10000.0;   /* pole column */
            args[4] = header[37] / 10000.0;   /* Cent. Long. */
            args[5] = header[36] / 1000.0;    /* Col. Inc in km */
/*
            printf("LambConfCone: %g %g %g %g %g %g\n", args[0], args[1],
                   args[2], args[3], args[4], args[5] );
*/
            info->Proj = new_projection( db, PROJ_LAMBERT, nr, nc, args );
            break;
         case 20:
            /* Cylindrical projection of cyl-equid. */
            args[0] = header[34] / 10000.0;  /* north bound */
            args[1] = header[35] / 10000.0;  /* west bound */
            args[2] = header[38] / M_PER_DEG;  /* row increment */
            args[3] = header[37] / M_PER_DEG;  /* column increment */
            /* header[36] = rotation, unsupported by Vis5D */
            info->Proj = new_projection( db, PROJ_CYLINDRICAL, nr, nc, args);
            break;
         case 21:
            /* Spherical projection of cyl-equid. */
            args[0] = header[34] / 10000.0;  /* north bound */
            args[1] = header[35] / 10000.0;  /* west bound */
            args[2] = header[38] / M_PER_DEG;  /* row increment */
            args[3] = header[37] / M_PER_DEG;  /* column increment */
            /* header[36] = rotation, unsupported by Vis5D */
            info->Proj = new_projection( db, PROJ_CYLINDRICAL, nr, nc, args);
            break;
         default:  /* Unknown projection! */
            printf("Warning: file %s, grid %d  has undefined projection: %d\n",
                   name, i, header[33] );
            args[0] = 0.0;   /* north */
            args[1] = 0.0;   /* west */
            args[2] = 1.0;   /* row inc */
            args[3] = 1.0;   /* column inc */
            info->Proj = new_projection( db, PROJ_GENERIC, nr, nc, args );
      }
         

      /*** Vertical Coordinate System ***/
/*
      printf("%05d %06d %-5s ",info->DateStamp,info->TimeStamp,info->VarName);
      printf("header[10] = %d  [11] = %d  [12] = %d\n",
             header[9], header[10], header[11] );
*/
      vertlevel = header[9];
      units = header[11];
      if (units==0x4d422020) {   /* "MB  " */
         units = 1;
      }
      else if (units==0x4b4d2020) {  /* "KM  " */
         units = 2;
      }
      if (vertlevel==1013) {
         /* Mean Sea Level */
         height = 0.0;    /* TODO: revisit */
      }
      else if (vertlevel==999) {
         /* Undefined */
         height = 0.0;    /* TODO: revisit */
         printf("Skipping grid %d, vertical level = undefined\n", i );
         continue;
      }
      else if (vertlevel==0) {
         /* Tropopause */
         height = 0.0;    /* TODO: fix */
         printf("Skipping grid %d, vertical level = TRO\n", i );
         continue;
      }
      else if (vertlevel==1001) {
         /* Surface */
         height = 0.0;
      }
      else if (vertlevel==0x4d415857) {
         /* "MAXW" */
         printf("Skipping grid %d, vertical level = MAXW\n", i );
         continue;
      }
      else {
         /* convert integer level to a float */
         height = vertlevel / pow( 10.0, (double) header[10] );
      }

      if (units==1 && height>0.0) {
         /* convert height from MB to KM */
         height = -7.2 * log( height / 1012.5 );
      }
      args[0] = height;
      args[1] = 1.0; /* WLH 7-12-96 */
      info->Vcs = new_vcs( db, VERT_EQUAL_KM, 1, 0, args );

      append_grid( info, db );
      grids++;


      if (height<minheight)  minheight = height;
      if (height>maxheight)  maxheight = height;

   }

/*
   printf("minheight=%g   maxheight=%g\n", minheight, maxheight );
*/

   return grids;
#else
   printf("Warning: can't read McIDAS GRID files on this system.\n");
   return 0;
#endif
}




/*
 * Get the grid data described by g.
 * Input:  g - pointer to a grid_info structure.  It tells us which grid
 *             in which file is needed.
 * Return:  pointer to grid data (can be free()'d when finished)
 *          OR return NULL if there's an error.
 */
float *get_grid_data( struct grid_info *g )
{
#ifdef HAVE_MCIDAS
   int *idata;
   float *data;
   int k, size, nr, nc, header[64];
   int pos;

   size = g->Nr * g->Nc;
   idata = (int *) MALLOC( size * sizeof(int) );

   pos = g->Position;
   k = IGGET( &g->FileNumber, &g->Position, &size, idata, &nr, &nc, header );

   if (k<0) {
      FREE( idata, 100 );
      return NULL;
   }

   data = (float *) MALLOC( size * sizeof(float) );

   /* convert and scale data from int to real */
   {
      int scale_exp;
      float scale;
      int i;
      scale_exp = header[7];
      scale = pow( 10.0, (double) scale_exp );
      scale = 1.0 / scale;
      for (i=0;i<size;i++) {
         if (idata[i]==0x80808080) {
            data[i] = MISSING;
         }
         else {
            data[i] = (float) idata[i] * scale;
         }
      }
   }

   if (Debug_i) {
      printf("IGGET: ");
      print_min_max( data, size );
   }

   FREE( idata, 101 );

   return data;
#else
   return NULL;
#endif
}

