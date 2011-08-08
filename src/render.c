/* render.c */



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

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#include "anim.h"
#include "api.h"
#include "box.h"
#include "matrix.h"
#include "globals.h"
#include "graphics.h"
#include "grid.h"
#include "labels.h"
#include "map.h"
#include "memory.h"
#include "misc.h"
#include "proj.h"
#include "queue.h"
#include "sounding.h"
#include "sync.h"
#include "topo.h"
#include "vis5d.h"
#include "volume.h"
#include "v5d.h"

extern int vis5d_verbose;

#define MAX(A,B)  ( (A) > (B) ? (A) : (B) )
#define MIN(A,B)  ( (A) < (B) ? (A) : (B) )

#define ABS(X)  ( (X) < 0.0f ? -(X) : (X) )


#define TICK_SIZE 0.05
#define CLOCK_SEGMENTS  36
#define VERT(Z) (ctx->VerticalSystem==VERT_NONEQUAL_MB ? height_to_pressure(Z) : (Z))
#define VERTPRIME(Z) (dtx->VerticalSystem==VERT_NONEQUAL_MB ? height_to_pressure(Z) : (Z))


/* Vertical spacing between rows of text: (in pixels) */
#define VSPACE 1


/*** float2string *****************************************************
   Convert a float into an ascii string.
**********************************************************************/
/* MJK 12.01.98 */
/* old now
void float2string( float f, char *str )
{
   if (f==0.0 || fabs(f)<0.01)
      strcpy( str, "0.0" );
   else if (f>=100.0 || f<=-100.0)
      sprintf(str, "%d", (int) f );
   else
      sprintf(str, "%4.2f", f );
}
*/

/* MJK 12.01.98 */
void float2string (Display_Context dtx, int icoord, float f, char *str)
{
   float        fmin, fmax, frng;

   if (fabs (f) < 0.01) f = 0.0;

   if (dtx->CoordFlag) {
      switch (icoord) {
         case 0:
            fmin = 1.0;
            fmax = dtx->Nc;
            break;
         case 1:
            fmin = 1.0;
            fmax = dtx->Nr;
            break;
         case 2:
            fmin = 1.0;
            fmax = dtx->MaxNl;
            break;
      }
   }
   else {
      switch (icoord) {
         case 0:
            fmin = dtx->WestBound;
            fmax = dtx->EastBound;
            break;
         case 1:
            fmin = dtx->SouthBound;
            fmax = dtx->NorthBound;
            break;
         case 2:
            fmin = VERTPRIME(dtx->BottomBound);
            fmax = VERTPRIME(dtx->TopBound);
            break;
      }
   }
   frng = (fmin < fmax) ? fmax - fmin : fmin - fmax;

   if (frng < 500.0)
      sprintf (str, "%.2f", f);
   else
      sprintf (str, "%.0f", f);
}

/*** check_view_side *****************************************************
   Determine if the plane of a clockwise series of points faces the camera.
   
   return:      -1      plane faces away from the camera
                 0      plane includes the camera
                 1      plane faces the camera
**********************************************************************/
int check_view_side (Context ctx, int type, int num)
{

    int         iside;
    float       xyz[3][3], xy[3][2], area;



    switch (type)
    {
/* TODO: Need to work on non-vertical slices */
        case VSLICE:
            xyz[0][0] = ctx->Variable[num]->VSliceRequest->X2;
            xyz[0][1] = ctx->Variable[num]->VSliceRequest->Y2;
            xyz[0][2] = ctx->dpy_ctx->Zmin;
            xyz[1][0] = ctx->Variable[num]->VSliceRequest->X1;
            xyz[1][1] = ctx->Variable[num]->VSliceRequest->Y1;
            xyz[1][2] = ctx->dpy_ctx->Zmin;
            xyz[2][0] = ctx->Variable[num]->VSliceRequest->X1;
            xyz[2][1] = ctx->Variable[num]->VSliceRequest->Y1;
            xyz[2][2] = ctx->dpy_ctx->Zmax;
            break;

        default:
            return 0;
    }


    project (&xyz[0][0], &xy[0][0], &xy[0][1]);
    project (&xyz[1][0], &xy[1][0], &xy[1][1]);
    project (&xyz[2][0], &xy[2][0], &xy[2][1]);

    area = ((xy[0][0] - xy[2][0]) * (xy[0][1] + xy[2][1])) +
           ((xy[1][0] - xy[0][0]) * (xy[1][1] + xy[0][1])) +
           ((xy[2][0] - xy[1][0]) * (xy[2][1] + xy[1][1]));

    iside = (area > 0.0) ? -1 : (area < 0.0) ? 1 : 0;


    return iside;
}


int flip_vslice_end_for_end (Context ctx, int time, int var)
{
   float		x;
   vslice_request	*vslice_req;

   vslice_req = ctx->Variable[var]->VSliceRequest;

   x              = vslice_req->R1;
   vslice_req->R1 = vslice_req->C1;
   vslice_req->C1 = x;

   x              = vslice_req->R2;
   vslice_req->R2 = vslice_req->C2;
   vslice_req->C2 = x;

   x              = vslice_req->X1;
   vslice_req->X1 = vslice_req->X2;
   vslice_req->X2 = x;

   x              = vslice_req->Y1;
   vslice_req->Y1 = vslice_req->Y2;
   vslice_req->Y2 = x;

   x                = vslice_req->Lat1;
   vslice_req->Lat1 = vslice_req->Lat2;
   vslice_req->Lat2 = x;

   x                = vslice_req->Lon1;
   vslice_req->Lon1 = vslice_req->Lon2;
   vslice_req->Lon2 = x;

   request_vslice (ctx, time, var, (time == ctx->CurTime));


   return 0;
}




/* same code as in uvwwidget.c */
char *return_var_plus_index( char *varname, int index )
{
   int yo;
/* WLH 20 Oct 98
   char whole[40];
*/
   /* WLH 20 Oct 98 */
   char* whole;

   char num[40];


   /* WLH 20 Oct 98 */
/* 01Feb06  Phil McDonald */
   whole = (char *) TMP_MALLOC (40);
/* end PM */

   if (index <0 || varname[0] == 0){
      whole[0] = 0;
      return whole;
   }
   for (yo = 0; yo < 17; yo++){
      if (varname[yo] == '\0' || varname[yo] == ' '){
         yo -=1;
         whole[yo+1] = '.';
         sprintf(num, "%d\n", index );
         if (index > 99 && yo < 15){
            whole[yo+2] = num[0];
            whole[yo+3] = num[1];
            whole[yo+4] = num[2];
            whole[yo+5] = '\0';
            return whole;
         }
         else if ( index > 9 && yo < 16){
            whole[yo+2] = num[0];
            whole[yo+3] = num[1];
            whole[yo+4] = '\0';
            return whole;
         }
         else{
            whole[yo+2] = num[0];
            whole[yo+3] = '\0';
            return whole;
         }
      }
      whole[yo] = varname[yo];
   }
   whole[yo] = '\0';
   return whole;
}


/*** plot_string  ******************************************************
   Plot (draw) a string in 3-D.  At this time, only strings of digits,
   periods, and dashes are implemented.
   Input: f - the string to plot.
          startx, y, z - the point in 3-D to start at.
          base - vector indicating text baseline.
          up - vector indicating upright direction for text.
          rjustify - non-zero value indicates right justify the text.
**********************************************************************/
void plot_string( char *str, float startx, float starty, float startz,
                  float base[], float up[], int rjustify )
{
   static float zero[] = { 0,0, 0,.8, .4,.8, .4,0, 0,0 },
      one[] = { 0,0, 0,.8 },
      two[] = { .4,0, 0,0, 0,.4, .4,.4, .4,.8, 0,.8 },
      three[] = { 0,0, .4,0, .4,.4, 0,.4, .4,.4, .4,.8, 0,.8 },
      four[] = { 0,.8, 0,.4, .4,.4, .4,.8, .4,0 },
      five[] = { 0,0, .4,0, .4,.4, 0,.4, 0,.8, .4,.8 },
      six[] = { .4,.8, 0,.8, 0,0, .4,0, .4,.4, 0,.4 },
      seven[] = { 0,.7, 0,.8, .4,.8, .4,0 },
      eight[] = { 0,0, 0,.8, .4,.8, .4,0, 0,0, 0,.4, .4,.4 },
      nine[] = { .4,.4, 0,.4, 0,.8, .4,.8, .4,0 },
      dash[] = { 0,.4, .4,.4 },
      dot[] = { 0,0, 0,.1, .1,.1, .1,0, 0,0 },
 /*MiB  03/2001 Longitudes*/
      west[] = {0.,0.8, 0.,0., 0.2,0.4, 0.4,0.0, 0.4,0.8},
      east[] = {0.4,0.8, 0.0,0.8, 0.0,0.4, 0.3,0.4, 0.0,0.4, 0.0,0.0, 0.4,0.0},
      north[] = {0.,0.0, 0.,0.8, 0.4,0.0, 0.4,0.8},
      south[] = {0.0,0.1, 0.1,0.0, 0.3,0.0, 0.4,0.1, 0.4,0.3, 
		 0.0,0.5, 0.0,0.7, 0.1,0.8, 0.3,0.8, 0.4,0.7};

   static float *index[16] = { zero, one, two, three, four, five, six,
                               seven, eight, nine, dash, dot, west, east,
				north, south};

   static float width[16] = { 0.6, 0.2, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,
      0.6, 0.6, 0.6, 0.3, 0.6, 0.6 , 0.6, 0.6};
   static int verts[16] = { 5, 2, 6, 7, 5, 6, 6, 4, 7, 5, 2, 5 ,5, 7, 4, 10};


   float *temp, plot[100][3];
   float cx, cy, cz;
   int i, j, k, len;

   cx = startx;  cy = starty;  cz = startz;
   len = strlen(str);

   if (rjustify) {
      /* draw right justified text */
      for (i=len-1; i>=0; i--) {
         if (str[i]=='-')
            k = 10;
         else if (str[i]=='.')
            k = 11;
/*MiB*/  else if (str[i]=='W')
            k = 12;
         else if (str[i]=='E')
            k = 13;
         else if (str[i]=='N')
            k = 14;
         else if (str[i]=='S')
/*MiB*/     k = 15;
         else if (str[i]>='0' && str[i]<='9')
            k = str[i] - '0';
         else
            continue;
         /* calculate position for this char */
         cx += width[k]*base[0];
         cy += width[k]*base[1];
         cz += width[k]*base[2];
         /* make the vertex array for this character */
         temp = index[k];
         for (j=0; j<verts[k]; j++) {
            float x, y;
            x = *temp++;
            y = *temp++;
            plot[j][0] = cx - x*base[0] + y*up[0];
            plot[j][1] = cy - x*base[1] + y*up[1];
            plot[j][2] = cz - x*base[2] + y*up[2];
         }
         polyline( plot, verts[k] );
      }

   }
   else {
      /* draw left justified text */
      for (i=0; i<len; i++) {
         if (str[i]=='-')
            k = 10;
         else if (str[i]=='.')
            k = 11;
/*MiB*/  else if (str[i]=='W')
            k = 12;
         else if (str[i]=='N')
            k = 14;
         else if (str[i]=='S')
            k = 15;
         else if (str[i]=='E')
/*MiB*/     k = 13;
         else if (str[i]>='0' && str[i]<='9')
            k = str[i] - '0';
         else
            continue;
         /* make the vertex array for this character */
         temp = index[k];
         for (j=0; j<verts[k]; j++) {
            float x, y;
            x = *temp++;
            y = *temp++;
            plot[j][0] = cx + x*base[0] + y*up[0];
            plot[j][1] = cy + x*base[1] + y*up[1];
            plot[j][2] = cz + x*base[2] + y*up[2];
         }
         polyline( plot, verts[k] );
         /* calculate position for next char */
         cx += width[k]*base[0];
         cy += width[k]*base[1];
         cz += width[k]*base[2];
      }
   }

}

/* TODO: this function should provide support for users
   to choose the policy by which a time is valid.
*/

int check_for_valid_time( Context ctx, int dtxcurtime)
{
   Display_Context dtx;
   int dtime, stime, ldtime, lstime;
   int ctxcurtime, yo, ctime, spandex;
   int ctxdtime, ctxstime;
   dtx = ctx->dpy_ctx;
   
   if (dtx->numofctxs == 1){
      return 1;
   }
   if (ctx->NumTimes == 1){
      return 1;
   }
   ctxcurtime = 
   ldtime = lstime = 0;


   for( yo = 0; yo < dtx->numofctxs; yo ++){
      spandex = dtx->TimeStep[dtxcurtime].owners[yo];
      ctime = dtx->TimeStep[dtxcurtime].ownerstimestep[yo];
      vis5d_get_ctx_time_stamp( spandex, ctime, &dtime, &stime);
      if (spandex == ctx->context_index){
         ctxcurtime = ctime;
         ctxdtime = dtime;
         ctxstime = stime;
      }
      else if (dtime > ldtime || (dtime == ldtime && stime > lstime)){
         ldtime = dtime;
         lstime = stime;
      }
   }
   if (ctxcurtime == 0 && (ctxdtime > ldtime || 
       (ctxdtime == ldtime && ctxstime > lstime))){
       return 0;
   }
   else if (ctxcurtime == ctx->NumTimes-1 && (ctxdtime < ldtime || 
       (ctxdtime == ldtime && ctxstime < lstime))){
       return 0;
   }
   else{
      return 1;
   }
}

/*
 * Draw the tick mark for a horizontal slice.
 * Input:  level - grid level
 *         z - graphics z coord
 *         height - geographic height coord
 */
static void draw_horizontal_slice_tick( Display_Context dtx, float level,
                                        float z, float height )
{
   float v[2][3];
   static float base[3] = { 0.035, -0.035, 0.0 };
   static float up[3] = { 0.0, 0.0, 0.07 };
   char str[1000];

   /* vertices for tick mark */
   v[0][0] = dtx->Xmax;
   v[0][1] = dtx->Ymin;
   v[0][2] = z;
   v[1][0] = dtx->Xmax + 0.05;
   v[1][1] = dtx->Ymin - 0.05;
   v[1][2] = z;
   polyline( v, 2 );

   /* the label */
   if (dtx->CoordFlag) {
      /* MJK 12.01.98 */
      float2string(dtx, 2,  level+1.0, str );
   }
   else {
      /* MJK 12.01.98 */
      float2string(dtx, 2, VERTPRIME(height), str );
   }
   plot_string( str, dtx->Xmax+0.07, dtx->Ymin-0.07, z, base,up, 0 );
}



/*
 * Draw a tick mark for a vertical slice.
 * Input:  row, col - position in grid coords
 *         x, y - position in graphics coords
 *         lat,lon - position in geographic coords
 */
static void draw_vertical_slice_tick( Display_Context dtx, float row, float col,
                                      float x, float y, float lat, float lon )
{
   float v[2][3];
   int cc, rr;
   /* base and up vectors for drawing 3-D text */
   static float b2[3] = { 0.05, 0.0, 0.0 }, u2[3] = { 0.0, 0.05, 0.0 };
   static float b3[3] = { -0.05, 0.0, 0.0 }, u3[3] = { 0.0, 0.05, 0.0 };
   char str[1000];

   cc = (int) (col);
   rr = (int) (row);
   
   if (cc <= 0) {
      /* draw on top-west edge */
      v[0][0] = x;
      v[0][1] = y;
      v[0][2] = dtx->Zmax;
      v[1][0] = x-0.05;
      v[1][1] = y;
      v[1][2] = dtx->Zmax;
      polyline( v, 2 );
      if (dtx->CoordFlag) {
         /* MJK 12.01.98 */
         float2string(dtx, 1, row+1, str );
      }
      else {
         /* MJK 12.01.98 */      
         float2string(dtx, 1, lat, str );
      }
      plot_string( str, x-0.07, y, dtx->Zmax, b3, u3, 1 );
   }
   else if (cc >= dtx->Nc-1) {
      /* draw on top-east edge */
      v[0][0] = x;
      v[0][1] = y;
      v[0][2] = dtx->Zmax;
      v[1][0] = x+0.05;
      v[1][1] = y;
      v[1][2] = dtx->Zmax;
      polyline( v, 2 );
      if (dtx->CoordFlag) {
         /* MJK 12.01.98 */      
         float2string(dtx, 1, row+1, str );
      }
      else {
         /* MJK 12.01.98 */      
         float2string(dtx, 1, lat, str );
      }
      plot_string( str, x+0.07, y, dtx->Zmax, b2, u2, 0 );
   }
   else if (rr <=0) {
      /* draw on top-north edge */
      v[0][0] = x;
      v[0][1] = y;
      v[0][2] = dtx->Zmax;
      v[1][0] = x;
      v[1][1] = y+0.05;
      v[1][2] = dtx->Zmax;
      polyline( v, 2 );
      if (dtx->CoordFlag) {
         /* MJK 12.01.98 */               
         float2string(dtx, 0, col+1.0, str );
      }
      else {
         /* MJK 12.01.98 */               
         float2string(dtx, 0, lon, str );
      }
      plot_string( str, x-0.07, y+0.07, dtx->Zmax, b2,u2, 0 );
   }
   else {
      /* draw on top-south edge */
      v[0][0] = x;
      v[0][1] = y;
      v[0][2] = dtx->Zmax;
      v[1][0] = x;
      v[1][1] = y-0.05;
      v[1][2] = dtx->Zmax;
      polyline( v, 2 );
      if (dtx->CoordFlag) {
         /* MJK 12.01.98 */                        
         float2string(dtx, 0, col+1.0, str );
      }
      else {
         /* MJK 12.01.98 */                        
         float2string(dtx, 0, lon, str );
      }
      plot_string( str, x-0.07, y-0.12, dtx->Zmax, b2,u2, 0 );
   }
}



/* 03Feb06  Phil McDonald	cleaned things up a little */
/*
 * Draw the tick mark for a horizontal slice.
 * Input:  level - grid level
 *         z - graphics z coord
 *         height - geographic height coord
 */
static void	draw_hslice_tick (Display_Context dtx, float level)
{
   static float	base[3] = { 0.035, -0.035, 0.0 };
   static float	up[3] = { 0.0, 0.0, 0.07 };
   float	v[2][3];
   float	z, hgt, xtmp;
   char		str[1000];


   clipping_off ();

   xtmp = 0.0;
   gridPRIME_to_xyzPRIME (dtx, -1, -1, 1, &xtmp, &xtmp, &level,
                          &xtmp, &xtmp, &z);
   xtmp = 0.0;
   xyzPRIME_to_geo (dtx, -1, -1, xtmp, xtmp, z, &xtmp, &xtmp, &hgt);

   /* vertices for tick mark */
   v[0][0] = dtx->Xmax;
   v[0][1] = dtx->Ymin;
   v[0][2] = z;
   v[1][0] = dtx->Xmax + 0.05;
   v[1][1] = dtx->Ymin - 0.05;
   v[1][2] = z;
   polyline( v, 2 );

   /* the label */
   if (dtx->CoordFlag) {
      /* MJK 12.01.98 */
      float2string(dtx, 2,  level+1.0, str );
   }
   else {
      /* MJK 12.01.98 */
      float2string(dtx, 2, VERTPRIME(hgt), str );
   }
   plot_string( str, dtx->Xmax+0.07, dtx->Ymin-0.07, z, base,up, 0 );

   clipping_on ();


   return;
}



/* 08Feb06  Phil McDonald	cleaned things up a little */
/*
 * Draw a tick mark for a vertical slice.
 * Input:  row1, col1, row2, col2 - positions in grid coords
 */
static void	draw_vslice_ticks (Display_Context dtx, float row1, float col1,
		                   float row2, float col2)
{

    float	x1, y1, x2, y2, lat, lon, xtmp;
    float	zbot, ztop, vert[4][3];


    clipping_off ();

    xtmp = 0.0;
    gridPRIME_to_xyzPRIME (dtx, -1, -1, 1, &row1, &col1, &xtmp,
                           &x1, &y1, &xtmp);
    xtmp = 0.0;
    gridPRIME_to_geo (dtx, -1, -1, 1, &row1, &col1, &xtmp, &lat, &lon, &xtmp);
    draw_vertical_slice_tick (dtx, row1, col1, x1, y1, lat, lon);

    xtmp = 0.0;
    gridPRIME_to_xyzPRIME (dtx, -1, -1, 1, &row2, &col2, &xtmp,
                           &x2, &y2, &xtmp);
    xtmp = 0.0;
    gridPRIME_to_geo (dtx, -1, -1, 1, &row2, &col2, &xtmp, &lat, &lon, &xtmp);
    draw_vertical_slice_tick (dtx, row2, col2, x2, y2, lat, lon);

    zbot = gridlevelPRIME_to_zPRIME (dtx, -1, -1, (float) dtx->LowLev);
    ztop = gridlevelPRIME_to_zPRIME (dtx, -1, -1,
                                     (float) (dtx->LowLev + dtx->Nl));

    vert[0][0] = vert[1][0] = vert[2][0] = vert[3][0] = (x1 + x2) * 0.5;
    vert[0][1] = vert[1][1] = vert[2][1] = vert[3][1] = (y1 + y2) * 0.5;
    vert[0][2] = ztop + TICK_SIZE;
    vert[1][2] = ztop;
    vert[2][2] = zbot;
    vert[3][2] = zbot - TICK_SIZE;

    set_line_width (5); /* WLH 3-5-96 */
    disjointpolyline (vert, 4);
    set_line_width (dtx->LineWidth); /* WLH 3-5-96 */

    clipping_on ();


    return;
}
/* end PM */


/*
 * Print the current cursor position.
 */
static void print_cursor_position( Display_Context dtx, int it )
{
   static float bx[3] = { 0.05, 0.0, 0.0 },      ux[3] = { 0.0, 0.05, 0.05 };
   static float by[3] = { -0.035, 0.0, -0.035 },  uy[3] = { 0.0, 0.07, 0.0 };
   static float bz[3] = { -0.035, -0.035, 0.0 }, uz[3] = { 0.0, 0.0, 0.07 };
   float v[6][3];
   float x, y, z, xx, yy;
   char str[100], xdir1[8],ydir1[8];

   /* MJK 12.01.98 */
   float lat, lon, hgt, row, col, lev;
   int ix;
   char fmt[] = {"%s: %9.3f %s  "};

/* MiB 03/2001   Longitudes in  [-180, 180]    Intitialize...*/
   xdir1[0] = ' ';
   xdir1[1] = '\0';
   ydir1[0] = ' ';
   ydir1[1] = '\0';

   /* MJK 12.01.98 begin*/
   if ((dtx->DisplayProbe) || (dtx->DisplaySound)){
      /* MJK 3.29.99 */
      if (dtx->Reversed){
         set_color (PACK_COLOR(0,0,0,255));
      }
      else{
         set_color (dtx->BoxColor);
      }
   }
   else{
      set_color( *dtx->CursorColor );
   }
   /* end MJK 12.01.98 */

   if (dtx->Projection==PROJ_LINEAR || dtx->Projection==PROJ_GENERIC
       /* ZLB */ || dtx->Projection==PROJ_GENERIC_NONEQUAL) {
      /* Rectangular box:  put labels along edge of box in 3-D */

      set_depthcue( dtx->DepthCue );

      /* draw tick marks */
      v[0][0] = v[1][0] = dtx->CursorX;
      v[0][1] = dtx->Ymin;  v[1][1] = dtx->Ymin-0.05;
      v[0][2] = dtx->Zmin;  v[1][2] = dtx->Zmin-0.05;
      v[2][0] = dtx->Xmin;  v[3][0] = dtx->Xmin-0.05;
      v[2][1] = v[3][1] = dtx->CursorY;
      v[2][2] = dtx->Zmin;  v[3][2] = dtx->Zmin-0.05;
      v[4][0] = dtx->Xmin;  v[5][0] = dtx->Xmin-0.05;
      v[4][1] = dtx->Ymin;  v[5][1] = dtx->Ymin-0.05;
      v[4][2] = v[5][2] = dtx->CursorZ;
      if (dtx->DisplaySound) {
         v[5][0] = 0;
         v[5][1] = 0;
         v[4][0] = 0;
         v[4][1] = 0;
      } 
      disjointpolyline( v, 6 );

      /* draw position labels */
      if (dtx->CoordFlag) {
         /* display cursor position in grid coordinates */
         xyzPRIME_to_gridPRIME( dtx, it, -1,
                      dtx->CursorX, dtx->CursorY, dtx->CursorZ, &y, &x, &z );
         x += 1.0;
         y += 1.0;
         z += 1.0;
      }
      else {
         /* display cursor position in geographic coordinates */
         xyzPRIME_to_geo( dtx, it, -1,
                     dtx->CursorX, dtx->CursorY, dtx->CursorZ, &y, &x, &z );
         z = VERTPRIME(z);
      }

/*MiB   03/2001 limit range to -180, 180 */
         xx = x;
         if (xx < -180.){
	     xx = 360. + xx;
         }
         if (xx > 180.){
	     xx = -360. + xx;
         }
/*MiB   03/2001 Define East/West  */
         if (xx > 0.){
	     	xdir1[0] = 'W';
         } else {
                xdir1[0] = 'E';
                xx = xx *(-1.);
                }
/*MiB   03/2001 Define North/South  */
	 yy = y;
         if (yy > 0.){
	     	ydir1[0] = 'N';
         } else {
                ydir1[0] = 'S';
                yy = yy *(-1.);
                }
  

      /* MJK 12.01.98 */                        
      float2string(dtx, 0, xx, str );
/*MiB*/  strcat(str,xdir1);
      plot_string( str, dtx->CursorX-0.04, dtx->Ymin-0.1,
                   dtx->Zmin-0.125, bx, ux, 0 );
      float2string(dtx, 1, yy, str );
/*MiB*/  strcat(str,ydir1);
      plot_string( str, dtx->Xmin-0.075, dtx->CursorY-0.02,
                   dtx->Zmin-0.075, by, uy, 1 );
      float2string(dtx, 2, z, str );


      if (!dtx->DisplaySound)
      plot_string( str, dtx->Xmin-0.07, dtx->Ymin-0.07,
                   dtx->CursorZ+0.005, bz, uz, 1 );

      set_depthcue( 0 );
   }


   /* MJK 12.01.98 */
   ix = dtx->Nr;
   if (dtx->Nc > ix) ix = dtx->Nc;
   if (dtx->MaxNl > ix) ix = dtx->MaxNl;
   x = ix;

   v[0][0] = dtx->Xmin, v[0][1] = dtx->Ymin, v[0][2] = dtx->Zmin;
   v[1][0] = dtx->Xmax, v[1][1] = dtx->Ymax, v[1][2] = dtx->Zmax;
   v[2][0] = dtx->Xmin, v[2][1] = dtx->Ymax, v[2][2] = dtx->Zmax;
   v[3][0] = dtx->Xmax, v[3][1] = dtx->Ymin, v[3][2] = dtx->Zmax;
   for (ix = 0; ix < 4; ix++)
   {
      xyzPRIME_to_geo (dtx, it, -1, v[ix][0], v[ix][1], v[ix][2],
                       &lat, &lon, &hgt);
      if (lat < 0.0) lat = -lat;
      if (lon < 0.0) lon = -lon;
      hgt = VERTPRIME(hgt);
      if (lat > x) x = lat;
      if (lon > x) x = lon;
      if (hgt > x) x = hgt;
   }

   sprintf (str, "%.3f", x);
   ix = strchr (str, '.') - str;
   if (ix < 1) ix = 1;
   fmt[5] = (ix + 4) + '0';
   sprintf (str, fmt, "XXX", x, "xx");
   ix = dtx->WinWidth - (text_width (dtx->gfx[WINDOW_3D_FONT]->font,str));

   if (dtx->CoordFlag) {
      /* display cursor position in grid coordinates */
      xyzPRIME_to_gridPRIME (dtx, it, -1,
                             dtx->CursorX, dtx->CursorY, dtx->CursorZ,
                             &row, &col, &lev);

      sprintf( str, "Row: %g", row+1.0 );
      draw_text( ix, (dtx->gfx[WINDOW_3D_FONT]->FontHeight+VSPACE), str );

      sprintf( str, "Col: %g", col+1.0 );
      draw_text( ix, 2*(dtx->gfx[WINDOW_3D_FONT]->FontHeight+VSPACE), str );

      if (!dtx->DisplaySound)
      {
         sprintf( str, "Lev: %g", lev+1.0 );
         draw_text( ix, 3*(dtx->gfx[WINDOW_3D_FONT]->FontHeight+VSPACE), str );
      }

   }
   else {
      /* display cursor position in geographic coordinates */
      xyzPRIME_to_geo (dtx, it, -1, dtx->CursorX, dtx->CursorY, dtx->CursorZ,
                       &lat, &lon, &hgt);
      if (dtx->Projection == PROJ_GENERIC /* ZLB ?????*/)
      {
         sprintf (str, fmt, "Row", lat, "");
         draw_text( ix, (dtx->gfx[WINDOW_3D_FONT]->FontHeight+VSPACE), str );

         sprintf (str, fmt, "Col", lon, "");
         draw_text( ix, 2*(dtx->gfx[WINDOW_3D_FONT]->FontHeight+VSPACE), str );
      }
      else
      {
         char hemi[] = {"X"};

         /* hemispheres are "Wisconsin-centric" ;-) */
         hemi[0] = 'N';
         if (lat < 0.0) lat = -lat, hemi[0] = 'S';
         sprintf (str, fmt, "Lat", lat, hemi);
         draw_text( ix, (dtx->gfx[WINDOW_3D_FONT]->FontHeight+VSPACE), str );

         hemi[0] = 'W';
         if (lon < 0.0) lon = -lon, hemi[0] = 'E';
/*MiB*/  if (lon > 180.) lon = -360. + lon, hemi[0] = 'W';
/*MiB*/  if (lon < 0.0) lon = -lon, hemi[0] = 'E';
         sprintf (str, fmt, "Lon", lon, hemi);
         draw_text( ix, 2*(dtx->gfx[WINDOW_3D_FONT]->FontHeight+VSPACE), str );
      }

      if (!dtx->DisplaySound)
      {
         if (dtx->VerticalSystem == VERT_NONEQUAL_MB)
            sprintf (str, fmt, "Hgt", VERTPRIME(hgt), "mb");
         else
            sprintf (str, fmt, "Hgt", hgt, "km");
         draw_text( ix, 3*(dtx->gfx[WINDOW_3D_FONT]->FontHeight+VSPACE), str );
      }
   }
   /* end MJK 12.01.98 */
}



/*
 * Examine the entries of a color lookup table to determine the alpha value.
 * If the alpha is the same for all entries return that value.  If the
 * alpha varies between entries return -1.
 * Input:  color_table - array of color values
 *         size - number of entries in table
 * Return:  -1 or constant alpha value
 */
int get_alpha( unsigned int color_table[], int size )
{
   int alpha, i;

   /* Check for variable vs. constant alpha */
   alpha = UNPACK_ALPHA( color_table[0] );
   for (i=0;i<size;i++) {
      if (UNPACK_ALPHA(color_table[i])!=alpha) {
         return -1;
      }
   }
   return alpha;
}



/* 03Feb06  Phil McDonald */
static int	render_get_read_lock (int animflag, int* p_lock)
{
    if (animflag) return cond_read_lock (p_lock);

    wait_read_lock (p_lock);

    return 1;
}



/* 09Feb06  Phil McDonald */
static unsigned int	render_get_color (Context ctx, int type, int var)
{
    int		icolor;

    icolor = (ctx->context_index * MAXVARS) + var;

    return ctx->dpy_ctx->Color[icolor][type];
}



/* 09Feb06  Phil McDonald */
static unsigned int*	render_get_colortable (Context ctx, int type, int var)
{
    int		icolor;

    icolor = (ctx->context_index * MAXVARS) + var;

    return ctx->dpy_ctx->ColorTable[type]->Colors[icolor];
}




/* 03Feb06  Phil McDonald */
static int	render_dump (char* str1, int ip, int it, char* str2,
		             int* p_lock, void* ptr, int num)
{
    fprintf (stderr, "%s: pid=%d   %d %d   %s   ",
             str1, getpid (), ip, it, str2);
    if (p_lock != NULL) fprintf (stderr, "lock = %d %p   ", *p_lock, p_lock);
    if (ptr != NULL) fprintf (stderr, "%p %d", ptr, num);
    fprintf (stderr, "\n");


    return 1;
}
/* end PM */




/* 06Feb08  Phil McDonald	cleaned things up a little */
/*
 * Render all isosurfaces selected for display.
 * Input:  ctx - the context
 *         time -  dtx time!! 
 *         transflag - transparency flag: 1=only draw opaque surfaces
 *                                 0=only draw transparent surfaces
 */
static void	render_isosurfaces (Context ctx, int dtxtime, int ctxtime,
		                    int transflag, int animflag)
{
    Display_Context	dtx;
    struct isosurface	*isosurf;
    int			time, var;



    dtx = ctx->dpy_ctx;

    for (var = 0; var < ctx->NumVars; var++)
    {
	if ((ctx->SameIsoColorVarOwner[var]) || (ctx->IsoColorVar[var] < 0))
        {
            time = ctxtime;
	}
	else
        {
            time = dtxtime;
	}
	
        isosurf = ctx->Variable[var]->SurfTable[time];

        if (ctx->DisplaySurf[var] && isosurf && isosurf->valid)
        {
/* 05Sep06  Phil McDonald */
            if (render_get_read_lock (animflag, &(isosurf->lock)))
            {
/* end PM */
                int		colorvar, colorown, alpha;
                unsigned int	*ctab;

                recent (ctx, ISOSURF, var);
                colorvar = isosurf->colorvar;
                colorown = isosurf->cvowner;

                /* Determine alpha for surface:  -1=variable, 0..255=constant */
                if (isosurf->colors)
                {
                    ctab = dtx->ColorTable[VIS5D_ISOSURF_CT]->
                           Colors[colorown*MAXVARS+colorvar];
                }
                else
                {
                    ctab = dtx->Color[ctx->context_index*MAXVARS+var];
                }
                alpha = get_alpha (ctab, 256);

                if (((transflag) && (alpha == 255)) ||
                    ((!transflag) && (alpha < 255)))
                {
                    int		num, fastdraw;
                    void	*verts, *norms, *colors;

                    num    = isosurf->numindex;
                    verts  = isosurf->verts;
                    norms  = isosurf->norms;
                    colors = isosurf->colors;
                    if (isosurf->colors)
                    {
                        vis5d_check_fastdraw (dtx->dpy_context_index,
                                              &fastdraw);
                        if ((fastdraw || animflag) && isosurf->deci_verts)
                        {
                             num    = isosurf->deci_numverts;
                             verts  = isosurf->deci_verts;
                             norms  = isosurf->deci_norms;
                             colors = isosurf->deci_colors;
                             if (isosurf->deci_colors)
                             {
                                 draw_colored_isosurface (num, isosurf->index,
                                                          verts, norms, 1,
                                                          colors, ctab, alpha);
                             }
                             else
                             {
                                 draw_isosurface (num, isosurf->index,
                                                  verts, norms, 1,
                                                  ctab[0], NULL, GL_COMPILE);
                             }
                        }
                        else
                        {
                            draw_colored_isosurface (num, isosurf->index,
                                                     verts, norms, 0,
                                                     colors, ctab, alpha);
                        }
                    }
                    else
                    {
                        draw_isosurface (num, isosurf->index,
                                         verts, norms, 0,
                                         ctab[0], NULL, 0);
                    }
                }
/* 05Sep06  Phil McDonald */
                done_read_lock (&isosurf->lock);
/* end PM */
            }
        }
    }


    return;
}
/* end PM */



static void render_hclips( Display_Context dtx, int animflag)
{
  int i;
   for (i = 0; i < 2; i++){
      if (dtx->HClipTable[i].highlight == 1){
         set_color( PACK_COLOR(100,25,240,255));
         set_line_width(4);
      }
      else{
         set_color( PACK_COLOR(50,200,75,255));
         set_line_width(1);
      }
      polyline( (void *) dtx->HClipTable[i].boxverts,
                      dtx->HClipTable[i].numboxverts);
      /* MJK 3.29.99 */
      if (dtx->Reversed){
         set_color( PACK_COLOR(0,0,0,255) );
      }
      else{
         set_color( dtx->BoxColor );
      }


      set_line_width(dtx->LineWidth);
   }
}

static void render_vclips( Display_Context dtx, int animflag)
{
   int i;
   for (i = 0; i < 4; i++){
      if (dtx->VClipTable[i].highlight == 1){
         set_color( PACK_COLOR(100,25,240,255));
         set_line_width(4);
      }
      else{
         set_color( PACK_COLOR(50,200,75,255));
         set_line_width(1);
      }
      polyline( (void *) dtx->VClipTable[i].boxverts,
                         dtx->VClipTable[i].numboxverts);
      if (dtx->VClipTable[i].highlight == 1){
         float vert[4][3];
         float zbot, ztop, x1, x2, y1, y2;
         float llev, hlev;
         llev = (float) dtx->LowLev;
         hlev = (float) (dtx->Nl-1+dtx->LowLev);
         gridPRIME_to_xyzPRIME( dtx, 0, 0, 1, &dtx->VClipTable[i].r1,
                                &dtx->VClipTable[i].c1, &llev,
                                &x1, &y1, &zbot);
         gridPRIME_to_xyzPRIME( dtx, 0, 0, 1, &dtx->VClipTable[i].r2,
                                &dtx->VClipTable[i].c2, &hlev,
                                &x2, &y2, &ztop);
         vert[0][0] = vert[1][0] = vert[2][0] = vert[3][0] = (x1 + x2)*0.5;
         vert[0][1] = vert[1][1] = vert[2][1] = vert[3][1] = (y1 +y2)*0.5;
         vert[0][2] = ztop+TICK_SIZE;
         vert[1][2] = ztop;
         vert[2][2] = zbot;
         vert[3][2] = zbot-TICK_SIZE;
         set_line_width(5);
         disjointpolyline( vert, 4 );
      }   
/* MJK 3.29.99 */
   if (dtx->Reversed){
      set_color( PACK_COLOR(0,0,0,255) );
   }
   else{
      set_color( dtx->BoxColor );
   }


      set_line_width(dtx->LineWidth);
   }
}


static void render_textplots( Irregular_Context itx, int time)
{
   int var;
   float a, b, c, d;

   vis5d_get_text_plot(itx->context_index, &var, &a, &b, &c, &d);
   if (itx->DisplayTextPlot && itx->TextPlotTable[time].valid){
      if (itx->TextPlotTable[time].colors){
         draw_colored_disjoint_lines(itx->TextPlotTable[time].numverts,
                               (void *) itx->TextPlotTable[time].verts,
                               (void *) itx->TextPlotTable[time].colors,
                               itx->dpy_ctx->ColorTable[VIS5D_TEXTPLOT_CT]->Colors[itx->context_index*
                               MAXVARS+var]);
      }
      else{
         draw_disjoint_lines( itx->TextPlotTable[time].numverts,
                              (void *) itx->TextPlotTable[time].verts,
                              itx->dpy_ctx->TextPlotColor[itx->context_index*
                             MAXVARS+var] , NULL, 0);
      }
   }
}

/* 09Feb06  Phil McDonald	cleaned things up a little */
/*
 * Render all horizontal contour slices selected for display.
 * Input:  ctx - the context
 *         time - the time step
 *         labels - draw labels flag.
 */
static void	render_hslices (Context ctx, int time, int labels, int animflag)
{
    Display_Context	dtx;
    struct hslice	*hslice;
    int			var, sfc, box, stipple, lwidth, show_hidden;
    unsigned int	color;

    dtx = ctx->dpy_ctx;

    for (var = 0; var < ctx->NumVars; var++)
    {
        hslice = ctx->Variable[var]->HSliceTable[time];

        if (ctx->DisplayHSlice[var] && hslice && hslice->valid) 
        {
            if (render_get_read_lock (animflag, &(hslice->lock)))
            {
                /* copy data that isn't protected by the lock */
                color   = render_get_color (ctx, HSLICE, var);
                stipple = ctx->Variable[var]->HSliceRequest->stipple;
                lwidth  = ctx->Variable[var]->HSliceRequest->linewidth;
                sfc     = ctx->DisplaySfcHSlice[var];
                box     = dtx->DisplayBox && !dtx->CurvedBox;

                recent (ctx, HSLICE, var);
                set_color (color);
                set_line_width (lwidth);
#ifdef HAVE_OPENGL
                if (stipple != VIS5D_SOLID_LINE)
                {
                    glEnable (GL_LINE_STIPPLE);
                    glLineStipple (1, (GLushort) stipple);
                }
#endif

                /* draw main contour lines */
                draw_disjoint_lines (hslice->num1, (void*) hslice->verts1,
                                     color, NULL, 0);

                show_hidden = !labels;
#ifdef USE_SYSTEM_FONTS
                show_hidden = 1;
#endif
                if (show_hidden)
                {
                    /* draw hidden contour lines */
                    draw_disjoint_lines (hslice->num2, (void*) hslice->verts2,
                                         color, NULL, 0);
                }
#ifdef HAVE_OPENGL
                glDisable (GL_LINE_STIPPLE);
#endif

                /* draw contour labels */
                if (labels)
                {
#ifdef USE_SYSTEM_FONTS
                    plot_strings (hslice->num3, hslice->labels,
                                  hslice->verts3, color, 
                                  dtx->gfx[CONTOUR_LABEL_FONT]->fontbase);
#else
                    draw_disjoint_lines (hslice->num3, (void*) hslice->verts3,
                                         color, NULL, 0);
#endif
                }

                set_line_width (dtx->LineWidth);

                if (box && !sfc)
                {
                    /* draw the bounding box */
                    polyline ((void*) hslice->boxverts, hslice->numboxverts);

                    /* draw position label */
                    draw_hslice_tick (dtx, hslice->level);
                }

                done_read_lock (&(hslice->lock));
            }
        }
    }


    return;
}
/* end PM */



/* 04Feb08  Phil McDonald	cleaned things up a little */
/*
 * Render all vertical contour slices selected for display.
 * Input:  ctx - the context
 *         time - the time step
 *         labels - draw labels flag.
 */
static void	render_vslices (Context ctx, int time, int labels, int animflag)
{
    Display_Context	dtx;
    struct vslice	*vslice;
    int			var, box, stipple, lwidth, show_hidden;
    unsigned int	color;

    dtx = ctx->dpy_ctx;

    for (var = 0; var < ctx->NumVars; var++)
    {
        vslice = ctx->Variable[var]->VSliceTable[time];

        if (ctx->DisplayVSlice[var] && vslice && vslice->valid) 
        {
            if (render_get_read_lock (animflag, &(vslice->lock)))
            {
                /* copy data that isn't protected by the lock */
                color   = render_get_color (ctx, VSLICE, var);
                stipple = ctx->Variable[var]->VSliceRequest->stipple;
                lwidth  = ctx->Variable[var]->VSliceRequest->linewidth;
                box     = dtx->DisplayBox && !dtx->CurvedBox;

#ifndef USE_SYSTEM_FONTS
                if (labels)
                {
                    if (check_view_side (ctx, VSLICE, var) < 0)
                    {
                        flip_vslice_end_for_end (ctx, time, var);
                    }
                }
#endif

                recent (ctx, VSLICE, var);
                set_color (color);
                set_line_width (lwidth);
#ifdef HAVE_OPENGL
                if (stipple != VIS5D_SOLID_LINE)
                {
                    glEnable (GL_LINE_STIPPLE);
                    glLineStipple (1, (GLushort) stipple);
                }
#endif

                /* draw main contour lines */
                draw_disjoint_lines (vslice->num1, (void*) vslice->verts1,
                                     color, NULL, 0);

                show_hidden = !labels;
#ifdef USE_SYSTEM_FONTS
                show_hidden = 1;
#endif
                if (show_hidden)
                {
                    /* draw hidden contour lines */
                    draw_disjoint_lines (vslice->num2, (void*) vslice->verts2,
                                         color, NULL, 0);
                }
#ifdef HAVE_OPENGL
                glDisable (GL_LINE_STIPPLE);
#endif

                /* draw contour labels */
                if (labels)
                {
#ifdef USE_SYSTEM_FONTS
                    plot_strings (vslice->num3, vslice->labels,
                                  vslice->verts3, color, 
                                  dtx->gfx[CONTOUR_LABEL_FONT]->fontbase);
#else
                    draw_disjoint_lines (vslice->num3, (void*) vslice->verts3,
                                         color, NULL, 0);
#endif
                }

                set_line_width (dtx->LineWidth);

                if (box)
                {
                    /* draw the bounding box */
                    polyline ((void*) vslice->boxverts, vslice->numboxverts);

                    /* draw position label */
                    draw_vslice_ticks (dtx, vslice->r1, vslice->c1,
                                       vslice->r2, vslice->c2);
                }

                done_read_lock (&(vslice->lock));
            }
        }
    }


    return;
}
/* end PM */




/* 09Feb06  Phil McDonald	cleaned things up a little */
/*
 * Render all horizontal colored slices selected for display.
 * Input:  ctx - the context
 *         time - the timestep
 *         trans - transparency flag: 1=only draw opaque slices
 *                                 0=only draw transparent slices
 *         animflag - 1=animating, 0=not animating
 */
static void	render_chslices (Context ctx, int time, int trans, int animflag)
{
    Display_Context	dtx;
    struct chslice	*chslice;
    int			var, box;
    unsigned int	color;
    unsigned int*	coltab;
  
    dtx = ctx->dpy_ctx;

    for (var = 0; var < ctx->NumVars; var++)
    {
        chslice = ctx->Variable[var]->CHSliceTable[time];

        if (ctx->DisplayCHSlice[var] && chslice && chslice->valid)
        {
            if (render_get_read_lock (animflag, &(chslice->lock)))
            {
                /* copy data that isn't protected by the lock */
                coltab = render_get_colortable (ctx, VIS5D_CHSLICE_CT, var);
                color  = render_get_color (ctx, CHSLICE, var);
                box    = dtx->DisplayBox && !dtx->CurvedBox;

                recent (ctx, CHSLICE, var);

                if (!trans)
                {
                    draw_color_quadmesh (chslice->rows, chslice->columns,
                                         (void*) chslice->verts,
                                         chslice->color_indexes,
                                         coltab, 0, NULL, 0);
                }
                else
                {
                    if (box)
                    {
                        /* draw position label */
                        set_color (color);
                        draw_hslice_tick (dtx, chslice->level);
                    }
                }

                done_read_lock (&(chslice->lock));
            }
        }
    }


    return;
}
/* end PM */





/* 05Feb08  Phil McDonald	cleaned things up a little */
/*
 * Render all horizontal colored slices selected for display.
 * Input:  ctx - the context
 *         time - the timestep
 *         trans - transparency flag: 1=only draw opaque slices
 *                                 0=only draw transparent slices
 *         animflag - 1=animating, 0=not animating
 */
static void	render_cvslices (Context ctx, int time, int trans, int animflag)
{
    Display_Context	dtx;
    struct cvslice	*cvslice;
    int			var, box;
    unsigned int	color;
    unsigned int*	coltab;
  
    dtx = ctx->dpy_ctx;

    for (var = 0; var < ctx->NumVars; var++)
    {
        cvslice = ctx->Variable[var]->CVSliceTable[time];

        if (ctx->DisplayCVSlice[var] && cvslice && cvslice->valid)
        {
            if (render_get_read_lock (animflag, &(cvslice->lock)))
            {
                /* copy data that isn't protected by the lock */
                coltab = render_get_colortable (ctx, VIS5D_CVSLICE_CT, var);
                color  = render_get_color (ctx, CVSLICE, var);
                box    = dtx->DisplayBox && !dtx->CurvedBox;

                recent (ctx, CVSLICE, var);

                if (!trans)
                {
                    draw_color_quadmesh (cvslice->rows, cvslice->columns,
                                         (void*) cvslice->verts,
                                         cvslice->color_indexes,
                                         coltab, 0, NULL, 0);
                }
                else
                {
                    if (box)
                    {
                        /* draw position label */
                        set_color (color);
                        draw_vslice_ticks (dtx, cvslice->r1, cvslice->c1,
                                           cvslice->r2, cvslice->c2);
                    }
                }

                done_read_lock (&(cvslice->lock));
            }
        }
    }


    return;
}
/* end PM */




/* 03Feb06  Phil McDonald	cleaned things up a little */
/*
 * Render all horizontal wind vector slices which are selected for display.
 * Input:  ctx - the context
 *         time - the timestep
 */
static void	render_hwind_slices (Context ctx, int time, int animflag)
{
    Display_Context	dtx;
    struct hwind	*hwind;
    int			w, sfc, box;
    unsigned int	color;

    dtx = ctx->dpy_ctx;

    for (w = 0; w < VIS5D_WIND_SLICES; w++)
    {
        hwind = &(dtx->HWindTable[w][time]);

        if (dtx->DisplayHWind[w] && hwind->valid &&
            (ctx->context_index == dtx->Uvarowner[w]))
        {
            if (render_get_read_lock (animflag, &(hwind->lock)))
            {
                /* copy data that isn't protected by the lock */
                color = dtx->HWindColor[w];
                sfc   = dtx->DisplaySfcHWind[w];
                box   = dtx->DisplayBox && !dtx->CurvedBox;

                recent (ctx, HWIND, w);
                set_color (color);

                /* draw wind vectors */
                if (hwind->barbs || sfc)
                {
                    draw_disjoint_lines (hwind->nvectors,
                                         (void*) hwind->verts, color, NULL, 0);
                }
                else
                {
                    draw_wind_lines ((hwind->nvectors / 4),
                                     (void*) hwind->verts, color);
                }

                if (box && !sfc)
                {
                    /* draw the bounding box */
                    polyline ((void*) hwind->boxverts, hwind->numboxverts);

                    /* draw position label */
                    draw_hslice_tick (dtx, hwind->level);
                }

                done_read_lock (&(hwind->lock));
            }
        }
    }


    return;
}
/* end PM */



/* 08Feb06  Phil McDonald	cleaned things up a little */
/*
 * Render all vertical wind vector slices which are selected for display.
 * Input:  ctx - the context
 *         time - the timestep
 */
static void	render_vwind_slices (Context ctx, int time, int animflag)
{
    Display_Context	dtx;
    struct vwind	*vwind;
    int			w, box;
    unsigned int	color;

    dtx = ctx->dpy_ctx;

    for (w = 0; w < VIS5D_WIND_SLICES; w++)
    {
        vwind = &(dtx->VWindTable[w][time]);

        if (dtx->DisplayVWind[w] && vwind->valid &&
            (ctx->context_index == dtx->Uvarowner[w]))
        {
            if (render_get_read_lock (animflag, &(vwind->lock)))
            {
                /* copy data that isn't protected by the lock */
                color = dtx->VWindColor[w];
                box   = dtx->DisplayBox && !dtx->CurvedBox;

                recent (ctx, VWIND, w);
                set_color (color);

                /* draw wind vectors */
                if (vwind->barbs)
                {
                    draw_disjoint_lines (vwind->nvectors,
                                         (void*) vwind->verts, color, NULL, 0);
                }
                else
                {
                    draw_wind_lines ((vwind->nvectors / 4),
                                     (void*) vwind->verts, color);
                }

                /* draw the bounding box */
                polyline ((void*) vwind->boxverts, vwind->numboxverts);

                if (box)
                {
                    /* draw position labels */
                    draw_vslice_ticks (dtx, vwind->r1, vwind->c1,
                                       vwind->r2, vwind->c2);
                }

                done_read_lock (&(vwind->lock));
            }
        }
    }


    return;
}
/* end PM */



/* 08Feb06  Phil McDonald	cleaned things up a little */
/*
 * Render all horizontal stream vector slices which are selected for display.
 * Input:  ctx - the context
 *         time - the timestep
 */
static void	render_hstream_slices (Context ctx, int time, int animflag)
{
    Display_Context	dtx;
    struct hstream	*hstream;
    int			w, sfc, box;
    unsigned int	color;

    dtx = ctx->dpy_ctx;

    for (w = 0; w < VIS5D_WIND_SLICES; w++)
    {
        hstream = &(dtx->HStreamTable[w][time]);

        if (dtx->DisplayHStream[w] && hstream->valid &&
            (ctx->context_index == dtx->Uvarowner[w]))
        {
            if (render_get_read_lock (animflag, &(hstream->lock)))
            {
                /* copy data that isn't protected by the lock */
                color = dtx->HStreamColor[w];
                sfc   = dtx->DisplaySfcHStream[w];
                box   = dtx->DisplayBox && !dtx->CurvedBox;

                recent (ctx, HSTREAM, w);
                set_color (color);

                /* draw the stream lines */
                draw_disjoint_lines (hstream->nlines, (void*) hstream->verts,
                                     color, NULL, 0);

                if (box && !sfc)
                {
                    /* draw the bounding box */
                    polyline ((void*) hstream->boxverts, hstream->numboxverts);

                    /* draw position label */
                    draw_hslice_tick (dtx, hstream->level);
                }

                done_read_lock (&(hstream->lock));
            }
        }
    }


    return;
}
/* end PM */



/* 08Feb06  Phil McDonald	cleaned things up a little */
/*
 * Render all vertical stream vector slices which are selected for display.
 * Input:  ctx - the context
 *         time - the timestep
 */
static void	render_vstream_slices (Context ctx, int time, int animflag)
{
    Display_Context	dtx;
    struct vstream	*vstream;
    int			w, box;
    unsigned int	color;

    dtx = ctx->dpy_ctx;

    for (w = 0; w < VIS5D_WIND_SLICES; w++)
    {
        vstream = &(dtx->VStreamTable[w][time]);

        if (dtx->DisplayVStream[w] && vstream->valid &&
            (ctx->context_index == dtx->Uvarowner[w]))
        {
            if (render_get_read_lock (animflag, &(vstream->lock)))
            {
                /* copy data that isn't protected by the lock */
                color = dtx->VStreamColor[w];
                box   = dtx->DisplayBox && !dtx->CurvedBox;

                recent (ctx, VSTREAM, w);
                set_color (color);

                /* draw the stream lines */
                draw_disjoint_lines (vstream->nlines,
                                     (void*) vstream->verts, color, NULL, 0);

                /* draw the bounding box */
                polyline ((void*) vstream->boxverts, vstream->numboxverts);

                if (box)
                {
                    /* draw position labels */
                    draw_vslice_ticks (dtx, vstream->r1, vstream->c1,
                                       vstream->r2, vstream->c2);
                }

                done_read_lock (&(vstream->lock));
            }
        }
    }


    return;
}
/* end PM */




static void render_trajectories( Context ctx, int it, int tf )
{
   int alpha, i, len, start;
   Display_Context dtx;

   dtx = ctx->dpy_ctx;
   for (i=0;i<dtx->NumTraj;i++) {
      struct traj *t = dtx->TrajTable[i];

      if (t->ctx_owner==ctx->context_index
          && dtx->DisplayTraj[t->group]
           && cond_read_lock(&t->lock)){

         assert( t->lock==1 );

         recent( ctx, TRAJ, t->group );
         alpha = UNPACK_ALPHA( dtx->TrajColor[t->group]);
         if ( (tf && alpha==255) || (tf==0 && alpha<255) ) {
            start = t->start[it];
            len = t->len[it];
            if (start!=0xffff && len>0) {
               if (t->kind==0) {
                  /* draw as line segments */
                  int colorvar = t->colorvar;
                  if (colorvar>=0) {
                     /* draw colored trajectory */
                     draw_colored_polylines( len,
                                     (void *) (t->verts + start*3),
                                     (void*)(t->colors + start),
                                     dtx->ColorTable[VIS5D_TRAJ_CT]->Colors[t->colorvarowner*MAXVARS+
                                     colorvar]);
                  }
                  else {
                     /* monocolored */
                     draw_polylines( len,
                                     (void *) (t->verts + start*3),
                                     dtx->TrajColor[t->group] );
                  }
               }
               else {
                  /* draw as triangle strip */
                  int colorvar = t->colorvar;
                  if (colorvar>=0) {
                     /* draw colored triangles */
                     draw_colored_triangle_strip( len,
                                       (void*)(t->verts + start*3),
                                       (void*)(t->norms + start*3),
                                       (void*)(t->colors + start),
                                       dtx->ColorTable[VIS5D_TRAJ_CT]->Colors[t->colorvarowner*MAXVARS+
                                       colorvar], alpha );
                  }
                  else {
                     /* monocolor */
                     draw_triangle_strip( len,
                                       (void*)(t->verts + start*3),
                                       (void*)(t->norms + start*3),
                                       dtx->TrajColor[t->group] );
                  }
               }
            }
         }
         done_read_lock( &t->lock );
      }
   }
}

int draw_clock_labels(int dtime, int stime, int CurTime,int NumTimes, 
							 int xpos, int ypos, XFontStruct *font,int JulianDate)
{
  char str[15];
  int i, width1, width2, width3;
  static char month[12][4] = {"Jan", "Feb", "Mar", "Apr", "May", "Jun",
										"Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};
  static int dds[24] = {31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365,
								31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366};

  i = stime;
  sprintf( str, "%02d:%02d:%02d", i/3600, (i/60)%60, i%60 );
  draw_text( xpos, ypos, str );
  
  width1 = text_width ( font, str);
  
  if (JulianDate) {
	 sprintf( str, "%7d", v5dDaysToYYDDD( dtime ));
  }
  else {
	 int iyyddd, iy, im, id, mon;
	 
	 iyyddd = v5dDaysToYYDDD( dtime );
	 iy = iyyddd / 1000;
	 id = iyyddd - (iy * 1000);
	 im = ((iy % 4) == 0 && ((iy % 100) != 0 || (iy % 400) == 0)) ? 12 : 0;
	 for (i=im; i<im+12; i++) {
		if (id <= dds[i]) {
		  mon = i-im;
		  if (mon > 0) id = id - dds[i-1];
		  break;
		}
	 }
	 sprintf(str, "%02d %s %02d", id, month[mon], iy);
  }
  draw_text( xpos, 2*ypos, str );

  width2 = text_width ( font, str);

  sprintf( str, "%d of %d", CurTime+1, NumTimes );
  draw_text( xpos, 3*ypos, str );

  width3 = text_width ( font, str);

  if(width1>width2 && width1>width3)
	 return width1;
  else if(width1>width2)
	 return width3;
  else if(width2>width3)
	 return width2;
  return width3;
}


/*
 * Draw the clock in the upper-left corner of the 3-D window.
 * Input:  ctx - the vis5d context
 *         c - the color to use.
 */
/* MJK 12.02.98 begin */
static void draw_clock( Display_Context dtx, unsigned int c )
{
   static char day[7][20] = {"Sunday", "Monday", "Tuesday", "Wednesday",
                             "Thursday", "Friday", "Saturday" };

   static float twopi = 2.0 * 3.141592;
   short pp[CLOCK_SEGMENTS+1][2];
   float ang, delta;
   float clk_size, clk_margin, clk_radius, clk_center_x, clk_center_y;
   char str[15];
   int i, time_str_width;
   int stime=0, stimeold, dtime=0, dtimeold;

   clk_size     = 4*(dtx->gfx[WINDOW_3D_FONT]->FontHeight+VSPACE);
   clk_margin   = clk_size / 16.0;
   clk_radius   = (clk_size / 2.0) - clk_margin;
   clk_center_x = clk_size / 2.0;
   clk_center_y = clk_size / 2.0;

   /* Draw the clock. */
   if (dtx->NumTimes)
      ang = twopi * (float) dtx->CurTime / (float) dtx->NumTimes;
   else
      ang = 0.0;

   pp[0][1] = clk_center_y;
   pp[0][0] = clk_center_x;
   pp[1][1] = clk_center_y - (clk_radius * cos(ang));
   pp[1][0] = clk_center_x + (clk_radius * sin(ang));
   pp[2][1] = pp[1][1] + 1;
   pp[2][0] = pp[1][0] + 1;
   pp[3][1] = pp[0][1] + 1;
   pp[3][0] = pp[0][0] + 1;
   pp[4][1] = pp[0][1] - 1;
   pp[4][0] = pp[0][0] + 1;
   pp[5][1] = pp[1][1] - 1;
   pp[5][0] = pp[1][0] + 1;

   set_color( c );
   polyline2d( pp, 6 );

   if (dtx->CircleClock){
      /* Draw a circle around the clock. */
      delta = twopi / ((float) CLOCK_SEGMENTS);
      ang   = 0.0;
      for (i = 0; i < CLOCK_SEGMENTS; i++)
      {
         pp[i][0] = clk_center_x + (clk_radius * sin (ang)) + 0.5;
         pp[i][1] = clk_center_y - (clk_radius * cos (ang)) + 0.5;
         ang += delta;
      }
      pp[i][0] = pp[0][0];
      pp[i][1] = pp[0][1];
      polyline2d (pp, CLOCK_SEGMENTS+1);
   }



   clk_size += clk_margin;

   dtimeold = -1;
   stimeold = -1;
   if(dtx->NumTimes)
	  vis5d_get_dtx_time_stamp( dtx->dpy_context_index,
										 dtx->CurTime,
										 &dtime, &stime);
	
	time_str_width = draw_clock_labels(dtime,stime, dtx->CurTime, 
												  dtx->NumTimes, clk_size, 
												  dtx->gfx[WINDOW_3D_FONT]->FontHeight+VSPACE,
												  dtx->gfx[WINDOW_3D_FONT]->font,
												  dtx->JulianDate);



   dtimeold = dtime;
   stimeold = stime;

   if (dtx->NumTimes == 1 ||
       ((dtx->Elapsed[dtx->NumTimes-1] - dtx->Elapsed[0])
         / (dtx->NumTimes - 1)) < 48*3600 ) {
     /* Print day of week */
     draw_text( clk_size, 4*(dtx->gfx[WINDOW_3D_FONT]->FontHeight+VSPACE),
                day[ (dtimeold+0) % 7 ] );
   }
   if (dtx->group_index > 0){
      sprintf( str, "   Group %d", dtx->group_index);
      draw_text( (clk_size + time_str_width), (dtx->gfx[WINDOW_3D_FONT]->FontHeight+VSPACE), str );
   }

	if(dtx->numofctxs+dtx->numofitxs>1){
	  int yo;
	  Context ctx;
	  Irregular_Context itx;
	  for(yo=0;yo<dtx->numofctxs;yo++){
		 ctx=dtx->ctxpointerarray[yo];

		 vis5d_get_ctx_time_stamp( ctx->context_index,
											ctx->CurTime,
											&dtime, &stime);

/* 05Apr06  Phil McDonald	TODO: maybe we still want this message */
/*
		 if((dtime!=dtimeold) || (stime!=stimeold))
			fprintf(stderr,"Time of ctx data %d differs from display time \n",yo);
*/
/* end PM */
	  }
	  for(yo=0;yo<dtx->numofitxs;yo++){
		 itx=dtx->itxpointerarray[yo];

		 vis5d_get_itx_time_stamp( itx->context_index,
											itx->CurTime,
											&dtime, &stime);

		 if((dtime!=dtimeold) || (stime!=stimeold)){
			time_str_width += draw_clock_labels(dtime,stime, itx->CurTime, 
														  itx->NumTimes, clk_size+time_str_width+clk_margin, 
														  dtx->gfx[WINDOW_3D_FONT]->FontHeight+VSPACE,
														  dtx->gfx[WINDOW_3D_FONT]->font,
														  dtx->JulianDate);

/* 05Apr06  Phil McDonald	TODO: maybe we still want this message */
/*
			fprintf(stderr,"Time of itx data %d differs from display time %d %d \n",yo,dtime,stime);
*/
/* end PM */
		 }
	  }

	}



}
/* MJK 12.02.98 end */



/*
 * Render all the 2-D text labels.
 */
static void render_text_labels( Display_Context dtx )
{
   struct label *lab;

  for (lab=dtx->FirstLabel; lab; lab=lab->next) {

	 if (dtx->Reversed  && lab->LabelColor==PACK_COLOR(255,255,255,255)){
      set_color( PACK_COLOR(0,0,0,255) );
	 }
	 else{
      set_color( lab->LabelColor );
	 }

    draw_text( lab->x, lab->y, lab->text );
    if (lab->state) {
       /* being edited -> draw cursor */
       short verts[4][2];
       verts[0][0] = lab->x2;     verts[0][1] = lab->y1;
       verts[1][0] = lab->x2;     verts[1][1] = lab->y2;
       verts[2][0] = lab->x2+1;   verts[2][1] = lab->y2;
       verts[3][0] = lab->x2+1;   verts[3][1] = lab->y1;
       polyline2d( verts, 4 );
    }
  }
}




static void draw_fake_pointer( Display_Context dtx )
{
   short pp[8][2];

   pp[0][0] = dtx->PointerX;       pp[0][1] = dtx->PointerY;
   pp[1][0] = dtx->PointerX+15;    pp[1][1] = dtx->PointerY+5;
   pp[2][0] = dtx->PointerX+5;     pp[2][1] = dtx->PointerY+15;
   pp[3][0] = dtx->PointerX;       pp[3][1] = dtx->PointerY;
   pp[4][0] = dtx->PointerX+20;    pp[4][1] = dtx->PointerY+20;

   polyline2d( pp, 5 );

}



/*
 * Print status info at bottom of window.
 */
static void print_info( Display_Context dtx )
{
   char str[1000];
   int size, waiters;
   size_t m;

   m = mem_used( dtx );
   get_queue_info( &size, &waiters );
   if (m>=0)
      sprintf(str, "Pending: %d   Memory Used: %d", size, m );
   else
      sprintf(str, "Pending: %d", size );

   draw_text( 10, dtx->WinHeight - dtx->gfx[WINDOW_3D_FONT]->FontHeight, str );
}



/*
 * Print the numeric value of each variable at the probe's current location.
 */
/* MJK 12.02.98 begin */
static void draw_probe( Display_Context dtx )
{
   float val;
   char str[1000];
   int y, var;
   int x;
   int yo;
   float rr,cc,ll;

   /* find widest parameter name, but only once */
   if (!dtx->do_not_recalc_probe_text_width){
      x = -1;
      for (yo = 0; yo < dtx->numofctxs; yo++){
         for (var=0;var<dtx->ctxpointerarray[yo]->NumVars; var++) {
            int w = text_width(dtx->gfx[WINDOW_3D_FONT]->font, dtx->ctxpointerarray[yo]->Variable[var]->VarName );
            int l = strlen( dtx->ctxpointerarray[yo]->Variable[var]->VarName );
            if (w < 1) w = 11 * l;
            if (w>x)
               x = w;
         }
      }
      if (dtx->numofctxs >0){
         x += 25;
      }
      dtx->do_not_recalc_probe_text_width = 1;
      dtx->probe_text_width = x;
   }
   x = dtx->probe_text_width;

   /* Draw from bottom of window upward */
   y = dtx->WinHeight - dtx->gfx[WINDOW_3D_FONT]->FontHeight;
   for (yo = 0; yo < dtx->numofctxs; yo++){
      Context ctx;
      int ipvar, npvar, lpvar;
      ctx = dtx->ctxpointerarray[yo];
      if (ctx->ProbeNumVars >= 0) {
         npvar = ctx->ProbeNumVars;
         lpvar = 1;
      }
      else {
         npvar = ctx->NumVars - 1;
         lpvar = 0;
      }
      for (ipvar = npvar; ipvar >= lpvar; ipvar--) {
         var = ipvar;
         if (ctx->ProbeNumVars > 0) {
            for (var=ctx->NumVars-1;var>=0;var--) {
               if (ctx->ProbeVar[var] == ipvar) break;
            }
         }
         if (var >= 0) {
            float r, c, l;
            xyzPRIME_to_grid( ctx, ctx->CurTime, var,
                         dtx->CursorX, dtx->CursorY, dtx->CursorZ, &rr, &cc, &ll );
            xyzPRIME_to_gridPRIME( ctx->dpy_ctx, dtx->CurTime, var, dtx->CursorX,
                                   dtx->CursorY, dtx->CursorZ, &r, &c, &l);
            if (ll < ctx->Variable[var]->LowLev || ll > ctx->Nl[var]-1 + ctx->Variable[var]->LowLev ||
                rr < 0 || rr > ctx->Nr-1 || cc < 0 || cc > ctx->Nc-1 ||
                !check_for_valid_time(ctx, dtx->CurTime)) {
               val = MISSING;
            }
            else if (dtx->CoordFlag==1) {
               /* discrete grid position */
               int row = (int) (r+0.01);
               int col = (int) (c+0.01);
               int lev = (int) (l+0.01);
               if (ctx->GridSameAsGridPRIME){
                  val = get_grid_value( ctx, ctx->CurTime, var, row, col, lev );
               }
               else{
                  vis5d_gridPRIME_to_grid(ctx->context_index, ctx->CurTime, var,
                                          (float)(row), (float)(col), (float)(lev),
                                           &rr, &cc, &ll);
                  if (ll < ctx->Variable[var]->LowLev || ll > ctx->Nl[var]-1 + ctx->Variable[var]->LowLev ||
                      rr < 0 || rr > ctx->Nr-1 || cc < 0 || cc > ctx->Nc-1) {
                     val = MISSING;
                  }
                  else{
                     val = get_grid_value( ctx, ctx->CurTime, var, rr, cc, ll);
                  }
               }
            }
            else {
               if (ctx->GridSameAsGridPRIME){
                  val = interpolate_grid_value( ctx, ctx->CurTime, var, r, c, l );
               }
               else{
                  val = interpolate_grid_value( ctx, ctx->CurTime, var, rr, cc, ll );
               }
            }
            sprintf( str, "%-4s", ctx->Variable[var]->VarName );


            if (dtx->numofctxs > 1){
               char yodle[40];
               strcpy(yodle, return_var_plus_index(str, ctx->context_index));
               draw_text( 10, y, yodle );
            }
            else{
               draw_text( 10, y, str);
            }
            if (IS_MISSING(val))
              sprintf( str, " = MISSING" );
            else
              sprintf( str, " = %.3g %s", val, ctx->Variable[var]->Units );
            draw_text( x+10, y, str );
            y -= (dtx->gfx[WINDOW_3D_FONT]->FontHeight+VSPACE);
         }
      }
   }
}
/* MJK 12.02.98 end */


/*
#define LEGEND_SPACE   20
#define LEGEND_WIDTH   25
#define LEGEND_HEIGHT 128
*/
#define TICK_LENGTH     4
#define NUM_TICKS       5


/*
 * Draw a color legend.
 * Input: var    = parameter index the color slice belongs to for which to
 *                draw the legend,
 *        type  = CHSLICE or CVSLICE,
 *        xleft = x-position for left side of color bar,
 *        ybot  = y-position for bottom of color bar.
 * Return: width of bar + numbers drawn, or height of bar + label.
 */
int draw_legend( Context ctx, int varowner, int var, int type, int xleft, int ybot )
{
   int   y, 
         lutindex,
         textwidth;
   int tick;
   short cline[2][2];
   uint  *lut;
   char  scrap[100], format[20];
   int legend_space, legend_width, legend_height;
   float label;
   Display_Context dtx;
   Context colorctx;

   /* MJK 12.02.98 */
   int bg_r, bg_g, bg_b, bg_a, fg_r, fg_g, fg_b, fg_a;

   dtx = ctx->dpy_ctx;
   colorctx = dtx->ctxpointerarray[return_ctx_index_pos(dtx, varowner)];

   legend_height = dtx->LegendSize;
   legend_width = (25 * dtx->LegendSize) / 128;
   legend_space = (20 * dtx->LegendSize) / 128;

   switch(type) {
      case VIS5D_ISOSURF:
         lut = dtx->ColorTable[VIS5D_ISOSURF_CT]->Colors[varowner*MAXVARS+var];
         break;
      case VIS5D_CHSLICE:
         lut = dtx->ColorTable[VIS5D_CHSLICE_CT]->Colors[varowner*MAXVARS+var];
         break;
      case VIS5D_CVSLICE:
         lut = dtx->ColorTable[VIS5D_CVSLICE_CT]->Colors[varowner*MAXVARS+var];
         break;
      case VIS5D_TRAJ:
         lut = dtx->ColorTable[VIS5D_TRAJ_CT]->Colors[varowner*MAXVARS+var];
         break;
      case VIS5D_VOLUME:
         lut = dtx->ColorTable[VIS5D_VOLUME_CT]->Colors[varowner*MAXVARS+var];
         break;
      case VIS5D_TOPO:
         lut = dtx->ColorTable[VIS5D_TOPO_CT]->Colors[varowner*MAXVARS+var];
         break;
      default:
         /* this should never happen */
         abort();
   }

   /* These line values never change */
   cline[0][0] = xleft;
   cline[1][0] = xleft + legend_width;

   /* Draw the colors */

   /* MJK 12.02.98 begin */

   bg_r = UNPACK_RED(dtx->BgColor);
   bg_g = UNPACK_GREEN(dtx->BgColor);
   bg_b = UNPACK_BLUE(dtx->BgColor);

   for (y = 0; y<legend_height; y++) {
      lutindex = (255 * y)/legend_height;
      cline[0][1] = cline[1][1] = ybot - y;

      fg_r = UNPACK_RED(lut[lutindex]);
      fg_g = UNPACK_GREEN(lut[lutindex]);
      fg_b = UNPACK_BLUE(lut[lutindex]);
      fg_a = UNPACK_ALPHA(lut[lutindex]);
      bg_a = 255 - fg_a;
      fg_r = ((fg_r * fg_a) + (bg_r * bg_a)) / 255;
      fg_g = ((fg_g * fg_a) + (bg_g * bg_a)) / 255;
      fg_b = ((fg_b * fg_a) + (bg_b * bg_a)) / 255;

      set_color (PACK_COLOR(fg_r, fg_g, fg_b, 255));
      polyline2d(cline, 2);

   }

   /* Draw a box around the legend colorbar */
      /* MJK 3.29.99 */
   if (dtx->Reversed) {
      set_color(PACK_COLOR( 0, 0, 0, 255));
   }
   else{
      set_color( dtx->BoxColor );
   }
   cline[0][0] = cline[1][0] = xleft;
   cline[0][1] = ybot;
   cline[1][1] = ybot - legend_height + 1;
   polyline2d(cline, 2);
   cline[0][0] = cline[1][0] = xleft + legend_width;
   polyline2d(cline, 2);
   cline[0][0] = xleft;
   cline[1][0] = xleft + legend_width;
   cline[0][1] = cline[1][1] = ybot;
   polyline2d(cline, 2);
   cline[0][1] = cline[1][1] = ybot - legend_height + 1;
   polyline2d(cline, 2);

   /* MJK 12.02.98 end */


   /* Determine largest value physical variable can have */
   label = ABS(colorctx->Variable[var]->MaxVal);
   if (ABS(colorctx->Variable[var]->MinVal) > label)
      label = ABS(colorctx->Variable[var]->MinVal);

   /* Create 'pretty' formatting string */
   sprintf(scrap, "% .0f", label); 
   sprintf(format, "%% %d.2f", (int)strlen(scrap)+3);

   /* Draw values and tick marks on the right hand side of the legend */
   textwidth = 0;

   /* Make sure we have a tick at the top of the legend @@ */

   cline[0][0] += TICK_LENGTH + legend_width;
      /* MJK 3.29.99 */
   if (dtx->Reversed) {
      set_color(PACK_COLOR( 0, 0, 0, 255));
   }
   else{
      set_color( dtx->BoxColor );
   }
XSync( GfxDpy, 0 );
   for (tick=0;tick<NUM_TICKS;tick++) {
      int ticky, texty;
      float value;

      ticky = ybot - tick * (legend_height-1) / (NUM_TICKS-1);
      texty = ybot - tick * (legend_height-dtx->gfx[WINDOW_3D_FONT]->FontHeight+dtx->gfx[WINDOW_3D_FONT]->FontDescent)
              / (NUM_TICKS-1);
      value = colorctx->Variable[var]->MinVal + (colorctx->Variable[var]->MaxVal-colorctx->Variable[var]->MinVal)*tick/4.0;

      cline[0][1] = cline[1][1] = ticky;
      polyline2d(cline, 2);
      
      sprintf(scrap, format, value);
      draw_text( xleft + legend_width + TICK_LENGTH + 2, texty, scrap );
XSync( GfxDpy, 0 );
      if (text_width(dtx->gfx[WINDOW_3D_FONT]->font,scrap) > textwidth)
         textwidth = text_width(dtx->gfx[WINDOW_3D_FONT]->font,scrap);
   }

   /* Print name of physical variable above legend */
   if (colorctx->Variable[var]->Units[0]) {
      sprintf( scrap, "%s (%s)", colorctx->Variable[var]->VarName, colorctx->Variable[var]->Units );
      draw_text( xleft, ybot - legend_height - dtx->gfx[WINDOW_3D_FONT]->FontDescent-2, scrap );
   }
   else {
      draw_text( xleft, ybot - legend_height - dtx->gfx[WINDOW_3D_FONT]->FontDescent-2,
                 colorctx->Variable[var]->VarName);
   }

   if (dtx->LegendPosition == VIS5D_BOTTOM ||
       dtx->LegendPosition == VIS5D_TOP) {
     return legend_width + TICK_LENGTH + 5 + textwidth + legend_space;
   }
   else {
     return legend_height + 5 + dtx->gfx[WINDOW_3D_FONT]->FontHeight + legend_space;
   }
}


/*
 * Draws color legends of activated color slices.
 * Since the space in the 3-D window is restricted only one row of legends
 * is drawn in the bottom of the 3-D window. The order of drawing is:
 * first the legends of the horizontal slices for parameter 0..NumVars-1
 * and then the vertical slices for parameter 0..NumVars-1.
 */
static void draw_color_legends( Display_Context dtx )
{
   int var, set;
   int left;      /* Left x position of current legend */
   int bottom;    /* Bottom y position of current legend */ 
   int inc, vert;
   int vindex;
   int ctxnum, cvar, cvowner;
   Context ctx;


   if (dtx->LegendPosition == VIS5D_BOTTOM) {
     left = 50+dtx->LegendMarginX;
     bottom  = dtx->WinHeight - 20 + dtx->LegendMarginY;
     vert = 0;
   }
   else if (dtx->LegendPosition == VIS5D_TOP) {
     left = 200+dtx->LegendMarginX;
     bottom  = dtx->LegendSize + 5 + dtx->gfx[WINDOW_3D_FONT]->FontHeight + 25 + dtx->LegendMarginY;
     vert = 0;
   }
   else if (dtx->LegendPosition == VIS5D_RIGHT) {
     left = dtx->WinWidth - ((35 * dtx->LegendSize) / 128) - 5 * dtx->gfx[WINDOW_3D_FONT]->FontHeight+dtx->LegendMarginX;
     bottom  = dtx->LegendSize + 5 + dtx->gfx[WINDOW_3D_FONT]->FontHeight + 50 + dtx->LegendMarginY;
     vert = 1;
   }
   else if (dtx->LegendPosition == VIS5D_LEFT) {
     left = 20+dtx->LegendMarginX;
     bottom  = dtx->LegendSize + 5 + dtx->gfx[WINDOW_3D_FONT]->FontHeight + 100 + dtx->LegendMarginY;
     vert = 1;
   }
   else {
     printf("draw_color_legends: bad LegendPosition\n");
     return;
   }

   /* Isosurface color tables */
   for (ctxnum=0; ctxnum<dtx->numofctxs; ctxnum++){
      ctx = dtx->ctxpointerarray[ctxnum];
      for (var=0; var<ctx->NumVars; var++) {
         cvar = ctx->IsoColorVar[var];
         cvowner = ctx->IsoColorVarOwner[var];
         if (ctx->DisplaySurf[var] && cvar>=0 ){
            /* Draw legend at position (xstart, ystart) = upper left corner */
            inc = draw_legend( ctx, cvowner, cvar, VIS5D_ISOSURF, left, bottom ); 
            if (vert) {
              bottom += inc;
              if (bottom > dtx->WinHeight - 50) return;
            }
            else {
              left += inc;
              if (left > dtx->WinWidth - 150) return;
            }
         }
      }
   }

   /* Find activated horizontal color slices */
   for (ctxnum=0; ctxnum<dtx->numofctxs; ctxnum++){
      ctx = dtx->ctxpointerarray[ctxnum];
      for (var=0; var<ctx->NumVars; var++) {
         if (ctx->DisplayCHSlice[var]) {
            /* Draw legend at position (xstart, ystart) = upper left corner */
            vindex = ctx->context_index;
            inc = draw_legend( ctx, vindex, var, VIS5D_CHSLICE, left, bottom ); 
            if (vert) {
              bottom += inc;
              if (bottom > dtx->WinHeight - 50) return;
            }
            else {
              left += inc;
              if (left > dtx->WinWidth - 150) return;
            }
         }
      }
   }

   /* Find activated vertical color slices */
   for (ctxnum=0; ctxnum<dtx->numofctxs; ctxnum++){
      ctx = dtx->ctxpointerarray[ctxnum];
      for (var=0; var<ctx->NumVars; var++) {
         if (ctx->DisplayCVSlice[var]) {
            /* Draw legend at position (xstart, ystart) = upper left corner */
            vindex = ctx->context_index;
            inc = draw_legend( ctx, vindex, var, VIS5D_CVSLICE, left, bottom ); 
            if (vert) {
              bottom += inc;
              if (bottom > dtx->WinHeight - 50) return;
            }
            else {
              left += inc;
              if (left > dtx->WinWidth - 150) return;
            }
         }
      }
   }

   /* Volume */
   if (dtx->CurrentVolume>=0  ) {
      /* Draw legend at position (xstart, ystart) = upper left corner */
      inc = draw_legend( ctx, dtx->CurrentVolumeOwner,
                         dtx->CurrentVolume, VIS5D_VOLUME, left, bottom ); 
      if (vert) {
        bottom += inc;
        if (bottom > dtx->WinHeight - 50) return;
      }
      else {
        left += inc;
        if (left > dtx->WinWidth - 150) return;
      }
   }

   /* Trajectory color tables */
   for (set=0; set<VIS5D_TRAJ_SETS; set++) {
      int cvar = dtx->TrajColorVar[set];
      int cvowner = dtx->TrajColorVarOwner[set];
      if (dtx->DisplayTraj[set] && cvar>=0 ){ 
         /* Draw legend at position (xstart, ystart) = upper left corner */
         inc = draw_legend( ctx, cvowner, cvar, VIS5D_TRAJ, left, bottom ); 
         if (vert) {
           bottom += inc;
           if (bottom > dtx->WinHeight - 50) return;
         }
         else {
           left += inc;
           if (left > dtx->WinWidth - 150) return;
         }
      }
   }

   /* Topo color table */
   if (dtx->topo->TopoColorVar>=0 && dtx->topo->DisplayTopo) {
      int cvar = dtx->topo->TopoColorVar;
      int cvowner = dtx->topo->TopoColorVarOwner;
      /* Draw legend at position (xstart, ystart) = upper left corner */
      inc = draw_legend( ctx, cvowner, cvar, VIS5D_TOPO, left, bottom ); 
      if (vert) {
        bottom += inc;
        if (bottom > dtx->WinHeight - 50) return;
      }
      else {
        left += inc;
        if (left > dtx->WinWidth - 150) return;
      }
   }

}




/*
 * Draw anything the user wants in 3D.
 * Drawing bounds are (Xmin,Ymin,Zmin) - (Xmax,Ymax,Zmax)
 */
static void draw_user_3d_graphics( Display_Context dtx )
{
}



/*
 * Draw anything the user wants in 2D.
 * Drawing bounds are (0,0) - (Width-1, Height-1), origin in upper-left corner.
 */
static void draw_user_2d_graphics( Display_Context dtx )
{
}


/*
 * Only draw the 3-D elements of the scene.  No matrix, viewport, etc
 * operations are done here.  This function is useful for the CAVE
 * since it controls the viewing parameters.
 * Input:  ctx - the context
 *         animflag - 1=animating, 0=not animating
 */
void render_3d_only( Display_Context dtx, int animflag )
{
   int yo,  labels, i;
   Context ctx;
   Irregular_Context itx;

/* 16Jan06  Phil McDonald */
/*
   set_3d (dtx->GfxProjection, dtx->FrntClip, dtx->Zoom, (float*) dtx->CTM);
*/

/* 23Jan06  Phil McDonald	TODO: clear a GL error from somewhere */
   set_current_window (dtx);
   glGetError ();
/* end PM */

   if (animflag){
      labels = !dtx->ContnumFlag;
   }
   else{
      labels = dtx->ContnumFlag;
   }

   /* Loop over antialiasing passes */
   for (i=0; i < (dtx->PrettyFlag ? AA_PASSES : 1); i++) {

      start_aa_pass(i);

      /*** Draw 3-D lines ***/

/* 16Jan06  Phil McDonald */
      if (dtx->DisplayBox)
      {
         draw_box (dtx, dtx->CurTime);
      }
/* end PM */

      clipping_off();

      if (dtx->DisplayCursor)
      {
         if (dtx->DisplayProbe)
         {
            if (dtx->Reversed){
               draw_cursor (dtx, 0, dtx->CursorX, dtx->CursorY, dtx->CursorZ,
                            PACK_COLOR(0,0,0,255));
            }
            else{
               draw_cursor (dtx, 0, dtx->CursorX, dtx->CursorY, dtx->CursorZ,
                            dtx->BoxColor);
            }
         }
         else if (dtx->DisplaySound)
         {
            if (dtx->Reversed){
               draw_cursor (dtx, 2, dtx->CursorX, dtx->CursorY, 0,
                            PACK_COLOR(0,0,0,255));
            }
            else{
               draw_cursor (dtx, 2, dtx->CursorX, dtx->CursorY, 0,
                            dtx->BoxColor);
            }

         }
         else
         {
            draw_cursor (dtx, dtx->RibbonFlag,
                         dtx->CursorX, dtx->CursorY, dtx->CursorZ,
                         *dtx->CursorColor);
         }

         if (dtx->DisplayBox)
         {
            print_cursor_position (dtx, dtx->CurTime);
         }
      }

      clipping_on();

      for (yo = 0; yo < dtx->numofitxs; yo++){
         itx = dtx->itxpointerarray[yo];
         render_textplots( itx, itx->CurTime);
      }

      for (yo= 0; yo < dtx->numofctxs; yo++){
         ctx = dtx->ctxpointerarray[yo];
         if (check_for_valid_time(ctx, dtx->CurTime)){
			  render_hslices( ctx, ctx->CurTime, labels, animflag );

			  render_vslices( ctx, ctx->CurTime, labels, animflag );

			  render_hwind_slices( ctx, ctx->CurTime, animflag );

			  render_vwind_slices( ctx, ctx->CurTime, animflag );

			  render_hstream_slices( ctx, ctx->CurTime, animflag );
				
			  render_vstream_slices( ctx, ctx->CurTime, animflag );
         }   
      }


      /* draw user graphics */
      draw_user_3d_graphics( dtx );



      /*** Draw opaque 3-D graphics ***/
      /* MJK 12.02.98 begin */
      set_depthcue( dtx->DepthCue );
      if (dtx->topo && dtx->topo->TopoFlag && dtx->topo->DisplayTopo) {
         set_depthcue(0);
         draw_topo( dtx, dtx->CurTime, dtx->DisplayTexture, 0 );
      }
      else if (dtx->DisplayTexture) {
         /* just draw flat textured image */
         set_depthcue(0);
         draw_topo( dtx, dtx->CurTime, 1, 1 );
      }

      if (dtx->MapFlag && dtx->DisplayMap) {
         if (dtx->DisplaySfcMap) {
            set_color( dtx->DarkMapColor );
            draw_map( dtx, dtx->CurTime, 0 );
         }
         else {
            set_color( dtx->LightMapColor );
            draw_map( dtx, dtx->CurTime, 1 );
         }
         if (dtx->Reversed){
            set_color( PACK_COLOR(0,0,0,255) );
         }
         else{
            set_color( dtx->BoxColor );
         }

      }
      set_depthcue(0);

      for (yo= 0; yo < dtx->numofctxs; yo++){      
         ctx = dtx->ctxpointerarray[yo];       
         if (check_for_valid_time(ctx, dtx->CurTime)){
			  
			  if(ctx->DataGridList){
				 glCallList(ctx->DataGridList);
			  }
			  
            render_trajectories( ctx, ctx->CurTime, 1);

            render_isosurfaces( ctx, dtx->CurTime, ctx->CurTime, 1, animflag );

            render_chslices( ctx, ctx->CurTime, 1, animflag );

            render_cvslices( ctx, ctx->CurTime, 1, animflag );
         }      
      }
      /*** Draw transparent 3-D objects ***/

      for (yo= 0; yo < dtx->numofctxs; yo++){
         ctx = dtx->ctxpointerarray[yo];
         if (check_for_valid_time(ctx, dtx->CurTime)){            
            render_trajectories( ctx, ctx->CurTime, 0);
            render_isosurfaces( ctx, dtx->CurTime, ctx->CurTime, 0, animflag );
            render_chslices( ctx, ctx->CurTime, 0, animflag );
            render_cvslices( ctx, ctx->CurTime, 0, animflag );
         }            
      }

      if (dtx->VolumeFlag==1 && dtx->CurrentVolume!=-1){
         ctx = dtx->ctxpointerarray[return_ctx_index_pos(dtx,
                                    dtx->CurrentVolumeOwner)];
         
         if (check_for_valid_time(ctx, dtx->CurTime)){                     
            draw_volume( ctx, ctx->CurTime, dtx->CurrentVolume,
                      dtx->ColorTable[VIS5D_VOLUME_CT]->Colors[ctx->context_index*MAXVARS+
                       dtx->CurrentVolume] );
         }
      }

      end_aa_pass(i);

   } /* aa passes */
}



/*
 * Only draw the 2-D elements of the scene.  No matrix, viewport, etc
 * operations are done here.
 * Input:  ctx - the context
 */
void render_2d_only( Display_Context dtx )
{

   if (dtx->DisplayClock) {
   /* MJK 3.29.99 */
      if (dtx->Reversed){
         draw_clock( dtx, PACK_COLOR(0,0,0,255) );
         draw_logo( dtx, PACK_COLOR(0,0,0,255) );
      }
      else{
         draw_clock( dtx, dtx->BoxColor );
         draw_logo( dtx, dtx->BoxColor );
      }
   }
   if (dtx->DisplayInfo) {
      print_info(dtx);
   }
   if (dtx->DisplayProbe) {
      draw_probe(dtx);
   }
   if (dtx->DisplayCursor && dtx->DisplayBox) {
      print_cursor_position( dtx, dtx->CurTime );
   }

   if (dtx->PointerX>=0 && dtx->PointerY>=0){
      draw_fake_pointer(dtx);  /* for remote widget mode */
   }

   /* MJK 3.29.99 */
   /*  Moved to render_text_labels  JPE 
   if (dtx->Reversed){
      set_color( PACK_COLOR(0,0,0,255) );
   }
   else{
      set_color( dtx->LabelColor );
   }
	*/

   render_text_labels(dtx);

   /* Draw color map legends of color slices (as much as fit in window) */
   if (dtx->DisplayLegends ) {
      draw_color_legends(dtx);
   }

   draw_user_2d_graphics( dtx );
}

/*
 *
 * draw everything that belongs inside of
 * the Sounding Graphics Window
 *
 *
 */
void render_sounding_only( Display_Context dtx, int pixmapflag )
{  
   if (dtx->DisplaySound) {
      if ( pixmapflag == 1 ){
         do_pixmap_art( dtx );
         /* MJK 12.02.98 */
         /*
         draw_sounding(dtx, dtx->CurTime);
         */
      }
      if ((dtx->CursorX != dtx->Sound.currentX || 
           dtx->CursorY != dtx->Sound.currentY) ||
           (dtx->CurTime != dtx->Sound.currentTime) ||
           /* MJK 12.02.98 */ (pixmapflag)){
         if (dtx->CurTime != dtx->Sound.currentTime){
            reload_sounding_data( dtx );
         }
         draw_sounding(dtx, dtx->CurTime); 
         dtx->Sound.currentX = dtx->CursorX;
         dtx->Sound.currentY = dtx->CursorY;
         dtx->Sound.currentTime = dtx->CurTime;
      } 
   }
}



/*
 * Redraw everything in the 3-D window but don't display it yet.  Call
 * swap_3d_window() to do that.
 * Input:  ctx - the vis5d context
 *         animflag - 1=animating, 0=notanimating
 * Return:  0 = ok, -1 = error.
 */
void render_everything( Display_Context dtx, int animflag )
{

   if (get_frame(dtx, dtx->CurTime)) {
      return;
   }

   /*** Draw 3-D Objects ***/
   set_3d( dtx->GfxProjection, dtx->FrntClip,
           dtx->Zoom, (float*) dtx->CTM);
   if(dtx->StereoOn){
      /*
       * The set_3d call should have set the model matrix and the
       * projection matrix correctly.  The stereo_set_3d_perspective
       * call just changes the projection matrix slightly.
       */

      /* left eye */
      stereo_set_3d_perspective(VIS5D_STEREO_LEFT,dtx->FrntClip);
      stereo_set_buff(VIS5D_STEREO_LEFT);
      clear_color (dtx->BgColor);
      clear_3d_window();

      if (dtx->DisplayBox){
         int i;
         for (i=0; i < (dtx->PrettyFlag ? AA_PASSES : 1); i++) {
            start_aa_pass(i);
            draw_box(dtx, dtx->CurTime);
            /* draw_tick_marks( dtx ); */
            end_aa_pass(i);
         }
      }

      clipping_on();
      render_3d_only( dtx, animflag );
      clipping_off();
   
      if (dtx->DisplayClips){
         render_vclips( dtx, animflag );
         render_hclips( dtx, animflag );
      }

      /* right eye */
      stereo_set_3d_perspective(VIS5D_STEREO_RIGHT,dtx->FrntClip);
      stereo_set_buff(VIS5D_STEREO_RIGHT);
      clear_color (dtx->BgColor);
      clear_3d_window();

      if (dtx->DisplayBox){
         int i;
         for (i=0; i < (dtx->PrettyFlag ? AA_PASSES : 1); i++) {
            start_aa_pass(i);
            draw_box(dtx, dtx->CurTime);
            /* draw_tick_marks( dtx ); */
            end_aa_pass(i);
         }
      }

      clipping_on();
      render_3d_only( dtx, animflag );
      clipping_off();
   
      if (dtx->DisplayClips){
         render_vclips( dtx, animflag );
         render_hclips( dtx, animflag );
      }

      /* reset to default draw buffer */
      stereo_set_buff(VIS5D_STEREO_BOTH);
   }else{

		 /*
		  * The clear_3d_window call was moved out of vis5d_draw_frame and moved
		  * down into this function so stereo and mono rendering loops can both
		  * work with the same higher level api function.
		  */
		 clear_3d_window();
		 clear_color (dtx->BgColor);
		 
		 if (dtx->DisplayBox){
			int i, listflag=0;
			  
			  for (i=0; i < (dtx->PrettyFlag ? AA_PASSES : 1); i++) {
				 start_aa_pass(i);
				 draw_box(dtx, dtx->CurTime); 
				 /* draw_tick_marks( dtx ); */
				 end_aa_pass(i);
			  }
		 }

		 clipping_on();
		 render_3d_only( dtx, animflag );

		 clipping_off();

		 if (dtx->DisplayClips){
			render_vclips( dtx, animflag );
			render_hclips( dtx, animflag );
		 }
	  }

/* 02Feb06  Phil McDonald */
   if (!dtx->CaveMode)
   {
/* end PM */
   /*** Draw 2-D objects ***/
   set_2d();
   render_2d_only( dtx );

   /*** Draw Sounding ***/
   render_sounding_only( dtx , 0); 
/* 17Jan06  Phil McDonald */
   }
/* end PM */
	
   if (dtx->AnimRecord) {
	  save_frame( dtx, dtx->CurTime ); 
   }
   finish_rendering();
}
