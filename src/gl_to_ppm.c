
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

#if defined(HAVE_OPENGL) || defined(HAVE_SGI_GL)

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "api.h"
#include "globals.h"
#include "graphics.h"
#include "gl_to_ppm.h"

#if HAVE_OPENGL
#  include <GL/gl.h>
#  include <GL/glu.h>
#  include <GL/glx.h>
#elif HAVE_SGI_GL
#  include <gl/gl.h>
#endif
/* 01Feb06  Phil McDonald */
#include "memory.h"
/* end PM */

static unsigned char *dataR;
static unsigned char *dataG;
static unsigned char *dataB;
static int current_x_offset;
static int current_y_offset;
static int big_x;
static int big_y;
static FILE *f;


static int write_ppm_val( int number )
{
   int val, valer;

   if (number > 999){
      valer = number / 1000;
      number -= valer*1000;
      val = valer + '0';
      if (fputc(val, f)==EOF){
         printf("Error: could not write to output file?\n");
         exit(0);
      }
      valer = number / 100;
      number -= valer*100;
      val = valer + '0';
      if (fputc(val, f)==EOF){
         printf("Error: could not write to output file?\n");
         exit(0);
      }
      valer = number / 10;
      number -= valer*10;
      val = valer + '0';
      if (fputc(val, f)==EOF){
         printf("Error: could not write to output file?\n");
         exit(0);
      }
      val = number + '0';
      if (fputc(val, f)==EOF){
         printf("Error: could not write to output file?\n");
         exit(0);
      }
      return 1;
   }
   if (number > 99){
      valer = number / 100;
      number -= valer*100;
      val = valer + '0';
      if (fputc(val, f)==EOF){
         printf("Error: could not write to output file?\n");
         exit(0);
      }
      valer = number / 10;
      number -= valer*10;
      val = valer + '0';
      if (fputc(val, f)==EOF){
         printf("Error: could not write to output file?\n");
         exit(0);
      }
      val = number + '0';
      if (fputc(val, f)==EOF){
         printf("Error: could not write to output file?\n");
         exit(0);
      }
      return 1;
   }
   else if (number > 9){
      valer = number / 10;
      number -= valer*10;
      val = valer + '0';
      if (fputc(val, f)==EOF){
         printf("Error: could not write to output file?\n");
         exit(0);
      }
      val = number + '0';
      if (fputc(val, f)==EOF){
         printf("Error: could not write to output file?\n");
         exit(0);
      }
      return 1;
   }
   else if (number > -1){
      val = number + '0';
      if (fputc(val, f)==EOF){
         printf("Error: could not write to output file?\n");
         exit(0);
      }
      return 1;
   }
   else{
      printf("Error: trying to write a negative number to a ppm file?\n");
      exit(0);
   }
}



static void free_pixel_data( void )
{
/* 01Feb06  Phil McDonald */
   TMP_FREE (dataR);
   TMP_FREE (dataG);
   TMP_FREE (dataB);
/* end PM */
}



int open_ppm_file( char *filename, int x, int y)
{
   /**************************/   
   /* alloc the pixel memory */
   /**************************/
/* 01Feb06  Phil McDonald */
   dataR = (unsigned char *) TMP_MALLOC (x*y*sizeof(char));
   if (!dataR){
      printf("Could not allocate memory to save ppm file\n");
      return 0;
   }
   dataG = (unsigned char *) TMP_MALLOC (x*y*sizeof(char));
   if (!dataG){
      printf("Could not allocate memory to save ppm file\n");
      TMP_FREE (dataR);
      return 0;
   }
   dataB = (unsigned char *) TMP_MALLOC (x*y*sizeof(char));
   if (!dataB){
      printf("Could not allocate memory to save ppm file\n");
      TMP_FREE (dataR);
      TMP_FREE (dataG);
      return 0;
   }
/* end PM */

   /****************************/   
   /*open the file for writing */
   /****************************/
   f = fopen( filename, "w");
   if (!f){
      printf("Could not open %s for writing\n", filename);
      return 0;
   }
   fseek(f, 0, SEEK_SET);


   /********************************/   
   /* write out simple header info */
   /********************************/
   if (fputc('P', f)==EOF){
      printf("Could not write to output file %s\n", filename);
      return 0;
   }
   if (fputc('6', f)==EOF){
      printf("Could not write to output file %s\n", filename);
      return 0;
   }
   if (fputc('\n', f)==EOF){
      printf("Could not write to output file %s\n", filename);
      return 0;
   }
   if (!write_ppm_val( x)){
      printf("Could not write to output file %s\n", filename);
      return 0;
   }
   if (fputc(' ', f)==EOF){
      printf("Could not write to output file %s\n", filename);
      return 0;
   }
   if (!write_ppm_val( y)){
      printf("Could not write to output file %s\n", filename);
      return 0;
   }
   if (fputc('\n', f)==EOF){
      printf("Could not write to output file %s\n", filename);
      return 0;
   }
   if (!write_ppm_val( 255)){
      printf("Could not write to output file %s\n", filename);
      return 0;
   }
   if (fputc('\n', f)==EOF){
      printf("Could not write to output file %s\n", filename);
      return 0;
   }
   
   /********************/   
   /* init some values */
   /********************/
   current_x_offset = 0;
   current_y_offset = 0;
   big_x = x;
   big_y = y;

   return 1;
}


   
int add_display_to_ppm_file(Display_Context dtx, int position)
{
   unsigned char *tempdataR;
   unsigned char *tempdataG;
   unsigned char *tempdataB;
   int x, y, xpos, ypos, p, q;
   int bx, by;

   set_current_window( dtx );
   glPixelStorei(GL_PACK_ALIGNMENT, 1);
   check_gl_error( "add_display_to_ppm_file (glPixelStore)" );
   
   /******************************************************/
   /*alloc the memory to retrieve data from frame buffer */
   /******************************************************/
/* 01Feb06  Phil McDonald */
   tempdataR = (unsigned char *) TMP_MALLOC
               (dtx->WinWidth * dtx->WinHeight*sizeof(unsigned char));
   if (!tempdataR){
      printf("Could not allocate enough memory to create ppm file\n");
      free_pixel_data();
      return 0;
   }
   tempdataG = (unsigned char *) TMP_MALLOC
               (dtx->WinWidth * dtx->WinHeight*sizeof(unsigned char));
   if (!tempdataG){
      printf("Could not allocate enough memory to create ppm file\n");
      free_pixel_data();
      TMP_FREE (tempdataR);
      return 0;
   }
   tempdataB = (unsigned char *) TMP_MALLOC
               (dtx->WinWidth * dtx->WinHeight*sizeof(unsigned char));
   if (!tempdataB){   
      printf("Could not allocate enough memory to create ppm file\n");   
      free_pixel_data();   
      TMP_FREE (tempdataR);
      TMP_FREE (tempdataG);
      return 0;   
   }   
/* end PM */


   glReadPixels( 0, 0, dtx->WinWidth, dtx->WinHeight, GL_RED, GL_UNSIGNED_BYTE, tempdataR);
   glReadPixels( 0, 0, dtx->WinWidth, dtx->WinHeight, GL_GREEN, GL_UNSIGNED_BYTE, tempdataG);
   glReadPixels( 0, 0, dtx->WinWidth, dtx->WinHeight, GL_BLUE, GL_UNSIGNED_BYTE, tempdataB);
   check_gl_error( "add_display_to_ppm_file (glReadPixels)" );

   /*************************************************/
   /* put data from frame buffer data into big data */
   /*************************************************/
   for (y = 0; y < dtx->WinHeight; y++){
      for (x = 0; x < dtx->WinWidth; x++){
         /* frame buffer data starts */
         /* from lower left corner */
         p = ((dtx->WinHeight-1-y)*dtx->WinWidth)+x;
         bx = current_x_offset + x;
         by = current_y_offset + y;
         q = by * big_x + bx;
         dataR[q] = tempdataR[p];
         dataG[q] = tempdataG[p];
         dataB[q] = tempdataB[p];
      }
   }
         

   /*********************/
   /* figure out offset */
   /**************************************/
   /* assuming raster order for displays!*/
   /**************************************/
   ypos = position/DisplayCols;
   xpos = position - (DisplayCols*ypos);
   if (xpos == DisplayCols-1){
      current_x_offset = 0;
      current_y_offset += dtx->WinHeight;
   }
   else{
      current_x_offset += dtx->WinWidth;
   }

/* 01Feb06  Phil McDonald */
   TMP_FREE (tempdataR);
   TMP_FREE (tempdataG);
   TMP_FREE (tempdataB);
/* end PM */
   tempdataR = NULL;
   tempdataG = NULL;
   tempdataB = NULL;

   return 1;
}

int close_ppm_file(void)
{
   int i;
 
   /*****************************************/
   /* tranfer data from dataR, dataG, dataB */
   /* to file                               */
   /*****************************************/
   for (i=0; i<big_x * big_y; i++){
      fputc(dataR[i], f);
      fputc(dataG[i], f);
      fputc(dataB[i], f);
   }
   fputc(EOF, f);
   if (fclose(f)==EOF){
      printf("Could not close output ppm file\n");
      free_pixel_data();
      return 0;
   }
   printf("Done writing ppm image file.\n");
   return 1;
}

#endif /* defined(HAVE_OPENGL) || defined(HAVE_SGI_GL) */
