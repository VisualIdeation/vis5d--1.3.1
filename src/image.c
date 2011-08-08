/*  image.c */


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
/* Texture/image mapping */



#include <assert.h>
#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#ifdef HAVE_SGI_GL
#  include <X11/Xlib.h>
#  include <gl/gl.h>
#endif
#ifdef HAVE_OPENGL
#  include <GL/gl.h>
#  include <GL/glu.h>
#endif
#include "globals.h"
#include "graphics.h"
#include "rgb.h"
/* 01Feb06  Phil McDonald */
#include "memory.h"
/* end PM */



#ifdef HAVE_SGI_GL
/*
 * Lighting and texturing properties:
 */
static float mat2[] = {
   AMBIENT, 0.5, 0.5, 0.5,   /* surface ambient color */
   DIFFUSE, 0.9, 0.9, 0.9,   /* surface diffuse color */
   LMNULL
};
static float lmodel2[] = {
   AMBIENT, 0.5, 0.5, 0.5,         /* ambient color */
   LOCALVIEWER, 0.0,               /* infinite viewpoint */
   TWOSIDE, 0,
   LMNULL
};
static float texprops[] = { TX_MINFILTER, TX_POINT,
                            TX_MAGFILTER, TX_POINT,
                            TX_NULL };
static float tevprops[] = { TV_MODULATE, TV_NULL };
static float subdiv_params[] = { 0.0, 0.0, 10.0 };
#endif /* HAVE_SGI_GL */



#ifdef HAVE_MCIDAS
/*
 * McIDAS stuff
 */
/* API - not to globals.c:
   only one call to init_mcidas_kludge for all Vis5D contexts */
int uc[1000];
int neguc[200];
int ttyfd[10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

#define OPNARA F77_FUNC(opnara,OPNARA)
#define READD F77_FUNC(readd,READD)
#define RHELP F77_FUNC(rhelp,RHELP)
#define CLSARA F77_FUNC(clsara,CLSARA)
#define REDARA F77_FUNC(redara,REDARA)
#define KLTWIN F77_FUNC(kltwin,KLTWIN)
#define LWINIT F77_FUNC(lwinit,LWINIT)

extern void OPNARA( int * );
extern void READD( int *, int * );
extern void RHELP( int *, int * );
extern void CLSARA( int * );
extern void REDARA( int *, int*, int *, int *, int *, void *);

#ifndef MCIDAS_SIDECAR
extern void KLTWIN( void );
extern void LWINIT( void );

static void init_mcidas_kludge( void ) /* API - call once for all contexts */
{
  int i;

  for (i=0; i<1000; i++) uc[i] = 0;
  for (i=0; i<200; i++) neguc[i] = 0;
  KLTWIN();
  LWINIT();
}
#endif /* end of ifndef MCIDAS_SIDECAR */


void F77_FUNC(zdest,ZDEST)( char *c, int *n, int len )
{
  int l;
  char s[100];

  l = (len < 100) ? len : 99;
  strncpy(s, c, l);
  s[l] = 0;
  printf("%s %d\n", s, *n);
}



/*
 * Load a mcidas area and store it in packed ABGR format.
 * Input:  area_num - number of mcidas area
 * Output:  width, height - size of image read
 *          data - pointer to image data
 * Return:  1 = ok, 0 = error.
 */
static int read_area( Display_Context dtx, int area_num,
                      int *width, int *height, unsigned char **data )
{
   int idir[64], bytes, e, band;
   unsigned char data1[80000];
   int x, y;
   int h, w;
   unsigned char *d;
   int ix;

#ifdef NAV
  float dx, dy;
  int qr, qc, i, j, status;
  float x, y, z;
  float lat, lon, hgt;
  float xline, xelem, xdum, aline, aelem, lpct, epct;
  static float zero = 0.0;
  int line, elem, lres, eres;

  /* get area navigation */
  if ((status = F77_FUNC(nvset,NVSET)("AREA        ", &area_num, 12)) != 0) {
    if (status == -1) printf("no area %d\n", area_num);
    else if (status == -2) printf("gen_nav: no navigation for area %d\n", area_num);
    else printf("gen_nav: navigation problem for area %d %d\n", area_num, status);
    return 0; /* **** or just don't navigate area **** */
  }
#endif

   /*printf("area_num %d\n", area_num);*/
   READD(&area_num, idir);
   if (idir[0] < 0) {
     /* unable to read area */
     printf("Error:  couldn't read AREA%04d, no areas read.\n",area_num);
     return 0;
   }

   OPNARA(&area_num);

#ifdef NAV
  line = idir[5] + line * idir[11]; /* 5 - area line, 11 - area line res */
  elem = idir[6] + elem * idir[12]; /* 6 - area elem, 12 - area elem res */

  if (lres == 0) {
    printf("gen_nav: bad line res %d\n", lres);
    return 0;
  }
  else if (lres > 0) lres = idir[11] * lres;
  else lres = idir[11] / (-lres);

  if (eres == 0) {
    printf("gen_nav: bad elem res %d\n", eres);
    return 0;
  }
  else if (eres > 0) eres = idir[12] * eres;
  else eres = idir[12] / (-eres);
#endif

   *height = h = idir[8];
/*   if (h > 768) *height = h = 768; */
   *width = w = idir[9];
   bytes = idir[10];
/* 01Feb06  Phil McDonald */
   d = (unsigned char *) mem_malloc (w * h * sizeof (unsigned char));
/* end PM */
   *data = d;

   e = 0;
   band = (bytes == 1) ? 0 : 8;
   /*printf("h %d w %d bytes %d band %d\n", h, w, bytes, band);*/

   /* new area stuff */
   bytes = 1;
   RHELP(&area_num, &bytes);

#ifdef NAV
  qr = dtx->qrows;
  qc = dtx->qcols;
  dx = (dtx->Xmax-dtx->Xmin) / (float) (qc-1);
  dy = (dtx->Ymax-dtx->Ymin) / (float) (qr-1);

  y = dtx->Ymax;
  for (i=0; i<qr; i++) {
    x = dtx->Xmin;
    for (j=0; j<qc; j++) {
      xyzPRIME_to_geo( dtx, -1, -1, x, y, 0.0, &lat, &lon, &hgt );

      if(F77_FUNC(nv1eas,NVLEAS)(&lat, &lon, &zero, &xline, &xelem, &xdum) == 0) {
        aline = (xline - line) / lres;
        aelem = (xelem - elem) / eres;
        GET min and max line and elem
        USE to calculate area sector
        TRANSFORM area coordinates to 0.0 to 1.0 coordinates
      }
      else {

      }
      x += dx;
    }
    y -= dy;
  }
#endif

   ix = 0;
   for (y=0;y<h;y++) {
      /* int yy = (h-y-1) * w; */
      if (bytes == 1) {
         REDARA(&area_num, &y, &e, width, &band, data1);
         for (x=0; x<w; x++) {
            d[ix++] = data1[x];
         }
      }
      else {
#ifdef LEAVEOUT
         REDARA(&area_num, &y, &e, width, &band, data2);
         for (x=0; x<w; x+=3) {
            d[ix++] = data2[x] / 2;
            d[ix++] = data2[x+1] / 2;
            d[ix++] = data2[x+2] / 2;
            d[ix++] = 0;;
         }
#endif
      }
   }

   CLSARA(&area_num);
   return 1;
}

#endif /* ifdef HAVE_MCIDAS */



/*
 * Find the value nearest to n which is also a power of two.
 */
static int round2( int n )
{
   int m;

   for (m=1; m<n; m*=2)
     ;

   /* m>=n */
   if (m-n <= n-m/2) {
      return m;
   }
   else {
      return m/2;
   }
}



/**********************************************************************/
/**********************************************************************/
/*****                    Public Functions                        *****/
/**********************************************************************/
/**********************************************************************/



/* 20Oct06  Phil McDonald */
void	*scale_image (int format, void *image, int *p_width, int *p_height)
{

    int			width_in, height_in, width_out, height_out, max;
    unsigned int	*image_in, *image_out;



#ifdef HAVE_OPENGL
    /* The texture size *MUST* be a power of two.  Rescale if needed. */
    width_in  = *p_width;
    height_in = *p_height;

    glGetIntegerv (GL_MAX_TEXTURE_SIZE, &max);
    width_out  = round2 (width_in);
    height_out = round2 (height_in);
    if (width_out  > max)  width_out  = max;
    if (height_out > max)  height_out = max;

    if ((width_out == width_in) && (height_out == height_in)) return image;

    fprintf (stderr, "Rescale texture from %d x %d to %d x %d. %d %s.\n",
             width_in, height_in, width_out, height_out,
             max, "is the GL max dim");
    image_in = image;
    image_out = (unsigned int *) mem_malloc (width_out * height_out * 4);
    gluScaleImage ((GLenum) format,
                   (GLint) width_in, (GLint) height_in,
                   GL_UNSIGNED_BYTE, (void *) image_in,
                   (GLint) width_out, (GLint) height_out,
                   GL_UNSIGNED_BYTE, (void *) image_out);

    mem_free (image_in);

    *p_width  = width_out;
    *p_height = height_out;
#endif  /* ifdef HAVE_OPENGL */


    return image_out;
}
/* end PM */



/*
 * Init global variables.
 */
void init_image( Display_Context dtx )
{
   int	i;
   void	*prev_image;

   prev_image = NULL;
   for (i=0;i<dtx->NumTimes;i++) {
      dtx->TexWidth[i] = dtx->TexHeight[i] = 0;
      dtx->TexComponents[i] = 0;
/* 01Feb06  Phil McDonald */
      if (dtx->TexImage[i] != NULL)
      {
         if (prev_image != dtx->TexImage[i]) mem_free (dtx->TexImage[i]);
         prev_image = dtx->TexImage[i];
      }
/* end PM */
      dtx->TexImage[i] = NULL;
      dtx->TexImageNew[i] = 1; 
   }
}



/*
 * Define the texture to be used for the given timestep.  The image data
 * is not copied; just the pointer to it is saved.
 */
void define_texture( Display_Context dtx, int time, int width, int height,
                     int components, void *image )
{
   assert( time>=0 && time<=dtx->NumTimes );

   dtx->TexWidth[time] = width;
   dtx->TexHeight[time] = height;
   dtx->TexComponents[time] = components;
/* 01Feb06  Phil McDonald */
   if (dtx->TexImage[time] != NULL) mem_free (dtx->TexImage[time]);
/* end PM */
   dtx->TexImage[time] = image;
}



/*
 * Read a texture map from the named file.  The image can be any SGI
 * .rgb file.
 * Input:  filename - name of image file
 * Return:  1 = ok, 0 = error.
 */
int read_texture_image( Display_Context dtx, char *filename )
{
   RGB_IMAGE *img;
   int width, height;
   int i;
   unsigned int *image;

   img = ReadRGB( filename );
   if (!img)
      return 0;

   width = img->sizeX;
   height = img->sizeY;
/* 10Oct06  Phil McDonald	why limit our texture size? */
#ifdef LEAVEOUT
   if (width>1024) {
      /* too wide */
      FreeRGB( img );
      return 0;
   }
#endif
/* end PM */

/* 23Oct06  Phil McDonald */
    image = (unsigned int *) img->data;

    image = scale_image (GL_RGBA, image, &width, &height);
/* end PM */

   /* use same texture for all timesteps */
   for (i=0;i<dtx->NumTimes;i++) {
      define_texture( dtx, i, width, height, 4, image );
   }

/* 23Oct06  Phil McDonald	img->data can't be freed here because it  */
/*				was either freed in scale_image or put in */
/*				dtx->TexImage, so just free img           */
    mem_free (img);
/* end PM */
   return 1;
}



/* WLH 11-25-94
 * Read a sequence of images from a file, to use for texture mapping.
 * Input:  name - name of file.
 * Return:  1 = ok, 0 = error.
 */
int read_texture_sequence( Display_Context dtx, char *name )
{
  int i, length, fd, head[3];
  int width, height, max;
  unsigned char *data;

  if ((fd = open(name, O_RDONLY, 0)) == -1) {
    return 0; /* cannot open file */
  }
  length = 3 * sizeof(int);
  if (read(fd, head, length) != length) {
    return 0; /* cannot read file header */
  }
  if (head[0] < dtx->NumTimes) {
    return 0; /* not enough time steps in file */
  }

  for (i=0;i<dtx->NumTimes;i++) {
    height = head[1];
    width = head[2];

    length = width * height * sizeof(unsigned char);
/* 01Feb06  Phil McDonald */
    data = (unsigned char *) mem_malloc (length);
/* end PM */

    if (read(fd, data, length) != length) {
      return 0; /* cannot read image from file */
    }

#ifdef HAVE_OPENGL
/* 23Oct06  Phil McDonald */
    data = scale_image (GL_LUMINANCE, data, &width, &height);
    check_gl_error( "read_texture_sequence" );
/* end PM */
#endif

    define_texture( dtx, i, width, height, 1, data );

  }
  return 1;
}



#ifdef HAVE_MCIDAS
/*
 * Read a sequence of McIDAS area files to use for texture mapping.
 * Input:  first - number of first AREA file in the sequence.
 * Return:  1 = ok, 0 = error.
 */
int read_texture_areas( Display_Context dtx, int first )
{
   int i;
   int width, height;
   unsigned int *image;

#ifndef MCIDAS_SIDECAR
   init_mcidas_kludge();
#endif

   for (i=0;i<dtx->NumTimes;i++) {
      if (!read_area( dtx, first+i, &width, &height,
                      (unsigned char **) &image )) {
         /* error */
         return 0;
      }
#ifdef HAVE_OPENGL
/* 23Oct06  Phil McDonald */
      image = scale_image (GL_LUMINANCE, image, &width, &height);
      check_gl_error( "read_texture_areas" );
/* end PM */
#endif
      define_texture( dtx, i, width, height, 1, image );
   }
  return 1;
}
#endif  /* ifdef HAVE_MCIDAS */




/*
 * Specify which texture to use.  If time==-1, disable texturing.
 */
int use_texture( Display_Context dtx, int time )
{

/* 25Oct06  Phil McDonald	cave and multi-proc need local data */
   static int		tex_id		= 0;
   static int		prev_time	= -1;
   static char		*new		= NULL;



   assert( time>=0 && time<dtx->NumTimes );


   /* Do one-time initializations */
/* 25Oct06  Phil McDonald */
   if (tex_id == 0)
   {
      if (new == NULL)
      {
          new = (char *) malloc (MAXTIMES * sizeof (char));
          memset (new, 1, MAXTIMES);
      }

#ifdef HAVE_SGI_GL
      /* define special material and lighting model for texturing */
      lmdef( DEFMATERIAL, 11, 0, mat2 );
      lmdef( DEFLMODEL, 31, 0, lmodel2 );
      Tex = 1;
#endif
#ifdef HAVE_OPENGL
      glGenTextures (1, &(tex_id));
      glEnable (GL_TEXTURE_2D);
      glTexEnvf (GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
      glBindTexture (GL_TEXTURE_2D, tex_id);
      glHint (GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
      glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
#endif
   } 
/* end PM */

   /* Give the texture to the graphics library / hardware */
/* 25Oct06  Phil McDonald */
   if (dtx->TexImage[time])
   {
#ifdef HAVE_SGI_GL
      if ((dtx->prev_time == -1) ||
          (dtx->TexImage[dtx->prev_time] != dtx->TexImage[time]) ||
          (dtx->TexImageNew[time] == 1))
      {
      if ((dtx->TexImageNew[time] == 1) || (dtx->prev_time == -1) ||
          (dtx->TexImage[dtx->prev_time] != dtx->TexImage[time]))
      {
         texdef2d (1, dtx->TexComponents[time],
                   dtx->TexWidth[time], dtx->TexHeight[time],
                   (unsigned int *) dtx->TexImage[time],
                   5, texprops);
         tevdef (1, 0, tevprops);
      }
#endif
#ifdef HAVE_OPENGL
      glEnable (GL_TEXTURE_2D);
      glBindTexture (GL_TEXTURE_2D, tex_id);

      if ((new[time] == 1) || (prev_time < 0) ||
          (dtx->TexImage[prev_time] != dtx->TexImage[time]))
      {
         glHint (GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
         glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
         glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
/* 08Nov06  Phil McDonald */
         if (dtx->TexComponents[time] == 1)
         {
            glTexImage2D (GL_TEXTURE_2D, 0,
                          1, dtx->TexWidth[time], dtx->TexHeight[time],
                          0, GL_LUMINANCE, GL_UNSIGNED_BYTE,
                          dtx->TexImage[time]);
         }
         else
         {
            glTexImage2D (GL_TEXTURE_2D, 0, 4,
                          dtx->TexWidth[time], dtx->TexHeight[time],
                          0, GL_RGBA, GL_UNSIGNED_BYTE, dtx->TexImage[time]);
         }
      }
/* end PM */
#endif
      dtx->TexImageNew[time] = 0;
      new[time] = 0;
   }
/* end PM */

   dtx->prev_time = time;
   prev_time = time;

   return 0;
}



/*
 * Enable or disable texture mapping.
 */
static void enable_texture( int enable )
{

   
   if (enable) {
#ifdef HAVE_SGI_GL
      texbind( TX_TEXTURE_0, 1 );
      tevbind( TV_ENV0, 1 );
      scrsubdivide( SS_DEPTH, subdiv_params );
      /* special lighting parameters */
      lmbind( MATERIAL, 11 );
      lmbind( LMODEL, 31 );
      /*lmcolor( LMC_AD );     Added on 4-11-95 by BEP */
      cpack( 0xffffffff );
#endif
#ifdef HAVE_OPENGL
      glEnable( GL_LIGHTING );
      glColor4f( 1.0, 1.0, 1.0, 1.0 );
      check_gl_error( "enable_texture (glEnable)" );
#endif
   }
   else {
      /* disable */
#ifdef HAVE_SGI_GL
      texbind( TX_TEXTURE_0, 0 );
      tevbind( TV_ENV0, 0 );
      scrsubdivide( SS_OFF, subdiv_params );
      /* regular lighting */
      lmbind( MATERIAL, 10 );
      lmbind( LMODEL, 30 );
#endif
#ifdef HAVE_OPENGL
      glDisable( GL_TEXTURE_2D );
      glDisable( GL_LIGHTING );
      check_gl_error( "enable_texture (glDisable)" );
#endif
   }
}



/* 24Oct07  Phil McDonald	optimized a bit */
/*
 * Draw a texture mapped polytriangle strips using the currently defined
 * texture.
 */
int draw_texture (Display_Context disp_ctx, float vert[][3], float norm[][3],
                  float texcoord[][2])
{
    int		i, j, nr, nc, nrc, row0, row1;



    nr  = disp_ctx->topo->qrows;
    nc  = disp_ctx->topo->qcols;
    nrc = nr * nc;

    /* turn on texture mapping */
    enable_texture(1);

    /* break mesh into strips */
#ifdef HAVE_SGI_GL
    for (i=0;i<nr-1;i++) {
        row0 = i * nc;
        row1 = (i+1) * nc;
        if (norm) {
            bgnqstrip();
            for (j=0;j<nc;j++) {

                t2f( texcoord[row0+j] );
                n3f( norm[row0+j] );
                v3f( vert[row0+j] );

                t2f( texcoord[row1+j] );
                n3f( norm[row1+j] );
                v3f( vert[row1+j] );
            }
            endqstrip();
        }
        else {
            /* no normals */
            float n[3];
            n[0] = 0.0;  n[1] = 0.0;  n[2] = 1.0;
            n3f(n);
            bgnqstrip();
            for (j=0;j<nc;j++) {
                t2f( texcoord[row0+j] );
                v3f( vert[row0+j] );
                t2f( texcoord[row1+j] );
                v3f( vert[row1+j] );
            }
            endqstrip();
        }
    }
#endif
#ifdef HAVE_OPENGL
    glPushMatrix ();
    if (norm != NULL)
    {
        /* Push the texture down a little so that data displayed on the */
        /* topo can be seen.  The texture will "peek through" anywhere  */
        /* topo data is transparent.                                    */

        glTranslatef (0.0, 0.0, -disp_ctx->topo->z_fuzz);

        for (row0 = 0, row1 = nc; row1 < nrc; row0 = row1, row1 += nc)
        {
            glBegin (GL_TRIANGLE_STRIP);
            for (j = 0; j < nc; j++)
            {
                glTexCoord2fv (texcoord[row1+j]);
                glNormal3fv (norm[row1+j]);
                glVertex3fv (vert[row1+j]);

                glTexCoord2fv (texcoord[row0+j]);
                glNormal3fv (norm[row0+j]);
                glVertex3fv (vert[row0+j]);
            }
            glEnd();
        }
    }
    else
    {
        /* no normals */
        for (row0 = 0, row1 = nc; row1 < nrc; row0 = row1, row1 += nc)
        {
            glBegin (GL_TRIANGLE_STRIP);
            for (j = 0; j < nc; j++)
            {
                glTexCoord2fv (texcoord[row1+j]);
                glVertex3fv (vert[row1+j]);

                glTexCoord2fv (texcoord[row0+j]);
                glVertex3fv (vert[row0+j]);
            }
            glEnd();
        }
    }
    glPopMatrix ();
#endif

    /* turn off texture mapping */
    enable_texture(0);


    return 0;
}
/* end PM */
