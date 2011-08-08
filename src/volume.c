/* volume.c */


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
 * Volume rendering.  Requires alpha blending capability.  This is
 * supported in on SGI hardware with VGX or higher graphics.  Otherwise,
 * it's done in software on some other systems.
 *
 * Volume rendering is also supported on everthing running OpenGL!
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "globals.h"
#include "graphics.h"
#include "grid.h"
#include "memory.h"
#include "proj.h"
#include "volume.h"

#if HAVE_OPENGL
#  include <GL/gl.h>
#elif HAVE_SGI_GL
#  include <gl/gl.h>
#endif

/* 25May06  Phil McDonald */
#include "linterp.h"
/* 11Aug06  Phil McDonald */
#include <string.h>
/* 11Sep06  Phil McDonald */
#include <sys/time.h>
/* 22Oct07  Phil McDonald */
#include <sys/types.h>
#include <unistd.h>
/* end PM */



#define ABS(A)  ( (A) < 0 ? -(A) : (A) )

/* 15Feb07  Phil McDonald */
#define VOLUME_X_DIR_FLAG			4
#define VOLUME_Y_DIR_FLAG			2
#define VOLUME_Z_DIR_FLAG			1
/* 28Mar06  Phil McDonald */
#define VOLUME_FUZZ				1.0e-3

/* 13Sep06  Phil McDonald */
#undef TIME_TEST
/* end PM */



/* 09Aug06  Phil McDonald */
struct VolData
{
/* 30Jun06  Phil McDonald	The following data are view-specific and,  */
/*				therefore, can not be shared in a CAVE     */
/*				environment.  In a non-CAVE environment    */
/*				it doesn't hurt to have them here instead  */
/*				of members of the volume data structure.   */

/* 27Jun06  Phil McDonald */
    int			dir;

/* 10Aug06  Phil McDonald */
    float		vu_vec[3];
    Context		data_ctx;

/* 26Feb07  Phil McDonald */
    int			nquad;
    unsigned int	*quad_ivert;

/* 27Feb07  Phil McDonald */
    int			vu_face;
};

struct ProcData
{
    int			proc_id;	/* process id */
    struct VolData	*vol_data;	/* volume data */
    struct ProcData	*next;		/* pointer to next in the list */
};
/* end PM */



/* 11Aug06  Phil McDonald */
struct VolData	*volume_data_for_current_process (Context data_ctx)
{
    int			pid;
    void		**p_tail;
    struct ProcData	*proc_data;
    struct VolData	*vol_data;
    struct volume	*v;



    v = data_ctx->Volume;

    pid       = getpid ();


    /*  See if volume data for the current process already exists. */

    wait_read_lock (&(v->lock));

    p_tail    = &(v->private_data);
    proc_data = v->private_data;
    vol_data  = NULL;
    while (proc_data != NULL)
    {
        if (pid == proc_data->proc_id)
        {
            vol_data = proc_data->vol_data;
            break;
        }
        p_tail    = (void **) &(proc_data->next);
        proc_data = proc_data->next;
    }

    done_read_lock (&(v->lock));

    if (vol_data != NULL) return vol_data;


    /*  The current process wasn't found in the list so create a hunk of */
    /*  volume data for this process and add it to the end of the list.  */

    wait_write_lock (&(v->lock));

    proc_data = (struct ProcData *) allocate (data_ctx,
                                              sizeof (struct ProcData));
    vol_data  = (struct VolData *)  allocate (data_ctx,
                                              sizeof (struct VolData));
    vol_data->dir        = -1;
    vol_data->data_ctx   = data_ctx;
    vol_data->quad_ivert = NULL;
    vol_data->nquad      = 0;

    proc_data->vol_data  = vol_data;
    proc_data->proc_id   = pid;
    proc_data->next      = NULL;

    *p_tail              = proc_data;

    done_write_lock (&(v->lock));


    return vol_data;
}



/* 11Aug06  Phil McDonald */
int	free_private_data (Context data_ctx)
{

    struct volume	*v;
    struct VolData	*vol_data;
    struct ProcData	*proc_data, *proc_next;



    v = data_ctx->Volume;

    proc_data = v->private_data;
    while (proc_data != NULL)
    {
        proc_next = proc_data->next;
        vol_data  = proc_data->vol_data;
        if (vol_data != NULL)
        {
            if ((vol_data->quad_ivert != NULL) && (vol_data->nquad > 0))
            {
                deallocate (data_ctx, vol_data->quad_ivert,
                            (vol_data->nquad * 4 * sizeof (unsigned int)));
            }
            deallocate (data_ctx, vol_data, sizeof (struct VolData));
        }
        deallocate (data_ctx, proc_data, sizeof (struct ProcData));
        proc_data = proc_next;
    }
    v->private_data = NULL;


    return 1;
}



/* 23May06  Phil McDonald */
void	free_volume_data (Context data_ctx)
{
    int			nxyz;
    struct volume	*v;



    if ((data_ctx == NULL) || (data_ctx->Volume == NULL)) return;

    v    = data_ctx->Volume;
    nxyz = v->nx * v->ny * v->nz;

    v->valid = 0;

    if (v->vertex != NULL)
    {
        deallocate (data_ctx, v->vertex, (nxyz * 3 * sizeof (float)));
        v->vertex = NULL;
    }

    if (v->index != NULL)
    {
        deallocate (data_ctx, v->index, (nxyz * sizeof (uint_1)));
        v->index = NULL;
    }

    free_private_data (data_ctx);


    return;
}



/* 23May06  Phil McDonald */
void	free_volume (Context data_ctx)
{
    if ((data_ctx == NULL) || (data_ctx->Volume == NULL)) return;

    free_volume_data (data_ctx);

    deallocate (data_ctx, data_ctx->Volume, sizeof (struct volume));
    data_ctx->Volume = NULL;


    return;
}



/* 16Feb07  Phil McDonald */
int	alloc_volume_data (Context data_ctx)
{
    size_t		nxyz, nvert, nind;
    struct volume	*v;



    if ((data_ctx == NULL) || (data_ctx->Volume == NULL)) return 0;


    v    = data_ctx->Volume;
    nxyz = v->nx * v->ny * v->nz;

    if (v->vertex != NULL) free_volume_data (data_ctx);


    nvert     = nxyz * 3 * sizeof (float);
    v->vertex = (float *) allocate (data_ctx, nvert);

    nind     = nxyz * sizeof (uint_1);
    v->index = (uint_1 *) allocate (data_ctx, nind);

    if ((v->vertex == NULL) || (v->index == NULL))
    {
        printf ("%s%d%s\n",
                "WARNING:  insufficient memory for volume rendering (",
                (nvert + nind), "bytes needed)");

        free_volume_data (data_ctx);
        data_ctx->dpy_ctx->VolRender = 0;

        return 0;
    }
    fprintf (stderr, "VOLUME ALLOC: %d x %d x %d = %d bytes\n",
             v->nx, v->ny, v->nz, (nvert + nind));

    v->valid = 0;
    v->ivar  = -1;
    v->itime = -1;


    return 1;
}




/*
 * Allocate a new volume structure and return a pointer to it.  If we
 * can't do volume rendering return NULL.
 * Input:  nr, nc, nl - dimensions of largest volume we'll encounter.
 * Return:  pointer to volume struct or NULL
 */
struct volume *alloc_volume( Context ctx, int nr, int nc, int nl )
{
   int volren = 0;
   struct volume *v = NULL;

   if (ctx->dpy_ctx->Projection == PROJ_CYLINDRICAL ||
       ctx->dpy_ctx->Projection == PROJ_SPHERICAL){
      ctx->dpy_ctx->VolRender = 0;
      return 0;
   }
   if (nl <= 1){
      ctx->dpy_ctx->VolRender = 0;
      return 0; /* no volume variables */
   }

#if defined (HAVE_SGI_GL) || defined (DENALI)
   alphabits = getgdesc( GD_BLEND ); 
   if (alphabits==0) {
      ctx->dpy_ctx->VolRender = 0; 
      /* no alpha support means no volume rendering */
      return NULL;
   }
   volren = 1;
#endif

#if defined(HAVE_OPENGL)
   volren = 1;
#endif

/* MJK 12.15.98 */
#ifdef HAVE_PEX
   volren = 1;
#endif

   
   if (volren) {
/* 11Aug06  Phil McDonald */
      v         = (struct volume *) allocate (ctx, sizeof (struct volume));
      v->itime  = -1;
      v->ivar   = -1;
      v->lock   = 0;
      v->valid  = 0;
      v->vertex = NULL;
      v->index  = NULL;

      v->private_data = NULL;

/* 15Feb07  Phil McDonald */
      v->nx     = 0;
      v->ny     = 0;
      v->nz     = 0;
/* end PM */

   }

   ctx->dpy_ctx->VolRender = (v != NULL);


   return v;
}



/* 15Feb07  Phil McDonald */
int	volume_dir (struct VolData *vol)
{

    int		dir, xdir, ydir, zdir, vu_face;
    float	xwin, ywin;
    GLdouble	x0, y0, z0, x1, y1, z1, dd;
    GLdouble	model_mat[16], proj_mat[16];
    GLint	vu_port[4];



    set_current_window (vol->data_ctx->dpy_ctx);

    glGetDoublev (GL_PROJECTION_MATRIX, (GLdouble *) proj_mat);
    glGetDoublev (GL_MODELVIEW_MATRIX, (GLdouble *) model_mat);
    glGetIntegerv (GL_VIEWPORT, (GLint *) vu_port);

    xwin = ((float) (vu_port[0] + vu_port[2])) / 2.0;
    ywin = ((float) (vu_port[1] + vu_port[3])) / 2.0;

    gluUnProject ((GLdouble) xwin, (GLdouble) ywin, (GLdouble) 0.0,
                  (GLdouble *) model_mat, (GLdouble *) proj_mat,
                  (GLint *) vu_port, &x0, &y0, &z0);

    gluUnProject ((GLdouble) xwin, (GLdouble) ywin, (GLdouble) 1.0,
                  (GLdouble *) model_mat, (GLdouble *) proj_mat,
                  (GLint *) vu_port, &x1, &y1, &z1);

    x1 -= x0;
    y1 -= y0;
    z1 -= z0;
    dd  = sqrt ((x1 * x1) + (y1 * y1) + (z1 * z1));

    vol->vu_vec[0] = x1 / dd;
    vol->vu_vec[1] = y1 / dd;
    vol->vu_vec[2] = z1 / dd;

    /* Volume quad direction is opposite the viewing direction */
    xdir     = (vol->vu_vec[0] <= 0.0) ? VOLUME_X_DIR_FLAG : 0;
    ydir     = (vol->vu_vec[1] <= 0.0) ? VOLUME_Y_DIR_FLAG : 0;
    zdir     = (vol->vu_vec[2] <= 0.0) ? VOLUME_Z_DIR_FLAG : 0;
    vol->dir = xdir | ydir | zdir;

    vu_face = 0;
    if (ABS (vol->vu_vec[1]) > ABS (vol->vu_vec[vu_face])) vu_face = 1;
    if (ABS (vol->vu_vec[2]) > ABS (vol->vu_vec[vu_face])) vu_face = 2;
    vol->vu_face = vu_face;


    return vol->dir;
}
/* end PM */



/* 13Feb07  Phil McDonald */
float	volume_interp_val (Context data_ctx, int itime, int ivar,
	                   float *data, float xyz[3])
{
    int		ir0, ir1, ic0, ic1, il0, il1, nrow, ncol, nlev;
    float	row, col, lev;
    float	dr0, dr1, dc0, dc1, dl0, dl1;
    float	d0, d1, d2, d3, d4, d5, d6, d7;
    float	val;



    nrow = data_ctx->Nr;
    ncol = data_ctx->Nc;
    nlev = data_ctx->Nl[ivar];

    xyzPRIME_to_grid (data_ctx, itime, ivar, xyz[0], xyz[1], xyz[2],
                      &row, &col, &lev);

    ir0 = ir1 = row;
    dr0 = row - ((float) ir0);
    if (dr0 > VOLUME_FUZZ) ir1++;
    if ((ir0 < 0) || (ir1 >= nrow)) return MISSING;
    dr1 = 1.0 - dr0;

    ic0 = ic1 = col;
    dc0 = col - ((float) ic0);
    if (dc0 > VOLUME_FUZZ) ic1++;
    if ((ic0 < 0) || (ic1 >= ncol)) return MISSING;
    dc1 = 1.0 - dc0;

    il0 = il1 = lev;
    dl0 = lev - ((float) il0);
    if (dl0 > VOLUME_FUZZ) il1++;
    if ((il0 < 0) || (il1 >= nlev)) return MISSING;
    dl1 = 1.0 - dl0;

    d0  = data[(((il0 * ncol) + ic0) * nrow) + ir0];
    d1  = data[(((il0 * ncol) + ic0) * nrow) + ir1];
    d2  = data[(((il0 * ncol) + ic1) * nrow) + ir0];
    d3  = data[(((il0 * ncol) + ic1) * nrow) + ir1];
    d4  = data[(((il1 * ncol) + ic0) * nrow) + ir0];
    d5  = data[(((il1 * ncol) + ic0) * nrow) + ir1];
    d6  = data[(((il1 * ncol) + ic1) * nrow) + ir0];
    d7  = data[(((il1 * ncol) + ic1) * nrow) + ir1];

    if (IS_MISSING (d0)) return MISSING;
    if (IS_MISSING (d1)) return MISSING;
    if (IS_MISSING (d2)) return MISSING;
    if (IS_MISSING (d3)) return MISSING;
    if (IS_MISSING (d4)) return MISSING;
    if (IS_MISSING (d5)) return MISSING;
    if (IS_MISSING (d6)) return MISSING;
    if (IS_MISSING (d7)) return MISSING;

    val = (dl1 * ((d0 * dr1 * dc1) +
                  (d1 * dr0 * dc1) +
                  (d2 * dr1 * dc0) +
                  (d3 * dr0 * dc0))) +
          (dl0 * ((d4 * dr1 * dc1) +
                  (d5 * dr0 * dc1) +
                  (d6 * dr1 * dc0) +
                  (d7 * dr0 * dc0)));


    return val;
}
/* end PM */



/* 16Feb07  Phil McDonald */
int	volume_colors (Context data_ctx, int itime, int ivar)
{
    int			ixyz, nxyz;
    float		*data, datmin, datmax, datfac, val, *p_xyz;
    uint_1		*p_ind;
    struct volume	*v;



    if ((data = get_grid (data_ctx, itime, ivar)) == NULL) return 0;

    datmin = data_ctx->Variable[ivar]->MinVal;
    datmax = data_ctx->Variable[ivar]->MaxVal;
    datfac = (datmin == datmax) ? 0.0 : 254.0 / (datmax - datmin);

    v     = data_ctx->Volume;
    nxyz  = v->nx * v->ny * v->nz;
    p_xyz = v->vertex;
    p_ind = v->index;
    for (ixyz = 0; ixyz < nxyz; ixyz++)
    {
        val = volume_interp_val (data_ctx, itime, ivar, data, p_xyz);
        if (IS_MISSING (val))
        {
            *p_ind = 255;
        }
        else if (val <= datmin)
        {
            *p_ind = 0;
        }
        else if (val >= datmax)
        {
            *p_ind = 254;
        }
        else
        {
            *p_ind = (val - datmin) * datfac;
        }
        p_ind++;
        p_xyz += 3;
    }
    release_grid (data_ctx, itime, ivar, data);


    return 1;
}
/* end PM */



/* 16Feb07  Phil McDonald */
int	volume_verts (Context data_ctx, int use_grid_z)
{
    int			i, n, ix, iy, iz, nx, ny, nz;
    float		x, xmin, xmax, dx;
    float		y, ymin, ymax, dy;
    float		z, zmin, zmax, dz;
    float		*p_xyz, r, c, l;
    struct volume	*v;



    v = data_ctx->Volume;


    r = data_ctx->Nr - 1;
    c = 0.0;
    l = 0.0;
    grid_to_xyzPRIME (data_ctx, -1, -1, 1, &r, &c, &l, &xmin, &ymin, &zmin);

    r = 0.0;
    c = data_ctx->Nc - 1;
    l = data_ctx->MaxNl - 1;
    grid_to_xyzPRIME (data_ctx, -1, -1, 1, &r, &c, &l, &xmax, &ymax, &zmax);

    nx = v->nx;
    ny = v->ny;
    nz = v->nz;
    dx = (xmax - xmin) / (float) (nx - 1);
    dy = (ymax - ymin) / (float) (ny - 1);
    dz = (zmax - zmin) / (float) (nz - 1);

    p_xyz = v->vertex;
    z     = zmin;
    for (iz = 0; iz < nz; iz++)
    {
        if (use_grid_z) z = gridlevel_to_z (data_ctx, -1, -1, iz);

        y = ymin;
        for (iy = 0; iy < ny; iy++)
        {
            x = xmin;
            for (ix = 0; ix < nx; ix++)
            {
                *(p_xyz++) = x;
                *(p_xyz++) = y;
                *(p_xyz++) = z;

                x += dx;
                if (x > xmax) x = xmax;
            }
            y += dy;
            if (y > ymax) y = ymax;
        }
        z += dz;
        if (z > zmax) z = zmax;
    }


    return 1;
}
/* end PM */



/* 25Oct07  Phil McDonald */
int	volume_vert_quads (int ix, int iy, int iz, int nx, int ny, int nz,
	                   int quad_off[], unsigned int **p_p_qiv)
{
    int			i, jx, jy, jz, iaxis, ivert;
    unsigned int	*p_qiv;



    p_qiv = *p_p_qiv;

    i = 0;
    for (iaxis = 0; iaxis < 3; iaxis++)
    {
        for (ivert = 0; ivert < 4; ivert++)
        {
            jx = ix - quad_off[i++];
            jy = iy - quad_off[i++];
            jz = iz - quad_off[i++];
            if ((jx < 0) || (jx >= nx) ||
                (jy < 0) || (jy >= ny) ||
                (jz < 0) || (jz >= nz))
            {
                /* The current quad extends beyond the limits of the     */
                /* volume.  Replace it with a null quad whose 4 vertices */
                /* are vertex[0].                                        */
                i          += (3 - ivert) * 3;
                p_qiv      -= ivert;
                *(p_qiv++)  = 0;
                *(p_qiv++)  = 0;
                *(p_qiv++)  = 0;
                *(p_qiv++)  = 0;
                break;
            }
            *(p_qiv++) = (((jz * ny) + jy) * nx) + jx;
        }
    }

    *p_p_qiv = p_qiv;


    return 1;
}
/* end PM */



/* 16Feb07  Phil McDonald */
int	volume_quads (struct VolData *vol)
{
    int		i, n, nx, ny, nz, nxy, quad_dir;

    int		quad_fwd[12][3] = { 0, 0, 0,
				    0, 1, 0,
				    0, 1, 1,
				    0, 0, 1,
		                    0, 0, 0, 
				    0, 0, 1,
				    1, 0, 1,
				    1, 0, 0,
				    0, 0, 0,
				    0, 1, 0,
				    1, 1, 0,
				    1, 0, 0 };
    int		quad_rev[12][3] = { 0, 0, 0,
				    0, 0, 1,
				    0, 1, 1,
				    0, 1, 0,
		                    0, 0, 0,
				    1, 0, 0,
				    1, 0, 1,
				    0, 0, 1,
				    0, 0, 0,
				    1, 0, 0,
				    1, 1, 0,
				    0, 1, 0 };
    int		*quad_off;

    int			ixstp, ixbeg, ixend, ix0, ix1;
    int			iystp, iybeg, iyend, iy0, iy1;
    int			izstp, izbeg, izend, iz0, iz1;
    unsigned int	*p_qiv;
    struct volume	*v;
#ifdef TIME_TEST
    double		usec0 = vis5d_time_usec (0.0);
#endif



    v   = vol->data_ctx->Volume;
    nx  = v->nx;
    ny  = v->ny;
    nz  = v->nz;
    nxy = nx * ny;


    if (vol->quad_ivert == NULL)
    {
        n               = nx * ny * nz * 3 * 4;
        vol->quad_ivert = allocate (vol->data_ctx, (n * sizeof (unsigned int)));
    }


    quad_dir = 0;

    ix0 = 0, ix1 = nx - 1;
    if (vol->dir & VOLUME_X_DIR_FLAG)
    {
        ixstp = 1;
        ixbeg = ix0;
        ixend = ix1;
        quad_dir++;
    }
    else
    {
        ixstp = -1;
        ixbeg = ix1;
        ixend = ix0;
    }

    iy0 = 0, iy1 = ny - 1;
    if (vol->dir & VOLUME_Y_DIR_FLAG)
    {
        iystp = 1;
        iybeg = iy0;
        iyend = iy1;
        quad_dir++;
    }
    else
    {
        iystp = -1;
        iybeg = iy1;
        iyend = iy0;
    }

    iz0 = 0, iz1 = nz - 1;
    if (vol->dir & VOLUME_Z_DIR_FLAG)
    {
        izstp = 1;
        izbeg = iz0;
        izend = iz1;
        quad_dir++;
    }
    else
    {
        izstp = -1;
        izbeg = iz1;
        izend = iz0;
    }
/* 30Jun06  Phil McDonald */
    /*  Backup a step to get the quads at the back faces */
    ixbeg -= ixstp;
    iybeg -= iystp;
    izbeg -= izstp;
/* end PM */

    quad_off = (quad_dir % 2) ? (int *) quad_rev : (int *) quad_fwd;
    for (i = 0; i < 36; i += 3)
    {
        if (quad_off[i+0] != 0) quad_off[i+0] = ixstp;
        if (quad_off[i+1] != 0) quad_off[i+1] = iystp;
        if (quad_off[i+2] != 0) quad_off[i+2] = izstp;
    }


    p_qiv = vol->quad_ivert;

    if (vol->vu_face == 0)
    {
        for (ix0 = ixbeg; ix0 != ixend; ix0 = ix1)
        {
            ix1 = ix0 + ixstp;
            for (iy0 = iybeg; iy0 != iyend; iy0 = iy1)
            {
                iy1 = iy0 + iystp;
                for (iz0 = izbeg; iz0 != izend; iz0 = iz1)
                {
                    iz1 = iz0 + izstp;

                    volume_vert_quads (ix1, iy1, iz1, nx, ny, nz,
                                       quad_off, &p_qiv);
                }
            }
        }
    }
    else if (vol->vu_face == 1)
    {
        for (iy0 = iybeg; iy0 != iyend; iy0 = iy1)
        {
            iy1 = iy0 + iystp;
            for (ix0 = ixbeg; ix0 != ixend; ix0 = ix1)
            {
                ix1 = ix0 + ixstp;
                for (iz0 = izbeg; iz0 != izend; iz0 = iz1)
                {
                    iz1 = iz0 + izstp;

                    volume_vert_quads (ix1, iy1, iz1, nx, ny, nz,
                                       quad_off, &p_qiv);
                }
            }
        }
    }
    else
    {
        for (iz0 = izbeg; iz0 != izend; iz0 = iz1)
        {
            iz1 = iz0 + izstp;
            for (iy0 = iybeg; iy0 != iyend; iy0 = iy1)
            {
                iy1 = iy0 + iystp;
                for (ix0 = ixbeg; ix0 != ixend; ix0 = ix1)
                {
                    ix1 = ix0 + ixstp;

                    volume_vert_quads (ix1, iy1, iz1, nx, ny, nz,
                                       quad_off, &p_qiv);
                }
            }
        }
    }
    vol->nquad = nx * ny * nz * 3;
#ifdef TIME_TEST
    vis5d_time_msg (usec0, "quads", "total");
#endif


    return 1;
}
/* end PM */



/* 25Oct07  Phil McDonald */
int	volume_new_dims (Context data_ctx, int new_dims[], int *p_use_grid_z)
{
    int		nx, ny, nz;



    nx = data_ctx->dpy_ctx->vol_nx;
    ny = data_ctx->dpy_ctx->vol_ny;
    nz = data_ctx->dpy_ctx->vol_nz;

    if ((nx < 0) || (ny < 0) || (nz < 0))
    {
        int     mx, my, mz, maxpts;
        float   x0, y0, z0, x1, y1, z1, xdim, ydim, zdim, res, xtmp;

        vis5d_grid_to_xyz (data_ctx->context_index, -1, -1,
                           0.0, 0.0, 0.0, &x0, &y0, &z0);
        vis5d_grid_to_xyz (data_ctx->context_index, -1, -1,
                           (float) (data_ctx->Nr - 1),
                           (float) (data_ctx->Nc - 1),
                           (float) (data_ctx->MaxNl - 1),
                           &x1, &y1, &z1);
        xdim = x1 - x0;
        ydim = y0 - y1;
        zdim = z1 - z0;

        mx = (nx != 0) ? nx : data_ctx->Nc;
        my = (ny != 0) ? ny : data_ctx->Nr;
        mz = (nz != 0) ? nz : data_ctx->MaxNl;

        maxpts = 0;
        if (nx < -1) maxpts = -nx;
        if (ny < -1) maxpts = -ny;
        if (nz < -1) maxpts = -nz;

        res = 0.0;
        if (maxpts > 0)
        {
            float       en_pts, en_dims;

            en_pts  = maxpts;
            en_dims = 3.0;
            if (mx > 0) en_pts /= (float) mx, en_dims -= 1.0, xdim = 1.0;
            if (my > 0) en_pts /= (float) my, en_dims -= 1.0, ydim = 1.0;
            if (mz > 0) en_pts /= (float) mz, en_dims -= 1.0, zdim = 1.0;
            xtmp = en_pts / (xdim * ydim * zdim);
            res  = pow ((double) xtmp, (double) (1.0 / en_dims));
        }
        else
        {
            if (mx > 0)
            {
                xtmp = (float) mx / xdim;
                if (xtmp > res) res = xtmp;
            }
            if (my > 0)
            {
                xtmp = (float) my / ydim;
                if (xtmp > res) res = xtmp;
            }
            if (nz > 0)
            {
                xtmp = (float) mz / zdim;
                if (xtmp > res) res = xtmp;
            }
        }

        if (res == 0.0)
        {
            fprintf (stderr, "%s: %s\n",
                     "vis5d_set_volume_dims",
                     "insufficient info to recompute dims");
            return 1;
        }

        if (nx < 0) nx = (xdim * res) + 0.5;
        if (ny < 0) ny = (ydim * res) + 0.5;
        if (nz < 0) nz = (zdim * res) + 0.5;
    }


    *p_use_grid_z = (nz == 0);

    if (nx == 0) nx = data_ctx->Nc;
    if (ny == 0) ny = data_ctx->Nr;
    if (nz == 0) nz = data_ctx->MaxNl;

    new_dims[0] = nx;
    new_dims[1] = ny;
    new_dims[2] = nz;

    if ((nx == data_ctx->Volume->nx) &&
        (ny == data_ctx->Volume->ny) &&
        (nz == data_ctx->Volume->nz)) return 0;


    return 1;
}
/* end PM */



/* 16Feb07  Phil McDonald */
int	make_volume (Context data_ctx, int itime, int ivar)
{
    int			new_dims[3], use_grid_z;
    struct volume	*v;



    if ((data_ctx == NULL) || (data_ctx->Volume == NULL)) return 0;

    v = data_ctx->Volume;


#ifdef TIME_TEST
    double	usec;
    usec = vis5d_time_usec (0.0);
#endif
    wait_write_lock (&(v->lock));
#ifdef TIME_TEST
    vis5d_time_msg (usec, "make_vol", "wait");
#endif


    /*
     *  Check to see if new volume grid dimensions have been specified.
     *  If so, get rid of the old vertex and quad arrays, allocate
     *  new ones, and find the verts for the new grid.
     */

    if (volume_new_dims (data_ctx, new_dims, &use_grid_z))
    {
        free_volume_data (data_ctx);

        v->nx = new_dims[0];
        v->ny = new_dims[1];
        v->nz = new_dims[2];
        alloc_volume_data (data_ctx);
    }


    if ((itime != v->itime) || (ivar != v->ivar) || (v->valid == 0))
    {
        if (ivar >= 0)
        {
            if (v->ivar < 0) volume_verts (data_ctx, use_grid_z);

            volume_colors (data_ctx, itime, ivar);
            v->valid = 1;
        }

        v->itime = itime;
        v->ivar  = ivar;
    }

    done_write_lock (&(v->lock));


#ifdef TIME_TEST
    vis5d_time_msg (usec, "make_vol", "total");
#endif
    return 1;
}
/* end PM */



/* 26Feb07  Phil McDonald	100 seems like an optimal number */
#define VOLUME_RENDER_MAX_QUADS         100
#define VOLUME_RENDER_MAX_VERTS         VOLUME_RENDER_MAX_QUADS*4

/* 11May06  Phil McDonald */
int	volume_render (struct VolData *vol, unsigned int color_tab[])
{
    int			i, iquad, nquad, iax, iqbeg, iqend, iv;
    int			ir, ig, ib, ia;
    unsigned int	color, *p_qiv;
    float		*p_xyz, wt[3];
    uint_1		*p_ind;
    struct volume	*v;
    float		qvert[VOLUME_RENDER_MAX_VERTS*3], *p_qv;
    unsigned int	qrgba[VOLUME_RENDER_MAX_VERTS], *p_qc;
#ifdef TIME_TEST
    double	usec0 = vis5d_time_usec (0.0);
#endif



    v = vol->data_ctx->Volume;

    if ((vol->nquad == 0) || (vol->quad_ivert == NULL)) return 0;


    wt[0] = ABS (vol->vu_vec[0]);
    wt[1] = ABS (vol->vu_vec[1]);
    wt[2] = ABS (vol->vu_vec[2]);


    glPushMatrix();
    glEnable (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDepthMask (GL_FALSE);


    nquad = vol->nquad;
    p_qiv = vol->quad_ivert;
    p_ind = v->index;
    p_xyz = v->vertex;


#define USE_GL_VERTEX_ARRAYS
#ifdef USE_GL_VERTEX_ARRAYS
    glEnableClientState (GL_COLOR_ARRAY);
    glColorPointer (4, GL_UNSIGNED_BYTE, 0, qrgba);
    glEnableClientState (GL_VERTEX_ARRAY);
    glVertexPointer (3, GL_FLOAT, 0, qvert);

    iqend = 0;
    while (iqend < nquad)
    {
        p_qv = qvert;
        p_qc = qrgba;

        iqbeg  = iqend;
        iqend += VOLUME_RENDER_MAX_QUADS;
        if (iqend > nquad) iqend = nquad;
        for (iquad = iqbeg; iquad < iqend; iquad++)
        {
            iax = iquad % 3;
            for (i = 0; i < 4; i++)
            {
                iv    = *(p_qiv++);

                color = color_tab[p_ind[iv]];
                ia    = UNPACK_ALPHA (color);
                if (ia != 255)
                {
                    ia   *= wt[iax];
#ifdef WORDS_BIGENDIAN
                    color = (color & 0xffffff00) | ia;
#else
                    color = (color & 0x00ffffff) | (ia << 24);
#endif
                }
                *(p_qc++) = color;

                *(p_qv++) = p_xyz[(iv*3)+0];
                *(p_qv++) = p_xyz[(iv*3)+1];
                *(p_qv++) = p_xyz[(iv*3)+2];
            }
        }
        glDrawArrays (GL_QUADS, 0, ((iqend-iqbeg)*4));
    }

    glDisableClientState (GL_VERTEX_ARRAY);
    glDisableClientState (GL_COLOR_ARRAY);
#else
    glBegin (GL_QUADS);

    for (iquad = 0; iquad < nquad; iquad++)
    {
        iax = iquad % 3;
        for (i = 0; i < 4; i++)
        {
            color  = color_tab[p_ind[*p_qiv]];
            ia     = UNPACK_ALPHA (color);
            ia    *= wt[iax];
            if (ia != 255)
            {
                ia   *= wt[iax];
#ifdef WORDS_BIGENDIAN
                color = (color & 0xffffff00) | ia;
#else
                color = (color & 0x00ffffff) | (ia << 24);
#endif
            }

            glColor4ubv ((GLubyte *) &color);
            glVertex3fv (&(p_xyz[(*p_qiv)*3]));
            p_qiv++;
        }
    }

    glEnd ();
#endif

    glDepthMask (GL_TRUE);
    glPopMatrix ();

    glDisable (GL_BLEND);
    glDisable (GL_ALPHA_TEST);
#ifdef TIME_TEST
    vis5d_time_msg (usec0, "render", "total");
#endif


    return 1;
}
/* end PM */



/* 18May06  Phil McDonald */
void	draw_volume (Context data_ctx, int itime, int ivar,
	             unsigned int *color_tab)
{

    int			cur_dir, cur_face;
    struct volume	*v;
    struct VolData	*vol;
#ifdef TIME_TEST
    static double	ntime			= 0.0;
    static double	time_tot		= 0.0;
    double		usec0 = vis5d_time_usec (0.0);
    char		time_str[32];
    vis5d_debug_msg2 ("draw_vol", "begin");
#endif



    if ((data_ctx == NULL) || (data_ctx->Volume == NULL)) return;

    v = data_ctx->Volume;

    if (v->valid == 0) return;


    /* Get the volume data for the current process.  If this is the */
    /* first time, the volume data will be created.  Since this     */
    /* data will not be shared with other processes (even though it */
    /* may be in shared memory) there is no need to lock it.        */

    vol           = volume_data_for_current_process (data_ctx);
    vol->data_ctx = data_ctx;

    cur_dir  = vol->dir;
    cur_face = vol->vu_face;
    volume_dir (vol);

    if ((vol->dir != cur_dir) || (vol->vu_face != cur_face))
    {
        volume_quads (vol);
    }
    volume_render (vol, color_tab);


#ifdef TIME_TEST
    ntime++;
    time_tot += vis5d_time_usec (usec0);
    vis5d_time_msg (usec0, "draw_vol", "total");
    sprintf (time_str, "average time = %.0lf usec", (time_tot / ntime));
    vis5d_debug_msg2 ("draw_vol", time_str);
#endif
    return;
}
