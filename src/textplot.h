/* textplot.h */

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

#ifndef TEXTPLOT_H
#define TEXTPLOT_H

#define MAX_TEXT_PLOT_VERTS 100000

extern int create_num_textplot( Irregular_Context itx, int time, float xs[],
                      float ys[], float zs[], double *numdata, int ploton[],
                      float vx[], 
                      float vy[], float vz[], int *numv);

extern int create_letter_textplot( Irregular_Context itx, int time, float xs[],
                      float ys[], float zs[], char *chardata, int ploton[], int var,
                      float vx[],
                      float vy[], float vz[], int *numv);

extern int create_color_num_textplot( Irregular_Context itx, int time, float xs[],
                     float ys[], float zs[], double *numdata,
                     int ploton[], float vx[],
                     float vy[], float vz[], int *numv,
                     uint_1 *color_indexes);

extern void space_plots( Irregular_Context itx, int time,       
                     int ploton[], float xs[],
                     float ys[], float zs[], int *numtouse);

#endif

