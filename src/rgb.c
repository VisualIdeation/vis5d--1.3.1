/* rgb.c */


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
 * Code to read SGI .rgb image files.  Borrowed from the SGI OpenGL TK
 * source code.
 */



#include <stdio.h>
#include <stdlib.h> 
#include <string.h>
#include "rgb.h"
/* 31Jan06  Phil McDonald */
#include "memory.h"
/* end PM */



#ifndef SEEK_SET
#  define SEEK_SET 0
#endif



typedef struct _rawImageRec {
    unsigned short imagic;
    unsigned short type;
    unsigned short dim;
    unsigned short sizeX, sizeY, sizeZ;
    unsigned long min, max;
    unsigned long wasteBytes;
    char name[80];
    unsigned long colorMap;
    FILE *file;
    unsigned char *tmp, *tmpR, *tmpG, *tmpB;
    unsigned long rleEnd;
    unsigned int *rowStart;
    int *rowSize;
} rawImageRec;



static void ConvertShort(unsigned short *array, long length)
{
    unsigned long b1, b2;
    unsigned char *ptr;

    ptr = (unsigned char *)array;
    while (length--) {
        b1 = *ptr++;
        b2 = *ptr++;
        *array++ = (b1 << 8) | (b2);
    }
}

static void ConvertLong( unsigned int *array, long length)
{
    unsigned long b1, b2, b3, b4;
    unsigned char *ptr;

    ptr = (unsigned char *)array;
    while (length--) {
        b1 = *ptr++;
        b2 = *ptr++;
        b3 = *ptr++;
        b4 = *ptr++;
        *array++ = (b1 << 24) | (b2 << 16) | (b3 << 8) | (b4);
    }
}

static void RawImageClose(rawImageRec *raw)
{

    fclose(raw->file);
/* 23Oct06  Phil McDonald */
    if (raw->tmp      != NULL) TMP_FREE (raw->tmp);
    if (raw->tmpR     != NULL) TMP_FREE (raw->tmpR);
    if (raw->tmpG     != NULL) TMP_FREE (raw->tmpG);
    if (raw->tmpB     != NULL) TMP_FREE (raw->tmpB);
    if (raw->rowStart != NULL) TMP_FREE (raw->rowStart);
    if (raw->rowSize  != NULL) TMP_FREE (raw->rowSize);
    TMP_FREE (raw);
/* end PM */
}

static rawImageRec *RawImageOpen(char *fileName)
{
    union {
        int testWord;
        char testByte[4];
    } endianTest;
    rawImageRec *raw;
    int swapFlag;
    int x;

    endianTest.testWord = 1;
    if (endianTest.testByte[0] == 1) {
        swapFlag = 1;
    } else {
        swapFlag = 0;
    }

/* 23Oct06  Phil McDonald */
    raw = (rawImageRec *) TMP_CALLOC (1, sizeof (rawImageRec));
/* end PM */
    if (raw == NULL) {
        fprintf(stderr, "Out of memory!\n");
        return NULL;
    }
    if ((raw->file = fopen(fileName, "rb")) == NULL) {
        perror(fileName);
        return NULL;
    }

    fread(raw, 1, 12, raw->file);

    if (swapFlag) {
        ConvertShort(&raw->imagic, 6);
    }

/* 23Oct06  Phil McDonald */
    raw->tmp  = (unsigned char *) TMP_MALLOC (raw->sizeX*256);
    raw->tmpR = (unsigned char *) TMP_MALLOC (raw->sizeX*256);
    raw->tmpG = (unsigned char *) TMP_MALLOC (raw->sizeX*256);
    raw->tmpB = (unsigned char *) TMP_MALLOC (raw->sizeX*256);
    if (raw->tmp == NULL || raw->tmpR == NULL || raw->tmpG == NULL ||
        raw->tmpB == NULL) {
        fprintf(stderr, "Out of memory!\n");
        RawImageClose (raw);
/* end PM */
        return NULL;
    }

    if ((raw->type & 0xFF00) == 0x0100) {
        x = raw->sizeY * raw->sizeZ * sizeof(unsigned int);
/* 23Oct06  Phil McDonald */
        raw->rowStart = (unsigned int *) TMP_MALLOC (x);
        raw->rowSize = (int *) TMP_MALLOC (x);
        if (raw->rowStart == NULL || raw->rowSize == NULL) {
            fprintf(stderr, "Out of memory!\n");
            RawImageClose (raw);
/* end PM */
            return NULL;
        }
        raw->rleEnd = 512 + (2 * x);
        fseek(raw->file, 512, SEEK_SET);
        fread(raw->rowStart, 1, x, raw->file);
        fread(raw->rowSize, 1, x, raw->file);
        if (swapFlag) {
            ConvertLong(raw->rowStart, x/sizeof(unsigned int));
            ConvertLong((unsigned int *)raw->rowSize, x/sizeof(int));
        }
    }
    return raw;
}

static void RawImageGetRow(rawImageRec *raw, unsigned char *buf, int y, int z)
{
    unsigned char *iPtr, *oPtr, pixel;
    int count;

    if ((raw->type & 0xFF00) == 0x0100) {
        fseek(raw->file, raw->rowStart[y+z*raw->sizeY], SEEK_SET);
        fread(raw->tmp, 1, (unsigned int)raw->rowSize[y+z*raw->sizeY],
              raw->file);

        iPtr = raw->tmp;
        oPtr = buf;
        while (1) {
            pixel = *iPtr++;
            count = (int)(pixel & 0x7F);
            if (!count) {
                return;
            }
            if (pixel & 0x80) {
                while (count--) {
                    *oPtr++ = *iPtr++;
                }
            } else {
                pixel = *iPtr++;
                while (count--) {
                    *oPtr++ = pixel;
                }
            }
        }
    } else {
        fseek(raw->file, 512+(y*raw->sizeX)+(z*raw->sizeX*raw->sizeY),
              SEEK_SET);
        fread(buf, 1, raw->sizeX, raw->file);
    }
}

static void RawImageGetData(rawImageRec *raw, RGB_IMAGE *final)
{
    unsigned char *ptr;
    int i, j;

/* 31Jan06  Phil McDonald */
    final->data = (unsigned char *) mem_malloc
                  ((raw->sizeX+1) * (raw->sizeY+1) * 4);
/* end PM */
    if (final->data == NULL) {
        fprintf(stderr, "Out of memory!\n");
        return;
    }

    ptr = final->data;
    for (i = 0; i < (int)(raw->sizeY); i++) {
       int ii = raw->sizeY - i - 1;    /* flip vertically */
        RawImageGetRow(raw, raw->tmpR, ii, 0);
        RawImageGetRow(raw, raw->tmpG, ii, 1);
        RawImageGetRow(raw, raw->tmpB, ii, 2);
        for (j = 0; j < (int)(raw->sizeX); j++) {
/* 16Jan07  Phil McDonald */
#ifdef LEAVEOUT
#ifdef HAVE_SGI_GL
            /* IRIS GL uses BGR order for textures */
            *ptr++ = *(raw->tmpR + j);
            *ptr++ = *(raw->tmpG + j);
            *ptr++ = *(raw->tmpB + j);
#else
            /* OpenGL, Mesa, etc use RGB order for textures */
            *ptr++ = *(raw->tmpB + j);
            *ptr++ = *(raw->tmpG + j);
            *ptr++ = *(raw->tmpR + j);
#endif
#else
#  ifdef WORDS_BIGENDIAN
            *ptr++ = *(raw->tmpB + j);
            *ptr++ = *(raw->tmpG + j);
            *ptr++ = *(raw->tmpR + j);
#  else
            *ptr++ = *(raw->tmpR + j);
            *ptr++ = *(raw->tmpG + j);
            *ptr++ = *(raw->tmpB + j);
#  endif
#endif
/* 08Nov06  Phil McDonald	set alpha to 0xff, not 0 */
            *ptr++ = 0xff;
/* end PM */
        }
    }
}

RGB_IMAGE *ReadRGB(char *fileName)
{
    rawImageRec *raw;
    RGB_IMAGE *final;

    raw = RawImageOpen(fileName);
    if (!raw) {
       return NULL;
    }
/* 31Jan06  Phil McDonald */
    final = (RGB_IMAGE *) mem_malloc (sizeof (RGB_IMAGE));
/* end PM */
    if (final == NULL) {
        fprintf(stderr, "Out of memory!\n");
        return NULL;
    }
    final->sizeX = raw->sizeX;
    final->sizeY = raw->sizeY;
    RawImageGetData(raw, final);
    RawImageClose(raw);
    return final;
}




void FreeRGB( RGB_IMAGE *image )
{
   if (image) {
/* 31Jan06  Phil McDonald */
      mem_free (image->data);
      mem_free (image);
/* end PM */
   }
}


