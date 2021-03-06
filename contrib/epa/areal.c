
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <fortran.h>

/* ... definitions ... */

#define UNKNOWN 0			/* UNKNOWN file type */
#define UAM_CONC 1			/* UAM concentration file type */
#define UAM_EMIS 2			/* UAM emissions file type */
#define UAM_BOUND 3			/* UAM boundary file type */

#define cw_EOD 017                      /* control word END-OF-DATA */
#define cw_EOF 016                      /* control word END-OF-FILE */
#define cw_EOR 010                      /* control word END-OF-RECORD */

#define SUCCESS 1                       /* successful function return */
#define FAILURE 0                       /* unsuccessful function return */
#define IOERROR -1                      /* unsuccessful i/o function return */

#define MAXBUF 4096                     /* max # bytes for UAM i/o buffer */
#define MAXSPEC 40			/* max # of species */

/* variables global to the UAM routines only */

static int arealncol;			/* # columns in data grid */
static int arealnrow;			/* # rows in data grid */
static int arealnlayer;			/* # levels (layers) in grid */
static int arealnrecord;			/* # records in file */
static int arealnspecies;			/* # species in the input file */
static int arealfilesize;			/* file size in bytes */
static int arealfiletype;			/* file type: CONC or other */
static int areallevellen;			/* # floats in a layer of data */
static int arealspeclen;			/* # floats in a species of data */
static int arealsteplen;			/* # bytes per time step */
static int arealntext;			/* # text header records */
static int arealreclen;			/* # bytes per UAM record */
static int arealdatablock1;		/* UAM starting data block # */
static int arealdatapos1;			/* UAM starting data word position */
static char arealspecname[MAXSPEC*5+1];	/* short name of all species */ 
static float sw_utmx;			/* lower left map x-coordinate */
static float sw_utmy;			/* lower left map y-coordinate */
static float ne_utmx;			/* upper right map x-coordinate */
static float ne_utmy;			/* upper right map y-coordinate */

/*******************************************************************/
/* areal_open                                                        */
/* Function: open input file and determine filesize and filetype   */
/* On error: return FAILURE and write error string into message    */
/* If no error: return SUCCESS                                     */
/*******************************************************************/
int areal_open(fd, filename, message)
int *fd;				/* INPUT:  file descriptor */
char *filename;				/* INPUT:  file name */
char *message;				/* OUTPUT: returned error message */
{
struct stat statbuf;			/* UNIX file status buffer */
register int istat;			/* I/O status from UNIX stat call */
register int status;			/* areal_open function return status */

status = FAILURE;			/* assume FAILURE */

if ((*fd = open(filename, O_RDONLY)) != IOERROR)
	{
	if ((istat = stat(filename, &statbuf)) != IOERROR)
		{
		arealfilesize = statbuf.st_size;
		if (arealfilesize > 0)
			{
			arealfiletype = get_areal_filetype(*fd, message);
			if ((arealfiletype == UAM_CONC) 
				|| (arealfiletype == UAM_EMIS))
				status = SUCCESS;
			else
				{
				sprintf(message,
					 "File is not recognized as UAM"); 
				close(*fd);
				}
			}
		else
			{
			sprintf(message, "File size of zero");
			close(*fd);
			}
		}
	else
		{
		sprintf(message, "Cannot get file status");
		close(*fd);
		}
	}
else
	sprintf(message, "Cannot open file %s", filename);
return(status);
}

/**********************************************************************/
/* get_areal_filetype                                                   */
/* Function: read 2 words from UAM file to determine filetype         */
/* On error: return UNKNOWN and write error string into message       */
/* If no error: return UAM_CONC or UAM_EMIS                           */
/**********************************************************************/
int get_areal_filetype(fd, message)
int fd;					/* INPUT:  file descriptor */
char *message;				/* OUTPUT: return error message */
{
static char buffer[32];
#define arealconc     "AVERAGE " 		/* UAM CRAY word 2 */
#define arealinst     "AIRQUALI" 		/* UAM CRAY word 2 */
#define arealemis     "EMISSION"		/* UAM CRAY word 2 */
#define arealbound    "B       O       "  /* UAM CRAY words 2 & 3 */

register int filetype;			/* return filetype */

filetype = UNKNOWN;		/* assume filetype cannot be determined */

if ((lseek(fd, (off_t) (1 * 8), SEEK_SET)) != IOERROR)
	{
	if ((read(fd, buffer, 8)) != IOERROR)
		{
		buffer[8] = '\0';
		if (!strcmp(buffer, arealconc))
			filetype = UAM_CONC;
		else if (!strcmp(buffer, arealinst))
			filetype = UAM_CONC;
		else if (!strcmp(buffer, arealemis))
			filetype = UAM_EMIS;
/*
		else if (!strcmp(buffer, arealbound))
			filetype = UAM_BOUND;
*/
		}
	else
		sprintf(message, "Cannot read filetype");
	}
else
	sprintf(message, "Cannot seek file position to determine filetype");
return(filetype);
}

int areal_inquire(fd, message, ncol, nrow, nlayer, nrecord, nspecies, specname)
int fd;
char *message;
int *ncol;
int *nrow;
int *nlayer;
int *nrecord;
int *nspecies;
char *specname;
{
if (areal_header(fd, message))
	{
	*ncol = arealncol;
	*nrow = arealnrow;
	*nlayer = arealnlayer;
	*nrecord = arealnrecord;
	*nspecies = arealnspecies;
	sprintf(specname, "%s", arealspecname);
	return(SUCCESS);
	}
return(FAILURE);
}

/*******************************************************************/
/* get_areal_offset                                                       */
/* Function: get the file offset for nwords input                  */
/*	     (nwords = CRAY whole words for low resolution files ) */
/* On error: return FAILURE and write error string into message    */
/* If no error: return SUCCESS                                     */
/*******************************************************************/
int get_areal_offset(nwords, offset, message)
int nwords;		/* INPUT:  # of words beyond header */ 
off_t *offset;		/* OUTPUT: file offset in bytes */
char *message;		/* OUTPUT: message buffer for error string */
{
register int status;

status = get_areal_unicos_offset(nwords, offset);
if (status == FAILURE)
	sprintf(message, "Error setting file offset");
return(status);
}

/*******************************************************************/
/* get_areal_buf                                                          */
/* Function: read nwords (CRAY words) into 'buffer',               */
/*           reading from file 'fd' of 'resolution' size at        */
/*           'offset' position.                                    */
/* On error: return FAILURE and write error string into message    */
/* If no error: return SUCCESS                                     */
/*******************************************************************/
int get_areal_buf(fd, buffer, nwords, offset, message)
int fd;					/* INPUT:  file descriptor */
char *buffer;				/* OUTPUT: character data buffer */
int nwords;				/* INPUT:  # words/halfwords */
off_t offset;				/* INPUT:  file offset in bytes */
char *message;				/* OUTPUT: buffer for error string */
{
register int status;

status = get_areal_unicos_buf(fd, buffer, nwords, offset);
if (status == FAILURE)
	sprintf(message, "Error reading data");
return(status);
}


/*******************************************************************/
/* get_areal_unicos_offset                                               */
/* Function: get the file offset for nwords input                  */
/*	     (nwords = CRAY whole words for low resolution files ) */
/* On error: return FAILURE                                        */
/* If no error: return SUCCESS                                     */
/*******************************************************************/
int get_areal_unicos_offset(nwords, offset)
int nwords;			/* INPUT:  # 8-byte CRAY words offset */
off_t *offset;			/* OUTPUT: file offset in bytes */
{
register int twords;		/* # words yet to count toward offset */
register int pos;		/* word offset position within offset block */
register int block;		/* 512-word block containing offset pos */
register int status;		/* return status */

status = SUCCESS;		/* assume SUCCESS */

/* begin offset search at EOR just beyond last header record */

block = arealdatablock1;
pos = arealdatapos1;

twords = nwords;		/* # words to count */
while (twords + pos > 511)	/* if # words to count span blocks ... */
	{
	++block;		/* . increment block # */
	twords -= (511 - pos);	/* . decrement # words by block remainder */
	pos = 0;		/* . reset position within new block */
	}
pos += twords;			/* set position to # words remaining */
*offset = (block * 512 + pos) * 8;	/* convert word offset to bytes */
if (*offset >= arealfilesize)	/* if offset is beyond EOF return with error */
	status = FAILURE;
return(status);
}

/*******************************************************************/
/* get_areal_unicos_buf                                                  */
/* Function: read nwords (CRAY words) into 'buffer',               */
/*           reading from file 'fd' of 'resolution' size at        */
/*           'offset' position.                                    */
/* On error: return FAILURE                                        */
/* If no error: return SUCCESS                                     */
/*******************************************************************/
int get_areal_unicos_buf(fd, buffer, nwords, offset)
int fd;				/* INPUT:  file descriptor */
char *buffer;			/* OUTPUT: character data buffer */
int nwords;			/* INPUT:  # words to store in buffer */
off_t offset;			/* INPUT:  beginning file position in bytes */
{
register int gotbytes;		/* # bytes read so far */
register int getbytes;		/* # bytes to read into buffer */
register int n;			/* buffer byte position to store next read */
register int status;		/* return status */

unsigned int cw;		/* control word */
int type;			/* control word type */

status = SUCCESS;		/* assume SUCCESS */

gotbytes = n = 0;		/* no bytes read, buffer position = 0 */

if ((lseek(fd, offset, SEEK_SET)) == IOERROR)
	status = FAILURE;

/* continue until # requested words read or until FATAL error encountered */

while ((gotbytes < nwords * 8)&& (status == SUCCESS))
	{
	if ((read(fd, &cw, 8)) == 8)	/* get control word (8 bytes) */
		{
		type = (cw) >> 60;	/* type of control word */	
		if ((type != cw_EOD) && (type != cw_EOF))
			{
			/* get # bytes to read in this record */
			getbytes = (cw & 0777) * 8;

			/* do not get more than requested */
			if (getbytes + n > nwords * 8)
				getbytes = nwords * 8 - n;

			/* read # bytes available/requested */
			if ((read(fd, buffer+n, getbytes)) == getbytes)
				{
				/* increment # bytes in buffer by # read */
				gotbytes += getbytes;

				/* set buffer position accordingly */
				n = gotbytes;
				}
			else
				status = FAILURE;
			}
		else
			status = FAILURE;
		}
	else
		status = FAILURE;
	}
return(status);
}

/*********************************************************************/
/* UAM  Header Routine                                               */
/*********************************************************************/
int areal_header(fd, message)
int fd;
char *message;		/* OUTPUT: return error string */
{
off_t offset;		/* file offset in bytes */
register int status;	/* routine return status */
register int i, j, k;	/* local loop variables */
static char buffer[MAXBUF];

int *ibuf = NULL;	/* integer value */
float *fbuf = NULL;	/* floating point value */

int utmzone;		/* UTM time zone */
int hour1;		/* starting hour # of file */
int hour2;		/* ending hour # of file */

float xgridsize;	/* UTM x grid cell size */
float ygridsize;	/* UTM y grid cell size */

status = FAILURE;	/* assume FAILURE */

/* 
 * File Description Header Record
 * word 1-10:	character	file name (10 characters, 1 character/word)
 * word 11-70:	character	file id (60 characters, 1 character/word)
 * word 71:	integer		number of segments (must be 1)
 * word 72:	integer		number of chemical species
 * word 73:	integer		Julian beginning date of the file
 * word 74:	float		beginning time of the file
 * word 75:	integer		Julian ending date of the file
 * word 76:	float		ending time of the file
 */

offset = 0;
if (get_areal_buf(fd, buffer, 15, offset, message))
	{
	ibuf = (int *) (buffer + 10 * 8);/* # of chem species */
	arealnspecies = *ibuf;

	/* beginning hour */
	fbuf = (float *) (buffer + 12 * 8);        
	hour1 = (int) *fbuf;
	
	/* ending hour */
	fbuf = (float *) (buffer + 14 * 8);
	hour2 = (int) *fbuf;

/* fix for 24 hours! needed if run starts at midnight! */

	if (hour2 == 0)
		hour2 = 23;

	arealnrecord = hour2 - hour1;
	--hour2;

/*
 * Region Description Header
 * word 1:	float		x-coordinate (UTM units)
 * word 2:	float		y-coordinate (UTM units)
 * word 3:	integer		UTM zone
 * word 4:	float		x-location (meters)
 * word 5:	float		y-location (meters)
 * word 6:	float		grid cell size in the x-direction (meters)
 * word 7:	float		grid cell size in the y-direction (meters)
 * word 8:	integer		number of grid cells in the x-direction
 * word 9:	integer		number of grid cells in the y-direction
 * word 10:	integer		number of grid cells in the z-direction
 * word 11:	integer		number of cells between surface layer and
 *					diffusion break
 * word 12:	integer		number of cells between diffusion break and
 *					top of region
 * word 13:	float		height of surface layer (meters)
 * word 14:	float		minimum height of cells between surface layer
 *					and diffusion break (meters)
 * word 15:	float		minimum height of cells between diffusion
 *					break and top of region (meters)
 */

/* Segment Description Header Record
 * word 1:	integer		x-location of the segment origin with
 *					respect to origin of the region
 * word 2:	integer		y-location of the segment origin with
 *					respect to origin of the region
 * word 3:	integer		number of grid cells in the segment in the
 *					x-direction
 * word 4:	integer		number of grid cells in the segment in the
 *					y-direction
 */

	offset = 16 * 8;	/* skip BCW + 15 words*/

	/* get region & segment description headers */

	if (get_areal_buf(fd, buffer, 19, offset, message))
		{

		/* # segment grid cells in x-dir */
		ibuf = (int *) (buffer + 17 * 8);
		arealncol = *ibuf;

		/* # segment grid cells in y-dir */
		ibuf = (int *) (buffer + 18 * 8);
		arealnrow = *ibuf;

		/* # grid cells in z-direction */
		if (arealfiletype == UAM_CONC)
			{
			ibuf = (int *) (buffer + 9 * 8);
			arealnlayer = *ibuf;
			}
		else if (arealfiletype == UAM_EMIS)
			arealnlayer = 1;		/* for emissions, 1 level */

		ibuf = (int *) (buffer + 2 * 8);	/* UTM zone */
		utmzone = *ibuf;

		fbuf = (float *) (buffer + 3 * 8);	/* UTM x sw corner */
		sw_utmx = *fbuf;

		fbuf = (float *) (buffer + 4 * 8);	/* UTM y sw corner */
		sw_utmy = *fbuf;
	
		fbuf = (float *) (buffer + 5 * 8);	/* x grid cell size */
		xgridsize = *fbuf;

		fbuf = (float *) (buffer + 6 * 8);	/* y grid cell size */
		ygridsize = *fbuf;

		ne_utmx = sw_utmx + arealncol * xgridsize;	/* set NE corner */
		ne_utmy = sw_utmy + arealnrow * ygridsize;

		areallevellen = 3 + (arealncol * arealnrow) + 1;
		arealspeclen = arealnlayer * areallevellen;
		arealsteplen = 4 + 1 + arealnspecies * arealspeclen;

/*
 * Species Description Header Record
 * word 1-10:	integer		species 1 (10 characters, 1 character/word) 
 * .........
 * .........
 * .........
 * word ..-10N:	integer		species N
 */

	/* skip BCW, 15 words + EOR, 15 words + EOR, 4 words */

		offset = 37 * 8;

	/* get species description header */

		if (get_areal_buf(fd, buffer, (arealnspecies*10)/8, offset, message))
			{
			k = 0;
			for (i = 0; i < arealnspecies; ++i)
				{
				for (j = 0; j < 4; ++j)
					arealspecname[k++] = buffer[i * 10 + j];
				if (i != arealnspecies - 1)
					{
					arealspecname[k] = ':';
					++k;
					}
				}
			arealspecname[k] = '\0';
	
			arealdatablock1 = 0;	/* if # species < MAXSPEC,
					 	headers stop in block 0 */

		/* word position of EOR at end of last header */

			arealdatapos1 = 37 + (arealnspecies * 10)/8 + 1;

	/* test filesize against # hours specified in header .. if descrepancy,
   	reduce hour2 and the number of steps:
	formula: (nwords - nBCW - headerwords)/arealsteplen */

			i = (arealfilesize/8 - arealfilesize/4096 - arealdatapos1)/arealsteplen;
			if (i < arealnrecord)
				{
				arealnrecord = i;
				hour2 = hour1 + arealnrecord;
				--hour2;
				}
			status = SUCCESS;
			}
		else
			sprintf(message,
				 "Error reading species description header.");
		}
	else
		sprintf(message,
			 "Error reading region & segment description headers.");
	}
else
	sprintf(message, "Error reading file description header.");
return(status);
}

int areal_get_date(fd, record_no, datestr, message)
int fd;
int record_no;
char *datestr;
char *message;
{
off_t offset;
int header_year;
int header_day;
int header_hour;
int *startdate = NULL;
float *starthour = NULL;
static char buffer[64];
register int i;
register int status;

/*
 * Time Step Header Record
 * word 1:	integer		Julian beginning date
 * word 2:	float		beginning time
 * word 3:	integer		Julian ending date
 * word 4:	float		ending time
 */

status = SUCCESS;
i = (record_no - 1) * arealsteplen;
if (!get_areal_offset(i, &offset, message))
	{
	sprintf(message, "Cannot get record offset for date.");
	status = FAILURE;
	}
else
	{
	if (get_areal_buf(fd, buffer, 2, offset, message))
		{
		startdate = (int *) (buffer);
		starthour = (float *) (buffer + 1 * 8);

                header_year = (*startdate)/1000;
                header_day = (*startdate) - (header_year * 1000);
		header_hour = (int) (*starthour + .5);

		julian2std(header_year + 1900, header_day,
			 header_hour, 0, 0, "EST", datestr);
		}
	else
		{
		sprintf(message, "Cannot get buffer for date.");
		status = FAILURE;
		}
	}
return(status);
}

int areal_get_data(fd, record_no, species_no, databuf, species, units, message)
int fd;
int record_no;
int species_no;
float *databuf;
char *species;
char *units;
char *message;
{

char *layerbuf = NULL;          /* character data layer buffer */
off_t offset;                   /* file offset in bytes */
float *fbuf = NULL;             /* floating point data value */
register int i, j, k;           /* local loop variables */
register int bindex;            /* index into input data layerbuf */
register int findex;            /* index into output field layerbuf */
register int ix, iy;            /* local loop variables for grid */
register int status;            /* routine return status */
 
status = SUCCESS;

i = (4 + arealncol * arealnrow) * 8;
if (	((layerbuf = ((char *) malloc(sizeof(char) * i))) != NULL)	)
	{
/*
 * Average Concentration Records
 *
 * ... within each time step 
 * ... ... within each layer
 * ... ... ... Average Concentration Record
 *
 * word 1:	integer		segment number (must be 1)
 * word 2-3:	character	species name (10 characters, 8 characters/word)
 * word 4-N:	float		concentrations averaged over time interval
 */

	i  = species_no - 1;
	k = 5 + (record_no - 1) * arealsteplen + i * arealspeclen;
	for (j = 0; ((j < arealnlayer) && (status == SUCCESS)); ++j) 
		{
		if (get_areal_offset(k, &offset, message))
			{
			if (get_areal_buf(fd, layerbuf, areallevellen - 1, offset,
					message))
				{
				k += areallevellen;
				bindex = 0;	/* layerbuf index */
				for (iy = 0; ((iy < arealnrow)
					 && (status == SUCCESS)); ++iy)
					{
					for (ix = 0; ((ix < arealncol)
						 && (status == SUCCESS)); ++ix)
						{
						fbuf = (float *) 
						(layerbuf + (3 + bindex) * 8);

/* calculate output data field index and store input value into the field */
 
			findex = j * arealnrow * arealncol + iy * arealncol + ix;
							databuf[findex] = *fbuf;
						++bindex;
						}
					 }
				 }
			else
				{
				sprintf(message, "Error reading data buffer.");
				status = FAILURE;
				}
			}
		else
			{
			sprintf(message, "Error getting data offset.");
			status = FAILURE;
			}
		}
	free ((char *) layerbuf);
	}
else
	{
	sprintf(message, "Error allocating data buffer.");
	status = FAILURE;
	}
sprintf(units, "%s", "PPM");

return(status);
}


