/* record.c */
 
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

#ifdef HAVE_LIBNETCDF
 
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "api.h"
#include "binio.h"
#include "record.h"
#include "globals.h"
#include "memory.h"
#include "imemory.h"
#include "proj.h"
#include "sync.h"
#include "irregular_v5d.h"
#include "v5d.h"
 
 
#ifndef M_PI
#  define M_PI 3.14159265
#endif
 
#define DEG2RAD (M_PI/180.0)
#define RAD2DEG (180.0/M_PI)
 


int init_record_cache( Irregular_Context itx, int maxbytes, float *ratio )
{
   int i, j;
   int numfloatvars, numcharvars, numberofchars; 
   int recordsize;
   int soundingsize;
   int totalnumrecs;

   /* dyamically allocate memory for RecordTable since */
   /* MAXTIMES*MAXRECS is way too large                */
   for (i = 0; i < itx->NumTimes; i++){
/* 01Feb06  Phil McDonald */
      itx->RecordTable[i] = (struct irreg_rec *) mem_calloc
                            (itx->NumRecs[i], sizeof (struct irreg_rec));
/* end PM */
   }

   ALLOC_LOCK(itx->Mutex);

   /*********************************/
   /* figure out size of one record */
   /*********************************/
   numfloatvars = 0;
   numcharvars = 0;
   numberofchars = 0;
   soundingsize = 0;
   for (i = 0; i < itx->NumVars; i++){
      if (itx->Variable[i]->VarType == NUMERICAL_VAR_1D ){
         numfloatvars++;
      }
      else if (itx->Variable[i]->VarType == NUMERICAL_VAR_2D ){
         soundingsize += itx->Levels;
      }
      else if (itx->Variable[i]->VarType == CHARACTER_VAR){
         numberofchars += itx->Variable[i]->CharVarLength;
         numcharvars++;
      }
      else{
         printf("Error in init_record_cache\n");
         return -1;
      }
   }
   if (itx->Type == SOUNDING){
      recordsize = numfloatvars * sizeof(double) +
                   soundingsize * sizeof(double) +
                   itx->Levels  * sizeof(float)  +
                   numberofchars  * sizeof(char);
   }
   else{
      recordsize = numfloatvars * sizeof(double) +
                   soundingsize * sizeof(double) +
                   numberofchars  * sizeof(char);
   }
   itx->CharArrayLength = numberofchars;


   /****************************************************/   
   /* figure out how many records can be put in memory */
   /****************************************************/
   itx->MaxCachedRecs = (int) maxbytes / recordsize;


   /**********************/
   /* get total num recs */
   /**********************/
   totalnumrecs = 0;
   for (i = 0; i < itx->NumTimes; i++){
      totalnumrecs += itx->NumRecs[i];
   }

   if (itx->MaxCachedRecs >= totalnumrecs){
      /* all records can be cached */
      itx->MaxCachedRecs = totalnumrecs;
      printf("Reading all records\n");
      *ratio = 1.0;
   }
   else{
      *ratio = ((float) itx->MaxCachedRecs)
             / ((float) totalnumrecs);
   }

   itx->NumCachedRecs = 0;

   printf("Cache size: %d records\n", itx->MaxCachedRecs);

   /****************************/
   /* allocate irregular cache */
   /****************************/
   itx->RecordCache = (struct cache_irreg_rec *) i_allocate( itx, itx->MaxCachedRecs
                                                 * sizeof(struct cache_irreg_rec) ); 

   /************/   
   /* check it */
   /************/
   if (!itx->RecordCache){
      printf("Error1: out of memory.  Couldn't allocate cache space.\n");
      return 0;
   }

   
   /******************************************************************/
   /* allocate memory for RecGeoPosition, this simply contains the   */
   /* lat, lon, alt for all the records.  This was removed from the  */
   /* cache_irreg_rec becuase if all the records can not fit into    */
   /* the cache, they will have to be loaded at some point in time   */
   /* into the cache in order to see whether or not to be used in a  */
   /* graphic, by loading all the lat, lon, and alts and keeping them*/
   /* this will eliminate this problem and the record cache will be  */
   /* more effective                                                 */
   /******************************************************************/

   for (i = 0; i < itx->NumTimes; i++){
      itx->RecGeoPosition[i] = (struct rec_geo_position *) i_allocate( itx,
                                                   itx->NumRecs[i] *
                                                   sizeof(struct rec_geo_position) );
      if (!itx->RecGeoPosition[i]){
         printf("Not enough memory to allocate for RecGeoPosition\n");
         return 0;
      }
   }


   itx->CacheClock = 1;

   for (i=0;i<itx->MaxCachedRecs;i++){
      /************************/
      /* alloc DataType array */
      /************************/
      itx->RecordCache[i].DataType =
                            i_allocate( itx, itx->NumVars*sizeof(int));
      /************/            
      /* check it */
      /************/
      if (!itx->RecordCache[i].DataType){
         printf("Error3: out of memory.  Couldn't allocate cache space.\n");
         return 0;
      }

      /***************/
      /* alloc Value */
      /***************/
      itx->RecordCache[i].Value =
                            i_allocate( itx, itx->NumVars*sizeof(double));
      /************/                        
      /* check it */
      /************/            
      if (!itx->RecordCache[i].Value){
         printf("Error4: out of memory.  Couldn't allocate cache space.\n");
         return 0;
      }

      /*************************/
      /* alloc Sounding Values */
      /*************************/
      if (soundingsize){
         itx->RecordCache[i].SoundingValue =
                           i_allocate( itx, soundingsize * sizeof(double));
         /************/
         /* check it */
         /************/
         if (!itx->RecordCache[i].SoundingValue){
            printf("Error5: out of memory.  Couldn't allocate cache space.\n");
            return 0;
         }
         itx->RecordCache[i].SoundingLevel = 
                           i_allocate( itx, itx->Levels * sizeof(float));
         if(!itx->RecordCache[i].SoundingLevel){
            printf("Error6: out of memory.  Couldn't allocate cache space.\n");
            return 0;
         }
      }

      /***********************/            
      /* Alloc char var data */
      /***********************/
      itx->RecordCache[i].CharData =
                            i_allocate( itx, numberofchars*sizeof(char ));

      /************/                                    
      /* check it */            
      /************/                        
      if (!itx->RecordCache[i].CharData){
         printf("Error7: out of memory.  Couldn't allocate cache space.\n");
         return 0;
      }


      itx->RecordCache[i].Locked = 0;
      itx->RecordCache[i].Timestep = 0;
   }


   for (i=0; i<itx->NumTimes; i++){
      for (j = 0; j < itx->NumRecs[i]; j++){
         itx->RecordTable[i][j].CachePos = -1;
         itx->RecordTable[i][j].DataType = NULL;
         itx->RecordTable[i][j].Value = NULL;
         itx->RecordTable[i][j].SoundingValue = NULL;
         itx->RecordTable[i][j].SoundingLevel = NULL;
         itx->RecordTable[i][j].CharData = NULL;
      }
   }
   return 1;
}

int get_empty_irreg_cache_pos(Irregular_Context itx)
{
   int g;

   /* fing g */
   if (itx->NumCachedRecs < itx->MaxCachedRecs ){
      /* There's an unused position. */
      g = itx->NumCachedRecs;
      itx->NumCachedRecs++;
   }
   else{
      int time, rec, minage, i, mini;
      minage = itx->CacheClock;
      for (i=0;i<itx->MaxCachedRecs;i++) {
         if (itx->RecordCache[i].Age<minage && itx->RecordCache[i].Locked==0) {
            minage = itx->RecordCache[i].Age;
            mini = i;
         }
      }
      g = mini;
      /* remove references to data being discarded */
      time = itx->RecordCache[g].Timestep;
      rec = itx->RecordCache[g].Rec;
      itx->RecordTable[time][rec].DataType = NULL;
      itx->RecordTable[time][rec].Value = NULL;
      itx->RecordTable[time][rec].SoundingValue = NULL;
      itx->RecordTable[time][rec].SoundingLevel = NULL;
      itx->RecordTable[time][rec].CharData = NULL;
      itx->RecordTable[time][rec].CachePos = -1;
   }
   itx->RecordCache[g].Locked = 1;
   return g;
}


static void load_record(Irregular_Context itx, int time, int record)
{
   int p;

   LOCK_ON( itx->Mutex );

   if (itx->RecordTable[time][record].CachePos > 0){
      /* already in the cache */
      p = itx->RecordTable[time][record].CachePos;
      if (p >= 0){
         itx->RecordCache[p].Locked = 1;
         itx->RecordCache[p].Age = itx->CacheClock++;
         itx->RecordCache[p].Locked = 0;
      }
      LOCK_OFF( itx->Mutex );
      return;
   }
   else{
      /* not in cache */
      int g;
      g = get_empty_irreg_cache_pos(itx);
      if (!irregular_v5dReadRecord( &itx->G, time, record,
          itx->RecordCache[g].Value,
          itx->RecordCache[g].SoundingValue,
          itx->RecordCache[g].CharData,
          itx->RecordCache[g].SoundingLevel)){
         printf("Error: unable to read record information\n");
         LOCK_OFF( itx->Mutex );
         itx->RecordCache[g].Locked = 0;
         return;
      }
      itx->RecordTable[time][record].DataType = itx->RecordCache[g].DataType;
      itx->RecordTable[time][record].Value = itx->RecordCache[g].Value;
      itx->RecordTable[time][record].SoundingValue = itx->RecordCache[g].SoundingValue;
      itx->RecordTable[time][record].SoundingLevel = itx->RecordCache[g].SoundingLevel;
      itx->RecordTable[time][record].CharData = itx->RecordCache[g].CharData;

      itx->RecordTable[time][record].CachePos = g;
      itx->RecordCache[g].Locked = 1;
      itx->RecordCache[g].Timestep = time;
      itx->RecordCache[g].Rec = record;
      itx->RecordCache[g].Age = itx->CacheClock++;
      LOCK_OFF(itx->Mutex);
      itx->RecordCache[g].Locked = 0;
      return;
   }
}

void get_all_record_numerical_data( Irregular_Context itx, int time,
                                int var, double *data)
{
   int i;

   if (itx->Variable[var]->VarType != NUMERICAL_VAR_1D){
      printf("wrong var type in get_record_numerical_data\n");
      return;
   }

   for (i = 0; i < itx->NumRecs[time]; i++){
      if (itx->RecordTable[time][i].CachePos < 0){
         load_record(itx, time, i);
      }
      data[i] = itx->RecordTable[time][i].Value[var];
   }
}

void get_some_record_numerical_data( Irregular_Context itx, int time, 
                                int var, int ploton[], double *data)
{
   int i, pcount = 0;

   if (itx->Variable[var]->VarType != NUMERICAL_VAR_1D){
      printf("wrong var type in get_record_numerical_data\n");
      return;
   }

   for (i = 0; i < itx->NumRecs[time]; i++){
      if (ploton[i]){
         if (itx->RecordTable[time][i].CachePos < 0){
            load_record(itx, time, i);
         }   
         data[pcount] = itx->RecordTable[time][i].Value[var];
         pcount++;
      }
   }
}

void get_all_record_char_data( Irregular_Context itx, int time,
                                int var, char *data)
{
   int i, j, count=0;
   
   if (itx->Variable[var]->VarType != CHARACTER_VAR){
      printf("wrong var type in get_record_char_data\n");
      return;
   }

   for (i = 0; i < itx->NumRecs[time]; i++){
      if (itx->RecordTable[time][i].CachePos < 0){
         load_record(itx, time, i);
      }
      for (j = itx->Variable[var]->CharPointer; j < itx->Variable[var]->CharPointer+
                                          itx->Variable[var]->CharVarLength; j++){
         data[count] = itx->RecordTable[time][i].CharData[j];
         count++;
      }
   }
}

            
void get_some_record_char_data( Irregular_Context itx, int time,
                                int var, int ploton[], char *data)
{
   int i, j, count=0;
   
   if (itx->Variable[var]->VarType != CHARACTER_VAR){
      printf("wrong var type in get_record_char_data\n");
      return;
   }

   for (i = 0; i < itx->NumRecs[time]; i++){
      if (ploton[i]){
         if (itx->RecordTable[time][i].CachePos < 0){
            load_record(itx, time, i);
         }
         for (j = itx->Variable[var]->CharPointer; j < itx->Variable[var]->CharPointer+
                                             itx->Variable[var]->CharVarLength; j++){
            data[count] = itx->RecordTable[time][i].CharData[j];
            count++;
         }
      }
   }
}

            

void get_record_location( Irregular_Context itx, int time, int rec,
                           float *lat, float *lon, float *alt)
{
   *lat = itx->RecGeoPosition[time][rec].Latitude;
   *lon = itx->RecGeoPosition[time][rec].Longitude*-1.0;
   *alt = itx->RecGeoPosition[time][rec].Altitude/1000.0;
}

void get_record_locations( Irregular_Context itx, int time, 
                           float *lat, float *lon, float *alt)
{
   int i;
   
   for (i = 0; i < itx->NumRecs[time]; i++){
      lat[i] = itx->RecGeoPosition[time][i].Latitude;
      lon[i] = itx->RecGeoPosition[time][i].Longitude*-1.0;
      alt[i] = itx->RecGeoPosition[time][i].Altitude/1000.0;
   }
}


void load_geo_data( Irregular_Context itx )
{
   int i, j;
   float lat, lon, alt;

   for (i = 0; i < itx->NumTimes; i++){
      for (j = 0; j < itx->NumRecs[i]; j++){
         if (!irregular_v5dReadRecordGeoData( &itx->G, i, j, &lat, &lon, &alt)){
            printf("Error in reading Geo Data\n");
            return;
         }
         itx->RecGeoPosition[i][j].Latitude = lat;
         itx->RecGeoPosition[i][j].Longitude = lon;
         itx->RecGeoPosition[i][j].Altitude = alt;
      }
   }
}


void preload_irregular_cache( Irregular_Context itx )
{
   int numrecs = 0; 
   int time = 0;
   int trec = 0;
	/*
   printf("Loading Records for time %d\n", time);
	*/
   while ( numrecs < itx->MaxCachedRecs){
      if (trec == itx->NumRecs[time]){
         trec = 0;
         time++;
			/*
         printf("Loading Records for time %d\n", time);   
			*/
      }
      load_record(itx, time, trec);
      numrecs++;
      trec++;
   }
}


int initially_open_recordfile( char filename[], irregular_v5dstruct *iv )
{
   char name[1000];

   strcpy( name, filename);

   if (!irregular_v5dOpenFile( name, iv)){
      printf("Error: datafile %s could not be loaded\n", filename);
      return 0;
   }
   return 1;
}




int open_recordfile(Irregular_Context itx, char filename[])
{
   int i, time, first;
 
   if (!initially_open_recordfile( filename, &itx->G)){
      return 0;
   }

  
   /* Copy header info from G to global variables */
   strcpy(itx->DataFile, filename);

   itx->Type = itx->G.Type;
   itx->Levels = itx->G.Levels; 
   itx->NumVars = itx->G.NumVars;
   itx->NumTimes = itx->G.NumTimes;
   itx->TopBound = itx->G.TopBound;
   itx->BottomBound = itx->G.BottomBound;
   itx->WestBound = itx->G.WestBound;
   itx->EastBound = itx->G.EastBound;
   itx->NorthBound = itx->G.NorthBound;
   itx->SouthBound = itx->G.SouthBound;
	
   for (i = 0; i < itx->NumVars; i++){
	  itx->Variable[i] = (vis5d_irregular_variable *) i_allocate(itx,sizeof(vis5d_irregular_variable));
	  strcpy(itx->Variable[i]->VarName, itx->G.VarName[i]);
	  itx->Variable[i]->VarType = itx->G.VarType[i];
	  itx->Variable[i]->CharVarLength = itx->G.CharVarLength[i];
	  itx->Variable[i]->CharPointer = itx->G.CharPointer[i];
	  itx->Variable[i]->SoundingPointer = itx->G.SoundingPointer[i];
	  itx->Variable[i]->MinVal = itx->G.VarMin[i];
	  itx->Variable[i]->MaxVal = itx->G.VarMax[i];

	}
  itx->TopBound = 10.0;
  itx->BottomBound = -0.1;
  if (itx->WestBound == itx->EastBound){
      itx->WestBound += 10.0;
      itx->EastBound -= 10.0;
   }
   if (itx->SouthBound == itx->NorthBound){
      itx->SouthBound -= 10.0;
      itx->NorthBound += 10.0;
   }

   /* do some checking */
   if (itx->NumVars > MAXVARS){
      printf("Error: Too many variables (%d) limit is %d\n",
             itx->NumVars, MAXVARS);
      return 0;
   }
   if (itx->NumTimes >MAXTIMES) {
      printf("Error: Too many time steps (%d) limit is %d\n", 
              itx->NumTimes, MAXTIMES);
      return 0;
   }

   for (time = 0; time < itx->NumTimes; time++){
      itx->TimeStamp[time] = v5dHHMMSStoSeconds( itx->G.TimeStamp[time] );
      itx->DayStamp[time] = v5dYYDDDtoDays( itx->G.DateStamp[time] );
      itx->NumRecs[time] = itx->G.NumRecs[time];
   }

   /* calculate elapsed time (in seconds) for each time since initial time */
   first = itx->DayStamp[0]* 24*60*60 + itx->TimeStamp[0];
   for (time=0;time<itx->NumTimes;time++) {
      itx->Elapsed[time] = itx->DayStamp[time] * 24*60*60
                         + itx->TimeStamp[time] - first;
   }

   return 1;
}

#endif /* HAVE_LIBNETCDF */
