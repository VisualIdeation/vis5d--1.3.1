/* sync.c */


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


/* A library of concurrency functions and macros usable in C on Stellar,
   SGI, and CRAY systems.

   Note that when declaring LOCKs or SEMAPHOREs, only one can be declared
   per line.  For example:

        LOCK a;      this is correct
        LOCK b;

        LOCK a, b;   this is incorrect
*/


#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "sync.h"


#ifdef HAVE_SGI_SPROC
#  include <unistd.h>
#  include <ulocks.h>
   usptr_t *Arena;
   static char *ArenaName;
#endif



/* 24Apr07  Phil McDonald */
#ifdef HAVE_PTHREADS
static int	(*sync_mutex_alloc_func) (pthread_mutex_t *p_mutex,
                                          const pthread_mutexattr_t *p_mutex_attr)
						= pthread_mutex_init;

static int	(*sync_mutex_free_func) (pthread_mutex_t *p_mutex)
						= pthread_mutex_destroy;

static int	(*sync_mutex_lock_func) (pthread_mutex_t *p_mutex)
						= pthread_mutex_lock;

static int	(*sync_mutex_unlock_func) (pthread_mutex_t *p_mutex)
						= pthread_mutex_unlock;



int	sync_mutex_alloc (pthread_mutex_t *p_mutex,
                          const pthread_mutexattr_t *p_mutex_attr)
{
    int	iret;
    iret = sync_mutex_alloc_func (p_mutex, p_mutex_attr);
    return iret;
}



int	sync_mutex_free (pthread_mutex_t *p_mutex)
{
    int	iret;
    iret = sync_mutex_free_func (p_mutex);
    return iret;
}



int	sync_mutex_lock (pthread_mutex_t *p_mutex)
{
    int	iret;
    iret = sync_mutex_lock_func (p_mutex);
    return iret;
}



int	sync_mutex_unlock (pthread_mutex_t *p_mutex)
{
    int	iret;
    iret = sync_mutex_unlock_func (p_mutex);
    return iret;
}



int	sync_mutex_set_funcs (
	    int (*alloc_func) (pthread_mutex_t *p_mutex,
	                       const pthread_mutexattr_t *p_mutex_attr),
	    int (*free_func) (pthread_mutex_t *p_mutex),
	    int (*lock_func) (pthread_mutex_t *p_mutex),
	    int (*unlock_func) (pthread_mutex_t *p_mutex)
	)
{
    if (alloc_func == NULL) alloc_func = pthread_mutex_init;
    sync_mutex_alloc_func = alloc_func;

    if (free_func == NULL) free_func = pthread_mutex_destroy;
    sync_mutex_free_func = free_func;

    if (lock_func == NULL) lock_func = pthread_mutex_lock;
    sync_mutex_lock_func = lock_func;

    if (unlock_func == NULL) unlock_func = pthread_mutex_unlock;
    sync_mutex_unlock_func = unlock_func;


    return 0;
}
#endif
/* end PM */



/*
 * The lock used to implement read/write locks below:
 */
static LOCK RWlock;



/*
 * This function should be called before using any of sychronization
 * primitives below.
 * Return:  1 = success, 0 = error;
 */
int init_sync( void )
{
#ifdef HAVE_SGI_SPROC
   ArenaName = tempnam( "", "LOCK" );
   if (!ArenaName) {
      fprintf(stderr,"Unable to allocate arena.\n");
      return 0;
   }
   /*
    * Increase CONF_INITUSERS to 20 to avoid the error:
    * New process pid 20029 could not join I/O arena:No space left on device
    * Need to raise CONF_INITUSERS? (usconfig(3P))
    */
   usconfig( CONF_INITUSERS, 20 );

   /*
    * This options causes the temporary lock file to be unlinked automatically.
    */
   usconfig( CONF_ARENATYPE, US_SHAREDONLY );

   /*
    * Initialize the shared arena
    */
   Arena = usinit(ArenaName);
   if (!Arena) {
      fprintf(stderr,"Unable to allocate arena.\n");
      return 0;
   }
#endif
   ALLOC_LOCK( RWlock );
   return 1;
}



/*** term_sync ********************************************************
   This function should be called prior to program termination to
   release resources associated with the synchronization primitives.
**********************************************************************/
void term_sync( void )
{
#ifdef HAVE_SGI_SPROC
   if (Arena) {
      unlink( ArenaName );
      Arena = NULL;
   }
#endif
}




/********************************/
/*** CRAY SEMAPHORE FUNCTIONS ***/
/********************************/

#ifdef cray

/* Semaphore value (n):
  if n<=0 then
       the semaphore is held; calls to wait will be blocked
  if n<0 then
       abs(n) is the number of waiters
  if n>0 then
       the semaphore is available; n calls to wait will not block.
*/



void alloc_sem( s, initval )
struct sem *s;
int initval;
{
   if (s) {
      LOCKASGN( &(s->mutex) );
      EVASGN( &(s->ev) );
      s->val = initval;
      s->evflag = 0;
      s->awaken = 0;
   }
}



void free_sem( s )
struct sem *s;
{
   if (s) {
      LOCKREL( &(s->mutex) );
      EVREL( &(s->ev) );
   }
}


/* changed to use awaken value: */


void wait_sem( s )
struct sem *s;
{
   LOCKON( &(s->mutex) );

   s->val--;

   while (s->val<0) {
      /* wait for a signal */
      LOCKOFF( &(s->mutex) );
      EVWAIT( &(s->ev) );   /* block until awaken>0 */
      LOCKON( &(s->mutex) );
      /* do EVCLEAR exactly once. */
      if (s->evflag) {
         EVCLEAR( &(s->ev) );
         s->evflag = 0;
      }

      if (s->awaken>0) {
         s->awaken--;
         LOCKOFF( &(s->mutex) );
         return;
      }
   }

   LOCKOFF( &(s->mutex) );
}



void signal_sem( s )
struct sem *s;
{
   LOCKON( &(s->mutex) );

   s->val++;
   if (s->val <= 0) {
      s->awaken++;
      if (s->evflag==0) {
         EVPOST( &(s->ev) );
         s->evflag = 1;
      }
   }

   LOCKOFF( &(s->mutex) );
}


#endif /*cray*/



/**********************************************/
/*** STELLAR / STARDENT SEMAPHORE FUNCTIONS ***/
/**********************************************/


#ifdef stellar

void alloc_sem( s, n )
struct sem *s;
int n;
{
   s->val = 0;
   s->mutex = 0;
   s->awaken = 0;
}



void wait_sem( s )
struct sem *s;
{
   LOCK_ON( s->mutex );

   s->val--;

   while (s->val<0) {
      /* wait for a signal */
      LOCK_OFF( s->mutex );
      wait_until_gt( &(s->awaken), 0 );  /* block until awaken>0 */
      LOCK_ON( s->mutex );
      if (s->awaken>0) {
         s->awaken--;
         LOCK_OFF( s->mutex );
         return;
      }
   }

   LOCK_OFF( s->mutex );
}



void signal_sem( s )
struct sem *s;
{
   LOCK_ON( s->mutex );

   s->val++;
   if (s->val <= 0)
      s->awaken++;  /* wake-up another blocked process */

   LOCK_OFF( s->mutex );
}

#endif /*stellar*/



/***********************************/
/*** SEQUENT SEMAPHORE FUNCTIONS ***/
/***********************************/


#ifdef sequent

#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>


int alloc_sem( initval )
int initval;
{
   key_t  key;
   int  s, val;

   key = 0;
   s = semget( key, 1, 0666 | IPC_CREAT );
   if (s<0) {
      perror("getting semaphore");
      return -1;
   }
   else {
      if (semctl( s, 0, SETVAL, initval ) < 0) {
         perror("initializing semaphore");
      }
      if (semctl( s, 0, GETVAL, 0 ) < 0) {
         perror("initializing semaphore");
      }
   }
   return s;
}



void free_sem( s )
int s;
{
   if (semctl( s, 0, IPC_RMID, 0 ) < 0) {
      perror("destroying semaphore");
   }
}




void wait_sem( s )
int s;
{
   struct sembuf sb;

   sb.sem_num = 0;
   sb.sem_op = -1;
   sb.sem_flg = SEM_UNDO;

   semop( s, &sb, 1 );
}



void signal_sem( s )
int s;
{
   struct sembuf sb;

   sb.sem_num = 0;
   sb.sem_op = 1;
   sb.sem_flg = SEM_UNDO;

   semop( s, &sb, 1 );
}

#endif /*sequent*/





/*
 * Read/write locking on an integer:
 *  if n==0 then lock is free
 *  if n>0 then there are (n) readers
 *  if n==-1 then there is a writer
 */


int cond_read_lock( int *lk )
{
   int result;
   LOCK_ON( RWlock );
   if (*lk>=0) {
      *lk = *lk + 1;
      result = 1;
   }
   else {
      result = 0;
   }
   LOCK_OFF( RWlock );
   return result;
}


void wait_read_lock( int *lk )
{
   int acquired = 0;
   do {
      LOCK_ON( RWlock );
      if (*lk>=0) {
         *lk = *lk + 1;
         acquired = 1;
      }
      LOCK_OFF( RWlock );
   } while (!acquired);
}


void wait_write_lock( int *lk )
{
   int acquired = 0;
   do {
      LOCK_ON( RWlock );
      if (*lk==0) {
         *lk = -1;
         acquired = 1;
      }
      LOCK_OFF( RWlock );
   } while (!acquired);
}


void done_read_lock( int *lk )
{
   LOCK_ON( RWlock );
/* 03Feb06  Phil McDonald */
/*
   assert( *lk > 0 );
*/
   if (*lk <= 0)
   {
      fprintf (stderr, "sync.c: pid=%d   %s    %d at %p   UNBALANCED\n",
               getpid(), "done_read_lock", *lk, lk);
      *lk = 1; /* TODO: this should never happen */
   }
/* end PM */
   *lk = *lk - 1;
   LOCK_OFF( RWlock );
}


void done_write_lock( int *lk )
{
   LOCK_ON( RWlock );
/* 03Feb06  Phil McDonald */
/*
   assert( *lk==-1 );
*/
   if (*lk != -1)
   {
      fprintf (stderr, "sync.c: pid=%d   %s    %d at %p   UNBALANCED\n",
               getpid(), "done_write_lock", *lk, lk);
   }
/* end PM */
   *lk = 0;
   LOCK_OFF( RWlock );
}
