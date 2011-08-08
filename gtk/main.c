/*
 * Initial main.c file generated by Glade. Edit as required.
 * Glade will not overwrite this file.
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif
#include <gtk/gtk.h>
#include <locale.h>
#include <mcheck.h>
#include <signal.h>

#include "window3D.h"
#include "support.h"


static void
enable (int sig)
{
  mtrace ();
  signal (SIGUSR1, enable);
}

static void
disable (int sig)
{
  muntrace ();
  signal (SIGUSR2, disable);
}




int
main (int argc, char *argv[])
{ 
  GtkWidget *window3D;

  /* MALLOC_TRACE for gnu systems */

  signal (SIGUSR1, enable);
  signal (SIGUSR2, disable);

  /*
  mtrace ();
  */
  gtk_set_locale ();

  /* This is necessary so that numeric text files are read in 
	  properly. 
  */

  setlocale(LC_NUMERIC,"C");

  /* TODO: need to look for gtkrc file, currently we just look in 
	  the source dir since we don't install this version anyway
  */
  gtk_rc_add_default_file(VIS5D_SOURCE_DIR "/gtk/gtkrc");


  gtk_init (&argc, &argv);

  
  window3D = new_window3D(NULL);

  /* when we handle command line options we will change this */
  on_open1_activate(NULL, (gpointer) window3D);

  gtk_main ();


  return 0;
}



