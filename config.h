/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.in by autoheader.  */

/* Define to installed location of Vis5d data files. */
#define DATA_PREFIX "${datarootdir}/vis5d+/"

/* Define to 1 if translation of program messages to the user's native
   language is requested. */
#define ENABLE_NLS 1

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef F77_DUMMY_MAIN */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
/* #undef F77_FUNC */

/* As F77_FUNC, but for C identifiers containing underscores. */
/* #undef F77_FUNC_ */

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* Define if the GNU dcgettext() function is already present or preinstalled.
   */
#define HAVE_DCGETTEXT 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you have the <fcntl.h> header file. */
#define HAVE_FCNTL_H 1

/* Define if the GNU gettext() function is already present or preinstalled. */
#define HAVE_GETTEXT 1

/* Define to 1 if you have the <GL/gl.h> header file. */
#define HAVE_GL_GL_H 1

/* Define if you have the iconv() function. */
/* #undef HAVE_ICONV */

/* Do we have Fortran idate function? */
/* #undef HAVE_IDATE */

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the <libintl.h> header file. */
#define HAVE_LIBINTL_H 1

/* Define to 1 if you have the `m' library (-lm). */
#define HAVE_LIBM 1

/* Define to 1 if you have the `mcidas_lib' library (-lmcidas_lib). */
/* #undef HAVE_LIBMCIDAS_LIB */

/* Do we have NetCDF library? */
#define HAVE_LIBNETCDF 1

/* Define to 1 if you have the `png' library (-lpng). */
#define HAVE_LIBPNG 1

/* Have Tcl library? */
/* #undef HAVE_LIBTCL */

/* Define to 1 if you have the `z' library (-lz). */
#define HAVE_LIBZ 1

/* Do we have McIDAS library? */
/* #undef HAVE_MCIDAS */

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Do we have mixkit library? */
/* #undef HAVE_MIXKIT */

/* Define to 1 if you have the <mixkit/mixio.h> header file. */
/* #undef HAVE_MIXKIT_MIXIO_H */

/* Define to 1 if you have the <netcdf.h> header file. */
#define HAVE_NETCDF_H 1

/* Define if we have OpenGL libraries. */
#define HAVE_OPENGL 1

/* Define if you have POSIX threads libraries and header files. */
/* #undef HAVE_PTHREAD */

/* Have POSIX threads? */
/* #undef HAVE_PTHREADS */

/* Define to 1 if you have the `setrlimit' function. */
#define HAVE_SETRLIMIT 1

/* Define if we have SGI GL libraries. */
/* #undef HAVE_SGI_GL */

/* Have SGI sproc? */
/* #undef HAVE_SGI_SPROC */

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the `strcasecmp' function. */
#define HAVE_STRCASECMP 1

/* Define to 1 if you have the `strdup' function. */
#define HAVE_STRDUP 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the `strncasecmp' function. */
#define HAVE_STRNCASECMP 1

/* Have SunOS threads? */
/* #undef HAVE_SUNOS_THREADS */

/* Define to 1 if you have the <sysmp.h> header file. */
/* #undef HAVE_SYSMP_H */

/* Define to 1 if you have the <sys/lock.h> header file. */
/* #undef HAVE_SYS_LOCK_H */

/* Define to 1 if you have the <sys/prctl.h> header file. */
#define HAVE_SYS_PRCTL_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/sysmp.h> header file. */
/* #undef HAVE_SYS_SYSMP_H */

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if you have the <X11/Xm/MwmUtil.h> header file. */
/* #undef HAVE_X11_XM_MWMUTIL_H */

/* Define to 1 if you have the `XMesaGetBackBuffer' function. */
/* #undef HAVE_XMESAGETBACKBUFFER */

/* path to Image Magick convert program */
#define IMCONVERT "/usr/bin/convert"

/* Have McIDAS sidecar? */
/* #undef MCIDAS_SIDECAR */

/* Name of package */
#define PACKAGE "vis5d+"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT ""

/* Define to the full name of this package. */
#define PACKAGE_NAME ""

/* Define to the full name and version of this package. */
#define PACKAGE_STRING ""

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME ""

/* Define to the version of this package. */
#define PACKAGE_VERSION ""

/* Define to necessary symbol if this constant uses a non-standard name on
   your system. */
/* #undef PTHREAD_CREATE_JOINABLE */

/* Define to disable multi-threading. */
#define SINGLE_TASK 1

/* The size of `float', as computed by sizeof. */
#define SIZEOF_FLOAT 4

/* The size of `int', as computed by sizeof. */
#define SIZEOF_INT 4

/* The size of `signed char', as computed by sizeof. */
#define SIZEOF_SIGNED_CHAR 1

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Version number of package */
#define VERSION "1.3.1"

/* directory containing message catalogs */
#define VIS5D_LOCALE_DIR "/usr/local//locale"

/* max. memory to use (MB), 0 for no maximum */
#define VIS5D_MAX_MEM 0

/* where the source was built */
#define VIS5D_SOURCE_DIR "/home/demo/workspace/vis5d--1.3.1"

/* Define to 1 if your processor stores words with the most significant byte
   first (like Motorola and SPARC, unlike Intel and VAX). */
/* #undef WORDS_BIGENDIAN */

/* Define to 1 if the X Window System is missing or not being used. */
/* #undef X_DISPLAY_MISSING */
