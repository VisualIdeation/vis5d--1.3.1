#include <stdio.h>
int bword (char *buf, int bloc, int bnum)
{
    int        i;
    int        word = 0;
    char       *byte = &word;



/*
	fprintf (stderr, "bword: %p %d %d: ", buf, bloc, bnum);
	for (i = 0; i < bnum; i++) fprintf (stderr, " %2.2x", (unsigned char) buf[bloc+i]);
*/
    buf += bloc;
#ifdef BIG_ENDIAN
    byte += bnum;
    for (i = 0; i < bnum; i++) *(--byte) = *(buf++);
#else
    byte += 4 - bnum;
    for (i = 0; i < bnum; i++) *(byte++) = *(buf++);
#endif
/*
	fprintf (stderr, ": %8.8x %d\n", word, word);
*/


    return word;
}
