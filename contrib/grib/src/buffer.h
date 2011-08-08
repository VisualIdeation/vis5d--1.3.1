#define BUFF 68000

#ifdef MAIN_R
int buf[BUFF];
float ubuf[BUFF];
#else
extern int buf[BUFF];
extern float ubuf[BUFF];
#endif
