/************************************************************************
*									*
*    Data Windowing							*
*									*
*					1996.1	N.Miyazaki		*
*									*
*	usage:								*
*		window [ infile ] [ options ] > outfile			*
*	options:							*
*		-l l	 :  input frame length		[256]		*
*		-L L	 :  output frame length		[l]		*
*		-n n	 :  type of normalization	[1]		*
*			n=0: none					*
*			n=1: normalize by power				*
*			n=2: normalize by magnitude			*
*		-w w	 :  type of window				*
*			w=0: blackman  window				*
*			w=1: hamming   window				*
*			w=2: hanning   window				*
*			w=3: bartlett  window				*
*			w=4: trapezoid window				*
*	infile:								*
*		stdin for default					*
*		input is assumed to be double				*
*	notice:								*
*		if L > l, (L-l)-zeros are padded			*
*									*
************************************************************************/
static char *rcs_id = "$Id: window-main.c,v 1.1 1996/01/25 12:47:30 nmiya Exp nmiya $";


/*  Standard C Libraries  */
#include <stdio.h>
#include <SPTK.h>
#include <string.h>


/*  Required Function */
double	window();


/*  Command Name  */
char*	cmnd;

int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, " %s - data windowing\n", cmnd);
	fprintf(stderr, "\n");
	fprintf(stderr, "  usage:\n"); 
	fprintf(stderr, "       %s [ options ] [ infile ] > outfile\n", cmnd); 
	fprintf(stderr, "  options:\n"); 
	fprintf(stderr, "       -l l  : frame length of input  [256]\n");
	fprintf(stderr, "       -L L  : frame length of output [l]\n");
	fprintf(stderr, "       -n n  : type of normalization  [1]\n");
	fprintf(stderr, "                 0 none\n");
	fprintf(stderr, "                 1 normalize by power\n");
	fprintf(stderr, "                 2 normalize by magnitude\n");
	fprintf(stderr, "       -w w  : type of window         [0]\n");
	fprintf(stderr, "                 0 (blackman)\n");
	fprintf(stderr, "                 1 (hamming)\n");
	fprintf(stderr, "                 2 (hanning)\n");
	fprintf(stderr, "                 3 (bartlett)\n");
	fprintf(stderr, "                 4 (trapezoid)\n");
	fprintf(stderr, "       -h    : print this message\n");
	fprintf(stderr, "  infile:\n"); 
	fprintf(stderr, "       data sequence (float)          [stdin]\n"); 
	fprintf(stderr, "  stdout:\n"); 
	fprintf(stderr, "       windowed data sequence (float)\n"); 
	fprintf(stderr, "\n");
	exit(1);
}


main(argc,argv)
int	argc;
char	*argv[];
{
	FILE	*fp = stdin;
	char	*s, c;
	int	l=256, outl = 0, zerol, normal=1, wintype;
	double  *x, *zero;	

	wintype = 0;

        if ((cmnd = strrchr(argv[0], '/')) == NULL)
	    cmnd = argv[0];
        else
	    cmnd++;

	while (--argc){
		if(*(s = *++argv) == '-') {
			c = *++s;
			if(*++s == '\0' && (c == 'n' || c == 'l' || c == 'w' || c == 'L')) {
				s = *++argv;
				--argc;
			}
			switch(c) {
			case 'w':
				wintype = atoi(s);
				break;
			case 'l':
				l = atoi(s);
				break;
			case 'L':
				outl = atoi(s);
				break;
			case 'n':
				normal = atoi(s);
				break;
			case 'h':
			default:
				usage();
			}
		}
		else
			 fp = getfp(*argv, "r");
	}

	if(outl < l)
	    outl = l;

	zerol = outl - l;
	
	x = dgetmem(l);

	if(zerol > 0){
	    zero = dgetmem(zerol);
	    
	    while(freadf(x, sizeof(*x), l, fp) == l) {
		window(windows[wintype], x, l, normal);
		fwritef(x, sizeof(*x), l, stdout);
		fwritef(zero, sizeof(*zero), zerol, stdout);
	    }
	}
	else
	    while(freadf(x, sizeof(*x), l, fp) == l) {
		window(windows[wintype], x, l, normal);
		fwritef(x, sizeof(*x), l, stdout);
	    }
	exit(0);
}
