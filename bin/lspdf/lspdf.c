/************************************************************************
*									*
*    LSP Digital Filter for Speech Synthesis				*
*									*
*					1996.9  K.Koishida		*
*									*
*	usage:								*
*		lspdf [ options ] [ infile ] > stdout			*
*	options:							*
*		-m m     :  order of coefficients	[24]		*
*		-p p     :  frame period		[100]		*
*		-i i     :  interpolation period	[1]		*
*		-s s	 :  sampling frequency		[10]		*
*		-o o	 :  input format		[0]		*
*	infile:								*
*		coefficients						*
*             	    input format      LSP                               *
*                     	   0    normalized frequency (0 ~ pi)           *
*                   	   1    normalized frequency (0 ~ 0.5)          *
*                   	   2    frequency (kHz)                         *
*                   	   3    frequency (Hz)                          *
*                   LSP                                                 *
*		           , K, f(1), ..., f(m),			*
*		excitation sequence					*
*		    , x(0), x(1), ..., 					*
*	stdout:								*
*		filtered sequence					*
*		    , y(0), y(1), ...,					*
*	require:							*
*		lspdf_even()						*
*									*
************************************************************************/

static char *rcs_id = "$Id: lspdf-main.c,v 1.1 1996/09/17 05:07:18 koishida Exp koishida $";


/*  Standard C Libraries  */
#include <stdio.h>
#include <string.h>
#include <SPTK.h>


/*  Required Functions  */
double	lspdf_even();


/*  Default Values  */
#define ORDER		24
#define	FPERIOD		100
#define	IPERIOD		1
#define ITYPE		0
#define SAMPLING	10


/*  Command Name  */
char	*cmnd;


void usage(int status)
{
    fprintf(stderr, "\n");
    fprintf(stderr, " %s - LSP digital filter for speech synthesis\n",cmnd);
    fprintf(stderr, "\n");
    fprintf(stderr, "  usage:\n");
    fprintf(stderr, "       %s [ options ] lspfile [ infile ] > stdout\n", cmnd);
    fprintf(stderr, "  options:\n");
    fprintf(stderr, "       -m m  : order of coefficients [%d]\n", ORDER);
    fprintf(stderr, "       -p p  : frame period          [%d]\n", FPERIOD);
    fprintf(stderr, "       -i i  : interpolation period  [%d]\n", IPERIOD);
    fprintf(stderr, "       -s s  : sampling frequency    [%d]\n", SAMPLING);
    fprintf(stderr, "       -o o  : input format          [%d]\n", ITYPE);
    fprintf(stderr, "                 0 (normalized frequency <0...pi>)\n");
    fprintf(stderr, "                 1 (normalized frequency <0...0.5>)\n");
    fprintf(stderr, "                 2 (frequency (kHz))\n");
    fprintf(stderr, "                 3 (frequency (Hz))\n");
    fprintf(stderr, "       -h    : print this message\n");
    fprintf(stderr, "  infile:\n");
    fprintf(stderr, "       filter input (float)          [stdin]\n");
    fprintf(stderr, "  stdout:\n");
    fprintf(stderr, "       filter output (float)\n");
    fprintf(stderr, "  lspfile:\n");
    fprintf(stderr, "       LSP coefficients (float)\n");
    fprintf(stderr, "\n");
    exit(status);
}


void main(int argc, char **argv)
{
    int		m = ORDER, fprd = FPERIOD, iprd = IPERIOD,
		itype = ITYPE, sampling = SAMPLING, i, j;
    FILE	*fp = stdin, *fpc = NULL;
    double	*c, *inc, *cc, *d, x;
    
    if ((cmnd = strrchr(argv[0], '/')) == NULL)
	cmnd = argv[0];
    else
	cmnd++;
    while (--argc)
	if (**++argv == '-') {
	    switch (*(*argv+1)) {
		case 'm':
		    m = atoi(*++argv);
		    --argc;
		    break;
		case 'p':
		    fprd = atoi(*++argv);
		    --argc;
		    break;
		case 'i':
		    iprd = atoi(*++argv);
		    --argc;
		    break;
		case 's':
		    sampling = atoi(*++argv);
		    --argc;
		    break;
		case 'o':
		    itype = atoi(*++argv);
		    --argc;
		    break;
		case 'h':
		    usage(0);
		default:
		    fprintf(stderr, "%s : Invalid option '%c' !\n", cmnd, *(*argv+1));
		    usage(1);
		}
	}
	else if (fpc == NULL)
	    fpc = getfp(*argv, "r");
	else
	    fp = getfp(*argv, "r");

    if (m % 2 != 0){
	fprintf(stderr,"%s : Order of LSP must be even!\n",cmnd);
	exit(1);
    }
    
    if(fpc == NULL){
	fprintf(stderr,"%s : Cannot open LSP file!\n",cmnd);
	exit(1);
    }

    c = dgetmem(5*m+4);
    cc  = c  + m + 1;
    inc = cc + m + 1;
    d   = inc+ m + 1;
    
    if(freadf(c, sizeof(*c), m+1, fpc) != m+1) exit(1);

    if(itype == 1)
        for(i=1; i<m+1; i++)
            c[i] *= PI2;
    else if (itype == 2 || itype == 3)
        for(i=1; i<m+1; i++)
            c[i] *= PI2 / sampling;
    
    if(itype == 3)
        for(i=1; i<m+1; i++)
            c[i] /= 1000;


    for(;;){
	if(freadf(cc, sizeof(*cc), m+1, fpc) != m+1) exit(0);

	if(itype == 1)
            for(i=1; i<m+1; i++)
                cc[i] *= PI2;
        else if (itype == 2 || itype == 3)
            for(i=1; i<m+1; i++)
                cc[i] *= PI2 / sampling ;
    
        if(itype == 3)
            for(i=1; i<m+1; i++)
                cc[i] /= 1000;
	
	for(i=0; i<=m; i++)
	    inc[i] = (cc[i] - c[i])*iprd / fprd;

	for(j=fprd, i=(iprd+1)/2; j--;){
	    if (freadf(&x, sizeof(x), 1, fp) != 1) exit(0);

	    x *= c[0];

	    x = lspdf_even(x, c, m, d);

	    fwritef(&x, sizeof(x), 1, stdout);
			
	    if (!--i){
		for (i=0; i<=m; i++) c[i] += inc[i];
		i = iprd;
	    }
	}
	movem(cc, c, sizeof(*cc), m+1);
    }
    exit(0);
}

