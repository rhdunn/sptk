/************************************************************************
*									*
*    LMA Digital Filter for Speech Synthesis				*
*									*
*					1986.6  K.Tokuda		*
*					1996.1  K.Koishida		*
*									*
*	usage:								*
*		lmadf [ options ] [ infile ] > stdout			*
*	options:							*
*		-m m     :  order of cepstrum		[25]		*
*		-p p     :  frame period		[100]		*
*		-i i     :  interpolation period	[1]		*
*		-P P	 :  order of pade approximation [4]		*
*	infile:								*
*		cepstral coefficients					*
*		    , c~(0), c~(1), ..., c~(M),				*
*		excitation sequence					*
*		    , x(0), x(1), ..., 					*
*	stdout:								*
*		filtered sequence					*
*		    , y(0), y(1), ...,					*
*	note:								*
*		P = 4 or 5						*
*	require:							*
*		lmadf()							*
*									*
************************************************************************/

static char *rcs_id = "$Id: lmadf-main.c,v 1.1 1996/04/01 09:24:25 koishida Exp koishida $";


/*  Standard C Libraries  */
#include <stdio.h>
#include <string.h>
#include <SPTK.h>


/*  Required Functions  */
double	lmadf(), exp();


/*  Default Values  */
#define ORDER		25
#define	FPERIOD		100
#define	IPERIOD		1
#define PADEORD 	4


/*  Command Name  */
char	*cmnd;


void usage(int status)
{
    fprintf(stderr, "\n");
    fprintf(stderr, " %s - LMA digital filter for speech synthesis\n",cmnd);
    fprintf(stderr, "\n");
    fprintf(stderr, "  usage:\n");
    fprintf(stderr, "       %s [ options ] cfile [ infile ] > stdout\n", cmnd);
    fprintf(stderr, "  options:\n");
    fprintf(stderr, "       -m m  : order of cepstrum           [%d]\n", ORDER);
    fprintf(stderr, "       -p p  : frame period                [%d]\n", FPERIOD);
    fprintf(stderr, "       -i i  : interpolation period        [%d]\n", IPERIOD);
    fprintf(stderr, "       -P P  : order of pade approximation [%d]\n", PADEORD);
    fprintf(stderr, "       -h    : print this message\n");
    fprintf(stderr, "  infile:\n");
    fprintf(stderr, "       filter input (float)                [stdin]\n");
    fprintf(stderr, "  stdout:\n");
    fprintf(stderr, "       filter output (float)\n");
    fprintf(stderr, "  cfile:\n");
    fprintf(stderr, "       cepstrum (float)\n");
    fprintf(stderr, "  note:\n");
    fprintf(stderr, "       P = 4 or 5\n");
    fprintf(stderr, "\n");
    exit(status);
}


void main(int argc, char **argv)
{
    int		m = ORDER, fprd = FPERIOD, iprd = IPERIOD, i, j, pd = PADEORD;
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
		case 'P':
		    pd = atoi(*++argv);
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

    if((pd < 4)||(pd > 5)){
	fprintf(stderr,"%s : Order of pade approximation is 4 or 5!\n",cmnd);
	exit(1);
    }

    if(fpc == NULL){
	fprintf(stderr,"%s : Cannot open cepstrum file!\n",cmnd);
	exit(1);
    }
	
    c = dgetmem(m+m+m+3+(m+1)*pd*2);
    cc  = c  + m + 1;
    inc = cc + m + 1;
    d   = inc+ m + 1;
    
    if(freadf(c, sizeof(*c), m+1, fpc) != m+1) exit(1);

    for(;;){
	if(freadf(cc, sizeof(*cc), m+1, fpc) != m+1) exit(0);
	
	for(i=0; i<=m; i++)
	    inc[i] = (cc[i] - c[i])*iprd / fprd;

	for(j=fprd, i=(iprd+1)/2; j--;){
	    if (freadf(&x, sizeof(x), 1, fp) != 1) exit(0);

	    x *= exp(c[0]);
	    x = lmadf(x, c, m, pd, d);
	    
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

