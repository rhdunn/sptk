/************************************************************************
*									*
*    Real Sequence to Phase   (filter command)				*
*									*
*					1989.6  K.Tokuda		*
*					1996.2  N.Isshiki Modified	*
*									*
*	usage:								*
*		phase [ options ] [ infile ] > stdout			*
*      	options:							*
*		-l l : FFT size	 			[256]		*
*		-p p : numerator cofficients file			*
*		-z z : denominator cofficients file			*
*		-m m : order of denominator polynomial	[L-1]		*
*		-n n : order of numerator polynomial    [L-1]		*
*	infile:								*
*		input sequense (float)					*
*	stdout:								*
*		P0, P1, ..., PL/2 (float)(pi radian)			*
*									*
************************************************************************/

static char *rcs_id = "$Id:";

/* Standard C Libraries */
#include <stdio.h>
#include <SPTK.h>
#include <string.h>


/* Required Functions */
void phase();


/* Default Values */
#define	FLMLNG	256


/* Command Name */
char 	*cmnd;


void usage(status)
int status;
{
	fprintf(stderr, "\n");
	fprintf(stderr, " %s - transform real seauence to phase\n", cmnd);
	fprintf(stderr, "\n");
	fprintf(stderr, "  usage:\n");
	fprintf(stderr, "       %s [ options ] [ infile ] > stdout\n", cmnd);
	fprintf(stderr, "  options:\n");
	fprintf(stderr, "       -l l  : flame length                    [%d]\n",FLMLNG);
	fprintf(stderr, "       -p p  : numerator cofficients file\n");
	fprintf(stderr, "       -z z  : denominator cofficients file\n");
	fprintf(stderr, "       -m m  : order of denominator polynomial [L-1]\n");
	fprintf(stderr, "       -n n  : order of numerator polynomial   [L-1]\n");
	fprintf(stderr, "       -h    : print this message\n");
	fprintf(stderr, "  infile:\n");
	fprintf(stderr, "       data sequence (float)                   [stdin]\n");
	fprintf(stderr, "  stdout:\n");
	fprintf(stderr, "       phase (float) [pi rad]\n");
	fprintf(stderr, "\n");
	exit(status);
}

main(argc, argv)
	int argc;
	char *argv[];
{

	int flng = FLMLNG, m = -1, n = -1;
	double *p, *z, *ph;
	char *file_z = "", *file_p = "";
	FILE *fp_z, *fp_p;
	int  no, i;

        if ((cmnd = strrchr(argv[0], '/')) == NULL)
	    cmnd = argv[0];
        else
	    cmnd++;

	while (--argc)
		if (**++argv == '-'){
			argc--;
			switch (*(*argv+1)) {
			case 'm':
				m = atoi(*++argv);
				break;
			case 'n':
				n = atoi(*++argv);
				break;
			case 'l':
				flng = atoi(*++argv);
				break;
			case 'z':
				file_z = *++argv;
				break;
			case 'p':
				file_p = *++argv;
				break;
			case 'h':
				usage(0);
			default:
				fprintf(stderr,
				"%s : not available option !\n", cmnd);
				usage(1);
			}
		} else {
			fprintf(stderr,
			"%s : input file is not available !\n", cmnd);
			usage(1);
		}
	if (m < 0) m = flng - 1;
	if (n < 0) n = flng - 1;
	no = flng / 2 + 1;

	p = dgetmem(2 * flng + n);
	z = &p[flng];
	ph = &z[flng];

	if (*file_z != '\0')
		fp_z = getfp(file_z, "r");

	if (*file_p != '\0')
		fp_p = getfp(file_p, "r");

	if (*file_z == '\0' && *file_p == '\0')
		fp_z = stdin;	
	for(;;){
		if (*file_z == '\0' && *file_p != '\0'){
			z[0] = 1;
			n = 0;
		}else{
			if (freadf(z, sizeof(*z), n+1, fp_z) == 0 )
				exit(0);
		}

		if (*file_p == '\0' ){
			p[0] = 1;
			m = 0;
		}else{
			if(freadf(p, sizeof(*p), m+1, fp_p) == 0 )
				exit(0);
		}
		phase(p, m, z, n, ph, flng);	
		fwritef(ph, sizeof(*ph), no, stdout);
	}
	exit(0);
}
