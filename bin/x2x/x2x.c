/************************************************************************
*									*
*    Data Type Transformation						*
*									*
*					1985.12 K.Tokuda		*
*					1996.5  K.Koishida		*
*									*
*	usage:								*
*		x2x [options] [infile] > stdout				*
*	options:							*
*		-r       :  specify rounding off when a real number 	*
*			    is substituted for a integer	[FALSE]	*
*		+type1   :  input data type 			[f]	*
*		+type2   :  output data type 			[type1]	*
*				c (char)     s (short)			*
*				i (int)      l (long)			*
*				f (float)    d (double)			*
*				a (ascii)				*
*		+a a     :  column number 			[1]	* 
*		%format  :  specify output format similar to 	[FALSE]	*
*                           "printf()" of C function, 			*
*                           if type2 is ascii formant.			* 
*									*
************************************************************************/

static char *rcs_id = "$Id:$";


/*  Standard C Libraries  */
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <SPTK.h>


typedef enum _Boolean {FA, TR} Boolean;
char *BOOL[] = {"FALSE", "TRUE"};


/*  Default Values  */
#define ROUND		FA
#define COL		1


/*  Command Name  */
char	*cmnd;


void usage(int status)
{
    fprintf(stderr, "\n");
    fprintf(stderr, " %s - data type transformation\n",cmnd);
    fprintf(stderr, "\n");
    fprintf(stderr, "  usage:\n", cmnd);
    fprintf(stderr, "       %s [ options ] [ infile ] > stdout\n", cmnd);
    fprintf(stderr, "  options:\n");
    fprintf(stderr, "       +type1  : input data type                             [f]\n");
    fprintf(stderr, "       +type2  : output data type                            [type1]\n");
    fprintf(stderr, "                 c (char)      s (short)\n");
    fprintf(stderr, "                 i (int)       l (long)\n");
    fprintf(stderr, "                 f (float)     d (double)\n");
    fprintf(stderr, "                 a (ascii)\n");
    fprintf(stderr, "       +a a    : column number                               [%d]\n",COL);
    fprintf(stderr, "       -r      : specify rounding off when a real number\n");
    fprintf(stderr, "                 is substituted for a integer                [%s]\n",BOOL[ROUND]);
    fprintf(stderr, "       %%format : specify output format similar to 'printf()', \n");
    fprintf(stderr, "                 if type2 is ascii.                          [N/A]\n");
    fprintf(stderr, "       -h      : print this message\n");
    fprintf(stderr, "  infile:\n");
    fprintf(stderr, "       data sequence                                    [stdin]\n");
    fprintf(stderr, "  stdout:\n");
    fprintf(stderr, "       transformed data sequence\n");
    fprintf(stderr, "\n");

    exit(status);
}

double r = 0.0;

void main(int argc, char **argv)
{
    char        c1, c2, *form = "%g";
    double      x;
    int         size1 = 0, size2 = 0, i = 1, col = COL, atoi();
    FILE	*fp = stdin;
    Boolean     round = ROUND;
    void        x2x();
    
    if ((cmnd = strrchr(argv[0], '/')) == NULL)
	cmnd = argv[0];
    else
	cmnd++;
    while (--argc)
	if (**++argv == '+') {
	    (*argv)++;
	    while(**argv != '\0'){
		    switch (**argv) {
			case 's':
			    if(size1 == 0){
				c1 = 's';
				size1 = sizeof(short);
			    } else {
				c2 = 's';
				size2 = sizeof(short);
			    }
			    break;
			case 'i':
			    if(size1 == 0){
				c1 = 'i';
				size1 = sizeof(int);
			    } else {
				c2 = 'i';
				size2 = sizeof(int);
			    }
				break;
			case 'l':
			    if(size1 == 0){
				c1 = 'l';
				size1 = sizeof(long);
			    } else {
				c2 = 'l';
				size2 = sizeof(long);
			    }
			    break;
			case 'f':
			    if(size1 == 0){
				c1 = 'f';
				size1 = sizeof(float);
			    } else {
				c2 = 'f';
				size2 = sizeof(float);
			    }
			    break;
			case 'd':
			    if(size1 == 0){
				c1 = 'd';
				size1 = sizeof(double);
			    } else {
				c2 = 'd';
				size2 = sizeof(double);
			    }
			    break;
			case 'c':
			    if(size1 == 0){
				c1 = 'c';
				size1 = sizeof(char);
			    } else {
				c2 = 'c';
				size2 = sizeof(char);
			    }
			    break;
			case 'a':
			    if(size1 == 0){
				c1 = 'a';
				size1 = -1;
			    } else {
				c2 = 'a';
				size2 = -1;
				if(*(argv+1) != '\0' && isdigit(**(argv+1))){
				    col = atoi(*++argv);
				    argc--;
				}

			    }
			    break;
			default:
			    fprintf(stderr, "%s : Invalid option '%c' !\n", cmnd, *(*argv+1));
			    usage(1);
			}
		(*argv)++;
		}
	}
	else if(**argv == '-'){
	    switch (*(*argv+1)) {
                case 'r':
		    round = 1 - round;
		    break;
		case 'h':
		    usage(0);
		default:
		    fprintf(stderr, "%s : Invalid option '%c' !\n", cmnd, *(*argv+1));
		    usage(1);
		}
	}
	else if(**argv == '%')
	    form = *argv++;
    	else 
	    fp = getfp(*argv, "r");

    if(round)
	r = 0.5;

    if(size1 == 0){
	size1 = sizeof(float);
	c1 = 'f';
    }
    if(size2 == 0){
	size2 = size1;
	c2 = c1;
    }

    if (c1 == 'a'){
	if (c2 == 'a')
	    while (fscanf(fp, "%le", &x) != EOF){
		printf(form, x);
		if (i == col){
		    i = 1;
		    printf("\n");
		}else{
		    i++;
		    printf("\t");
		}
	    }
	else
	    while (fscanf(fp, "%le", &x) != EOF){
		x2x(&x, &x, 'd', c2);
		fwrite(&x, size2, 1, stdout);
	    }
    }
    else{
	if (c2 == 'a'){
	    while (fread(&x, size1, 1, fp) == 1){
		x2x(&x, &x, c1, 'd');
		printf(form, x);
		if (i == col){
		    i = 1;
		    printf("\n");
		}else{
		    i++;
		    printf("\t");
		}
	    }
	}
	else
	    while (fread(&x, size1, 1, fp) == 1){
		x2x(&x, &x, c1, c2);
		fwrite(&x, size2, 1, stdout);
	    }
    }
    
    exit(0);
    
}

void x2x(x1, x2, c1, c2)
char *x1, *x2;
char c1, c2;
{
    double x;

    switch (c1){
        case 's' :
	    x = *(short *)x1;
	    break;
	case 'i' :
	    x = *(int *)x1;
	    break;
	case 'l' :
	    x = *(long *)x1;
	    break;
	case 'f' :
	    x = *(float *)x1;
	    break;
	case 'd' :
	    x = *(double *)x1;
	    break;
	case 'c' :
	    x = *(char *)x1;
	    break;
    }

    switch (c2){
	case 's' :
	    *(short *)x2 = x + r;
	    break;
	case 'i' :
	    *(int *)x2 = x + r;
	    break;
	case 'l' :
	    *(long *)x2 = x + r;
	    break;
	case 'f' :
	    *(float *)x2 = x;
	    break;
	case 'd' :
	    *(double *)x2 = x;
	    break;
	case 'c' :
	    *(char *)x2 = x + r;
	    break;
    }
}

