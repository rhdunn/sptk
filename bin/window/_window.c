/****************************************************************
$Id: window.c,v 1.1 1996/01/25 12:47:34 nmiya Exp nmiya $

			Window function
			---------------

	double  window( name, x, size, pnflg );

		char	*name;	window name

			blackman, hamming,
			hanning,  bartlett, trapezoid

		real	*x;	1 frame data
		int	size;	window(frame) size
		int	nflg;	normalizing flag

			nflg = 0 : don't normalize
			       1 : normalize by power
			       2 : normalize by magnitude

    set windowed value to "*x" and return "normalizing gain".
*****************************************************************/
#include	<stdio.h>
/*#include	<SPTK.h>*/

#define	M_2PI   (2 * 3.14159265358979323846) 

double	cos();
double	sqrt();
double	*dgetmem();

double*	blackman();
double*	hamming();
double*	hanning();
double*	bartlett();
double*	trapezoid();

char	*windows[] = {
			"blackman",
			"hamming",
			"hanning",
			"bartlett",
			"trapezoid",
	};

double window( name, x, size, nflg )
	char	*name;
	double	*x;
	int	size, nflg;
{
	register int	i, n, k;
	double	g;
	static double	*w=NULL;
	static int	size_w;

	for (i=0, n=-1; i<5; i++)
		if (!strcmp(name, windows[i])) n=i;

	if ( n==-1 ) {
		fprintf(stderr, "window : unknown window name\n");
		exit(1);
	}


	if (w == NULL) {
		w = dgetmem(size);
		size_w = size;
	}
	if ( size > size_w) {
		free(w);
		w = dgetmem(size);
		size_w = size;
	}
	switch ( n )  {
	    case 1:	hamming(w, size);		break;
	    case 2:	hanning(w, size);		break;
	    case 3:	bartlett(w, size);		break;
	    case 4:	trapezoid(w, size);		break;
	    case 0:
	    default:	blackman(w, size);		break;
	}

	switch ( nflg ) {
	    case 1:
		g = 0.0;
		for (n=0; n<size; n++)
			g += w[n] * w[n];
		g = sqrt(g);
		break;
	    case 2:
		g = 0.0;
		for (n=0; n<size; n++)
			g += w[n];
		break;
	    case 0:
	    default:
		g = 1.0;
	}

	for (n=0; n<size; n++)
		x[n] = x[n] * w[n] / g;

	return ( g );
}
/************************************************
		Bartlett window

	double  *bartlett(w, leng)

	double	*w;	window values
	int	leng;	window length
************************************************/

double  *bartlett(w, leng)
	double	*w;
	int	leng;
{
	register int	k, l;
	register double	*p;

	l = leng - 1;
	p = w;

	for (k=0; k<l/2; k++)
		*p++ = (2.0 * k) / l;
	for ( ; k<=l; k++)
		*p++ = 2.0 - (2.0 * k) / l;

	return ( w );
}

/************************************************
		Blackman window

	double  *blackman(w, leng)

	double	*w;	window values
	int	leng;	window length
************************************************/

double  *blackman(w, leng)
	double	*w;
	int	leng;
{
	register int	i;
	double		arg, x;
	register double	*p;

	arg = M_2PI / (leng - 1);
	for (p=w, i=0; i<leng; i++)  {
		x = arg * i;
		*p++ = 0.42 - 0.50 * cos(x) + 0.08 * cos(x+x);
	}
	return( w );
}


/************************************************
		Hamming window

	double  *hamming(w, leng)
	double	*w;	window values
	int	leng;	window length
************************************************/

double  *hamming(w, leng)
	double	*w;
	int	leng;
{
	register int	i;
	double		arg;
	register double	*p;

	arg = M_2PI / (leng - 1);
	for (p=w, i=0; i<leng; i++)
		*p++ = 0.54 - 0.46 * cos(i * arg);

	return ( w );
}


/************************************************
		Hanning window

	double  *hanning(w, leng)
	double	*w;	window values
	int	leng;	window length
************************************************/

double  *hanning(w, leng)
	double	*w;
	int	leng;
{
	register int	i;
	double		arg;
	register double	*p;

	arg = M_2PI / (leng - 1);
	for (p=w, i=0; i<leng; i++)
		*p++ = 0.5 * (1 - cos(i * arg));

	return ( w );
}

/************************************************
		trapezoid window

	double  *trapezoid(w, leng)
	double	*w;	window values
	int	leng;	window length
************************************************/

double  *trapezoid(w, leng)
	double	*w;
	int	leng;
{
	register int	k;
	int		l, m1, m2;
	register double	*p;

	p = w;
	l = leng - 1;
	m1 = l / 4;
	m2 = 3 * m1;
	for (k=0; k<m1; k++)
		*p++ = (4.0 * k) / l;
	for (k=m1; k<m2; k++)
		*p++ = 1.0;
	for (k=m2; k<=l; k++)
		*p++ = 4.0 - (4.0 * k) / l;

	return ( w );
}
