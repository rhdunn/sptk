/***********************************************************
  $Id: acorr.c,v 1.2 1996/02/02 08:26:08 isshiki Exp $
	obtain autocorrelation sequence
		acorr(x,l,r,np)

		double	*x 	: input sequence
		int	l	: flame length
		double	*r	: autocorrelation sequence
		int 	np	: order of sequence

			Naohiro Isshiki,	Dec.1995
***********************************************************/

void acorr(x,l,r,np)
double 	*x,*r;
int	l,np;
{
	double	d;
	register int	k, i;

	for(k = 0; k <= np; ++k) {
		for(d = i = 0; i < l - k; ++i)
			d += x[i] * x[i + k];
		r[k] = d;
	}
}
