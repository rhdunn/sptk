/***************************************************************
    $Id: c2sp.c,v 1.1 1996/04/08 06:42:02 koishida Exp koishida $

    Transform cepstrum to spectrum

	void	c2sp(c, m, sp, y, l)

	double	*c  : cepstrum
	int	m   : order of cepstrum
	double	*x  : real part of spectrum
	double	*y  : imaginary part of spectrum
	double	l   : FFT length

***************************************************************/

#include <stdio.h>
#include <SPTK.h>

void c2sp(c, m, x, y, l)
double	*c, *x, *y;
int	m, l;
{
    int            m1;

    m1 = m + 1;
    
    movem(c, x, sizeof(*c), m1);
    fillz(x+m1, sizeof(*x), l-m1);

    fftr(x, y, l);
}
	
