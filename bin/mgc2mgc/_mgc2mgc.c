/***********************************************************************

    $Id: mgc2mgc.c,v 1.1 1996/04/08 07:17:36 koishida Exp koishida $

    Frequency and Generalized Cepstral Transformation

	void mgc2mgc(c1, m1, a1, g1, c2, m2, a2, g2)

	double	*c1   : minimum phase sequence (input)
	int	m1    : order of c1
	double	a1    : alpha of c1
	double	g1    : gamma of c1
	double	*c2   : transformed sequence (output)
	int	m2    : order of c2
	double	a2    : alpha of c2
	double	g2    : gamma of c2

***********************************************************************/

#include <stdio.h>
#include <SPTK.h>

void mgc2mgc(c1, m1, a1, g1, c2, m2, a2, g2)
double *c1, *c2, a1, a2, g1, g2;
int m1, m2;
{
    double 	  a;
    static double *ca = NULL;
    static int    size;

    if(ca == NULL){
	ca = dgetmem(m2+1);
	size = m2;
    }
    if(m2 > size){
	free(ca);
	ca = dgetmem(m2+1);
	size = m2;
    }

    a = (a2 - a1) / (1 - a1*a2);

    if(a == 0){
	movem(c1, ca, sizeof(*c1), m1+1);
	fillz(ca+m1+1, sizeof(*ca), m2-m1);
    }
    else
	freqt(c1, m1, ca, m2, a);
    
    gnorm(ca, ca, m2, g1);
    gc2gc(ca, m2, g1, c2, m2, g2);
    ignorm(c2, c2, m2, g2);
}

