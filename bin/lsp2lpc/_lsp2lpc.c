/****************************************************************

    $Id: lsp2lpc.c,v 1.1 1996/03/30 07:53:09 koishida Exp koishida $

    Transformation LSP to LPC

	void lsp2lpc(lsp, a, m)

	double  *lsp : LSP
	double  *a   : LPC
	int    	m    : order of LPC

*****************************************************************/

#include <stdio.h>
#include <SPTK.h>

void lsp2lpc(lsp, a, m)
double *lsp, *a;
int m;
{
    int		  i, k, mh;
    double	  xx, xf, cos();
    static double *f = NULL, *p, *q, *a0, *a1, *a2, *b0, *b1, *b2;
    static int	  size;
    
    mh = m / 2;

    if(f == NULL){
	f = dgetmem(5*m+6);
	p  = f  + m;      q  = p  + mh;
	a0 = q  + mh;     a1 = a0 + (mh+1);
	a2 = a1 + (mh+1); b0 = a2 + (mh+1);
	b1 = b0 + (mh+1); b2 = b1 + (mh+1);
	size = m;
    }
    if(m > size){
	free(f);
	f = dgetmem(5*m+6);
	p  = f  + m;      q  = p  + mh;
	a0 = q  + mh;     a1 = a0 + (mh+1);
	a2 = a1 + (mh+1); b0 = a2 + (mh+1);
	b1 = b0 + (mh+1); b2 = b1 + (mh+1);
	size = m;
    }
    
    movem(lsp, f, sizeof(*lsp), m);

    fillz(a0, sizeof(*a0), mh+1); fillz(a1, sizeof(*a1), mh+1);
    fillz(a2, sizeof(*a2), mh+1); fillz(b0, sizeof(*b0), mh+1);
    fillz(b1, sizeof(*b1), mh+1); fillz(b2, sizeof(*b2), mh+1);

    /* lsp filter parameters */
    for (i=0,k=0; i<mh; i++,k+=2){
	p[i] = -2.0 * cos(PI2 * f[k]);
	q[i] = -2.0 * cos(PI2 * f[k+1]);
    }

    /* impulse response of analysis filter */
    xf = 0.0;
    for (k=0; k<=m; k++){
	xx = 0.0;
	if (k == 0)
	    xx = 1.0;
	a0[0] = xx + xf;
	b0[0] = xx - xf;
	xf = xx;
	
	for (i=0; i<mh; i++){
	    a0[i+1] = a0[i] + p[i] * a1[i] + a2[i];
	    b0[i+1] = b0[i] + q[i] * b1[i] + b2[i];
	    a2[i] = a1[i]; a1[i] = a0[i];
	    b2[i] = b1[i]; b1[i] = b0[i];
	}
	
	if (k != 0)
	    a[k-1] = -0.5 * (a0[mh] + b0[mh]);
    }

    for (i=m-1; i>=0; i--)
	a[i+1] = -a[i];
    a[0] = 1.0;
}
