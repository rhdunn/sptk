/****************************************************************

    $Id: lpc2lsp.c,v 1.1 1996/03/30 07:39:37 koishida Exp koishida $

    Transformation LPC to LSP

	void	lpc2lsp(a, lsp, m, n, maxp, maxq, eps)

	double *a    : LPC
	double *lsp  : LSP
	int    m     : order of LPC
	int    n     : split number of unit circle
	int    maxp  : maximum number of interpolation for P(z)
	int    maxq  : maximum number of interpolation for Q(z)
	double eps   : end condition of interpolation

*****************************************************************/

#include <stdio.h>
#include <SPTK.h>

void lpc2lsp(a, lsp, m, n, maxp, maxq, eps)
double *a, *lsp, eps;
int m, n, maxp, maxq;
{
    int 	  i, j, mp, mh, nf, mb, jc;
    double	  fr, pxr, tpxr, tfr, pxm, pxl, fl, qxl, tqxr, 
    		  ang, fm, qxm, qxr, tqxl, cos(), sin(), fabs();
    static double *p = NULL, *q, *f;
    static int    size;
    
    if(p == NULL){
	p = dgetmem(m+m+m+1);
	q = p + m; f = q + m;
	size = m;
    }
    if(m > size){
	free(p);
	p = dgetmem(m+m);
	q = p + m; f = q + m;
	size = m;
    }
	
    mp = m + 1;
    mh = m / 2;

    if(n  < -1) n  = 128;
    if(maxp < -1) maxp = 4;
    if(maxq < -1) maxq = 15;

    if(eps < -1.0) eps = 10e-6;
    
    /* generate p and q polynomials */
    for (i=0; i<mh; i++){
	p[i] = a[i+1] + a[m-i];
	q[i] = a[i+1] - a[m-i];
    }
    
    /* compute p at f=0.0 */
    fl = 0.0;
    for (pxl=1.0,j=0; j<mh; j++) pxl += p[j];

    /* search for zeros of p */
    nf = 0;
    for (i=1; i<=n; pxl=tpxr, fl=tfr, i++){
	mb = 0;
	fr = i * (0.5 / n);
	pxr = cos(mp * PI * fr);
	for (j=0; j<mh; j++){
	    jc = mp - (j+1)*2;
	    ang = jc * PI * fr;
	    pxr += cos(ang) * p[j];
	}
	tpxr = pxr;
	tfr = fr;
	
	if (pxl*pxr > 0.0) continue;

    	do{
	    mb++;
	    fm = fl + (fr-fl) / (pxl-pxr) * pxl;
	    pxm = cos(mp * PI * fm);
    
	    for (j=0; j<mh; j++){
		jc = mp - (j+1) * 2;
		ang = jc * PI * fm;
		pxm += cos(ang) * p[j];
	    }
	    
	    (pxm*pxl > 0.0) ? (pxl=pxm, fl=fm) : (pxr=pxm, fr=fm);

	} while ((fabs(pxm) > eps) && (mb < maxp));

	f[nf] = fl + (fr-fl) / (pxl-pxr) * pxl;
	nf += 2;
	if (nf > m-2) break;
    }

    /* search for the zeros of q(z) */
    f[m] = 0.5;
    fl = f[0];
    qxl = sin(mp * PI * fl);
    for (j=0; j<mh; j++){
	jc = mp - (j+1) * 2;
	ang = jc * PI * fl;
	qxl += sin(ang) * q[j];
    }

    for (i=2; i<mp; qxl=tqxr, fl=tfr, i+=2){
	mb = 0;
	fr = f[i];
	qxr = sin(mp * PI * fr);
	for (j=0; j<mh; j++){
	    jc = mp - (j+1) * 2;
	    ang = jc * PI * fr;
	    qxr += sin(ang) * q[j];
	}

	tqxl = qxl;
	tfr = fr;
	tqxr = qxr;
    
	do{
	    mb++;
	    fm = (fl+fr) * 0.5;
	    qxm = sin(mp * PI * fm);

	    for (j=0; j<mh; j++){
		jc = mp - (j+1) * 2;
		ang = jc * PI * fm;
		qxm += sin(ang) * q[j];
	    }
	    
	    (qxm*qxl > 0.0) ? (qxl=qxm, fl=fm) : (qxr=qxm, fr=fm);
	    
	} while ((fabs(qxm) > eps) && (mb < maxq));

	f[i-1] = fl + (fr-fl) / (qxl-qxr) * qxl;
    }
    movem(f, lsp, sizeof(*f), m);
}
