/****************************************************************

    $Id: lspdf.c,v 1.1 1996/09/17 05:07:18 koishida Exp koishida $

    LSP Speech Synthesis Digital Filter

	double	lspdf_even(x, a, m, d)

	double	x   : input
	double	*f  : LSP coefficients
	int	m   : order of coefficients
	double  *d  : delay

	return value : filtered data

*****************************************************************/
#include <stdio.h>

double	lspdf_even(x, f, m, d)
double	x, *f, *d;
int	m;
{
    double  	  *d1, *d2, *lsp, x1, x2, cos();
    register int  i;

    d1 = d + 1;
    d2 = d1 + m;
    lsp = f + 1;
    x1 = x2 = d[0];

    for(i=0; i<m; i+=2){
	d1[i] -= 2.0 * x1 * cos(lsp[i]);
	d2[i] -= 2.0 * x2 * cos(lsp[i+1]);
	d1[i+1] += x1;
	d2[i+1] += x2;
	x += d1[i] + d2[i];
	x1 = d1[i+1];
	x2 = d2[i+1];
    }
    
    x -= d2[m-1] - d1[m-1];

    for(i=m-1; i>0; i--){
	d1[i] = d1[i-1];
	d2[i] = d2[i-1];
    }
    d1[0] = d2[0] = d[0];
    d[0] = -0.5 * x;

    return(x);
}

