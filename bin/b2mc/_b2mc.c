/***************************************************************
    $Id: b2mc.c,v 1.1 1996/02/28 04:48:33 koishida Exp koishida $

    Transform MLSA Digital Filter Coefficients to Mel Cepstrum 

	void	b2mc(b, mc, m, a)

	double	*b  : MLSA digital filter coefficients
	double	*mc : mel cepstral coefficients
	int	m   : order of mel cepstrum
	double	a   : all-pass constant

***************************************************************/

void b2mc(b, mc, m, a)
double *b, *mc, a;
int m;
{
    double d, o;
    
    d = mc[m] = b[m];
    for(m--; m>=0; m--){
	o = b[m] + a * d;
	d = b[m];
	mc[m] = o;
    }
}
