/***************************************************************
    $Id: gc2gc.c,v 1.1 1996/04/03 07:41:30 koishida Exp koishida $

    Generalized Cepstral Transformation	

	void gc2gc(c1, m1, g1, c2, m2, g2)

	double	*c1	: normalized generalized cepstrum (input)
	int	m1	: order of gc1
	double	g1	: gamma of gc1
	double	*c2	: normalized generalized cepstrum (output)
	int	m2	: order of gc2
	double	g2	: gamma of gc2

*****************************************************************/

void gc2gc(c1, m1, g1, c2, m2, g2)
double 	*c1, *c2, g1, g2;
int 	m1, m2;
{
    double 	ss1, ss2, cc;
    int 	i, min, k, mk;
    
    c2[0] = c1[0];
    for (i = 1; i <= m2; i++){
	ss1 = ss2 = 0.0;
	min = m1 < i ? m1 : i - 1;
	for (k = 1; k <= min; k++){
	    mk = i - k;
	    cc = c1[k] * c2[mk];
	    ss2 += k * cc;
	    ss1 += mk * cc;
	}

	if (i <= m1)
	    c2[i] = c1[i] + (g2*ss2 - g1*ss1)/i;
	else
	    c2[i] = (g2*ss2 - g1*ss1)/i;
    }
}
