/***************************************************************
    $Id:$

    Generalized Exponential Function (real argument)

	double gexp(r, x)

	double	r : gamma
	double  x : argument

*****************************************************************/

double gexp(r, x)
double r, x;
{
    double exp(), pow();

    if (r == 0.0)
	return(exp(x));
    else
	if (r < 0.)
	    return(pow(1/(1 + r*x), -1/r));
	else
	    return(pow(1 + r*x, 1/r));
}
