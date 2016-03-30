#define _USE_MATH_DEFINES
#include <math.h>
#include <mex.h>
#include <random>

double rand_exp(const double uii1, const double lambda) {
	return log(1 - uii1) / (-lambda);
}

double rand_levy(const double uii1, const double uii2, const double alpha, const double c) {
	double u, v, t, s;

	u = M_PI * (uii1 - 0.5);
	v = rand_exp(uii2, 1);

	if (alpha == 1)               /* cauchy case */
	{
		t = tan(u);
		return c * t;
	}

	if (alpha == 2)             /* gaussian case */
	{
		t = 2 * sin(u) * sqrt(v);
		return c * t;
	}

	/* general case */

	t = sin(alpha * u) / pow(cos(u), 1 / alpha);
	s = pow(cos((1 - alpha) * u) / v, (1 - alpha) / alpha);

	return c * t * s;
}

double rand_levy_skew(const double uii1, const double uii2, const double alpha, const double beta, const double c) {
	double V, W, X;

	V = M_PI * (uii1 - 0.5);
	W = rand_exp(uii2, 1);

	if (beta == 0)  /* symmetric case */
	{
		return rand_levy(uii1, uii2, c, alpha);
	}

	if (alpha == 1)
	{
		X = ( (M_PI_2 + beta * V) * tan(V) - beta * log( M_PI_2 * W * cos(V) / (M_PI_2 + beta * V) ) ) / M_PI_2;
		return c * (X + beta * log(c) / M_PI_2);
	}
	else
	{
		double t = beta * tan(M_PI_2 * alpha);
		double B = atan(t) / alpha;
		double S = pow(1 + t * t, 1 / (2 * alpha));

		X = S * sin(alpha * (V + B)) / pow(cos(V), 1 / alpha) * pow(cos(V - alpha * (V + B)) / W, (1 - alpha) / alpha);
		return c * X;
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	const size_t numOfRands = mxGetScalar(prhs[0]);
	const double beta = mxGetScalar(prhs[1]);
	const double alpha = mxGetScalar(prhs[2]);
	const double c = mxGetScalar(prhs[3]);
	
	
	plhs[0] = mxCreateDoubleMatrix(1, numOfRands, mxREAL);
	double *rands = (double*)mxGetData(plhs[0]);


	std::random_device rd;

	std::uniform_real_distribution<double> uii(0, 1);
	std::mt19937 mt(rd());
	
	for (size_t i = 0; i < numOfRands; i++) {
		double uii1 = 1 - uii(mt);
		double uii2 = 1 - uii(mt);
		rands[i] = rand_levy_skew(uii1, uii2, alpha, beta, c);
	}
}