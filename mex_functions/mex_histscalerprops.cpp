#include <mex.h>
#include <functional>
#include <algorithm>
#include <cstdint>
#include <vector>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	/*
	input:
		var1 [array] = data to be binned.
		var2 [array] = scaling data.
		var3 [array] = [minBin binSize maxBin].

	output:
		var1 [array] = scaler mean;
		var2 [array] = scaler root mean square;
		var3 [array] = scaler variances;
		var4 [array] = bins;
	*/

	bool reterror = false;
	int errnum = 0;

	size_t M(0), N(0);

	if (nrhs >= 2) {
		M = std::max(mxGetDimensions(prhs[0])[0], mxGetDimensions(prhs[0])[1]);
		N = std::max(mxGetDimensions(prhs[1])[0], mxGetDimensions(prhs[1])[1]);

		if (M == N) {
			double* data1 = (double*)mxGetData(prhs[0]);
			double* data2 = (double*)mxGetData(prhs[1]);

			double minBin(0), stepBin(0), maxBin(0), numOfBins(100);
			const size_t B = std::max(mxGetDimensions(prhs[2])[0], mxGetDimensions(prhs[2])[1]);
			if (B == 3) {
				double* binprops = (double*)mxGetData(prhs[2]);
				minBin	= binprops[0];
				stepBin	= binprops[1];
				maxBin	= binprops[2];
				numOfBins = (maxBin - minBin) / stepBin;
			}
			else {
				minBin	= *std::min(data1, data1 + M);
				maxBin	= *std::max(data1, data1 + M);
				stepBin	= (maxBin - minBin) / numOfBins;
			}
			//numOfBins += 2;

			//mexPrintf("[%f, %f (%f), %f]\n", minBin, stepBin, numOfBins, maxBin);
			
			std::vector<double> bins(numOfBins);
			std::vector<double> hist(numOfBins);
			std::vector<double> means(numOfBins);
			std::vector<double> rms(numOfBins);
			std::vector<double> vars(numOfBins);

			for (size_t i = 0; i < M; i++) {
				size_t index =(size_t) ((data1[i] - minBin) / stepBin);

				if (index < 0) {
					index = 0;
				}

				if (index >= numOfBins) {
					index = numOfBins - 1;
				}

				means[index] += data2[i];
				rms[index] += data2[i] * data2[i];
				hist[index]++;
			}

			for (size_t i = 0; i < means.size(); i++) {
				bins[i]		 = i*stepBin;
				if (hist[i] > 0) {
					means[i] /= hist[i];
					rms[i] /= hist[i];
					rms[i] = sqrt(rms[i]);
				}
			}

			for (size_t i = 0; i < M; i++) {
				size_t index = (size_t)((data1[i] - minBin) / stepBin);

				if (index < 0) {
				index = 0;
				}

				if (index >= numOfBins) {
				index = numOfBins - 1;
				}

				vars[index] += std::pow(data2[i] - means[index], 2.);
			}

			for (size_t i = 0; i < means.size(); i++) {
				if (hist[i] > 0) {
					vars[i] /= hist[i];
				}
			}

			plhs[0] = mxCreateDoubleMatrix(1, means.size(),	mxREAL);
			plhs[1] = mxCreateDoubleMatrix(1, rms.size(),	mxREAL);
			plhs[2] = mxCreateDoubleMatrix(1, vars.size(),	mxREAL);
			plhs[3] = mxCreateDoubleMatrix(1, bins.size(),	mxREAL);

			double* omeans	= (double*)mxGetData(plhs[0]);
			double* orms	= (double*)mxGetData(plhs[1]);
			double* ovars	= (double*)mxGetData(plhs[2]);
			double* obins	= (double*)mxGetData(plhs[3]);

			memcpy(omeans,	&means[0],	sizeof(double)*means.size());
			memcpy(orms,	&rms[0],	sizeof(double)*rms.size());
			memcpy(ovars,	&vars[0],	sizeof(double)*vars.size());
			memcpy(obins,	&bins[0],	sizeof(double)*bins.size());
		}
		else {
			reterror = true;
			errnum = 1;
		}
	}
	else {
		reterror = true;
		errnum = 2;
	}

	if (reterror) {
		plhs[0] = mxCreateDoubleScalar(errnum);
	}
}