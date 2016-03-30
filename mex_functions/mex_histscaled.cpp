#include <mex.h>
#include <functional>
#include <algorithm>
#include <cstdint>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	/*
	input:
		var1 = data to be binned.
		var2 = scaling data.
		var3 = [minBin binSize maxBin].

	output:
		out1 = frequencies.
		out2 = bin centres.
	*/

	bool reterror = false;
	int errnum = 0;

	size_t M(0), N(0);

	if (nrhs > 1) {
		M = std::max(mxGetDimensions(prhs[0])[0], mxGetDimensions(prhs[0])[1]);
		N = std::max(mxGetDimensions(prhs[1])[0], mxGetDimensions(prhs[1])[1]);
	}
	else {
		reterror = true;
		errnum = 1;
	}

	double *binrange;
	size_t numOfBinPars = 0;
	if (nrhs == 3) {
		numOfBinPars = std::max(mxGetDimensions(prhs[2])[0], mxGetDimensions(prhs[2])[1]);
		if (numOfBinPars == 3) {
			binrange = (double*)mxGetData(prhs[2]);
		}
		else {
			reterror = true;
			errnum = 2;
		}
	}

	if (M > 0) {
		if (M == N) {
			double *data = (double*)mxGetData(prhs[0]);
			double *scaler = (double*)mxGetData(prhs[1]);

			double maxElem = *std::max_element(data, data + M);
			double minElem = *std::min_element(data, data + M);

			double binStart(0), binFinish(maxElem), binSize(1.);

			if (numOfBinPars == 3) {
				binStart = binrange[0];
				binSize = binrange[1];
				binFinish = binrange[2];
			}

			size_t numOfBins = ((binFinish - binStart) / binSize);

			plhs[0] = mxCreateNumericMatrix(1, numOfBins + 1, mxDOUBLE_CLASS, mxREAL);
			plhs[1] = mxCreateNumericMatrix(1, numOfBins + 1, mxDOUBLE_CLASS, mxREAL);

			double *hist = (double*) mxGetData(plhs[0]);
			double *bincentres = (double*)mxGetData(plhs[1]);

			for (size_t i = 0; i < numOfBins + 1; i++) {
				bincentres[i] = (double) i * binSize;
			}

			for (size_t i = 0; i < M; i++) {
				double x = data[i];
				size_t y = (x - binStart) / binSize;
				hist[y] += scaler[i];
			}
		}
		else {
			reterror = true;
			errnum = 4;
		}
	}
	else {
		reterror = true;
		errnum = 5;
	}


	if (reterror) {
		plhs[0] = mxCreateDoubleScalar(errnum);
	}
}