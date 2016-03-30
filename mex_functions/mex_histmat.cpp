#include <mex.h>
#include <functional>
#include <algorithm>
#include <cstdint>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	/*
		input:
			-[array] Data 1.
			-[array] [MinValue StepValue MaxValue].
			-[array] Data 2.
			-[array] [MinValue StepValue MaxValue].
			-[array] (Optional) Scaler.

		output:
			-[double mat] histmat.
			-[double array] Bins 1.
			-[Double array] Bins 2.
	*/


	bool reterror = false;
	int errnum = 0;

	size_t M(0), N(0), sclrlength(0), numOfBinPars1(0), numOfBinPars2(0);

	if (nrhs == 4 || nrhs == 5) {
		M = std::max(mxGetDimensions(prhs[0])[0], mxGetDimensions(prhs[0])[1]);
		N = std::max(mxGetDimensions(prhs[2])[0], mxGetDimensions(prhs[2])[1]);

		numOfBinPars1 = std::max(mxGetDimensions(prhs[1])[0], mxGetDimensions(prhs[1])[1]);
		numOfBinPars2 = std::max(mxGetDimensions(prhs[3])[0], mxGetDimensions(prhs[3])[1]);
	}

	if (nrhs == 5) {
		sclrlength = std::max(mxGetDimensions(prhs[4])[0], mxGetDimensions(prhs[4])[1]);
	}

	//mexPrintf("[%d, %d, %d]\n", M, N, sclrlength);

	if (M > 0) {
		if (M == N) {
			double *data1 = (double*)mxGetData(prhs[0]);
			double *data2 = (double*)mxGetData(prhs[2]);
			double *sclr;

			if (sclrlength == M) {
				sclr = (double*)mxGetData(prhs[4]);
			}

			size_t numOfBins1(100), numOfBins2(100);
			double binminv1, binstep1, binmaxv1;
			double binminv2, binstep2, binmaxv2;

			if (numOfBinPars1 == 3) {
				binminv1 = (double)((double*)mxGetData(prhs[1]))[0];
				binstep1 = (double)((double*)mxGetData(prhs[1]))[1];
				binmaxv1 = (double)((double*)mxGetData(prhs[1]))[2];

				numOfBins1 = (binmaxv1 - binminv1) / binstep1;
			}
			else {
				binminv1 = *std::min_element(data1, data1+M);
				binmaxv1 = *std::max_element(data1, data1 + M);
				binstep1 = (binmaxv1 - binminv1) / numOfBins1;
			}

			if (numOfBinPars2 == 3) {
				binminv2 = (double)((double*)mxGetData(prhs[3]))[0];
				binstep2 = (double)((double*)mxGetData(prhs[3]))[1];
				binmaxv2 = (double)((double*)mxGetData(prhs[3]))[2];

				numOfBins2 = (binmaxv2 - binminv2) / binstep2;
			}
			else {
				binminv2 = *std::min_element(data2, data1 + M);
				binmaxv2 = *std::max_element(data2, data1 + M);
				binstep2 = (binmaxv2 - binminv2) / numOfBins2;
			}

			//mexPrintf("1: [%f %f %f] = %d\n", binminv1, binstep1, binmaxv1, numOfBins1);
			//mexPrintf("2: [%f %f %f] = %d\n", binminv2, binstep2, binmaxv2, numOfBins2);

			size_t numOfCols = numOfBins1 + 1;
			size_t numOfRows = numOfBins2 + 1;

			plhs[0] = mxCreateNumericMatrix(numOfRows, numOfCols, mxDOUBLE_CLASS, mxREAL);

			plhs[1] = mxCreateNumericMatrix(1, numOfCols, mxDOUBLE_CLASS, mxREAL);
			plhs[2] = mxCreateNumericMatrix(1, numOfRows, mxDOUBLE_CLASS, mxREAL);

			double *histmat = (double*)mxGetData(plhs[0]);
			double *dom1 = (double*)mxGetData(plhs[1]);
			double *dom2 = (double*)mxGetData(plhs[2]);

			size_t length = numOfRows*numOfCols;

			size_t numOfBinsMax = std::max(numOfRows, numOfCols);
			for (size_t i = 0; i < numOfBinsMax; i++) {
				if (i < numOfCols) {
					dom1[i] = binminv1 + i*binstep1;
				}

				if (i < numOfRows) {
					dom2[i] = binminv2 + i*binstep2;
				}
			}
			
			for (size_t i = 0; i < M; i++) {
				int x = (size_t) ( ( (double) data1[i] - binminv1 ) / binstep1 );
				int y = (size_t) ( ( (double) data2[i] - binminv2 ) / binstep2 );

				//Ensure that x and y are within the bin limits. The first and last bins will collect
				//frequencies which are out of bounds.
				if (x < 0) {
					x = 0;
				}
				else if (x > numOfBins1) {
					x = numOfBins1;
				}
				
				if (y < 0) {
					y = 0;
				}
				else if (y > numOfBins2) {
					y = numOfBins2;
				}

				if (sclrlength == M) {
					histmat[x*numOfRows + y] += sclr[i];
				} else {
					histmat[x*numOfRows + y]++;
				}
			}
		} 
		else {
			reterror = true;
			errnum = 2;
		}
	}
	else {
		reterror = true;
		errnum = 1;
	}

	if (reterror) {
		plhs[0] = mxCreateDoubleScalar(errnum);
	}
}