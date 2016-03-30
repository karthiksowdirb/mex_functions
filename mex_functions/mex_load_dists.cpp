#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <string>
#include <cstdint>
#include <sys/stat.h>
#include <algorithm>

inline bool fexists(const std::string& name) {
	struct stat buffer;
	return (stat(name.c_str(), &buffer) == 0);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	if (nrhs == 3) {
		std::uint32_t S(0), R(0);

		const mwSize *dims = mxGetDimensions(prhs[0]);
		size_t strlen = dims[1] + 1;
		char *dir = (char*)malloc(sizeof(char)*strlen);
		mxGetString(prhs[0], dir, strlen);
		std::string dirstr(dir);
		free(dir);

		size_t LR = (size_t)mxGetScalar(prhs[1]);
		size_t UR = (size_t)mxGetScalar(prhs[2]);

		mexPrintf("[%d, %d]\n", LR, UR);

		//mexPrintf("%s\n", dirstr.c_str());

		std::string fname = dirstr + "/segment_dist.dat";
		std::uint16_t L;

		if (fexists(fname)) {
			//mexPrintf("File Exists\n");

			FILE *file = fopen(fname.c_str(), "rb");
			fread(&S, sizeof(std::uint32_t), 1, file);
			fread(&R, sizeof(std::uint32_t), 1, file);

			mexPrintf("[%d %d]\n", S, R);

			size_t bestSize = (size_t) std::min<size_t>((size_t) (UR - LR), (size_t) S);

			plhs[0] = mxCreateNumericMatrix(1, bestSize, mxDOUBLE_CLASS, mxREAL);
			plhs[1] = mxCreateNumericMatrix(1, bestSize, mxDOUBLE_CLASS, mxREAL);
			plhs[2] = mxCreateNumericMatrix(1, bestSize, mxDOUBLE_CLASS, mxREAL);
			plhs[3] = mxCreateNumericMatrix(1, bestSize, mxDOUBLE_CLASS, mxREAL);
			double* meandist = (double*)mxGetData(plhs[0]);
			double* vardist = (double*)mxGetData(plhs[1]);
			double* maxdist = (double*)mxGetData(plhs[2]);
			double* lifetimes = (double*)mxGetData(plhs[3]);

			std::fill(maxdist, maxdist + bestSize, 0);

			size_t ss = 0;
			for (size_t s = 0; s < S; s++) {
				if (s >= UR) {
					break;
				}

				if (s >= LR) {
					fread(&L, sizeof(std::uint16_t), 1, file);

					lifetimes[ss] = (double)L;

					float *data = new float[L]();
					fread(data, sizeof(float), L, file);

					for (size_t l = 0; l < L; l++) {
						meandist[ss] += (double) data[l];

						if (maxdist[ss] < (double) data[l]) {
							maxdist[ss] = (double) data[l];
						}
					}
					meandist[ss] /= (double) L;

					for (size_t l = 0; l < L; l++) {
						vardist[ss] += std::pow((double)data[l] - meandist[ss], 2.);
					}
					vardist[ss] /= (double) L;

					delete[] data;
	
					ss++;
				}
			}

			fclose(file);

			mexPrintf("Load Complete\n");
		}
		else {
			plhs[0] = mxCreateLogicalMatrix(1, 1);
			mexPrintf("File Does Not Exists.\n");
		}
	}
}