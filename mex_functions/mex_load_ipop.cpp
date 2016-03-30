#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <string>
#include <cstdint>
#include <sys/stat.h>
#include <algorithm>
#include <vector>

inline bool fexists(const std::string& name) {
	struct stat buffer;
	return (stat(name.c_str(), &buffer) == 0);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/*
	input:
		-[string]		directory,

	output:
		-[double array] crowd population sizes.
		-[double array] crowd population lifetimes.
	*/


	if (nrhs == 1) {
		std::uint32_t S(0), R(0);

		const mwSize *dims = mxGetDimensions(prhs[0]);
		size_t strlen = dims[1] + 1;
		char *dir = (char*)malloc(sizeof(char)*strlen);
		mxGetString(prhs[0], dir, strlen);
		std::string dirstr(dir);
		free(dir);

		//mexPrintf("[%d, %d]\n", LR, UR);

		//mexPrintf("%s\n", dirstr.c_str());

		std::string fname = dirstr + "/segment_pop.dat";
		mexPrintf("%s\n", fname.c_str());
		std::uint16_t L;

		std::vector<double> pops;
		std::vector<double> lifetimes;

		if (fexists(fname)) {
			mexPrintf("Loading...");

			FILE *file = fopen(fname.c_str(), "rb");
			fread(&S, sizeof(std::uint32_t), 1, file);
			fread(&R, sizeof(std::uint32_t), 1, file);

			size_t ss = 0;
			for (size_t s = 0; s < S; s++) {
				if (pops.size() == pops.max_size() - 1) {
					break;
				}

				L = 0;
				fread(&L, sizeof(std::uint16_t), 1, file);

				std::uint16_t *data = new std::uint16_t[L]();
				fread(data, sizeof(std::uint16_t), L, file);
				bool useable = true;
				for (size_t lftm = 1; lftm < L; lftm++) {
					double sim = std::min(data[lftm], data[lftm - 1]) / sqrt(data[lftm] * data[lftm - 1]);
					if (sim < 0.7071){
						useable = false;
					}
				}

				if (useable) {
					double maxpop = 0;
					pops.push_back(data[0]);
					lifetimes.push_back(L);
				}

				delete[] data;
			}

			fclose(file);

			plhs[0] = mxCreateNumericMatrix(1, pops.size(), mxDOUBLE_CLASS, mxREAL);
			plhs[1] = mxCreateNumericMatrix(1, lifetimes.size(), mxDOUBLE_CLASS, mxREAL);

			double* popout = (double*)mxGetData(plhs[0]);
			memcpy(popout, &pops[0], sizeof(double)*pops.size());

			double* poplftms = (double*)mxGetData(plhs[1]);
			memcpy(poplftms, &lifetimes[0], sizeof(double)*lifetimes.size());

			mexPrintf("Load Complete\n");
		}
		else {
			plhs[0] = mxCreateLogicalScalar(0);
			mexPrintf("File Does Not Exists.\n");
		}
	}
	else {
		mexPrintf("Wrong Number of Input Arguments.\n");
		plhs[0] = mxCreateLogicalScalar(0);
	}
}