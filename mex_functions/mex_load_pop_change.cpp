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
	intput:
	-[string]		basedir.

	output:
	-[double array] lifetime.
	-[double array] change in population size.
	*/

	if (nrhs == 1) {
		std::uint32_t S(0), R(0);

		const mwSize *dims = mxGetDimensions(prhs[0]);
		size_t strlen = dims[1] + 1;
		char *dir = (char*)malloc(sizeof(char)*strlen);
		mxGetString(prhs[0], dir, strlen);
		std::string dirstr(dir);
		free(dir);

		std::string fname = dirstr + "/segment_pop.dat";
		mexPrintf("%s\n", fname.c_str());
		std::uint16_t L = 0;

		if (fexists(fname)) {
			mexPrintf("Loading...");

			FILE *file = fopen(fname.c_str(), "rb");

			if (file != NULL) {
				fread(&S, sizeof(std::uint32_t), 1, file);
				fread(&R, sizeof(std::uint32_t), 1, file);

				mexPrintf("[%u, %u]\n", S, R);

				//std::vector<double> inits;
				std::vector<double> lifetimes;
				std::vector<double> deltas;

				for (size_t s = 0; s < S; s++) {
					if (lifetimes.size() == lifetimes.max_size() - 1) {
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
							for (size_t lftm = 1; lftm < L; lftm++) {
								for (size_t il = 0; il + lftm < L; il++) {
									double cdel = ((double)data[il+lftm] - (double)data[il]);
									deltas.push_back(cdel);
									lifetimes.push_back(lftm);
								}
							}
					}


					delete[] data;
				}

				fclose(file);

				//plhs[0] = mxCreateNumericMatrix(1, inits.size(), mxDOUBLE_CLASS, mxREAL);
				plhs[0] = mxCreateNumericMatrix(1, lifetimes.size(), mxDOUBLE_CLASS, mxREAL);
				plhs[1] = mxCreateNumericMatrix(1, deltas.size(), mxDOUBLE_CLASS, mxREAL);

				//double* initpop = (double*)mxGetData(plhs[0]);
				double* lifetime = (double*)mxGetData(plhs[0]);
				double* deltapop = (double*)mxGetData(plhs[1]);

				//memcpy(initpop, &inits[0], sizeof(double)*inits.size());
				memcpy(lifetime, &lifetimes[0], sizeof(double)*lifetimes.size());
				memcpy(deltapop, &deltas[0], sizeof(double)*deltas.size());


				mexPrintf("Load Complete\n");
			}
			else {
				plhs[0] = mxCreateLogicalScalar(0);
				mexPrintf("File Not Loaded!\n");
			}
		}
		else {
			plhs[0] = mxCreateLogicalScalar(0);
			mexPrintf("File Does Not Exists.\n");
		}
	}
	else {
		mexPrintf("Wrong Number of Input Arguments.\n USAGE \n\t[ipop, dpop] = mex_load_pop_change(basedir); clear mex_load_pop_change;\n");
		plhs[0] = mxCreateLogicalScalar(0);
	}
}