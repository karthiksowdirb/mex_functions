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
	-[string]		Directory,
	output:
	-[double array] Diameter (Maximum Path Lengths).
	-[double array] Population.
	*/

	if (nrhs == 1) {
		std::uint32_t S(0), R(0);

		const mwSize *dims = mxGetDimensions(prhs[0]);
		size_t strlen = dims[1] + 1;
		char *dir = (char*)malloc(sizeof(char)*strlen);
		mxGetString(prhs[0], dir, strlen);
		std::string dirstr(dir);
		free(dir);

		mexPrintf("%s\n", dirstr.c_str());

		std::string fname = dirstr + "/segment_pop.dat";
		std::string afname = dirstr + "/segment_mpl.dat";

		std::uint16_t L;

		if (fexists(afname) & fexists(fname)) {
			//mexPrintf("File Exists\n");

			FILE *file = fopen(fname.c_str(), "rb");
			fread(&S, sizeof(std::uint32_t), 1, file);
			fread(&R, sizeof(std::uint32_t), 1, file);

			std::vector<double> popvect;
			std::vector<bool>	usevect(S, true);

			for (size_t s = 0; s < S; s++) {
				if (popvect.size() == popvect.max_size() - 1) {
					break;
				}

				L = 0;
				fread(&L, sizeof(std::uint16_t), 1, file);

				std::uint16_t *data = new std::uint16_t[L]();
				fread(data, sizeof(std::uint16_t), L, file);
				for (size_t lftm = 1; lftm < L; lftm++) {
					double sim = std::min(data[lftm], data[lftm - 1]) / sqrt(data[lftm] * data[lftm - 1]);
					if (sim < 0.7071){
						usevect[s] = false;
					}
				}
				if (usevect[s]) {
					usevect.push_back(true);
					for (size_t lftm = 0; lftm < L; lftm++) {
						popvect.push_back((double)data[lftm]);
					}
				}

				delete[] data;
			}

			fclose(file);

			///////////////////////////////////////////////////////////

			std::vector<double> aplvect;

			FILE *afile = fopen(afname.c_str(), "rb");
			fread(&S, sizeof(std::uint32_t), 1, afile);
			fread(&R, sizeof(std::uint32_t), 1, afile);

			for (size_t s = 0; s < S; s++) {
				if (aplvect.size() == aplvect.max_size() - 1) {
					break;
				}

				L = 0;
				fread(&L, sizeof(std::uint16_t), 1, afile);

				float *data = new float[L]();
				fread(data, sizeof(float), L, afile);

				if (usevect[s] == true) {
					for (size_t lftm = 0; lftm < L; lftm++) {
						aplvect.push_back((double)data[lftm]);
					}
				}

				delete[] data;
			}

			fclose(afile);

			//////////////////////////////////////////////////////////

			plhs[0] = mxCreateNumericMatrix(1, popvect.size(), mxDOUBLE_CLASS, mxREAL);
			double* pops = (double*)mxGetData(plhs[0]);
			memcpy(pops, &popvect[0], sizeof(double)*popvect.size());

			plhs[1] = mxCreateNumericMatrix(1, aplvect.size(), mxDOUBLE_CLASS, mxREAL);
			double* apls = (double*)mxGetData(plhs[1]);
			memcpy(apls, &aplvect[0], sizeof(double)*aplvect.size());



			mexPrintf("Load Complete\n");
		}
		else {
			plhs[0] = mxCreateLogicalMatrix(1, 1);
			mexPrintf("File Does Not Exists.\n");
		}
	}
}