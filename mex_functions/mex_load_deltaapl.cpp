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
	-[double array] Delta Average Path Lengths.
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

		std::string afname = dirstr + "/segment_apl.dat";
		std::string pfname = dirstr + "/segment_pop.dat";

		std::uint16_t L;

		if (fexists(afname) & fexists(pfname)) {
			//mexPrintf("File Exists\n");

			//Open the population file for reading.
			FILE *pfile = fopen(pfname.c_str(), "rb");
			fread(&S, sizeof(std::uint32_t), 1, pfile);
			fread(&R, sizeof(std::uint32_t), 1, pfile);

			std::vector<bool> useable;
			std::vector<double> daplvect;
			std::vector<double> lifetimes;

			for (size_t s = 0; s < S; s++) {
				//Get lifetimes from file. (Should be the same).
				fread(&L, sizeof(std::uint16_t), 1, pfile);

				//Read the population data into a buffer and check to make sure the crowd data is valid.
				std::uint16_t *pdata = new std::uint16_t[L]();
				fread(pdata, sizeof(std::uint16_t), L, pfile);
				bool tmpuseable = true;
				for (size_t lftm = 1; lftm < L; lftm++) {
					double sim = std::min(pdata[lftm], pdata[lftm - 1]) / sqrt(pdata[lftm] * pdata[lftm - 1]);
					if (sim < 0.7071){
						tmpuseable = false;
					}
				}

				//If the data is valid then read the apls into a vector.
				if (tmpuseable) {
					useable.push_back(true);
				}
				else {
					useable.push_back(false);
				}

				//Delete the temporary buffers.
				delete[] pdata;
			}

			fclose(pfile);
			////////////////////////////////

			//Open the APL file for reading.
			FILE *afile = fopen(afname.c_str(), "rb");
			fread(&S, sizeof(std::uint32_t), 1, afile);
			fread(&R, sizeof(std::uint32_t), 1, afile);

			for (size_t s = 0; s < S; s++) {
				//Get lifetimes from file. (Should be the same).
				fread(&L, sizeof(std::uint16_t), 1, afile);

				//Load the APL data for the crowd into a buffer.
				float *adata = new float[L]();
				fread(adata, sizeof(float), L, afile);

				if (useable[s]) {
					double del = adata[0] - adata[L - 1];
					daplvect.push_back(del);
					lifetimes.push_back(L);
				}

				delete[] adata;
			}

			fclose(afile);

			///////////////////////////////////////////////////////////

			//Output the valid APL data to MatLab.
			plhs[0] = mxCreateNumericMatrix(1, daplvect.size(), mxDOUBLE_CLASS, mxREAL);
			plhs[1] = mxCreateNumericMatrix(1, lifetimes.size(), mxDOUBLE_CLASS, mxREAL);

			double* apls = (double*)mxGetData(plhs[0]);
			memcpy(apls, &daplvect[0], sizeof(double)*daplvect.size());

			double* lftms = (double*)mxGetData(plhs[1]);
			memcpy(lftms, &lifetimes[0], sizeof(double)*lifetimes.size());

			mexPrintf("Load Complete\n");
		}
		else {
			plhs[0] = mxCreateLogicalMatrix(1, 1);
			mexPrintf("File Does Not Exists.\n");
		}
	}
}