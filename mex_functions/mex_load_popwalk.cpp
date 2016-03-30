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
	-[scalar]		Lifetime.

	output:
	-[double matrix] population at each time step for each crowd.
	*/

	if (nrhs == 2) {
		std::uint32_t S(0), R(0);

		const mwSize *dims = mxGetDimensions(prhs[0]);
		size_t strlen = dims[1] + 1;
		char *dir = (char*)malloc(sizeof(char)*strlen);
		mxGetString(prhs[0], dir, strlen);
		std::string dirstr(dir);
		free(dir);

		size_t LFTM = mxGetScalar(prhs[1]);

		std::string fnamev1 = dirstr + "/segment_pop.dat";
		std::string fnamev2 = dirstr + "/segment_pop.v2.dat";
		std::string fname;
		size_t vers = 0;
		if (fexists(fnamev1)) {
			vers = 1;
			fname = fnamev1;
		}
		else if (fexists(fnamev2)) {
			vers = 2;
			fname = fnamev2;
		}

		mexPrintf("%s\n", fname.c_str());
		std::uint16_t L = 0;

		if (fexists(fname)) {
			mexPrintf("Loading...");

			FILE *file = fopen(fname.c_str(), "rb");

			//std::vector<double> inits;
			std::vector<std::vector<double>> walks;

			if (file != NULL) {
				if (vers == 1) {
					fread(&S, sizeof(std::uint32_t), 1, file);
					fread(&R, sizeof(std::uint32_t), 1, file);

					

					for (size_t s = 0; s < S; s++) {
						if (walks.size() == walks.max_size() - 1) {
							break;
						}

						L = 0;
						fread(&L, sizeof(std::uint16_t), 1, file);

						std::uint16_t *data = new std::uint16_t[L]();
						fread(data, sizeof(std::uint16_t), L, file);

						if (L == LFTM) {
							bool useable = true;
							/*
							for (size_t lftm = 1; lftm < L; lftm++) {
							double sim = std::min(data[lftm], data[lftm - 1]) / sqrt(data[lftm] * data[lftm - 1]);
							if (sim < 0.7071){
							useable = false;
							}
							}
							*/

							if (useable) {
								std::vector<double> walk;
								for (size_t lftm = 0; lftm < L; lftm++) {
									walk.push_back((double)data[lftm]);
								}
								walks.push_back(walk);
							}
						}

						delete[] data;
					}
				}
				else if (vers == 2) {
					fread(&R, sizeof(std::uint32_t), 1, file);

					//std::vector<double> inits;
					std::vector<std::vector<double>> walks;
					for (size_t r = 1; r < R; r++) {
						fread(&S, sizeof(std::uint32_t), 1, file);
						for (size_t s = 0; s < S; s++) {
							if (walks.size() == walks.max_size() - 1) {
								break;
							}

							L = 0;
							fread(&L, sizeof(std::uint16_t), 1, file);

							std::uint16_t *data = new std::uint16_t[L]();
							fread(data, sizeof(std::uint16_t), L, file);

							if (L == LFTM) {
								bool useable = true;
								/*
								for (size_t lftm = 1; lftm < L; lftm++) {
								double sim = std::min(data[lftm], data[lftm - 1]) / sqrt(data[lftm] * data[lftm - 1]);
								if (sim < 0.7071){
								useable = false;
								}
								}
								*/

								if (useable) {
									std::vector<double> walk;
									for (size_t lftm = 0; lftm < L; lftm++) {
										walk.push_back((double)data[lftm]);
									}
									walks.push_back(walk);
								}
							}

							delete[] data;
						}
					}
				}

				fclose(file);

				//plhs[0] = mxCreateNumericMatrix(1, inits.size(), mxDOUBLE_CLASS, mxREAL);
				plhs[0] = mxCreateNumericMatrix(LFTM, walks.size(), mxDOUBLE_CLASS, mxREAL);

				//double* initpop = (double*)mxGetData(plhs[0]);
				double* cwalks = (double*)mxGetData(plhs[0]);

				for (size_t c = 0; c < walks.size(); c++) {
					cwalks[c*LFTM] = 1;
					for (size_t l = 1; l < walks[c].size(); l++) {
						cwalks[c*LFTM + l] = walks[c][l] / walks[c][l - 1];
					}
				}

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