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
	-[double array] rms pop.
	-[double array] mean pop.
	-[double array[ total change in pop.
	-[double array] crowd lifetime.
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

		std::uint16_t L;

		std::vector<double> rmspops;
		std::vector<double> meanpops;
		std::vector<double> totcpops;
		std::vector<double> lifetimes;

		if (fexists(fname)) {
			mexPrintf("Loading...");

			FILE *file = fopen(fname.c_str(), "rb");

			if (vers == 1) {
				fread(&S, sizeof(std::uint32_t), 1, file);
				fread(&R, sizeof(std::uint32_t), 1, file);

				size_t ss = 0;
				for (size_t s = 0; s < S; s++) {
					if (rmspops.size() == rmspops.max_size() - 1) {
						break;
					}

					L = 0;
					fread(&L, sizeof(std::uint16_t), 1, file);

					std::uint16_t *data = new std::uint16_t[L]();
					fread(data, sizeof(std::uint16_t), L, file);
					bool useable = true;

					double rms = 0;
					double tot = 0;
					for (size_t lftm = 1; lftm < L; lftm++) {
						double sim = std::min(data[lftm], data[lftm - 1]) / sqrt(data[lftm] * data[lftm - 1]);
						rms += pow((double)(data[lftm] - data[lftm-1]), 2.);
						tot += (double) (data[lftm] - data[0]);

						if (sim < 0.7071){
							useable = false;
						}
					}


					if (useable) {
						rms /= (double) L;
						tot /= (double)(L - 1);

						rmspops.push_back(sqrt(L));
						meanpops.push_back(tot);
						totcpops.push_back((double)(data[L - 1] - data[0]));
						lifetimes.push_back(L);
					}

					delete[] data;
				}
			}
			else if (vers == 2) {

				fread(&R, sizeof(std::uint32_t), 1, file);

				for (size_t r = 0; r < R; r++) {
					fread(&S, sizeof(std::uint32_t), 1, file);
					for (size_t s = 0; s < S; s++) {
						if (rmspops.size() == rmspops.max_size() - 1) {
							break;
						}

						L = 0;
						fread(&L, sizeof(std::uint16_t), 1, file);

						std::uint16_t *data = new std::uint16_t[L]();
						fread(data, sizeof(std::uint16_t), L, file);
						bool useable = true;

						double rms = 0;
						double tot = 0;
						for (size_t lftm = 1; lftm < L; lftm++) {
							double sim = std::min(data[lftm], data[lftm - 1]) / sqrt(data[lftm] * data[lftm - 1]);
							rms += pow((double)(data[lftm] - data[0]), 2.);
							tot += (double)(data[lftm] - data[lftm-1]);
							if (sim < 0.7071){
								useable = false;
							}
						}


						if (useable) {
							rms /= (double) L;
							tot /= (double)(L - 1);
							rmspops.push_back(sqrt(rms));
							meanpops.push_back(tot);
							totcpops.push_back((double)(data[L - 1] - data[0]));
							lifetimes.push_back(L);
						}

						delete[] data;
					}
				}
			}

			fclose(file);

			plhs[0] = mxCreateNumericMatrix(1, rmspops.size(), mxDOUBLE_CLASS, mxREAL);
			double* popout = (double*)mxGetData(plhs[0]);
			memcpy(popout, &rmspops[0], sizeof(double)*rmspops.size());

			plhs[1] = mxCreateNumericMatrix(1, meanpops.size(), mxDOUBLE_CLASS, mxREAL);
			double* mpopout = (double*)mxGetData(plhs[1]);
			memcpy(mpopout, &meanpops[0], sizeof(double)*meanpops.size());

			plhs[2] = mxCreateNumericMatrix(1, totcpops.size(), mxDOUBLE_CLASS, mxREAL);
			double* tcpopout = (double*)mxGetData(plhs[2]);
			memcpy(tcpopout, &totcpops[0], sizeof(double)*totcpops.size());

			plhs[3] = mxCreateNumericMatrix(1, lifetimes.size(), mxDOUBLE_CLASS, mxREAL);
			double* lftmout = (double*)mxGetData(plhs[3]);
			memcpy(lftmout, &lifetimes[0], sizeof(double)*lifetimes.size());

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