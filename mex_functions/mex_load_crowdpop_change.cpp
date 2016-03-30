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
		-[double array] initial population size.
		-[double array] change in population size.
		-[double array] lifetimes.
		-[double]		Number of Realisations.
	*/

	if (nrhs == 1) {
		std::uint32_t S(0), R(0);

		const mwSize *dims = mxGetDimensions(prhs[0]);
		size_t strlen = dims[1] + 1;
		char *dir = (char*)malloc(sizeof(char)*strlen);
		mxGetString(prhs[0], dir, strlen);
		std::string dirstr(dir);
		free(dir);

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
			std::vector<double> ipops;
			std::vector<double> maxpops;
			std::vector<double> minpops;
			std::vector<double> avepops;
			std::vector<double> cdeltas;
			std::vector<double> lifetimes;

			if (file != NULL) {
				if (vers == 2) {
					fread(&R, sizeof(std::uint32_t), 1, file);

					mexPrintf("[%u, %u]\n", S, R);

					for (size_t r = 0; r < R; r++) {
						fread(&S, sizeof(std::uint32_t), 1, file);
						for (size_t s = 0; s < S; s++) {
							if (ipops.size() == ipops.max_size() - 1) {
								break;
							}

							L = 0;
							fread(&L, sizeof(std::uint16_t), 1, file);

							std::uint16_t *data = new std::uint16_t[L]();
							fread(data, sizeof(std::uint16_t), L, file);

							bool useable = true;

							double tot = (double)data[0];
							for (size_t lftm = 1; lftm < L; lftm++) {
								tot += data[lftm];
								double sim = std::min(data[lftm], data[lftm - 1]) / sqrt(data[lftm] * data[lftm - 1]);
								if (sim < 0.7071){
									useable = false;
								}
							}


							if (useable) {

								/*
								for (size_t lftm = 1; lftm < L; lftm++) {
									ipops.push_back(data[0]);
									double cdel = ((double)data[lftm] - (double)data[0]);
									cdeltas.push_back(cdel);
									lifetimes.push_back(lftm);
								}
								*/
								

								
								ipops.push_back(data[0]);
								maxpops.push_back(*std::max_element(data, data + L));
								minpops.push_back(*std::min_element(data, data + L));
								avepops.push_back(tot / L);
								double cdel = ((double)data[L - 1] - (double)data[0]);
								cdeltas.push_back(cdel);
								lifetimes.push_back(L);
							}


							delete[] data;
						}
					}
				}
				else if (vers == 1) {
					fread(&S, sizeof(std::uint32_t), 1, file);
					fread(&R, sizeof(std::uint32_t), 1, file);

					mexPrintf("[%u, %u]\n", S, R);
					
					for (size_t s = 0; s < S; s++) {
						if (ipops.size() == ipops.max_size() - 1) {
							break;
						}

						L = 0;
						fread(&L, sizeof(std::uint16_t), 1, file);

						std::uint16_t *data = new std::uint16_t[L]();
						fread(data, sizeof(std::uint16_t), L, file);

						bool useable = true;

						double tot = (double) data[0];
						for (size_t lftm = 1; lftm < L; lftm++) {
							tot += (double) data[lftm];
							double sim = std::min(data[lftm], data[lftm - 1]) / sqrt(data[lftm] * data[lftm - 1]);
							if (sim < 0.7071){
								useable = false;
							}
						}
						

						if (useable) {
							/*
							for (size_t lftm = 1; lftm < L; lftm++) {
								ipops.push_back(data[0]);
								double cdel = ((double)data[lftm] - (double)data[0]);
								cdeltas.push_back(cdel);
								lifetimes.push_back(lftm);
							}
							*/
							
							ipops.push_back(data[0]);
							maxpops.push_back(*std::max_element(data, data + L));
							minpops.push_back(*std::min_element(data, data + L));
							avepops.push_back(tot / L);
							double cdel = ((double)data[L - 1] - (double)data[0]);
							cdeltas.push_back(cdel);
							lifetimes.push_back(L);
						}


						delete[] data;
					}
				}

				fclose(file);

				//plhs[0] = mxCreateNumericMatrix(1, inits.size(), mxDOUBLE_CLASS, mxREAL);
				plhs[0] = mxCreateNumericMatrix(1, ipops.size(), mxDOUBLE_CLASS, mxREAL);
				plhs[1] = mxCreateNumericMatrix(1, cdeltas.size(), mxDOUBLE_CLASS, mxREAL);
				plhs[2] = mxCreateNumericMatrix(1, lifetimes.size(), mxDOUBLE_CLASS, mxREAL);
				plhs[3] = mxCreateNumericMatrix(1, maxpops.size(), mxDOUBLE_CLASS, mxREAL);
				plhs[4] = mxCreateNumericMatrix(1, minpops.size(), mxDOUBLE_CLASS, mxREAL);
				plhs[5] = mxCreateNumericMatrix(1, avepops.size(), mxDOUBLE_CLASS, mxREAL);

				//double* initpop = (double*)mxGetData(plhs[0]);
				double* ipop = (double*)mxGetData(plhs[0]);
				double* cdeltapop = (double*)mxGetData(plhs[1]);
				double* lifepop = (double*)mxGetData(plhs[2]);
				double* maxpop = (double*)mxGetData(plhs[3]);
				double* minpop = (double*)mxGetData(plhs[4]);
				double* avepop = (double*)mxGetData(plhs[5]);

				//memcpy(initpop, &inits[0], sizeof(double)*inits.size());
				memcpy(ipop, &ipops[0], sizeof(double)*ipops.size());
				memcpy(cdeltapop, &cdeltas[0], sizeof(double)*cdeltas.size());
				memcpy(lifepop, &lifetimes[0], sizeof(double)*lifetimes.size());
				memcpy(maxpop, &maxpops[0], sizeof(double)*maxpops.size());
				memcpy(minpop, &minpops[0], sizeof(double)*minpops.size());
				memcpy(avepop, &avepops[0], sizeof(double)*avepops.size());

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