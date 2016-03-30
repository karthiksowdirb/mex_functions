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

inline double getHurst(const std::vector<double> x, const std::vector<double> y) {
	double mx(0), my(0), mxy(0), mx2(0);

	for (size_t i = 0; i < x.size(); i++) {
		mx += x[i];
		mx2 += x[i] * x[i];
		my += y[i];
		mxy += x[i] * y[i];
	}
	mx /= (double)x.size();
	mx2 /= (double)x.size();
	my /= (double)y.size();
	mxy /= (double)x.size();

	return (mxy - mx*my) / (mx2 - mx*mx); 
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/*
	* Rescaled Range analysis.
	input:
	-[string]		directory,
	output:
	-[double array] population sizes.
	*/


	if (nrhs == 2) {
		std::uint32_t S(0), R(0);

		const mwSize *dims = mxGetDimensions(prhs[0]);
		size_t strlen = dims[1] + 1;
		char *dir = (char*)malloc(sizeof(char)*strlen);
		mxGetString(prhs[0], dir, strlen);
		std::string dirstr(dir);
		free(dir);

		const size_t lftmThreshold = (size_t) mxGetScalar(prhs[1]);


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
		std::uint16_t L;

		std::vector<double> hurst;
		std::vector<double> lifetimes;

		if (fexists(fname)) {
			mexPrintf("Loading...");

			FILE *file = fopen(fname.c_str(), "rb");

			if (vers == 2) {
				fread(&R, sizeof(std::uint32_t), 1, file);
				for (size_t r = 0; r < R; r++) {
					fread(&S, sizeof(std::uint32_t), 1, file);
					for (size_t s = 0; s < S; s++) {
						if (hurst.size() == hurst.max_size() - 1) {
							break;
						}

						L = 0;
						fread(&L, sizeof(std::uint16_t), 1, file);

						std::uint16_t *data = new std::uint16_t[L]();
						fread(data, sizeof(std::uint16_t), L, file);

						if (L >= lftmThreshold) {
							bool useable = true;
							for (size_t lftm = 1; lftm < L; lftm++) {
								double sim = std::min(data[lftm], data[lftm - 1]) / sqrt(data[lftm] * data[lftm - 1]);
								if (sim < 0.7071){
									useable = false;
								}
							}

							if (useable) {
								const size_t N = log((double)L) / log(2) - 2;
								size_t *lengths = new size_t[N]();
								double *means = new double[L*N]();
								double *stdev = new double[L*N]();
								double *Y = new double[L*N]();		//The mean adjusted series.
								double *Z = new double[L*N]();		//The cumulative deviate series.
								double *minZ = new double[L*N]();	//minimum Z cumulative deviate series;
								double *maxZ = new double[L*N]();	//maximum Z cumulative deviate series;
								double *R = new double[L*N]();		//Range of the cumulative deviate series;

								//Calculate the lengths of each chunk.
								for (size_t c = 0; c < N; c++) {
									lengths[c] = (size_t)((double)L / pow(2, c));
								}

								//Calculate the means of each chunk.
								for (size_t lftm = 0; lftm < L; lftm++) {
									double rolsum = 0;
									size_t ct = lftm;
									for (size_t c = 0; c < N; c++) {
										size_t index = lftm / lengths[c];
										means[L*c + index] += (double)data[lftm];
									}
								}

								for (size_t c = 0; c < N; c++) {
									size_t maxindex = pow(2, c);
									for (size_t index = 0; index < maxindex; index++) {
										means[L*c + index] /= (double)lengths[c];
									}
								}

								//Calculate the stdev, Ys and Zs.
								for (size_t lftm = 0; lftm < L; lftm++) {
									for (size_t c = 0; c < N; c++) {
										size_t index = lftm / lengths[c];
										Y[L*c + lftm] = (double)data[lftm] - means[L*c + index];
										stdev[L*c + index] += Y[L*c + lftm] * Y[L*c + lftm];
										Z[L*c + lftm] += Y[L*c + lftm];

										if (lftm % lengths[c] == 0) {
											Z[L*c + lftm] = Y[L*c + lftm];
											minZ[L*c + index] = Z[L*c + lftm];
											maxZ[L*c + index] = Z[L*c + lftm];
										}

										//Get the minimum Z for each chunk.
										if (minZ[L*c + index] > Z[L*c + lftm]) {
											minZ[L*c + index] = Z[L*c + lftm];
										}

										//Get the maximum Z for each chunk.
										if (maxZ[L*c + index] < Z[L*c + lftm]) {
											maxZ[L*c + index] = Z[L*c + lftm];
										}
									}
								}

								for (size_t c = 0; c < N; c++) {
									for (size_t index = 0; index < pow(2, c); index++) {
										R[L*c + index] = maxZ[L*c + index] - minZ[L*c + index];
									}
								}

								std::vector<double> log_averR;
								std::vector<double> log_chunklength;

								//Calculate the average rescaled range at each level.
								for (size_t c = 0; c < N; c++) {
									size_t maxindex = pow(2, c);

									double rR = 0;
									for (size_t index = 0; index < maxindex; index++) {
										stdev[L*c + index] /= (double)lengths[c];
										stdev[L*c + index] = sqrt(stdev[L*c + index]);
										if (stdev[L*c + index] > 0) {
											rR += R[L*c + index] / stdev[L*c + index];
										}
									}
									rR /= (double)maxindex;

									log_averR.push_back(std::log(rR));
									log_chunklength.push_back(std::log(lengths[c]));
								}

								double H = getHurst(log_chunklength, log_averR);
								hurst.push_back(H);
								lifetimes.push_back(L);

								delete[] lengths;
								delete[] means;
								delete[] stdev;
								delete[] Y;
								delete[] Z;
								delete minZ;
								delete maxZ;
								delete R;
							}
						}
						delete[] data;
					}
				}
			}
			else {
				fread(&S, sizeof(std::uint32_t), 1, file);
				fread(&R, sizeof(std::uint32_t), 1, file);

				for (size_t s = 0; s < S; s++) {
					if (hurst.size() == hurst.max_size() - 1) {
						break;
					}

					L = 0;
					fread(&L, sizeof(std::uint16_t), 1, file);

					std::uint16_t *data = new std::uint16_t[L]();
					fread(data, sizeof(std::uint16_t), L, file);

					if (L >= lftmThreshold) {
						bool useable = true;
						for (size_t lftm = 1; lftm < L; lftm++) {
							double dp = abs(data[lftm] - data[lftm - 1]);
							if (dp < data[lftm-1] && dp > (double) data[lftm-1]/2){
								useable = false;
							}
						}

						if (useable) {
							const size_t N = log((double)L) / log(2) - 2;
							size_t *lengths = new size_t[N]();
							double *means = new double[L*N]();
							double *stdev = new double[L*N]();
							double *Y = new double[L*N]();		//The mean adjusted series.
							double *Z = new double[L*N]();		//The cumulative deviate series.
							double *minZ = new double[L*N]();	//minimum Z cumulative deviate series;
							double *maxZ = new double[L*N]();	//maximum Z cumulative deviate series;
							double *R = new double[L*N]();		//Range of the cumulative deviate series;

							//Calculate the lengths of each chunk.
							for (size_t c = 0; c < N; c++) {
								lengths[c] = (size_t)((double)L / pow(2, c));
							}

							//Calculate the means of each chunk.
							for (size_t lftm = 0; lftm < L; lftm++) {
								double rolsum = 0;
								size_t ct = lftm;
								for (size_t c = 0; c < N; c++) {
									size_t index = lftm / lengths[c];
									means[L*c + index] += (double)data[lftm];
								}
							}

							for (size_t c = 0; c < N; c++) {
								size_t maxindex = pow(2, c);
								for (size_t index = 0; index < maxindex; index++) {
									means[L*c + index] /= (double)lengths[c];
								}
							}

							//Calculate the stdev, Ys and Zs.
							for (size_t lftm = 0; lftm < L; lftm++) {
								for (size_t c = 0; c < N; c++) {
									size_t index = lftm / lengths[c];
									Y[L*c + lftm] = (double)data[lftm] - means[L*c + index];
									stdev[L*c + index] += Y[L*c + lftm] * Y[L*c + lftm];
									Z[L*c + lftm] += Y[L*c + lftm];

									if (lftm % lengths[c] == 0) {
										Z[L*c + lftm] = Y[L*c + lftm];
										minZ[L*c + index] = Z[L*c + lftm];
										maxZ[L*c + index] = Z[L*c + lftm];
									}

									//Get the minimum Z for each chunk.
									if (minZ[L*c + index] > Z[L*c + lftm]) {
										minZ[L*c + index] = Z[L*c + lftm];
									}

									//Get the maximum Z for each chunk.
									if (maxZ[L*c + index] < Z[L*c + lftm]) {
										maxZ[L*c + index] = Z[L*c + lftm];
									}
								}
							}

							for (size_t c = 0; c < N; c++) {
								for (size_t index = 0; index < pow(2, c); index++) {
									R[L*c + index] = maxZ[L*c + index] - minZ[L*c + index];
								}
							}

							std::vector<double> log_averR;
							std::vector<double> log_chunklength;

							//Calculate the average rescaled range at each level.
							for (size_t c = 0; c < N; c++) {
								size_t maxindex = pow(2, c);

								double rR = 0;
								for (size_t index = 0; index < maxindex; index++) {
									stdev[L*c + index] /= (double)lengths[c];
									stdev[L*c + index] = sqrt(stdev[L*c + index]);
									if (stdev[L*c + index] > 0) {
										rR += R[L*c + index] / stdev[L*c + index];
									}
								}
								rR /= (double)maxindex;

								log_averR.push_back(std::log(rR));
								log_chunklength.push_back(std::log(lengths[c]));
							}

							double H = getHurst(log_chunklength, log_averR);
							hurst.push_back(H);
							lifetimes.push_back(L);

							delete[] lengths;
							delete[] means;
							delete[] stdev;
							delete[] Y;
							delete[] Z;
							delete minZ;
							delete maxZ;
							delete R;
						}
					}
					delete[] data;
				}
			}

			fclose(file);

			plhs[0] = mxCreateNumericMatrix(1, hurst.size(), mxDOUBLE_CLASS, mxREAL);
			double* hurstout = (double*)mxGetData(plhs[0]);
			memcpy(hurstout, &hurst[0], sizeof(double)*hurst.size());

			plhs[1] = mxCreateNumericMatrix(1, lifetimes.size(), mxDOUBLE_CLASS, mxREAL);
			double* lftmout = (double*)mxGetData(plhs[1]);
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