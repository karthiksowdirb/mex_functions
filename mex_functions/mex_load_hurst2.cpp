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
	-[string]		directory.
	output:
	-[double array] widths.
	-[double array] chunk sizes.

	*/


	if (nrhs == 3) {
		std::uint32_t S(0), R(0);

		const mwSize *dims = mxGetDimensions(prhs[0]);
		size_t strlen = dims[1] + 1;
		char *dir = (char*)malloc(sizeof(char)*strlen);
		mxGetString(prhs[0], dir, strlen);
		std::string dirstr(dir);
		free(dir);
		const size_t L_lt = (size_t) mxGetScalar(prhs[1]);
		const size_t L_ut = (size_t) mxGetScalar(prhs[2]);;

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
		mexEvalString("drawnow;");
		std::uint16_t L;	

		std::vector< std::vector<double> > chunks;
		std::vector<double> widths;
		std::vector<double> chunksize;

		if (fexists(fname)) {
			mexPrintf("Loading...");
			mexEvalString("drawnow;");

			FILE *file = fopen(fname.c_str(), "rb");

			if (vers == 2) {
				fread(&R, sizeof(std::uint32_t), 1, file);
				for (size_t r = 0; r < R; r++) {
					fread(&S, sizeof(std::uint32_t), 1, file);
					for (size_t s = 0; s < S; s++) {
						
						if (chunks.size() == chunks.max_size() - 1) {
							break;
						}
						

						L = 0;
						fread(&L, sizeof(std::uint16_t), 1, file);

						std::uint16_t *data = new std::uint16_t[L]();
						fread(data, sizeof(std::uint16_t), L, file);

						if (L >= L_lt && L <= L_ut) {
							bool useable = true;
							for (size_t lftm = 1; lftm < L; lftm++) {
								double sim = std::min(data[lftm], data[lftm - 1]) / sqrt(data[lftm] * data[lftm - 1]);
								if (sim < 0.7071){
									useable = false;
								}
							}

							if (useable) {

								if (L > chunks.size()) {
									chunks.resize(L);
								}


								for (size_t win = 2; win < L; win++) {
									for (size_t i = 0; i < L - win; i+=win) {
										double sum = 0;
										for (size_t j = i; j < i+win; j++) {
											sum += (double)data[j];
										}
										double mean = sum / (double)win;

										mean = sum / (double)win;
										double dev = 0;
										for (size_t j = i; j < i + win; j++) {
											dev += std::pow((double)data[j] - mean, 2.);
										}
										dev /= (double)win;
										dev = sqrt(dev);
										chunks[win].push_back(dev);
									}
								}
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
					/*
					if (chunks.size() == chunks.max_size() - 1) {
						break;
					}
					*/

					L = 0;
					fread(&L, sizeof(std::uint16_t), 1, file);

					std::uint16_t *data = new std::uint16_t[L]();
					fread(data, sizeof(std::uint16_t), L, file);

					if (L >= L_lt && L <= L_ut) {
						bool useable = true;
						for (size_t lftm = 1; lftm < L; lftm++) {
							double sim = std::min(data[lftm], data[lftm - 1]) / sqrt(data[lftm] * data[lftm - 1]);
							if (sim < 0.7071){
								useable = false;
							}
						}

						if (useable) {

							if (L > chunks.size()) {
								chunks.resize(L);
							}


							for (size_t win = 2; win < L; win++) {
								double sum = 0;
								for (size_t i = 0; i < win; i++) {
									sum += (double)data[i];
								}
								double mean = sum / (double)win;

								double dev = 0;
								for (size_t i = 0; i < win; i++) {
									dev += std::pow((double)data[i] - mean, 2.);
								}
								dev /= (double)win;
								dev = sqrt(dev);
								chunks[win].push_back(dev);

								for (size_t win = 2; win < L; win++) {
									for (size_t i = 0; i < L - win; i += win) {
										double sum = 0;
										for (size_t j = i; j < i + win; j++) {
											sum += (double)data[j];
										}
										double mean = sum / (double)win;

										mean = sum / (double)win;
										double dev = 0;
										for (size_t j = i; j < i + win; j++) {
											dev += std::pow((double)data[j] - mean, 2.);
										}
										dev /= (double)win;
										dev = sqrt(dev);
										chunks[win].push_back(dev);
									}
								}
							}
						}
					}
					delete[] data;
				}
			}

			fclose(file);

			
			widths.resize(chunks.size());
			chunksize.resize(chunks.size());
			for (size_t c = 0; c < chunks.size(); c++) {
				double sum = 0;
				for (size_t i = 0; i < chunks[c].size(); i++) {
					sum += chunks[c][i];
				}

				double mean = 0.;
				if (chunks[c].size() > 0) {
					mean = sum / (double)chunks[c].size();
				}

				widths[c] = mean;
				chunksize[c] = c;
			}
			

			plhs[0] = mxCreateNumericMatrix(1, widths.size(), mxDOUBLE_CLASS, mxREAL);
			double* widthsout = (double*)mxGetData(plhs[0]);
			memcpy(widthsout, &widths[0], sizeof(double)*widths.size());

			plhs[1] = mxCreateNumericMatrix(1, chunksize.size(), mxDOUBLE_CLASS, mxREAL);
			double* chunksizeout = (double*)mxGetData(plhs[1]);
			memcpy(chunksizeout, &chunksize[0], sizeof(double)*chunksize.size());

			mexPrintf("Load Complete\n");
			mexEvalString("drawnow;");
		}
		else {
			plhs[0] = mxCreateLogicalScalar(0);
			mexPrintf("File Does Not Exists.\n");
			mexEvalString("drawnow;");
		}
	}
	else {
		mexPrintf("Wrong Number of Input Arguments.\n");
		mexEvalString("drawnow;");
		plhs[0] = mxCreateLogicalScalar(0);
	}
}