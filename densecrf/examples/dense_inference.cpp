
#include <cstdio>
#include <cmath>
#include <iostream>
#include <io.h>
#include <fstream>
#include "densecrf.h"
#include "common.h"
using namespace std;
int class_num = 4; // The class number
string data_path = ""; // "The path of training data"
string parameter_path = "";

bool load_parameters(VectorXf& unary, VectorXf& pairwise,
	VectorXf& kernel, const string parameter_path)
{
	ifstream infile(parameter_path);
	if (!infile) return false;

	char str[512];
	vector<float> vec_unary, vec_pairwise, vec_kernel;
	while (infile.getline(str, 512))
	{
		char* pch;
		float tmp;
		pch = strtok(str, ": ");
		if (0 == strcmp(pch, "Unary"))
		{
			pch = strtok(NULL, ": "); // eat the word "parameters
			pch = strtok(NULL, ": ");
			while (pch != NULL)
			{
				tmp = atof(pch);
				vec_unary.push_back(tmp);
			}
		}
		else if (0 == strcmp(pch, "Pairwise"))
		{
			pch = strtok(NULL, ": ");
			pch = strtok(NULL, ": ");
			while (pch != NULL)
			{
				tmp = atof(pch);
				vec_pairwise.push_back(tmp);
				pch = strtok(NULL, ": ");

			}
		}
		else if (0 == strcmp(pch, "Kernel"))
		{
			pch = strtok(NULL, ": ");
			pch = strtok(NULL, ": ");
			while (pch != NULL)
			{
				tmp = atof(pch);
				vec_kernel.push_back(tmp);
				pch = strtok(NULL, ": ");
			}
		}
	}


	int n = vec_unary.size();
	unary.resize(n);
	memcpy(unary.data(), vec_unary.data(), n * sizeof(float));

	n = vec_pairwise.size();
	pairwise.resize(n);
	memcpy(pairwise.data(), vec_pairwise.data(), n * sizeof(float));

	n = vec_kernel.size();
	kernel.resize(n);
	memcpy(kernel.data(), vec_kernel.data(), n * sizeof(float));

	infile.close();
	return true;
}

void CRFInference()
{
	const int NIT = 5;

	vector<string> all_files;
	get_all_filenames(all_files, data_path + "\\*.points");

	VectorXf unary, pairwise, kernel;
	bool has_parameter = load_parameters(unary, pairwise, kernel, parameter_path);

	int n = all_files.size();
	cout << "Inference ... ";
#pragma omp parallel for
	for (int i = 0; i < n; ++i)
	{
		cout << "Processing: " << all_files[i] + "\n";
		VectorXs label;
		MatrixXf pt, normal, prob;
		load_pointcloud(pt, normal, label, all_files[i]);

		string filename_prob = all_files[i];
		filename_prob.replace(filename_prob.rfind('.') + 1, string::npos, "probo.txt");
		load_prob(prob, filename_prob);
		float* ptr = prob.data();
		for (int j = 0; j < prob.size(); ++j)
		{
			ptr[j] = -log(ptr[j] + 0.01f);
		}

		// build crf
		int num = prob.cols();
		int ch = prob.rows();
		PointCloudCRF pcCRF(num, ch);
		pcCRF.setUnaryEnergy(prob);
		pcCRF.addPairwiseGaussian(2, 2, 2, pt.data(), new PottsCompatibility(3));
		pcCRF.addPairwiseBilateral(4, 4, 4, 0.3, 0.3, 0.3, pt.data(), normal.data(),
			new MatrixCompatibility(MatrixXf::Identity(ch, ch)));

		// set paramters
		if (unary.size()) pcCRF.setUnaryParameters(unary);
		if (pairwise.size()) pcCRF.setLabelCompatibilityParameters(pairwise);
		if (kernel.size()) pcCRF.setKernelParameters(kernel);

		// inference
		VectorXs map = pcCRF.map(NIT);
		string filename = all_files[i] + ".crf.txt";
		save_label(map, filename);
	}

	cout << "Done!" << endl;
}

int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		cout << "Usage: deanse_learning.exe <data_path> <class_num> [parameter file]";
		return 0;
	}

	data_path = argv[1];
	class_num = atoi(argv[2]);
	if (argc > 3) parameter_path = argv[3];

	CRFInference();

	return 0;
}
