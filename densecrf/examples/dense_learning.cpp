/*
    Copyright (c) 2013, Philipp Krähenbühl
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
        * Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.
        * Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.
        * Neither the name of the Stanford University nor the
        names of its contributors may be used to endorse or promote products
        derived from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY Philipp Krähenbühl ''AS IS'' AND ANY
    EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL Philipp Krähenbühl BE LIABLE FOR ANY
    DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "densecrf.h"
#include "optimization.h"
#include <cstdio>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <vector>
#include "ppm.h"
#include "common.h"
using namespace std;
int class_num = 4; // The class number
string data_path = ""; // "The path of training data"


class CRFEnergy : public EnergyFunction {
protected:
	VectorXf initial_u_param_, initial_lbl_param_, initial_knl_param_;
	DenseCRF & crf_;
	const ObjectiveFunction& objective_;
	int NIT_;
	bool unary_, pairwise_, kernel_;
	float l2_norm_;
public:
	CRFEnergy(DenseCRF & crf, const ObjectiveFunction & objective, int NIT, bool unary = 1,
		bool pairwise = 1, bool kernel = 1) :crf_(crf), objective_(objective), NIT_(NIT), unary_(unary),
		pairwise_(pairwise), kernel_(kernel), l2_norm_(0.f) {
		initial_u_param_ = crf_.unaryParameters();
		initial_lbl_param_ = crf_.labelCompatibilityParameters();
		initial_knl_param_ = crf_.kernelParameters();
	}
	void setL2Norm(float norm) {
		l2_norm_ = norm;
	}
	virtual VectorXf initialValue() {
		VectorXf p(unary_*initial_u_param_.rows() + pairwise_*initial_lbl_param_.rows() + kernel_*initial_knl_param_.rows());
		p << (unary_ ? initial_u_param_ : VectorXf()), (pairwise_ ? initial_lbl_param_ : VectorXf()), (kernel_ ? initial_knl_param_ : VectorXf());
		return p;
	}
	virtual double gradient(const VectorXf & x, VectorXf & dx) {
		int p = 0;
		if (unary_) {
			crf_.setUnaryParameters(x.segment(p, initial_u_param_.rows()));
			p += initial_u_param_.rows();
		}
		if (pairwise_) {
			crf_.setLabelCompatibilityParameters(x.segment(p, initial_lbl_param_.rows()));
			p += initial_lbl_param_.rows();
		}
		if (kernel_)
			crf_.setKernelParameters(x.segment(p, initial_knl_param_.rows()));

		VectorXf du = 0 * initial_u_param_, dl = 0 * initial_u_param_, dk = 0 * initial_knl_param_;
		double r = crf_.gradient(NIT_, objective_, unary_ ? &du : NULL, pairwise_ ? &dl : NULL, kernel_ ? &dk : NULL);
		dx.resize(unary_*du.rows() + pairwise_*dl.rows() + kernel_*dk.rows());
		dx << -(unary_ ? du : VectorXf()), -(pairwise_ ? dl : VectorXf()), -(kernel_ ? dk : VectorXf());
		r = -r;
		if (l2_norm_ > 0) {
			dx += l2_norm_ * x;
			r += 0.5*l2_norm_ * (x.dot(x));
		}

		return r;
	}
};

class CRFEnergys : public EnergyFunction
{
public:
	vector<CRFEnergy*> pEnergys;
	~CRFEnergys()
	{
		cleanCRFEnerges();
	}
	void setCRFEnerges(vector<PointCloudCRF*>& crfs,
		vector<IntersectionOverUnion*>& objectives, int NIT,
		bool unary = 1, bool pairwise = 1, bool kernel = 1)
	{
		cleanCRFEnerges();
		int n = crfs.size();
		pEnergys.resize(n);
		for (int i = 0; i < n; ++i)
		{
			pEnergys[i] = new CRFEnergy(*crfs[i], *objectives[i],
				NIT, unary, pairwise, kernel);
		}
	}

	void cleanCRFEnerges()
	{
		for (int i = 0; i < pEnergys.size(); ++i)
		{
			delete pEnergys[i];
		}
	}

	void setL2Norm(float Norm)
	{
		for (int i = 0; i < pEnergys.size(); ++i)
		{
			pEnergys[i]->setL2Norm(Norm);
		}
	}
	virtual VectorXf initialValue()
	{
		return pEnergys[0]->initialValue();
	}

	virtual double gradient(const VectorXf & x, VectorXf & dx)
	{
		int n = pEnergys.size();
		vector<double> rs(n);
		vector<VectorXf> dxs(n);
		#pragma omp parallel for
		for (int i = 0; i < n; ++i)
		{
			rs[i] = pEnergys[i]->gradient(x, dxs[i]);
		}

		dx = dxs[0];
		double r = rs[0];
		for (int i = 1; i < n; ++i)
		{
			dx += dxs[i];
			r += rs[i];
		}

		dx /= (double)n;
		return r / double(n);
	}
};

void CRFLearning()
{
	vector<string> all_files;
	get_all_filenames(all_files, data_path + "\\*.points");

	int n = all_files.size();
	vector<PointCloudCRF*> pcCRFs(n);
	vector<IntersectionOverUnion*> objectives(n);
	cout << "Load data ...";
#pragma omp parallel for
	for (int i = 0; i < n; ++i)
	{
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

		int num = prob.cols();
		int ch = prob.rows();
		pcCRFs[i] = new PointCloudCRF(num, ch);
		pcCRFs[i]->setUnaryEnergy(prob);
		pcCRFs[i]->addPairwiseGaussian(4, 4, 4, pt.data(), new PottsCompatibility(3));
		pcCRFs[i]->addPairwiseBilateral(8, 8, 8, 0.3, 0.3, 0.3, pt.data(), normal.data(),
			new MatrixCompatibility(MatrixXf::Identity(ch, ch)));

		objectives[i] = new IntersectionOverUnion(label);
	}

	int NIT = 5;
	const bool verbose = true;

	MatrixXf learning_params(2, 3);
	// Optimize the CRF in 3 phases:
	//  * pairwise
	//  * Full CRF
	learning_params << 1, 1, 0, 1, 1, 1;
	cout << "Optimize ...";
	for (int i = 0; i < learning_params.rows(); i++)
	{
		cout << "Phase : " << i << endl;
		// Setup the energy
		CRFEnergys energy;
		energy.setCRFEnerges(pcCRFs, objectives, NIT,
			learning_params(i, 0), learning_params(i, 1), learning_params(i, 2));
		energy.setL2Norm(1e-3);

		// Minimize the energy
		VectorXf p = minimizeLBFGS(energy, 2, true);

		// Save the values
		for (int k = 0; k < pcCRFs.size(); ++k)
		{
			int id = 0;
			PointCloudCRF& crf = *pcCRFs[k];
			if (learning_params(i, 0))
			{
				crf.setUnaryParameters(p.segment(id, crf.unaryParameters().rows()));
				id += crf.unaryParameters().rows();
			}

			if (learning_params(i, 1))
			{
				crf.setLabelCompatibilityParameters(p.segment(id, crf.labelCompatibilityParameters().rows()));
				id += crf.labelCompatibilityParameters().rows();
			}

			if (learning_params(i, 2))
			{
				crf.setKernelParameters(p.segment(id, crf.kernelParameters().rows()));
			}
		}
	}
	// Return the parameters
	PointCloudCRF& crf = *pcCRFs[0];
	std::cout << "Unary parameters: " << crf.unaryParameters().transpose() << std::endl;
	std::cout << "Pairwise parameters: " << crf.labelCompatibilityParameters().transpose() << std::endl;
	std::cout << "Kernel parameters: " << crf.kernelParameters().transpose() << std::endl;

	// Do map inference
	for (int i = 0; i < n; ++i)
	{
		VectorXs map = pcCRFs[i]->map(NIT);
		string filename = all_files[i] + ".crf.txt";
		save_label(map, filename);
	}

	for (int i = 0; i < pcCRFs.size(); ++i)
	{
		delete pcCRFs[i];
		delete objectives[i];
	}
}

int main(int argc, char* argv[])
{

	if (argc < 3)
	{
		cout << "Usage: deanse_learning.exe <data_path> <class_num>";
		return 0;
	}

	data_path = argv[1];
	class_num = atoi(argv[2]);

	CRFLearning();

	return 0;
}
