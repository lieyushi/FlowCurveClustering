/*
 * evrot.h
 *
 *  https://github.com/pthimon/clustering
 *  Class to compute the gradient of the eigenvectors
 *  alignment quality
 *
 *  Lihi Zelnik (Caltech) March.2005
 *
 *  Created on: 02-Mar-2009
 *      Author: sbutler
 *
 */

#ifndef EVROT_H_
#define EVROT_H_

#include <iostream>
#include <vector>
#include <math.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <exception>

using namespace Eigen;
using namespace std;

#define EPS 2.2204e-8

class Evrot {

public:
	Evrot(const Eigen::MatrixXf& X, int method);
	virtual ~Evrot();
	float getQuality() { return mQuality; }
	std::vector<std::vector<int> > getClusters() { return mClusters; }
	Eigen::MatrixXf& getRotatedEigenVectors() { return mXrot; }

protected:
	void evrot();
	void cluster_assign();
	float evqual(const Eigen::MatrixXf& X);
	float evqualitygrad(const Eigen::VectorXf& theta, const int& angle_index);
	Eigen::MatrixXf rotate_givens(const Eigen::VectorXf& theta);
	Eigen::MatrixXf build_Uab(const Eigen::VectorXf& theta, const int& a, const int& b);
	Eigen::MatrixXf gradU(const Eigen::VectorXf& theta, const int& k);

	//constants

	int mMethod;
	const int mNumDims;
	const int mNumData;
	int mNumAngles;
	Eigen::VectorXi ik;
	Eigen::VectorXi jk;

	//current iteration
	Eigen::MatrixXf mX;
	Eigen::MatrixXf mXrot;
	float mQuality;

	std::vector<std::vector<int> > mClusters;
};

#endif /* EVROT_H_ */
