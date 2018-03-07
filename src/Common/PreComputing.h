#ifndef PRECOMPUTING_H
#define PRECOMPUTING_H

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/SVD>
#include <climits>
#include <float.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>
#include <sys/time.h>
using namespace std;
using namespace Eigen;


struct MultiVariate
{	
	Matrix3f covariance;
	Vector3f meanVec;
	MultiVariate()
	{}
	~MultiVariate()
	{}
};

void getSequence(const VectorXf& array, 
				 const int& size, 
				 std::vector<float>& rowSequence);

const float getRotation(const Eigen::VectorXf& array, 
						const int& size);

void getNormalMultivariate(const VectorXf& array, 
				 	 	   const int& size, 
				 	 	   MultiVariate& rowSequence);

void getEachFixedSequence(const VectorXf& array, 
				 		  const int& size, 
				 		  std::vector<float>& rowSequence);

void getUnnormalizedMultivariate(const VectorXf& array, 
				 	 	  		 const int& size, 
				 	 	  		 MultiVariate& rowSequence);

void getUnitDirection_byEach(const VectorXf& array, 
							 const int& pointNum, 
							 VectorXf& direction);

/* get the bin-based histogram for signature */
void getSignatureHist(const Eigen::VectorXf& array,
					  const int& binSize,
					  std::vector<float>& histogram);


void getPairWise_byEach(const VectorXf& data,
						const int& size,
					 	std::vector<float>& wiseVec,
					 	std::vector<float>& wiseNorm);

/* compute pseudo-inverse for trivial matrix for B-Metric clustering */
template<typename _Matrix_Type_>
_Matrix_Type_ pseudoInverse(const _Matrix_Type_ &a, 
							double epsilon = std::numeric_limits<double>::epsilon())
{
    Eigen::JacobiSVD< _Matrix_Type_ > svd(a ,Eigen::ComputeThinU | Eigen::ComputeThinV);
    double tolerance = epsilon * std::max(a.cols(), a.rows()) 
    				  *svd.singularValues().array().abs()(0);
    return svd.matrixV() *  (svd.singularValues().array().abs() > tolerance).select
           (svd.singularValues().array().inverse(), 0).matrix().asDiagonal() 
           * svd.matrixU().adjoint();
}

#endif
