/*
 * @brief The class contains some preliminary functions to calculate before the distance matrix starts
 * @author Lieyu Shi
 */


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
#include <set>
#include <queue>
#include <cmath>

using namespace std;
using namespace Eigen;


/*
 * @brief The class to store the variate normal distribution, e.g., the covariance, meanVec for the 3-variate variables
 */
struct MultiVariate
{	
	Matrix3f covariance;
	Vector3f meanVec;
	MultiVariate()
	{}
	~MultiVariate()
	{}
};


/*
 * @brief The class that contains the index and the discrete curvature
 */
struct CurvatureObject
{
	float curvature;
	int index;

	CurvatureObject(const float& curvature, const int& i): curvature(curvature), index(i)
	{}

	CurvatureObject()
	{}
};


/*
 * @brief Define the compare function for calling the heap operation
 */
class CompareFunc
{
public:
	bool operator()(const CurvatureObject& first, const CurvatureObject& second)
	{
		return first.curvature < second.curvature;
	}
};


/*
 * @brief calculate the sequence values for a given streamline coordinate
 *
 * @param[in] array The coordinate of a streamline
 * @param[in] size The size of coordinate
 * @param[out] rowSequence The rowSequence vector
 */
void getSequence(const VectorXf& array, 
				 const int& size, 
				 std::vector<float>& rowSequence);


/*
 * @brief Calculate the average discrete curvatures on the streamline coordinate
 *
 * @param[in] array The coordinate of the streamline
 * @param[in] size The size of coordinate
 * @return The float value of the average rotation
 */
const float getRotation(const Eigen::VectorXf& array, 
						const int& size);


/*
 * @brief Get the normalized multivariate w.r.t. segments of the streamline
 *
 * @param[in] array The coordinate of the streamline
 * @param[in] size The size of coordinate
 * @param[out] rowSequence The row sequence to be updated
 */
void getNormalMultivariate(const VectorXf& array, 
				 	 	   const int& size, 
				 	 	   MultiVariate& rowSequence);


/*
 * @brief Get the fixed sequence for the streamline
 *
 * @param[in] array The coordinate of the streamline
 * @param[in] The size of coordinate
 * @param[out] rowSequence The row sequence to be updated
 */
void getEachFixedSequence(const VectorXf& array, 
				 		  const int& size, 
				 		  std::vector<float>& rowSequence);


/*
 * @brief Get the unnormalized multivariate w.r.t. segments of the streamline coordinates
 *
 * @param[in] array The coordinate of the streamline
 * @param[in] The size of coordinate
 * @param[out] rowSequence The row sequence to be updated
 */
void getUnnormalizedMultivariate(const VectorXf& array, 
				 	 	  		 const int& size, 
				 	 	  		 MultiVariate& rowSequence);


/*
 * @brief Get the unit direction for each streamline
 *
 * @param[in] array The coordinate of a streamline
 * @param[in] pointNum The size of point in the line
 * @param[out] direction The direction vector to be updated
 */
void getUnitDirection_byEach(const VectorXf& array, 
							 const int& pointNum, 
							 VectorXf& direction);


/*
 * @brief The difference compared to the former function is that, this will resample on high-curvature points
 *
 * @param[in] array The coordinate of a streamline
 * @param[in] binNum The number of points on the line
 * @param[out] histogram The vector values making up the histogram to be updated
 */
void getSignatureHist(const Eigen::VectorXf& array,
					  const int& binSize,
					  std::vector<float>& histogram);


/*
 * @brief Get the bin-based histogram for signature, the difference is that we should try to get maximal binNum-1 points
 *
 * @param[in] array The coordinate of the streamline
 * @param[in] binNum The number of bins on the histogram
 * @param[out] histogram The histogram of values as an output of vector
 */
void getSignatureHistSampled(const Eigen::VectorXf& array,
					  	  	 const int& binSize,
							 std::vector<float>& histogram);


/*
 * @brief Get the linear and angular entropy
 *
 * @param[in] array The coordinate of the streamline
 * @param[in] bundleSize The size of bundles for the attributes generated
 * @param[out] histogram The histogram of values as output
 */
void getLinearAngularEntropy(const Eigen::VectorXf& array,
 	  	 	 	 	 	 	 const int& bundleSize,
							 std::vector<float>& histogram);


/*
 * @brief Get the pairwise each attribute for the streamline
 *
 * @param[in] array The coordinate of a streamline
 * @param[in] size The size of the coordinate
 * @param[out] wiseVec The pairwise intersection angle of each segment
 * @param[out] wiseNorm The pairwise norm vector of each segment
 */
void getPairWise_byEach(const VectorXf& data,
						const int& size,
					 	std::vector<float>& wiseVec,
					 	std::vector<float>& wiseNorm);


/*
 * @brief Compute pseudo-inverse for trivial matrix for B-Metric clustering
 *
 * @param[in] a The input matrix
 * @param[in] epsilon The epsilon as the accuracy for the pseudo-inverse computation
 */
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
