/*
 * @brief This is a class to prepare the functions for computing the distance for integral curves
 * @author Lieyu Shi
 */


#include "Metric.h"


/*
 * @brief the bin size for calculating the d_S (14) similarity measure
 */
const int& BIN_SIZE = 20;


/*
 * @brief the number of bundles for calculating the entropy-based similarity (16)
 */
const int& BUNDLE_SIZE = 20;


/*
 * @brief Calculate the mean rotation for all the streamlines
 *
 * @param[in] data The matrix coordinates of the streamlines
 * @param[in] Row The row size
 * @param[in] Column The column size
 * @param[out] rotation The rotation vector
 */
void computeMeanRotation(const Eigen::MatrixXf& data, 
						 const int& Row, 
						 const int& Column, 
						 std::vector<float>& rotation)
{
	rotation = std::vector<float>(Row, 0.0);
	const int& pointNum = Column/3-2;
#pragma omp parallel for schedule(static) num_threads(8)
	for (int i = 0; i < Row; ++i)
	{
		rotation[i] = getRotation(data.row(i), pointNum);
	}
}


/*
 * @brief Calculate the rotation sequence for all the streamlines
 *
 * @param[in] data The matrix coordinates
 * @param[in] Row The row size
 * @param[in] Column The column size
 * @param[out] rotationSequence The rotation sequence vector to be updated
 */
void getRotationSequence(const Eigen::MatrixXf& data, 
						 const int& Row, 
						 const int& Column, 
						 std::vector<std::vector<float> >& rotationSequence)
{
	const int& pointNum = Column/3-2;
#pragma omp parallel for schedule(static) num_threads(8)
	for (int i = 0; i < Row; ++i)
	{
		getSequence(data.row(i), pointNum, rotationSequence[i]);
	}
}


/*
 * @brief Calculate the multivariate normal sequence for all the streamlines
 *
 * @param[in] data The matrix coordinates
 * @param[in] Row The row size
 * @param[in] Column The column size
 * @param[out] normalMultivariate The multivariate normal sequence vector to be updated
 */
void getNormalSequence(const Eigen::MatrixXf& data, 
					   const int& Row, 
					   const int& Column, 
					   std::vector<MultiVariate>& normalMultivariate)
{
	const int& pointNum = Column/3-1;
#pragma omp parallel for schedule(static) num_threads(8)
	for (int i = 0; i < Row; ++i)
	{
		getNormalMultivariate(data.row(i), pointNum, normalMultivariate[i]);
	}
}


/*
 * @brief Calculate the values of segment to a fixed direction for all the streamlines
 *
 * @param[in] data The matrix coordinates
 * @param[in] Row The row size
 * @param[in] Column The column size
 * @param[out] rotationSequence The rotation sequence vector to be updated
 */
void getFixedSequence(const Eigen::MatrixXf& data, 
					  const int& Row, 
					  const int& Column, 
					  std::vector<std::vector<float> >& rotationSequence)
{
	const int& pointNum = Column/3-1;
#pragma omp parallel for schedule(static) num_threads(8)
	for (int i = 0; i < Row; ++i)
	{
		getEachFixedSequence(data.row(i), pointNum, rotationSequence[i]);
	}
}


/*
 * @brief Get the unnormalized seuqnce of the multivariate segments for all the streamlines
 *
 * @param[in] data The matrix coordinates of the streamlines
 * @param[in] Row The row size
 * @param[in] Column The column size
 * @param[out] normalMultivariate The multivariate normal sequence vector to be updated
 */
void getUnnormalizedSequence(const Eigen::MatrixXf& data, 
					   		 const int& Row, 
					  		 const int& Column, 
					  		 std::vector<MultiVariate>& normalMultivariate)
{
	const int& pointNum = Column/3-1;
#pragma omp parallel for schedule(static) num_threads(8)
	for (int i = 0; i < Row; ++i)
	{
		getUnnormalizedMultivariate(data.row(i), pointNum, normalMultivariate[i]);
	}
}


/*
 * @brief Get the unit directions of segments for all the streamlines
 *
 * @param[in] data The matrix coordinates of streamlines
 * @param[in] Row The row size
 * @param[in] Column The column size
 * @param[out] unitLength The unit directions of segments for all the streamlines
 */
void getUnitDirection(const Eigen::MatrixXf& data, 
					  const int& Row, 
					  const int& Column, 
					  std::vector<VectorXf>& unitLength)
{
	const int& pointNum = Column/3;
#pragma omp parallel for schedule(static) num_threads(8)
	for (int i = 0; i < Row; ++i)
	{
		getUnitDirection_byEach(data.row(i), pointNum, unitLength[i]);
	}
}


/*
 * @brief Calculate the pari-wise distance values
 *
 * @param[in] data The matrix coordinates of streamlines
 * @param[in] Row The row size
 * @param[in] Column The column size
 * @param[out] pairwise The pairwise vector to be updated
 * @param[out] pairwiseNorm The pairwise norm vector to be updated
 */
void computePairWise(const Eigen::MatrixXf& data, 
					 const int& Row, 
					 const int& Column, 
					 std::vector<std::vector<float> >& pairwise,
					 std::vector<std::vector<float> >& pairwiseNorm)
{
	const int& pointNum = Column/3; // how many line segments
#pragma omp parallel for schedule(static) num_threads(8)
	for (int i = 0; i < Row; ++i)
	{
		getPairWise_byEach(data.row(i), pointNum, pairwise[i], pairwiseNorm[i]);
	}
}


/*
 * @brief Get signature-based bins for histogram for similarity computation
 *
 * @param[in] data The matrix coordinates of the streamlines
 * @param[in] Row The row size
 * @param[in] Column The column size
 * @param[out] pairwise The pairwise vector to be updated
 */
void getSignatureBin(const Eigen::MatrixXf& data,
					 const int& Row,
					 const int& Column,
					 std::vector<std::vector<float> >& pairwise)
{
#pragma omp parallel for schedule(static) num_threads(8)
	for (int i=0;i<Row;++i)
	{
		//getSignatureHist(data.row(i),BIN_SIZE,pairwise[i]);
		getSignatureHistSampled(data.row(i),BIN_SIZE,pairwise[i]);
	}
}


/*
 * @brief Get streamline linear entropy and angular entropy in http://vis.cs.ucdavis.edu/papers/pg2011paper.pdf
 *
 * @param[in] data The matrix coordinates of streamlines
 * @param[in] Row The row size
 * @param[in] Column The column size
 * @param[out] pairwise The pairwise vector to be updated that has linear entropy and angular entropy values
 */
void getBundleEntropy(const Eigen::MatrixXf& data,
		 	 	 	  const int& Row,
					  const int& Column,
					  std::vector<std::vector<float> >& pairwise)
{
#pragma omp parallel for schedule(static) num_threads(8)
	for (int i=0;i<Row;++i)
	{
		getLinearAngularEntropy(data.row(i),BUNDLE_SIZE,pairwise[i]);
	}
}


/*
 * @brief Get the calculus of int_0^1(sqrt(a+2bt+ct^2))d_t
 *
 * @param[in] a The coefficient in the formula sqrt(a+2bt+ct^2)
 * @param[in] b The coefficient in the formula sqrt(a+2bt+ct^2)
 * @param[in] c The coefficient in the formula sqrt(a+2bt+ct^2)
 * @return A float value from the calculus computation
 */
const float get_calculus(const float& a, const float& b, const float& c)
{
	typedef boost::multiprecision::cpp_dec_float_50 mp_type;
	float result = integral(0.0F, 1.0F, 0.00001F, cyl_bessel_j_integral_rep<float>(a,b,c));
	return result;
}

