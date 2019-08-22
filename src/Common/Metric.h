/*
 * @brief This is a class to prepare the functions for computing the distance for integral curves
 * @author Lieyu Shi
 */


#ifndef _METRIC_H
#define _METRIC_H

#include "PreComputing.h"
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>


/*
 * @brief the bin size for calculating the d_S (14) similarity measure
 */
extern const int& BIN_SIZE;

/*
 * @brief the number of bundles for calculating the entropy-based similarity (16)
 */
extern const int& BUNDLE_SIZE;


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
						 std::vector<float>& rotation);


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
						 std::vector<std::vector<float> >&rotationSequence);


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
					   std::vector<MultiVariate>& normalMultivariate);


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
					  std::vector<std::vector<float> >& rotationSequence);


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
					  		 std::vector<MultiVariate>& normalMultivariate);


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
					  std::vector<VectorXf >& unitLength);


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
					 std::vector<std::vector<float> >& pairwiseNorm);


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
					 std::vector<std::vector<float> >& pairwise);


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
					  std::vector<std::vector<float> >& pairwise);


/*
 * @brief The class for pre-calculating and storing preliminary variables for distance computation, e.g., entropy,
 * signature-based histograms, e.t.c..
 *
 * @details
 * 	To make it versatile and robust towards different similarity measures, we assign several member variables for
 * 	potential or possible assignment.
 */
struct MetricPreparation
{
	std::vector<float> rotation;
	std::vector<std::vector<float> > rotationSequence;
	std::vector<MultiVariate> normalMultivariate;
	std::vector<VectorXf > unitLength;
	std::vector<std::vector<float> > pairwise;
	std::vector<std::vector<float> > pairwiseNorm;

	int row, column;


	/*
	 * @brief The constructor to assign the size from the input
	 *
	 * @param[in] Row The row size
	 * @param[in] Column The column size
	 */
	MetricPreparation(const int& Row,
					  const int& Column)
	{
		/* don't need to allocate redundant memory since the dataset size would be huge */
		row = Row;
		column = Column;
	}

	/*
	 * @brief The default constructor
	 */
	MetricPreparation()
	{
		row = column = 0;
	}

	/*
	 * @brief The destructor
	 */
	~MetricPreparation()
	{
		row = column = -1;
	}

	/*
	 * @brief Pre-calculate and store some necessary attributes related to norm option
	 *
	 * @param[in] data The matrix coordinates
	 * @param[in] Row The row size
	 * @param[in] Column The column size
	 * @param[in] normOption The norm option
	 */
	void preprocessing(const Eigen::MatrixXf& data,
					   const int& Row,
					   const int& Column,
					   const int& normOption)
	{
		switch(normOption)
		{
			case 2:
			case 5:
			case 8:
				{
					pairwise = std::vector<std::vector<float> >(Row, std::vector<float>(Column-3));
					pairwiseNorm = std::vector<std::vector<float> >(Row, std::vector<float>((Column-3)/3));
					computePairWise(data, Row, Column-3, pairwise, pairwiseNorm);
				}
				break;

			case 4:
				{
					/*  if rotation used for judge similarity difference, has to use pre-defined cache */
					computeMeanRotation(data, Row, Column, rotation);
				}
				break;

			/*  pre-defined cache for sequence mean and standard deviation */
			case 3:
				{
					rotationSequence = std::vector<std::vector<float> >(Row,std::vector<float>(2));
					getRotationSequence(data, Row, Column, rotationSequence);
				}
				break;

			case 6:
				{
					normalMultivariate = std::vector<MultiVariate>(Row, MultiVariate());
					getNormalSequence(data, Row, Column, normalMultivariate);
				}
				break;

			case 7:
				{
					rotationSequence = std::vector<std::vector<float> >(Row,std::vector<float>(2));
					getFixedSequence(data, Row, Column, rotationSequence);
				}
				break;

			case 9:
				{
					normalMultivariate = std::vector<MultiVariate>(Row, MultiVariate());
					getUnnormalizedSequence(data, Row, Column, normalMultivariate);
				}
				break;

			case 10:
				{
					unitLength = std::vector<VectorXf >(Row, VectorXf(Column));
					getUnitDirection(data, Row, Column, unitLength);
				}
				break;

			case 14:
				{
					pairwise = std::vector<std::vector<float> >(Row, std::vector<float>(BIN_SIZE));
					getSignatureBin(data, Row, Column, pairwise);
				}
				break;

			case 16:
				{
					pairwise = std::vector<std::vector<float> >(Row, std::vector<float>(2));
					getBundleEntropy(data, Row, Column, pairwise);
				}

			default:
				break;
		}
	}
};


/*
 * @brief The template for computing numerical integration
 * @details
 * 	The implementation by boost C++ library can be seed at
 * 	https://www.boost.org/doc/libs/1_63_0/libs/multiprecision/doc/html/boost_multiprecision/tut/floats/fp_eg/gi.html
 *
 * @param[in] a The value type of a
 * @param[in] b The value type of b
 * @param[in] tol The step size
 * @param[in] func The function type
 */
template<typename value_type, typename function_type>
inline value_type integral(const value_type a,
                           const value_type b,
                           const value_type tol,
                           function_type func)
{
	unsigned n = 1U;

	value_type h = (b - a);
	value_type I = (func(a) + func(b)) * (h / 2);

	for(unsigned k = 0U; k < 8U; k++)
	{
		h /= 2;

		value_type sum(0);
		for(unsigned j = 1U; j <= n; j++)
		{
			sum += func(a + (value_type((j * 2) - 1) * h));
		}

		const value_type I0 = I;
		I = (I / 2) + (h * sum);

		const value_type ratio     = I0 / I;
		const value_type delta     = ratio - 1;
		const value_type delta_abs = ((delta < 0) ? -delta : delta);

		if((k > 1U) && (delta_abs < tol))
		{
			break;
		}

		n *= 2U;
	}

	return I;
}


/*
 * @brief The class for creating the integrated function
 */
template<typename value_type>
class cyl_bessel_j_integral_rep
{
public:
	cyl_bessel_j_integral_rep(const value_type& a, const value_type& b, const value_type& c) : a(a), b(b), c(c)
	{}

	value_type operator()(const value_type& t) const
	{
		// pi * Jn(x) = Int_0^pi [cos(x * sin(t) - n*t) dt]
		// return cos(x * sin(t) - (n * t));
		return sqrt(a+2.0*b*t+c*t*t);
	}

private:
	const value_type a, b, c;
};


/*
 * @brief Get the calculus of int_0^1(sqrt(a+2bt+ct^2))d_t
 *
 * @param[in] a The coefficient in the formula sqrt(a+2bt+ct^2)
 * @param[in] b The coefficient in the formula sqrt(a+2bt+ct^2)
 * @param[in] c The coefficient in the formula sqrt(a+2bt+ct^2)
 * @return A float value from the calculus computation
 */
const float get_calculus(const float& a, const float& b, const float& c);


#endif
