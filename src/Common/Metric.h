#ifndef _METRIC_H
#define _METRIC_H

#include "PreComputing.h"
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>


// the number of histograms for calculating the signature-based similarity measure (14)
extern const int& BIN_SIZE;

// the number of entropy segments for calculating entropy-based similarity measure (16)
extern const int& BUNDLE_SIZE;


void computeMeanRotation(const Eigen::MatrixXf& data, 
						 const int& Row, 
						 const int& Column, 
						 std::vector<float>& rotation);


void getRotationSequence(const Eigen::MatrixXf& data, 
						 const int& Row, 
						 const int& Column, 
						 std::vector<std::vector<float> >&rotationSequence);


void getNormalSequence(const Eigen::MatrixXf& data, 
					   const int& Row, 
					   const int& Column, 
					   std::vector<MultiVariate>& normalMultivariate);


void getFixedSequence(const Eigen::MatrixXf& data, 
					  const int& Row, 
					  const int& Column, 
					  std::vector<std::vector<float> >& rotationSequence);


void getUnnormalizedSequence(const Eigen::MatrixXf& data, 
					   		 const int& Row, 
					  		 const int& Column, 
					  		 std::vector<MultiVariate>& normalMultivariate);


void getUnitDirection(const Eigen::MatrixXf& data, 
					  const int& Row, 
					  const int& Column, 
					  std::vector<VectorXf >& unitLength);


void computePairWise(const Eigen::MatrixXf& data, 
					 const int& Row, 
					 const int& Column, 
					 std::vector<std::vector<float> >& pairwise,
					 std::vector<std::vector<float> >& pairwiseNorm);

/* get signature-based bin for histogram for dissimilarity computation */
void getSignatureBin(const Eigen::MatrixXf& data,
					 const int& Row,
					 const int& Column,
					 std::vector<std::vector<float> >& pairwise);


/* get streamline linear entropy and angular entropy in http://vis.cs.ucdavis.edu/papers/pg2011paper.pdf */
void getBundleEntropy(const Eigen::MatrixXf& data,
		 	 	 	  const int& Row,
					  const int& Column,
					  std::vector<std::vector<float> >& pairwise);


struct MetricPreparation
{
	std::vector<float> rotation;
	std::vector<std::vector<float> > rotationSequence;
	std::vector<MultiVariate> normalMultivariate;
	std::vector<VectorXf > unitLength;
	std::vector<std::vector<float> > pairwise;
	std::vector<std::vector<float> > pairwiseNorm;

	int row, column;

	MetricPreparation(const int& Row,
					  const int& Column)
	{
		/* don't need to allocate redundant memory since the dataset size would be huge */
//		rotationSequence = std::vector<std::vector<float> >(Row,std::vector<float>(2));
//		normalMultivariate = std::vector<MultiVariate>(Row, MultiVariate());
//		unitLength = std::vector<VectorXf >(Row, VectorXf(Column));
//		pairwise = std::vector<std::vector<float> >(Row, std::vector<float>(Column-3));
//		pairwiseNorm = std::vector<std::vector<float> >(Row, std::vector<float>((Column-3)/3));

		row = Row;
		column = Column;
	}

	MetricPreparation()
	{
		row = column = 0;
	}

	~MetricPreparation()
	{
		row = column = -1;
	}

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


/* the template for computing numerical integration from
 * https://www.boost.org/doc/libs/1_63_0/libs/multiprecision/doc/html/boost_multiprecision/tut/floats/fp_eg/gi.html
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

/* get the calculus of int_0^1(sqrt(a+2bt+ct^2))d_t */
const float get_calculus(const float& a, const float& b, const float& c);


#endif
