#include "Metric.h"

const int& BIN_SIZE = 25;

const int& BUNDLE_SIZE = 20;

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


void getRotationSequence(const Eigen::MatrixXf& data, 
						 const int& Row, 
						 const int& Column, 
						 std::vector<std::vector<float> >&rotationSequence)
{
	const int& pointNum = Column/3-2;
#pragma omp parallel for schedule(static) num_threads(8)
	for (int i = 0; i < Row; ++i)
	{
		getSequence(data.row(i), pointNum, rotationSequence[i]);
	}
}


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

void getFixedSequence(const Eigen::MatrixXf& data, 
					  const int& Row, 
					  const int& Column, 
					  std::vector<std::vector<float> >&rotationSequence)
{
	const int& pointNum = Column/3-1;
#pragma omp parallel for schedule(static) num_threads(8)
	for (int i = 0; i < Row; ++i)
	{
		getEachFixedSequence(data.row(i), pointNum, rotationSequence[i]);
	}
}


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


/* get signature-based bin for histogram for dissimilarity computation */
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


/* get streamline linear entropy and angular entropy in http://vis.cs.ucdavis.edu/papers/pg2011paper.pdf */
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


/* get the calculus of int_0^1(sqrt(a+2bt+ct^2))d_t */
const float get_calculus(const float& a, const float& b, const float& c)
{
	typedef boost::multiprecision::cpp_dec_float_50 mp_type;
	float result = integral(0.0F, 1.0F, 0.00001F, cyl_bessel_j_integral_rep<float>(a,b,c));
	return result;
}

