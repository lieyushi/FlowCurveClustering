#include "Metric.h"


// the bin size for calculating the d_S (14) similarity measure
const int& BIN_SIZE = 20;

// the number of bundles for calculating the entropy-based similarity (16)
const int& BUNDLE_SIZE = 20;


// calculate the mean rotation for all the streamlines
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


// calculate the roataion sequence for all the streamlines
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


// calculate the multivariate normal sequence for all the streamlines
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


// calculate values of segment to a fixed direction for all the streamlines
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


// get the unnormalized seuqnce of the multivariate segments for all the streamlines
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


// get the unit directions of segments for all the streamlines
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


// calculate the pair-wise distance values
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


/* get signature-based bins for histogram for dissimilarity computation */
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

