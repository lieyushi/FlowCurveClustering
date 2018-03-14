#include "Metric.h"

const int& BIN_SIZE = 15;

void computeMeanRotation(const Eigen::MatrixXf& data, 
						 const int& Row, 
						 const int& Column, 
						 std::vector<float>& rotation)
{
	rotation = std::vector<float>(Row, 0.0);
	const int& pointNum = Column/3-2;
#pragma omp parallel for schedule(dynamic) num_threads(8)
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
#pragma omp parallel for schedule(dynamic) num_threads(8)
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
#pragma omp parallel for schedule(dynamic) num_threads(8)
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
#pragma omp parallel for schedule(dynamic) num_threads(8)
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
#pragma omp parallel for schedule(dynamic) num_threads(8)
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
#pragma omp parallel for schedule(dynamic) num_threads(8)
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
#pragma omp parallel for schedule(dynamic) num_threads(8)
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
#pragma omp parallel for schedule(dynamic) num_threads(8)
	for (int i=0;i<Row;++i)
	{
		//getSignatureHist(data.row(i),BIN_SIZE,pairwise[i]);
		getSignatureHistSampled(data.row(i),BIN_SIZE,pairwise[i]);
	}
}

