#ifndef _METRIC_H
#define _METRIC_H

#include "PreComputing.h"

extern const int& BIN_SIZE;


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

			default:
				break;
		}
	}
};



#endif
