#ifndef _METRIC_H
#define _METRIC_H

#include "PreComputing.h"


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

struct MetricPreparation
{
	std::vector<float> rotation;
	std::vector<std::vector<float> > rotationSequence;
	std::vector<MultiVariate> normalMultivariate;
	std::vector<VectorXf > unitLength;
	std::vector<std::vector<float> > pairwise;
	std::vector<std::vector<float> > pairwiseNorm;

	MetricPreparation(const int& Row,
					  const int& Column)
	{
		rotationSequence = std::vector<std::vector<float> >(Row,std::vector<float>(2));
		normalMultivariate = std::vector<MultiVariate>(Row, MultiVariate());
		unitLength = std::vector<VectorXf >(Row, VectorXf(Column));
		pairwise = std::vector<std::vector<float> >(Row, std::vector<float>(Column-3));
		pairwiseNorm = std::vector<std::vector<float> >(Row, std::vector<float>((Column-3)/3));
	}

	MetricPreparation()
	{}

	~MetricPreparation()
	{}

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
				computePairWise(data, Row, Column-3, pairwise, pairwiseNorm);
				break;

			case 4:
			/*  if rotation used for judge similarity difference, has to use pre-defined cache */
				computeMeanRotation(data, Row, Column, rotation);
				break;

			/*  pre-defined cache for sequence mean and standard deviation */
			case 3:
				getRotationSequence(data, Row, Column, rotationSequence);
				break;

			case 6:
				getNormalSequence(data, Row, Column, normalMultivariate);
				break;

			case 7:
				getFixedSequence(data, Row, Column, rotationSequence);
				break;

			case 9:
				getUnnormalizedSequence(data, Row, Column, normalMultivariate);
				break;

			case 10:		
				getUnitDirection(data, Row, Column, unitLength);
				break;

			default:
				break;
		}
	}
};



#endif
