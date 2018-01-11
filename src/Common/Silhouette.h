#ifndef SILHOUETTE_H
#define SILHOUETTE_H

#include "Distance.h"

class Silhouette
{
public:
	/* used to store Silhouette value for all streamlines */
	std::vector<float> sData;
	/* used to store Silhouette value for all clusters */
	std::vector<float> sCluster;
	/* used to store Silhouette value for average datasets */
	float sAverage;


	Silhouette();

	~Silhouette();

	void computeValue(const int& normOption,
					  const MatrixXf& array, 
					  const int& Row, 
					  const int& Column,
					  const std::vector<int>& group, 
					  const MetricPreparation& object,
					  const int& groupNumber);

	void computeValue(const int& normOption,
					  const MatrixXf& array,
					  const int& Row,
					  const int& Column,
					  const std::vector<int>& group,
					  const MetricPreparation& object,
					  const int& groupNumber,
					  const std::vector<vector<int> >& storage);

	void reset();

private:
	const float getDist(const int& first,
						const int& second,
						const MetricPreparation& object,
						const MatrixXf& array,
						const int& normOption);

	const float getA_i(const std::vector<std::vector<int> >& storage,
					   const std::vector<int>& group,
					   const MatrixXf& array,
					   const int& index,
					   const MetricPreparation& object,
					   const int& normOption);

	const float getB_i(const std::vector<std::vector<int> >& storage,
					   const std::vector<int>& group,
					   const MatrixXf& array,
					   const int& index,
					   const MetricPreparation& object,
					   const int& normOption);

};

#endif
