/*
 * The class is to compute the Silhouette, Gamma statistics and DB index for the clustering results
 */


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

	/* DB index for evaluating clustering result */
	float dbIndex;

	/* Hubert's Gamma Statistic */
	float gammaStatistic = -1.0;


	Silhouette();

	~Silhouette();

	/* non-PBF based evaluation */
	void computeValue(const int& normOption,
					  const MatrixXf& array, 
					  const int& Row, 
					  const int& Column,
					  const std::vector<int>& group, 
					  const MetricPreparation& object,
					  const int& groupNumber,
					  const bool& isPBF);

	void computeValue(const int& normOption,
					  const MatrixXf& array,
					  const int& Row,
					  const int& Column,
					  const std::vector<int>& group,
					  const MetricPreparation& object,
					  const int& groupNumber,
					  const bool& isPBF,
					  const std::vector<vector<int> >& storage);

	/* this is for PCA silhouette computing where only Euclidean space in projected-space is considered */
	void computeValue(const Eigen::MatrixXf& cArray,
					  const std::vector<int>& group,
					  const int& groupNo,
					  const bool& isPBF);


	void reset();

private:

	// general norm option of get the distance between two lines
	const float getDist(const int& first,
						const int& second,
						const MetricPreparation& object,
						const MatrixXf& array,
						const int& normOption);

	// get the a_i with general norm option
	const float getA_i(const std::vector<std::vector<int> >& storage,
					   const std::vector<int>& group,
					   const MatrixXf& array,
					   const int& index,
					   const MetricPreparation& object,
					   const int& normOption);

	// get the b_i with general norm option
	const float getB_i(const std::vector<std::vector<int> >& storage,
					   const std::vector<int>& group,
					   const MatrixXf& array,
					   const int& index,
					   const MetricPreparation& object,
					   const int& normOption);

	/* this is for PCA silhouette computing where only Euclidean space in projected-space is considered */
	void getMatrixM(const Eigen::MatrixXf& cArray,
			  	  	const std::vector<int>& group,
					const std::vector<std::vector<int> >& storage,
					Eigen::MatrixXf& distM,
					Eigen::MatrixXf& idealDistM);


	/* this is for the other clustering algorithm for non-PBF datasets which could compuote Gamma statistics */
	void getMatrixM(const Eigen::MatrixXf& cArray,
			  	  	const std::vector<int>& group,
					const std::vector<std::vector<int> >& storage,
					Eigen::MatrixXf& idealDistM);


	/* get A_i in low-dimension projected-space */
	const float getA_i(const std::vector<std::vector<int> >& storage,
					   const std::vector<int>& group,
					   const Eigen::MatrixXf& array,
					   const int& index,
					   const bool& isPBF,
					   const Eigen::MatrixXf& distM);


	/* get B_i in low-dimension projected-space */
	const float getB_i(const std::vector<std::vector<int> >& storage,
					   const std::vector<int>& group,
					   const Eigen::MatrixXf& array,
					   const int& index,
					   const bool& isPBF,
					   const Eigen::MatrixXf& distM);


	/* this function is used to compute silhouette value for clustering */
	void computeSilhouette(const Eigen::MatrixXf& array,
						   const std::vector<int>& group,
						   const bool& isPBF,
						   const std::vector<std::vector<int> >& storage,
						   const Eigen::MatrixXf& distM);


	/* this function is used to compute silhouette value for clustering */
	void computeSilhouette(const Eigen::MatrixXf& array,
						   const std::vector<int>& group,
						   const std::vector<std::vector<int> >& storage,
						   const MetricPreparation& object,
						   const int& normOption);


	/* compute DB index */
	void computeDBIndex(const Eigen::MatrixXf& array,
						const std::vector<int>& group,
						const std::vector<std::vector<int> >& storage);


	/* compute DB index with normOption as input */
	void computeDBIndex(const Eigen::MatrixXf& array,
						const std::vector<int>& group,
						const std::vector<std::vector<int> >& storage,
						const MetricPreparation& object,
						const int& normOption);


	/* for computing gamma statistics for PCA only which requires two input */
	void computeGammaStatistic(const Eigen::MatrixXf& distM,
							   const Eigen::MatrixXf& idealDistM);


	/* for computing gamma statistics for non-PBF datasets */
	void computeGammaStatistic(const Eigen::MatrixXf& idealDistM);

};

#endif
