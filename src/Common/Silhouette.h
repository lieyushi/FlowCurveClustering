/*
 * @brief This class is to calculate the silhouette, Gamma statistics and DB index for clustering evaluation
 * @author Lieyu Shi
 */


#ifndef SILHOUETTE_H
#define SILHOUETTE_H

#include "Distance.h"


/*
 * @brief The class to compute the clustering evaluation metrics, silhouette, DB-index and the Gamma statistics
 */
class Silhouette
{
public:

	/*
	 * @brief used to store Silhouette value for all streamlines
	 */
	std::vector<float> sData;

	/*
	 * @brief used to store Silhouette value for all clusters
	 */
	std::vector<float> sCluster;

	/*
	 * @brief used to store Silhouette value for average data sets
	 */
	float sAverage;

	/*
	 * @brief DB index for evaluating clustering result
	 */
	float dbIndex;

	/*
	 * @brief Hubert's Gamma Statistic
	 */
	float gammaStatistic = -1.0;


	/*
	 * @brief Default constructor to create storage for the silhouette value
	 */
	Silhouette();


	/*
	 * @brief Destructor to empty the vector memory
	 */
	~Silhouette();


	/*
	 * @brief Compute the silhouette value with given input information
	 *
	 * @param[in] normOption The norm type
	 * @param[in] array The coordinate matrix of the streamlines
	 * @param[in] Row The size of the row
	 * @param[in] Column The column size
	 * @param[in] group The labels for all the streamlines
	 * @param[in] object A MetricPerparation object for distance matrix computation
	 * @param[in] groupNumber how many clusters are as input
	 * @param[in] isPBF A bool variable to tell whether the coordinates are from PBF fluid simulation or not
	 */
	void computeValue(const int& normOption,
					  const MatrixXf& array, 
					  const int& Row, 
					  const int& Column,
					  const std::vector<int>& group, 
					  const MetricPreparation& object,
					  const int& groupNumber,
					  const bool& isPBF);


	/*
	 * @brief Compute the three clustering evaluation metrics for general norm input
	 *
	 * @param[in] normOption The norm option
	 * @param[in] array The matrix coordinates of the streamlines
	 * @param[in] Row The row size
	 * @param[in] Column The column size
	 * @param[in] group The labels for all the streamlines
	 * @param[in] object The MetricPreparation object for distance matrix computation
	 * @param[in] groupNumber number of clusters as input
	 * @param[in] isPBF The bool to say whether it is from PBF or not
	 * @param[in] storage The candidates for each cluster
	 */
	void computeValue(const int& normOption,
					  const MatrixXf& array,
					  const int& Row,
					  const int& Column,
					  const std::vector<int>& group,
					  const MetricPreparation& object,
					  const int& groupNumber,
					  const bool& isPBF,
					  const std::vector<vector<int> >& storage);


	/*
	 * @brief Compute the silhouette value for PCA based input, so it will be calculated using Euclidean distance
	 *
	 * @param[in] array The matrix coordinates of the dimensionality reduced space
	 * @param[in] group The labels for different streamlines
	 * @param[in] groupNumber The number of clusters that are formed
	 * @param[in] isPBF The bool variable to tell whether it's PBF fluid simulation or not
	 */
	void computeValue(const Eigen::MatrixXf& cArray,
					  const std::vector<int>& group,
					  const int& groupNo,
					  const bool& isPBF);


	/*
	 * @brief Reset the member variables in the silhouette class
	 */
	void reset();

private:

	/*
	 * @brief Compute the similarity distance between two streamlines for a general norm option
	 *
	 * @param[in] first The index of first streamline
	 * @param[in] second The index of second streamline
	 * @param[in] object The MetricPreparation class object
	 * @param[in] array The matrix coordinates
	 * @param[in] normOption The norm option
	 */
	const float getDist(const int& first,
						const int& second,
						const MetricPreparation& object,
						const MatrixXf& array,
						const int& normOption);


	/*
	 * @brief Compute the A_i for silhouette value calculation for a general norm option
	 *
	 * @param[in] storage The candidates included for all the clusters
	 * @param[in] group The cluster labels of streamlines
	 * @param[in] array The coordinates of the streamlines
	 * @param[in] index The target streamline
	 * @param[in] object The MetricPreparation object
	 * @param[in] normOption The norm option
	 */
	const float getA_i(const std::vector<std::vector<int> >& storage,
					   const std::vector<int>& group,
					   const MatrixXf& array,
					   const int& index,
					   const MetricPreparation& object,
					   const int& normOption);


	/*
	 * @brief Compute the B_i for silhouette value calculation
	 *
	 * @param[in] storage The candidates included for all the clusters
	 * @param[in] group The cluster labels of streamlines
	 * @param[in] array The coordinates of the streamlines
	 * @param[in] index The target streamline
	 * @param[in] object The MetricPreparation object
	 * @param[in] normOption The norm option
	 */
	const float getB_i(const std::vector<std::vector<int> >& storage,
					   const std::vector<int>& group,
					   const MatrixXf& array,
					   const int& index,
					   const MetricPreparation& object,
					   const int& normOption);


	/*
	 * @brief Compute the distM (distance matrix) and idealDistM (ideal distance matrix) for PCA only
	 *
	 * @param[in] cArray The matrix coordinates of the streamlines
	 * @param[in] group The labels for all the streamlines
	 * @param[in] storage The candidates inside all the clusters
	 * @param[out] distM The distance matrix to be assigned the value
	 * @param[out] idealDist The ideal distance matrix with 0 and 1
	 */
	void getMatrixM(const Eigen::MatrixXf& cArray,
			  	  	const std::vector<int>& group,
					const std::vector<std::vector<int> >& storage,
					Eigen::MatrixXf& distM,
					Eigen::MatrixXf& idealDistM);


	/*
	 * @brief Compute the distM idealDistM (ideal distance matrix) for non-PBF case and distM already been stored
	 *
	 * @param[in] cArray The matrix coordinates of the streamlines
	 * @param[in] group The labels for all the streamlines
	 * @param[in] storage The candidates inside all the clusters
	 * @param[out] idealDist The ideal distance matrix with 0 and 1
	 */
	void getMatrixM(const Eigen::MatrixXf& cArray,
			  	  	const std::vector<int>& group,
					const std::vector<std::vector<int> >& storage,
					Eigen::MatrixXf& idealDistM);


	/*
	 * @brief Compute the A_i for silhouette value calculation for PCA case (use Euclidean distance)
	 *
	 * @param[in] storage The candidates included for all the clusters
	 * @param[in] group The cluster labels of streamlines
	 * @param[in] array The coordinates of the streamlines
	 * @param[in] index The target streamline
	 * @param[in] isPBF A bool tag whether the data is from PBF or not
	 * @param[in] distM The calculated distance matrix that has been calculated before
	 */
	const float getA_i(const std::vector<std::vector<int> >& storage,
					   const std::vector<int>& group,
					   const Eigen::MatrixXf& array,
					   const int& index,
					   const bool& isPBF,
					   const Eigen::MatrixXf& distM);


	/*
	 * @brief Compute the B_i for silhouette value calculation for PCA case (use Euclidean distance)
	 *
	 * @param[in] storage The candidates included for all the clusters
	 * @param[in] group The cluster labels of streamlines
	 * @param[in] array The coordinates of the streamlines
	 * @param[in] index The target streamline
	 * @param[in] isPBF A bool tag whether the data is from PBF or not
	 * @param[in] distM The calculated distance matrix that has been calculated before
	 */
	const float getB_i(const std::vector<std::vector<int> >& storage,
					   const std::vector<int>& group,
					   const Eigen::MatrixXf& array,
					   const int& index,
					   const bool& isPBF,
					   const Eigen::MatrixXf& distM);


	/*
	 * @brief Compute silhouette value for PCA-based only (no norm option)
	 *
	 * @param[in] array The matrix coordinates of streamlines
	 * @param[in] group The labels of all the streamlines
	 * @param[in] isPBF The bool tag whether is PBF or not
	 * @param[in] storage The candidates in all the clusters
	 * @param[in] distM The distance matrix as input
	 */
	void computeSilhouette(const Eigen::MatrixXf& array,
						   const std::vector<int>& group,
						   const bool& isPBF,
						   const std::vector<std::vector<int> >& storage,
						   const Eigen::MatrixXf& distM);


	/*
	 * @brief Compute silhouette value for general norm input
	 *
	 * @param[in] array The matrix coordinates of streamlines
	 * @param[in] group The labels of all the streamlines
	 * @param[in] storage The candidates in all the clusters
	 * @param[in] object The MetricPreparation class object
	 * @param[in] normOption The norm option
	 */
	void computeSilhouette(const Eigen::MatrixXf& array,
						   const std::vector<int>& group,
						   const std::vector<std::vector<int> >& storage,
						   const MetricPreparation& object,
						   const int& normOption);


	/*
	 * @brief Compute DB index for PCA case (only Euclidean distance used)
	 *
	 * @param[in] array The matrix coordinates of streamlines
	 * @param[in] group The labels of all the streamlines
	 * @param[in] storage The candidates in all the clusters
	 */
	void computeDBIndex(const Eigen::MatrixXf& array,
						const std::vector<int>& group,
						const std::vector<std::vector<int> >& storage);


	/*
	 * @brief Compute DB index for general input of norm option
	 *
	 * @param[in] array The matrix coordinates of streamlines
	 * @param[in] group The labels of all the streamlines
	 * @param[in] storage The candidates in all the clusters
	 * @param[in] object The MetricPreparation class object
	 * @param[in] normOption The norm option
	 */
	void computeDBIndex(const Eigen::MatrixXf& array,
						const std::vector<int>& group,
						const std::vector<std::vector<int> >& storage,
						const MetricPreparation& object,
						const int& normOption);


	/*
	 * @brief Compute the Gamma statistics between two matrices
	 *
	 * @param[in] distM The distance matrix
	 * @param[in] idealDistM The ideal distance matrix with only value 1 and 0
	 */
	void computeGammaStatistic(const Eigen::MatrixXf& distM,
							   const Eigen::MatrixXf& idealDistM);


	/*
	 * @brief Compute the Gamma statistics for the general input of norm options
	 *
	 * @param[in] idealDistM The ideal distance matrix with 1 and 0 only
	 */
	void computeGammaStatistic(const Eigen::MatrixXf& idealDistM);

};

#endif
