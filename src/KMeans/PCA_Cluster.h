/*
 * @brief The class to perform PCA-based clustering and k-means clustering on the input streamlines/pathlines
 * @author Lieyu Shi
 */


#ifndef _PCA_CLUSTER_H
#define _PCA_CLUSTER_H

#include "IOHandler.h"
#include "Initialization.h"
#include "Silhouette.h"
#include "ValidityMeasurement.h"
#include "Predefined.h"


/*
 * @brief Build a struct to record the old and new index of the cluster when sorted by the cluster size
 */
struct Ensemble
{
	int number;
	int newIndex;
	int oldIndex;
	Ensemble(const int& number, const int& index):number(number),newIndex(-1), oldIndex(index)
	{}
	bool operator<(const Ensemble& object) const
	{
		return number<object.number;
	}
};


/*
 * @brief Record event list and time list
 */
struct TimeRecorder
{
	std::vector<string> eventList;
	std::vector<string> timeList;
};


/*
 * @brief The class to perform PCA and k-means clustering algorithm
 */
class PCA_Cluster
{
public: 

	/*
	 * @brief Perform PCA-based clustering with input data
	 * @details
	 * 	It will first perform the PCA clustering, then perform either AHC-average or k-means clustering on the
	 * 	dimensionality reduced space
	 *
	 * @param[in] data The matrix coordinates
	 * @param[in] Row The row size
	 * @param[in] Column The column size
	 * @param[out] massCenter The centroid coordinates of the clusters
	 * @param[out] group The labels of all streamlines
	 * @param[out] totalNum The size of different clusters
	 * @param[out] closest The closest streamline representatives of all the clusters
	 * @param[out] furthest The furthest streamline representatives of the clusters
	 * @param[out] tr The TimeRecorder object to record the time
	 * @param[out] sil The Silhouette class for the clustering evaluation
	 */
	static void performPCA_Clustering(const Eigen::MatrixXf& data, 
									  const int& Row, 
									  const int& Column, 
									  std::vector<MeanLine>& massCenter,
		 							  std::vector<int>& group, 
		 							  std::vector<int>& totalNum, 
		 							  std::vector<ExtractedLine>& closest,
		 							  std::vector<ExtractedLine>& furthest,
									  TimeRecorder& tr,
									  Silhouette& sil);


	/*
	 * @brief Perform PCA-based clustering with input data and number of clusters as input
	 * @details
	 * 	It will first perform the PCA clustering, then perform either AHC-average or k-means clustering on the
	 * 	dimensionality reduced space
	 *
	 * @param[in] data The matrix coordinates
	 * @param[in] Row The row size
	 * @param[in] Column The column size
	 * @param[out] massCenter The centroid coordinates of the clusters
	 * @param[out] group The labels of all streamlines
	 * @param[out] totalNum The size of different clusters
	 * @param[out] closest The closest streamline representatives of all the clusters
	 * @param[out] furthest The furthest streamline representatives of the clusters
	 * @param[in] Cluster Number of clusters as input
	 * @param[out] tr The TimeRecorder object to record the time
	 * @param[out] sil The Silhouette class for the clustering evaluation
	 */
	static void performPCA_Clustering(const Eigen::MatrixXf& data, 
									  const int& Row, 
									  const int& Column, 
									  std::vector<MeanLine>& massCenter, 
									  std::vector<int>& group, 
									  std::vector<int>& totalNum, 
									  std::vector<ExtractedLine>& closest, 
									  std::vector<ExtractedLine>& furthest, 
									  const int& Cluster,
									  TimeRecorder& tr,
									  Silhouette& sil);


	/*
	 * @brief Perform k-means clustering with input data
	 *
	 * @param[in] data The matrix coordinates
	 * @param[in] Row The row size
	 * @param[in] Column The column size
	 * @param[out] massCenter The centroid coordinates of the clusters
	 * @param[out] group The labels of all streamlines
	 * @param[out] totalNum The size of different clusters
	 * @param[out] closest The closest streamline representatives of all the clusters
	 * @param[out] furthest The furthest streamline representatives of the clusters
	 * @param[in] normOption The norm option as input
	 * @param[out] tr The TimeRecorder object to record the time
	 * @param[out] sil The Silhouette class for the clustering evaluation
	 */
	static void performDirectK_Means(const Eigen::MatrixXf& data, 
									 const int& Row, 
									 const int& Column, 
									 std::vector<MeanLine>& massCenter,
									 std::vector<int>& group, 
									 std::vector<int>& totalNum, 
									 std::vector<ExtractedLine>& closest,
									 std::vector<ExtractedLine>& furthest, 
									 const int& normOption,
									 TimeRecorder& tr,
									 Silhouette& sil);


	/*
	 * @brief Perform k-means clustering with input data and the number of clusters as input
	 *
	 * @param[in] data The matrix coordinates
	 * @param[in] Row The row size
	 * @param[in] Column The column size
	 * @param[out] massCenter The centroid coordinates of the clusters
	 * @param[out] group The labels of all streamlines
	 * @param[out] totalNum The size of different clusters
	 * @param[out] closest The closest streamline representatives of all the clusters
	 * @param[out] furthest The furthest streamline representatives of the clusters
	 * @param[in] Cluster The number of clusters as input
	 * @param[in] normOption The norm option as input
	 * @param[out] tr The TimeRecorder object to record the time
	 * @param[out] sil The Silhouette class for the clustering evaluation
	 */
	static void performDirectK_Means(const Eigen::MatrixXf& data, 
									 const int& Row, 
									 const int& Column, 
									 std::vector<MeanLine>& massCenter,
									 std::vector<int>& group, 
									 std::vector<int>& totalNum, 
									 std::vector<ExtractedLine>& closest, 
									 std::vector<ExtractedLine>& furthest, 
									 const int& Cluster, 
									 const int& normOption,
									 TimeRecorder& tr,
									 Silhouette& sil);

private:


	/*
	 * @brief Perform the SVD for the input of matrix coordinates
	 * @details
	 * 	After SVD decomposition, it will select number of dimensions that can add up to 99.9% of the total variance,
	 * 	which usually results in 3 or 4 dimensions.
	 *
	 * @param[out] cArray The reduced-dimension matrix of the coordinates
	 * @param[in] data The matrix coordinates of the streamlines
	 * @param[in] Row The row size
	 * @param[in] Column The column size
	 * @param[out] PC_Number The number of PCs as output
	 * @param[out] SingVec The singular vectors
	 * @param[out] meanTrajectory The mean coordinate of the trajectory
	 * @param[out] tr The TimeRecorder class object
	 */
	static void performSVD(MatrixXf& cArray, 
						   const Eigen::MatrixXf& data, 
						   const int& Row, 
						   const int& Column, 
						   int& PC_Number, 
						   MatrixXf& SingVec,
	                       VectorXf& meanTrajectory,
						   TimeRecorder& tr);


	/*
	 * @brief Perform the k-means clustering algorithm on the PCs
	 *
	 * @param[in] cArray The dimension-reduced matrix coordinates
	 * @param[in] Row The row size
	 * @param[in] Column The column size
	 * @param[in] PC_Number The number of PCs
	 * @param[in] SingVec The singular vector matrix
	 * @param[in] meanTrajectory The mean coordinate of the trajectory
	 * @param[in] Cluster The number of clusters as input
	 * @param[out] group The labels for all the streamlines
	 * @param[out] totalNum The size of all the clusters
	 * @param[out] closest The coordinates of the closest extracted lines
	 * @param[out] furthest The coordinates of the furthest extracted lines
	 * @param[in] data The matrix coordinates of all the streamlines
	 * @param[out] tr The TimeRecorder class object
	 * @param[out] sil The Silhouette class object
	 */
	static void performPC_KMeans(const MatrixXf& cArray, 
								 const int& Row,
								 const int& Column, 
								 const int& PC_Number, 
								 const MatrixXf& SingVec, 
								 const VectorXf& meanTrajectory, 
								 std::vector<MeanLine>& massCenter, 
								 const int& Cluster, 
								 std::vector<int>& group,
								 std::vector<int>& totalNum, 
								 std::vector<ExtractedLine>& closest, 
								 std::vector<ExtractedLine>& furthest, 
								 const Eigen::MatrixXf& data,
								 TimeRecorder& tr,
								 Silhouette& sil);


	/*
	 * @brief Perform the k-means directly on similarity measures
	 *
	 * @param[in] data The matrix coordinates of the streamlines
	 * @param[in] Row The row size
	 * @param[in] Column The column size
	 * @param[out] massCenter The centroid coordinates of the clusters
	 * @param[in] Cluster The number of clusters
	 * @param[in] totalNum The size of clusters
	 * @param[out] closest The closest extracted lines of the clusters
	 * @param[out] furthest The furthest extracted lines of the clusters
	 * @param[in] normOption The norm option
	 * @param[out] tr The TimeRecorder object
	 * @param[out] sil The Silhouette object
	 */
	static void performFullK_MeansByClusters(const Eigen::MatrixXf& data, 
											 const int& Row, 
											 const int& Column, 
											 std::vector<MeanLine>& massCenter,
						   					 const int& Cluster, 
						   					 std::vector<int>& group, 
						   					 std::vector<int>& totalNum, 
						   					 std::vector<ExtractedLine>& closest, 
						   					 std::vector<ExtractedLine>& furthest, 
						   					 const int& normOption,
											 TimeRecorder& tr,
											 Silhouette& sil);


	/*
	 * @brief Perform the AHC-average on the dimensionality-reduced space coordinates
	 *
	 * @param[in] cArray The dimension-reduced matrix coordinates
	 * @param[in] PC_Number The number of PCs
	 * @param[in] SingVec The singular vector matrix
	 * @param[in] meanTrajectory The mean coordinate of the trajectory
	 * @param[out] massCenter The centroid coordinates of the clusters
	 * @param[in] Cluster The number of clusters as input
	 * @param[out] group The labels for all the streamlines
	 * @param[out] totalNum The size of all the clusters
	 * @param[out] closest The coordinates of the closest extracted lines
	 * @param[out] furthest The coordinates of the furthest extracted lines
	 * @param[in] data The matrix coordinates of all the streamlines
	 * @param[out] tr The TimeRecorder class object
	 * @param[out] sil The Silhouette class object
	 */
	static void perform_AHC(const MatrixXf& cArray,
							const int& PC_Number,
							const MatrixXf& SingVec,
							const VectorXf& meanTrajectory,
							std::vector<MeanLine>& massCenter,
							const int& Cluster,
							std::vector<int>& group,
							std::vector<int>& totalNum,
							std::vector<ExtractedLine>& closest,
							std::vector<ExtractedLine>& furthest,
							const Eigen::MatrixXf& data,
							TimeRecorder& tr,
							Silhouette& sil);

	/*
	 * @brief Perform AHC merging by given an input number of clusters
	 *
	 * @param[out] nodeMap The hash map for nodes
	 * @param[out] dNodeVec The node vector for nodes
	 * @param[out] nodeVec The vector of AHC hierarchical clustering node
	 * @param[in] reduced_dist_matrix The distance matrix of the dimensionality reduced space coordinates
	 * @param[in] cArray The dimensionality reduce coordinates
	 * @param[in] numberOfClusters The number of clusters as input
	 * @param[out] tr The TimeRecorder object
	 */
	static void hierarchicalMerging(std::unordered_map<int, AHC_node>& nodeMap, std::vector<DistNode>& dNodeVec,
			std::vector<AHC_node>& nodeVec, const Eigen::MatrixXf& reduced_dist_matrix, const Eigen::MatrixXf& cArray,
			const int& numberOfClusters, TimeRecorder& tr);


	/*
	 * @brief Get the distance between two nodes by a given linkage type
	 *
	 * @param[in] firstList The first node that contains the candidates
	 * @param[in] secondList The second node that contains the candidates
	 * @param[in] reduced_dist_matrix The distance matrix of the dimensionality reduced coordinates
	 * @return The float value between two nodes
	 */
	static float getDistAtNodes(const vector<int>& firstList, const vector<int>& secondList,
			const Eigen::MatrixXf& reduced_dist_matrix);


	/*
	 * @brief Set the nodes and perform necessary merges for nodes before the start of AHC clustering
	 *
	 * @param[out] dNodeVec The vector of nodes with distance
	 * @param[in] reduced_data The matrix coordinates of the dimensionality reduce coordinates
	 * @param[in] reduced_dist_matrix The distance matrix of the dimensionality reduced coordinates
	 */
	static void setValue(std::vector<DistNode>& dNodeVec, const Eigen::MatrixXf& reduced_data,
					     const Eigen::MatrixXf& reduced_dist_matrix);


	/*
	 * @brief Set the labels for streamlines from the clustering results
	 *
	 * @param[in] nodeVec The vector of AHC nodes for the AHC clustering results
	 * @param[out] neighborVec The candidates that belongs to the clusters
	 * @param[out] storage The size of clusters
	 * @param[out] centroid The centroid streamline coordinates of all the clusters
	 * @param[in] cArray The coordinates of the dimensionality reduced space
	 * @param[out] recorder The recorder vector
	 */
	static void setLabel(const std::vector<AHC_node>& nodeVec, vector<vector<int> >& neighborVec,
			vector<int>& storage, Eigen::MatrixXf& centroid, const Eigen::MatrixXf& cArray, std::vector<int>& recorder);

};

#endif

