/* @brief Standard agglomerative hierarchical clustering methods.
 * @details
 * 	Assume input object number is N, then form a N*N distance matrix to be stored.
 * 	Basic procedures are
 * 		1. store distance matrix
 * 		2. min-heap to sort mutual distance
 * 		3. every time merge two nodes with smallest distance and update the min-heap
 * 		4. merge until only one cluster or given cluster number is obtained
 * 	Performance and memory usage analysis
 * 		1. Performance
 * 			N*N (distance matrix) + 2N*N logN (build min-heap) + 2N*N*log(N*N) (min-heap update)
 * 		2. Memory usage
 * 			N*N (distance matrix) + N*N(min-heap)
 *	 	Would expect this algorithm to be super time-consuming and memory-wasting
 *
 * @author Lieyu Shi
 */

#ifndef _AHC_H_
#define _AHC_H_

#include <unordered_set>
#include <unordered_map>
#include <map>
#include <string>
#include <algorithm>

#include "Predefined.h"
#include "ValidityMeasurement.h"
#include "DetermClusterNum.h"


/*
 * @brief The class to perform agglomerative hierarchical clustering algorithms and related clustering analysis
 */
class AHC
{

public:


	/*
	 * @brief The default constructor
	 */
	AHC();


	/*
	 * @brief A constructor with parameters
	 * @details
	 * 	To set the data set and perform some reading operation into the member variables from the file
	 *
	 * @param[in] argc The count of argv
	 * @param[in] argv The argument string line
	 */
	AHC(const int& argc, char **argv);


	/*
	 * @brief The destructor of the class AHC
	 * @details
	 * 	Delete the global pointer distance matrix
	 *
	 */
	~AHC();


	/*
	 * @brief Perform the clustering for selected similarity measure labels w.r.t. user input and data set type
	 * @details
	 *	If the number of clusters are pre-stored in the "cluster_number" file, the code will read the numbers first.
	 *	Then for different similarity measures, it will decide whether the L-method is activated or not. If L-method
	 *	is activated, the number of clusters is set to be 1, otherwise it will be set as user input. Hierarchical merging
	 *	operation for the tree is called after parameter setting is finished.
	 */
	void performClustering();


private:


	/*
	 * @brief Extract the features and calculate the evaluation metrics for clustering results
	 * @details
	 * 	Based on the clustering result, the cluster representatives will be extracted first for each cluster based on the
	 * 	closest/furthest candidate to the cluster centroid.
	 * 	The clustering evaluation metrics will be computed for the quantitative analysis of the clustering result.
	 * 	All the information (cluster representatives stored in .vtk file, clustering evaluation metrics stored in readme)
	 * 	will be recorded and stored in the designated folders for further batch processing.
	 *
	 * @param[in] storage The size of each cluster as input as input
	 * @param[in] neighborVec The candidate vector for each cluster as input
	 * @param[in] centroid The centroid for each cluster as input
	 */
	void extractFeatures(const std::vector<int>& storage, const std::vector<std::vector<int> >& neighborVec,
            			 const Eigen::MatrixXf& centroid);


	/*
	 * @brief metric preparation object to be stored ahead of time
	 */
	MetricPreparation object;

	/*
	 * @brief input norm option
	 */
	int normOption;

	/*
	 * @brief the tag to tell whether it's a PBF or not
	 */
	bool isPBF;

	/*
	 * @brief group information for integral curves
	 */
	std::vector<int> group;

	/*
	 * @brief activityList vector to store event
	 */
	std::vector<string> activityList;

	/*
	 * @brief timeList vector to store time information
	 */
	std::vector<string> timeList;

	/*
	 * @brief store dataset information
	 */
	DataSet ds;

	/*
	 * @brief how many clusters to be needed
	 */
	int numberOfClusters;

	/*
	 * @brief k-means initialization option
	 */
	int initializationOption;

	/*
	 * @brief linkage choice
	 */
	int linkageOption;

	/*
	 * @brief whether used L-method to detect optimal number of clusters
	 */
	bool lMethod;

	/*
	 * @brief whether read cluster by input or read from txt
	 */
	bool readCluster;

	/*
	 * @brief whether is pathline or not
	 */
	bool isPathlines;

	/*
	 * @brief used to test the curve of some function
	 */
	std::vector<float> curveValue[4];


	/*
	 * @brief Set the data set and perform necessary operations with user parameter input
	 * @details
	 * 	The function will read in the coordinates of the streamlines/pathlines from the given argument.
	 * 	Then it will provide necessary sampling strategy based on user input and data set type
	 * 	Furthmore, parameter input will be enforced from the console.
	 *
	 * @param[in] argc The count of argument string
	 * @param[in] argv The char* array for argument
	 */
	void setDataset(const int& argc, char **argv);


	/*
	 * @brief Get the distance between two nodes with a given linkage type
	 *
	 * @param[in] firstList The first node with candidates
	 * @param[in] secondList The second node with candidates
	 * @param[in] Linkage The linkage type, 0 for single, 1 for complete and 2 for average
	 * @return The distance value between two nodes in AHC clustering
	 */
	const float getDistAtNodes(const vector<int>& firstList, const vector<int>& secondList, const int& Linkage);


	/*
	 * @brief Perform hierarchical merge for the trees by a given required cluster number
	 * @details
	 * 	Hiarachically merge the nodes until the number of cluster is reached. Then based on whether L-method is activated
	 * 	or not, the posterior operation will be called on either finding the clustering information or finding the optimal
	 * 	number of clusters
	 *
	 * @param[out] node_map The initial node with each node representing one streamline/pathline
	 * @param[out] dNodeVec The DistNode vector which has the indices of two nodes and their distance
	 * @param[out] nodeVec The vector of Ensemble which has candidate index of the cluster
	 */
	void hierarchicalMerging(std::unordered_map<int, Ensemble>& node_map, std::vector<DistNode>& dNodeVec,
			  std::vector<Ensemble>& nodeVec);


	/*
	 * @brief Set the labels and compute the centroid and cluster related information
	 * @details
	 *	With generated node information, the cluster size, cluster centroids and candidates belonging to the same cluster
	 *	will be determined for further clustering evaluation metric calculation.
	 *
	 * @param[in] nodeVec The remained node vector
	 * @param[out] neighborVec The candidate vector of each cluster to be updated
	 * @param[out] storage The size of each cluster to be updated
	 * @param[out] centroid The centroid coordinates to be updated
	 */
	void setLabel(const std::vector<Ensemble>& nodeVec, vector<vector<int> >& neighborVec,
			      vector<int>& storage, Eigen::MatrixXf& centroid);


	/*
	 * @brief Get the string for linkage type
	 * @return A string type for linkage
	 */
	string getLinkageStr();	


	/*
	 * @brief Calculate the normalized entropy
	 *
	 * @param[in] storage The size of different clusters
	 * @param[out] EntropyRatio The normalized entropy to be updated
	 */
	void getEntropyRatio(const std::vector<int>& storage, float& EntropyRatio);


	/*
	 * @brief Get the string type of input similarity measure
	 * @return A string type
	 */
	string getNormStr();


	/*
	 * @brief Get the string type for entropy value
	 *
	 * @param[out] EntropyRatio The normalized entropy value
	 * @return The string of the float value
	 */
	string getEntropyStr(const float& EntropyRatio);	


	/*
	 * @brief Set the merged nodes and perform necessary merge operations before the starting of AHC
	 *
	 * @param[out] dNodeVec The node vector to be updated
	 * @param[out] node_map The map to record the index and node
	 */
	void setValue_merge(std::vector<DistNode>& dNodeVec, std::unordered_map<int, Ensemble>& node_map);


	/*
	 * @brief Set value for the dNodeVec and node_map as initialization of the AHC procedure
	 *
	 * @param[out] dNodeVec The vector of nodes to be updated
	 * @param[out] node_map The map for recording index and candidate streamlines in the cluster
	 */
	void setValue(std::vector<DistNode>& dNodeVec, std::unordered_map<int, Ensemble>& node_map);


	/*
	 * @brief The function to perform AHC clustering by a given norm
	 *
	 * @details
	 * 	The function will first judge whether the local storage of distance matrix exists or not. If it exists, the program
	 * 	will read in the distance matrix from the local file, otherwise it will calculate the distance matrix. In this format
	 * 	the time for calculating distance matrix can be saved for different clustering techniques.
	 * 	Then it will read in number of clusters as input for the clustering operation.
	 */
	void performClustering_by_norm();

};

#endif
