/*
 * @brief Agglomerative hierarchical clustering implementation with single thread.
 * @details
 * 	It's very hard to achieve multi-threaded speedup because it involves many
 * 	merge operation.
 * 	This implementation seems not to be a pure and standard agglomerative hierarchical
 * 	clustering method. Instead, it's more like a Birch  method which would merge all relavant
 * 	objects within a given threshold. Hence, it could not get any type of required cluster
 * 	number.
 * @author Lieyu Shi
 */


#ifndef _AHC_H_
#define _AHC_H_


#include "Predefined.h"
#include "ValidityMeasurement.h"
#include <unordered_set>
#include <map>
#include <string>


/*
 * @brief This is a class for AHC with distance threshold input and the hierarchical merging is continuously operating
 * as long as the merged distance is below the distance threshold
 */
class AHC
{

public:

	/*
	 * @brief The default constructor
	 */
	AHC();


	/*
	 * @brief The constructor with parameters
	 * @details
	 *	It firstly sets up the data set from the argument string and reads in from the local file
	 *	then create the distance matrix for AHC clustering
	 *
	 * @param[in] argc Count of argument
	 * @param[in] argv Char* array of argument
	 */
	AHC(const int& argc, char **argv);


	/*
	 * @brief The destructor for the class
	 */
	~AHC();


	/*
	 * @brief Perform AHC clustering with distance input
	 * @details
	 * 	It will select the type of AHC clustering, either by input of a cluster number, or a distance threshold. Then
	 * 	perform the hiararchical clustering by bottom-up merge of the tree. Then posterior calculation on feature
	 * 	extraction and clustering evaluation is performed for quantitative and visual analysis
	 *
	 */
	void performClustering();

private:

	/*
	 * @brief Extract the features and compute the evaluation metrics
	 *
	 * @param[in] storage size of clusters as input
	 * @param[in] neighborVec candidate vectors for each cluster
	 * @param[in] centroid The centroids for each cluster
	 */
	void extractFeatures(const std::vector<int>& storage, const std::vector<std::vector<int> >& neighborVec,
            			 const Eigen::MatrixXf& centroid);

	/*
	 * @brief Metric preparation object to be stored ahead of time for distance matrix computation
	 */
	MetricPreparation object;

	/*
	 * @brief whether the input dataset is PBF or not
	 */
	bool isPBF;

	/*
	 * @brief Input norm option
	 */
	int normOption;

	/*
	 * @brief Group information for different integral curves
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
	 * @brief distanc threshold
	 */
	float distanceThreshold;

	/*
	 * @brief store dataset information
	 */
	DataSet ds;

	/*
	 * @brief How many clusters to be needed
	 */
	int numberOfClusters;

	/*
	 * @brief expected cluster number as input
	 */
	int expectedClusters;

	/*
	 * @brief distance range recorded
	 */
	vector<float> distRange;		

	/*
	 * @brief linkage choice
	 */
	int linkageOption;

	/*
	 * @brief set dataset from user command
	 */
	void setDataset(const int& argc, char **argv);

	/*
	 * @brief set norm option, must be within 0-17
	 */
	void setNormOption();


	/*
	 * @brief Calculate the range of distance matrix
	 */
	void getDistRange();	


	/*
	 * @brief Get distance between nodes by linkage type
	 *
	 * @param[in] firstList The first node
	 * @param[in] secondList The second node
	 * @param[in] Linkage The linkage type
	 * @return A float value for the distance
	 */
	const float getDistAtNodes(const vector<int>& firstList, const vector<int>& secondList, const int& Linkage);


	/*
	 * @brief Perform hierarchical merging for AHC if the merged distance is lower than the distance threshold
	 *
	 * @param[out] nodeVec The Ensemble vector to be updated
	 */
	void hierarchicalMerging(std::vector<Ensemble>& nodeVec);


	/*
	 * @brief Set label for streamlines from the hierarchical tree
	 *
	 * @param[in] nodeVec The Ensemble vector as input
	 * @param[out] neighborVec The candidate vector for each cluster
	 * @param[out] storage The size of clusters to be updated
	 * @param[out] centroid The centroids of streamlines to be updated
	 */
	void setLabel(const std::vector<Ensemble>& nodeVec, vector<vector<int> >& neighborVec,
			      vector<int>& storage, Eigen::MatrixXf& centroid);


	/*
	 * @breif Perform bottom-up clustering with cluster number as input
	 * @details
	 * 	Iteratively merge the hierarchical tree until the tree size is approximately close to the required number
	 *
	 * @param[out] nodeVec A vector for Ensemble objects
	 */
	void bottomUp_byGroup(std::vector<Ensemble>& nodeVec);


	/*
	 * @brief Perform the bottom-up clustering by distance threshold of the linkage type
	 *
	 * @param[out] nodeVec The Ensemble vector to be updated
	 */
	void bottomUp_byThreshold(std::vector<Ensemble>& nodeVec);		


	/*
	 * @brief Get string for the linkage type
	 * @return String of linkage type
	 */
	string getLinkageStr();	


	/*
	 * @brief Get entropy ratio
	 * @param[in] storage Size of clusters
	 * @param[out] EntropyRatio The normalized entropy to be updated
	 */
	void getEntropyRatio(const std::vector<int>& storage, float& EntropyRatio);


	/*
	 * @brief Get string for norm
	 * @return String of norm
	 */
	string getNormStr();


	/*
	 * @brief Get the string for entropy value
	 * @param[in] EntropyRatio The entropy value
	 * @return String type
	 */
	string getEntropyStr(const float& EntropyRatio);	

};

#endif
