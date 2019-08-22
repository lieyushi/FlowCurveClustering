/*
 * @brief The source cpp for implementing the clustering algorithm DBSCAN for flow visualization
 * @author Lieyu Shi
 */


#ifndef _DENSITYCLUSTERING_H
#define _DENSITYCLUSTERING_H


#include "IOHandler.h"
#include "Initialization.h"
#include "Silhouette.h"
#include "ValidityMeasurement.h"
#include <queue>


/*
 * @brief point type in DBSCAN
 */
enum PointType
{
	CORE = 0,
	BORDER,
	NOISE
};


/*
 * @brief point node struct, it has point type (PointType), visited (bool), group (cluser label)
 */
struct PointNode
{
	int type;
	bool visited;
	int group;
	PointNode():type(-1), visited(false), group(-1)
	{}

	~PointNode()
	{}
};


/*
 * @brief data set struct to record the coordinates, matrix representation, max elements, vertex count and dimension
 */
struct DataSet
{
	vector<vector<float> > dataVec;
	Eigen::MatrixXf dataMatrix;
	int maxElements;
	int vertexCount;
	int dimension;

	string strName;
	string fullName;

	DataSet()
	{}

	~DataSet()
	{}
};


/*
 * @brief DBSCAN clustering class object
 */
class DensityClustering
{
public:


	/*
	 * @brief DBSCAN constructor with parameters
	 *
	 * @param[in] argc The count of arguments
	 * @param[in] argv The char* array with data sets
	 */
	DensityClustering(const int& argc, char **argv);


	/*
	 * @brief Destructor
	 */
	~DensityClustering();


	/*
	 * @brief Perform DBSCAN clustering
	 */
	void performClustering();


private:

	/*
	 * @brief The point node vector for representating all the candidates
	 */
	vector<PointNode> nodeVec;

	/*
	 * @brief The MetricPreparation object for computing the distance matrix
	 */
	MetricPreparation object;

	/*
	 * @brief The norm option
	 */
	int normOption;

	/*
	 * @brief The data set object to store geometric information
	 */
	DataSet ds;

	/*
	 * @brief whether it is a PBF dataset or not
	 */
	bool isPBF;

	/*
	 * @brief Whether the data set is pathline or streamline
	 */
	bool isPathlines;


	/*
	 * @brief Set the data set from the arguments
	 *
	 * @param[in] argc Count of arguments
	 * @param[in] argv The argument string type
	 */
	void setDataset(const int& argc, char **argv);	//get dataset from file


	/*
	 * @brief Read the norm option from user input
	 */
	void setNormOption();


	/*
	 * @brief Perform the DBSCAN clustering with given parameters
	 *
	 * @param[in] radius_eps The radius to check the neighboring information
	 * @param[in] minPts The minPts for the DBSCAN clustering
	 */
	void DBSCAN(const float& radius_eps, const int& minPts);


	/*
	 * @brief Expand the cluster with candidates that lie within range of the target
	 *
	 * @param[in] index The index of candidate streamlines
	 * @param[in] neighbor The neighborhood candidates found
	 * @param[in] cluster_id The label for the clusters
	 * @param[in] radius_eps The radius for searching around the neighborhood
	 * @param[in] minPts The min number of points for the DBSCAN clustering
	 */
	void expandCluster(const int& index,
					   vector<int>& neighbor,
					   const int& cluster_id,
					   const float& radius_eps,
					   const int& minPts);


	/*
	 * @brief Perform the region-based query for the target streamline within a given radius
	 *
	 * @param[in] index The index of the target streamlines
	 * @param[in] radius_eps The radius of neighborhood checking
	 * @return A vector<int> object that contains the region candidates
	 */
	const vector<int> regionQuery(const int& index,
								  const float& radius_eps);


	/*
	 * @brief Calculate the minimal and maximal distance range for the user input of radius
	 *
	 * @param[in] minDist The minimal distance to be updated
	 * @param[in] maxDist The maximal distance to be updated
	 */
	void getDistRange(float& minDist, 
					  float& maxDist);


	/*
	 * @brief Set the minPts parameter, default value 6 is already enough for our paper
	 */
	const int setMinPts();


	/*
	 * @brief Select a ratio between [0,1] for the radius input
	 *
	 * @param[in] minDist The minimal distance value as input
	 * @param[in] maxDist The maximal distance value as input
	 * @return A radius for the DBSCAN clustering after the user input
	 */
	const float setTimesMin(const float& minDist, 
					  		const float& maxDist);


	/*
	 * @brief Extract features for the clustering results and calculate the clustering evaluation metrics
	 *
	 * @param[in] radius_eps The radius for DBSCAN clustering
	 * @param[in] minPts The minPts parameter for DBSCAN
	 */
	void extractFeatures(const float& radius_eps,
							   const int& minPts);

	/*
	 * @brief Get the distance threshold
	 *
	 * @param[in] minPts The min points for the DBSCAN
	 */
	const float getDistThreshold(const int& minPts);


	/*
	 * @brief Compute the minPts-th dist for all candidates
	 *
	 * @param[in] minPts The min point numbers for the DBSCAN clustering
	 */
	const float getAverageDist(const int& minPts);
};

#endif
