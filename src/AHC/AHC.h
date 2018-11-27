/* Standard agglomerative hierarchical clustering methods.
 * Assume input object number is N, then form a N*N distance matrix to be stored.
 * Basic procedures are
 * 	1. store distance matrix
 * 	2. min-heap to sort mutual distance
 * 	3. every time merge two nodes with smallest distance and update the min-heap
 * 	4. merge until only one cluster or given cluster number is obtained
 */

/* Performance and memory usage analysis
 * 1. Performance
 * 		N*N (distance matrix) + 2N*N logN (build min-heap) + 2N*N*log(N*N) (min-heap update)
 * 2. Memory usage
 * 		N*N (distance matrix) + N*N(min-heap)
 *
 * Would expect this algorithm to be super time-consuming and memory-wasting
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

class AHC
{

public:

/* default constructor */
	AHC();

/* argument constructor with argc and argv */
	AHC(const int& argc, char **argv);

/* destructor */
	~AHC();

/* perform clustering function */
	void performClustering();

private:
/* extract features from datasets as representative curves */
	void extractFeatures(const std::vector<int>& storage, const std::vector<std::vector<int> >& neighborVec,
            			 const Eigen::MatrixXf& centroid);

/* metric preparation object to be stored ahead of time */
	MetricPreparation object;

/* input norm option */
	int normOption;

/* the tag to tell whether it's a PBF or not */
	bool isPBF;

/* group information */
	std::vector<int> group;

/* activityList vector to store event */
	std::vector<string> activityList;

/* timeList vector to store time information */
	std::vector<string> timeList;

/* store dataset information */
	DataSet ds;

/* how many clusters to be needed */
	int numberOfClusters;

/* k-means initialization option */
	int initializationOption;

/* linkage choice */
	int linkageOption;

/* whether used L-method to detect optimal number of clusters */
	bool lMethod;

/* whether read cluster by input or read from txt */
	bool readCluster;

/* whether is pathline or not */
	bool isPathlines;

/* used to test the curve of some function */
	std::vector<float> curveValue[4];

/* set dataset from user command */
	void setDataset(const int& argc, char **argv);

/* compute distance between two clusters based on likage type */
	const float getDistAtNodes(const vector<int>& firstList, const vector<int>& secondList, const int& Linkage);

/* perform AHC merging by given a distance threshold */
	void hierarchicalMerging(std::unordered_map<int, Ensemble>& node_map, std::vector<DistNode>& dNodeVec,
			  std::vector<Ensemble>& nodeVec);

/* perform group-labeling information */
	void setLabel(const std::vector<Ensemble>& nodeVec, vector<vector<int> >& neighborVec,
			      vector<int>& storage, Eigen::MatrixXf& centroid);

/* get string for linkage type */
	string getLinkageStr();	

/* get entropy ratio */
	void getEntropyRatio(const std::vector<int>& storage, float& EntropyRatio);

/* get norm string */
	string getNormStr();

/* get entropy ratio string */ 
	string getEntropyStr(const float& EntropyRatio);	

/* set a vector for min-heap */
	void setValue_merge(std::vector<DistNode>& dNodeVec, std::unordered_map<int, Ensemble>& node_map);

/* set a vector for min-heap */
	void setValue(std::vector<DistNode>& dNodeVec, std::unordered_map<int, Ensemble>& node_map);

/* perform clustering on normOption */
	void performClustering_by_norm();

};

#endif
