/* Agglomerative hierarchical clustering implementation with single thread.
 * It's very hard to achieve multi-threaded speedup because it involves many
 * merge operation.
 */

/* This implementation seems not to be a pure and standard agglomerative hierarchical
 * clustering method. Instead, it's more like a Birch  method which would merge all relavant
 * objects within a given threshold. Hence, it could not get any type of required cluster
 * number.
 */

#ifndef _AHC_H_
#define _AHC_H_

#include "Predefined.h"
#include "ValidityMeasurement.h"
#include <unordered_set>
#include <map>
#include <string>

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

/* whether the input dataset is PBF or not */
	bool isPBF;

/* input norm option */
	int normOption;

/* group information */
	std::vector<int> group;

/* activityList vector to store event */
	std::vector<string> activityList;

/* timeList vector to store time information */
	std::vector<string> timeList;

/* distanc threshold */
	float distanceThreshold;

/* store dataset information */
	DataSet ds;

/* how many clusters to be needed */
	int numberOfClusters;

/* expected cluster number as input */
	int expectedClusters;

/* distance range recorded */
	vector<float> distRange;		

/* k-means initialization option */
	int initializationOption;

/* linkage choice */
	int linkageOption;

/* set dataset from user command */
	void setDataset(const int& argc, char **argv);

/* set norm option, must be within 0-12 */
	void setNormOption();

/* set threshold for AHC function */
	void getDistRange();	

/* get distance between two treeNode nodes */
	const float getDistAtNodes(const vector<int>& firstList, const vector<int>& secondList, const int& Linkage);

/* perform AHC merging by given a distance threshold */
	void hierarchicalMerging(std::vector<Ensemble>& nodeVec);

/* perform group-labeling information */
	void setLabel(const std::vector<Ensemble>& nodeVec, vector<vector<int> >& neighborVec,
			      vector<int>& storage, Eigen::MatrixXf& centroid);

/* perform hierarchical clustering by given a group */
	void bottomUp_byGroup(std::vector<Ensemble>& nodeVec);

/* perform hierarchical clustering by given a threshold */
	void bottomUp_byThreshold(std::vector<Ensemble>& nodeVec);		

/* get string for linkage type */
	string getLinkageStr();	

/* get entropy ratio */
	void getEntropyRatio(const std::vector<int>& storage, float& EntropyRatio);

/* get norm string */
	string getNormStr();

/* get entropy ratio string */ 
	string getEntropyStr(const float& EntropyRatio);	

};

#endif
