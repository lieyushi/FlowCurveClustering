/* Agglomerative hierarchical clustering implementation with single thread.
 * It's very hard to achieve multi-threaded speedup because it involves many
 * merge operation.
 */


#ifndef _AHC_H_
#define _AHC_H_

#include "Predefined.h"
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

/* k-means initialization option */
	int initializationOption;

/* linkage choice */
	int linkageOption;

/* set dataset from user command */
	void setDataset(const int& argc, char **argv);

/* set norm option, must be within 0-12 */
	void setNormOption();

/* set threshold for AHC function */
	void setThreshold();

/* get distance between two treeNode nodes */
	const float getDistAtNodes(const vector<int>& firstList, const vector<int>& secondList, const int& Linkage);

};

#endif
