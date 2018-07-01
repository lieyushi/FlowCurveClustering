/* Affinity propagation is a newly emerging clustering techniques based on "message passing" between data points
 * The wikipage could be referenced at https://en.wikipedia.org/wiki/Affinity_propagation
 * A sample C++ code could be seen at https://github.com/jincheng9/AffinityPropagation/blob/master/affinity_propagation.cpp,
 * or https://github.com/nojima/affinity-propagation-sparse/blob/master/ap.cpp
 */

#ifndef _AFFINITY_PROPAGATION_H_
#define _AFFINITY_PROPAGATION_H_

#include "Predefined.h"
#include "ValidityMeasurement.h"
#include <unordered_set>
#include <map>
#include <string>


#define LAMBDA 0.9


struct Para
{

	/* 1.directly filling with last vertex; 2. uniform sampling, 3. equal-arc sampling */
	int sampled;

	/* extraction option, 1. centroid, closest and furthest, 2. median, 3. statistical representation */
	int extractOption;

	/* max iteration for AP clustering */
	int maxIteration;
};



class AffinityPropagation
{

public:

/* default constructor */
	AffinityPropagation();

/* argument constructor with argc and argv */
	AffinityPropagation(const int& argc, char **argv, const Para& p, bool& automatic);

/* destructor */
	~AffinityPropagation();

/* perform clustering function */
	void performClustering();

private:

/**********************************************************************************************************
 **************************************   Private member variables   **************************************
 **********************************************************************************************************/

/* metric preparation object to be stored ahead of time */
	MetricPreparation object;

/* input norm option */
	int normOption = -1;

/* group information */
	std::vector<int> group;

/* activityList vector to store event */
	std::vector<string> activityList;

/* timeList vector to store time information */
	std::vector<string> timeList;

/* store dataset information */
	DataSet ds;

/* how many clusters to be needed */
	int numberOfClusters = -1;

/* extraction option, 1. centroid, closest and furthest, 2. median, 3. statistical representation */
	int extractOption = -1;

/* max iteration */
	int maxIteration = -1;

/* tell whether it is a tag for PBF dataset */
	bool isPBF;

/**********************************************************************************************************
 **************************************   Private member functions   **************************************
 **********************************************************************************************************/

/* extract features from datasets as representative curves */
	void extractFeatures(const std::vector<int>& storage, const std::vector<std::vector<int> >& neighborVec,
            			 const Eigen::MatrixXf& centroid);

/* set dataset from user command */
	void setDataset(const int& argc, char **argv);

/* set parameter */
	void getParameterUserInput();

/* set automatic parameter */
	void setParameterAutomatic(const Para& p);

/* run clustering based on different norm */
	void clusterByNorm(const int& norm);

/* perform group-labeling information */
	void setLabel(vector<vector<int> >& neighborVec, vector<int>& storage, Eigen::MatrixXf& centroid);

/* get entropy ratio */
	void getEntropyRatio(const std::vector<int>& storage, float& EntropyRatio);

/**********************************************************************************************************
 **************************************   Affinity Propagation Steps   *************************************
 **********************************************************************************************************/

/* get matrix S from distance matrix */
	void getMatrixS(Eigen::MatrixXf& matrixS);

/* initialize matrix S, R, A */
	void initializeMatrices(Eigen::MatrixXf& matrixS, Eigen::MatrixXf& matrixR, Eigen::MatrixXf& matrixA);

/* update responsibility matrix R */
	void updateResponsibility(Eigen::MatrixXf& matrixR, const Eigen::MatrixXf& matrixA,
							  const Eigen::MatrixXf& matrixS);

/* update availability matrix */
	void updateAvailability(Eigen::MatrixXf& matrixA, const Eigen::MatrixXf& matrixR);

/* get assignment by three matrices */
	void getGroupAssignment(const Eigen::MatrixXf& matrixR, const Eigen::MatrixXf& matrixA,
			  	  	  	  	const Eigen::MatrixXf& matrixS, std::vector<std::vector<int> >& neighborVec,
							std::vector<int>& storage);

};


#endif
