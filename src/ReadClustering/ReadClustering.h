/*
 * ReadClustering.h
 *
 *  Created on: Mar 13, 2018
 *      Author: lieyu
 */

#ifndef SRC_READCLUSTERING_READCLUSTERING_H_
#define SRC_READCLUSTERING_READCLUSTERING_H_

#include "IOHandler.h"
#include "Initialization.h"
#include "Silhouette.h"
#include <unordered_map>
#include "ValidityMeasurement.h"
#include <sstream>

using namespace std;

struct Dataset
{
	/* original coordinates */
	std::vector<std::vector<float> > dataVec;

	/* sampled array */
	Eigen::MatrixXf array;

	/* label information */
	unordered_map<string, std::vector<int> > groupAggregate;

	/* cluster number */
	unordered_map<string, int> maxGroup;

	/* number of elements inside */
	int numOfElements;

	std::vector<std::vector<int> > neighborVec;
};


class ReadClustering {
public:
	ReadClustering();
	virtual ~ReadClustering();

/* the public function called by main.cpp */
	void getEvaluation(const char* fileName);

private:

/* activityList vector to store event */
	std::vector<string> activityList;

/* timeList vector to store time information */
	std::vector<string> timeList;

/* Dataset object */
	Dataset ds;

/* judge whether it is a PBF or not */
	bool isPBF;

/* read data from file and store it into Dataset ds */
	void readData(const char* fileName);

/* compute four analysis evaluation measures */
	void computeEvaluation();

/* write those analysis framework */
	void writeAnalysis();

/* compute evaluation based on norm option */
	void computeEvaluation(std::unordered_map<string, std::vector<int> >::const_iterator& iter);

};

#endif /* SRC_READCLUSTERING_READCLUSTERING_H_ */
