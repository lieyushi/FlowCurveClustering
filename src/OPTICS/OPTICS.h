#ifndef _OPTICS_H_
#define _OPTICS_H_

#include "Predefined.h"
#include "ValidityMeasurement.h"
#include <string>
#include <queue>

/* algorithm implemented based on https://en.wikipedia.org/wiki/OPTICS_algorithm */


class DensityClustering
{
public:
	DensityClustering(const int& argc,
					  char **argv);

	~DensityClustering();

	void performClustering();


private:
	/* final output vector or ordered list */
	vector<int> orderedList;

	vector<PointNode> nodeVec;
	MetricPreparation object;
	int normOption;
	DataSet ds;

	bool isPBF;

	void setDataset(const int& argc,
				    char **argv);	//get dataset from file

	void setNormOption();	//get normOption from user input

	void OPTICS(const float& radius_eps,
				const int& minPts);

	void update(const int& index,
				const vector<int>& neighbor,
				LinkedList& seeds,
				const float& radius_eps,
				const int& minPts);

	const vector<int> regionQuery(const int& index,
								  const float& radius_eps);

	void getDistRange(float& minDist, 
					  float& maxDist);

	const int setMinPts();

	const float setTimesMin(const float& minDist, 
					  		const float& maxDist);

	const float getReachability(const int& first,
								const int& target,
								const int& minPts);

	void extractFeatures(const float& radius_eps,
							   const int& minPts);

	void computeCoredDistance(const float& radius_eps,
							  const int& minPts);

	/* how to get group information based on reachability-plot */
	void getGroup(const float& radius_eps);

	void writeReachability();


	/* set the eps as averaged minPt-th dist */
	const float getMinPt_thDist(const int& minPts);
};

#endif
