/*
 * DetermClusterNum.h
 *
 *  Created on: Aug 22, 2018
 *      Author: lieyu
 *      This is the implementation for the paper, Determing the Number of Clusters/Segments in Hierarchical
 *      Clustering/Segmentation Algorithms
 *      link: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1374239
 */

#ifndef SRC_COMMON_DETERMCLUSTERNUM_H_
#define SRC_COMMON_DETERMCLUSTERNUM_H_

#include <eigen3/Eigen/Dense>
#include <float.h>
#include <iostream>
#include <fstream>
#include <map>

using namespace std;
using namespace Eigen;

class DetermClusterNum {
public:
	DetermClusterNum();
	virtual ~DetermClusterNum();

	/* public accessor for get num of clusters */
	const int& getFinalNumOfClusters()
	{
		return finalNumOfClusters;
	}

	/* use iterative refinement of knee to get optimal number for hierarchical clustering */
	void iterativeRefinement(std::map<int, float>& eval_graph);

	/* write the number in a file */
	void recordLMethodResult(const int& normOption);

private:
	int finalNumOfClusters;

	/* return a knee by L method given a cutoff */
	const int LMethod(const std::map<int, float>& eval_graph, const int& cutoff);

	/* remove extremely dissimilar cluster merges */
	void removeExtreme(std::map<int, float>& eval_graph);

};

#endif /* SRC_COMMON_DETERMCLUSTERNUM_H_ */
