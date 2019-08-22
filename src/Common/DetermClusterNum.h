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
#include <vector>


using namespace std;
using namespace Eigen;


/*
 * @brief The class of hierarchical L-method to find the optimal number of clusters on the input of
 * number-merged_distance
 */
class DetermClusterNum {

public:
	/*
	 * @brief The default constructor
	 */
	DetermClusterNum();


	/*
	 * @brief Destructor
	 */
	virtual ~DetermClusterNum();


	/*
	 * @brief public accessor for get num of clusters
	 */
	const int& getFinalNumOfClusters()
	{
		return finalNumOfClusters;
	}


	/*
	 * @brief Use iterative refinement of knee to get optimal number for hierarchical clustering
	 * @param[out] eval_graph The map that contains the cluster number and its merged distance
	 */
	void iterativeRefinement(std::map<int, float>& eval_graph);


	/*
	 * @brief Record the L-method result in the local file
	 * @param[in] normOption The norm option as input
	 */
	void recordLMethodResult(const int& normOption);

private:

	/*
	 * @brief The final number of clusters as the optimal value from the L-method
	 */
	int finalNumOfClusters;


	/*
	 * @brief Find the knee value by the L-method with a given cutoff value
	 * @param[in] eval_graph The map with cluster numbers and their relative merged distance
	 * @param[in] cutoff The cutoff index point
	 * @return An index found by the L-method which is related to the knee
	 */
	const int LMethod(const std::map<int, float>& eval_graph, const int& cutoff);


	/*
	 * @brief Remove extremely dissimilarity mcluster merges
	 * @param[out] eval_graph The map including the cluster numbers and their merged distance
	 */
	void removeExtreme(std::map<int, float>& eval_graph);

};

#endif /* SRC_COMMON_DETERMCLUSTERNUM_H_ */
