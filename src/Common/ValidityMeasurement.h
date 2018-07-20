/*
 * ValidityMeasurement.h
 *
 *  Created on: Jun 24, 2018
 *      Author: lieyu
 */

#ifndef SRC_COMMON_VALIDITYMEASUREMENT_H_
#define SRC_COMMON_VALIDITYMEASUREMENT_H_

// A C++ implementation for the paper https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=4761242
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <assert.h>
#include <tuple>
#include "Distance.h"

using namespace std;

class ValidityMeasurement {
public:

	// f_c=h(DDc)*g(Sc)
	float f_c;

	ValidityMeasurement();
	virtual ~ValidityMeasurement();

	// function API for computing the validity measurement for general cases
	void computeValue(const int& normOption, const MatrixXf& array, const std::vector<int>& group,
			const MetricPreparation& object, const bool& isPBF);

	// function API for computing the validity measurement for PCA only
	void computeValue(const MatrixXf& array, const std::vector<int>& group);

private:
	// min and max of S_c
	float min_Sc, max_Sc;

	// get MST for each cluster given index and pair-wise distance for general cases
	void getMST_Parent_Node(std::tuple<float, float, float>& values,
							const std::vector<int>& clusterNode,
							const MetricPreparation& object,
							const int& normOption,
							const MatrixXf& array,
							const bool& isPBF);

	// get MST for each cluster given index and pair-wise distance, PCA case only
	void getMST_Parent_Node(std::tuple<float, float, float>& values,
				const std::vector<int>& clusterNode, const MatrixXf& array);

	// compute the Sc by input range value for general cases
	const float get_Sc_by_range(const bool& isPBF, const Eigen::MatrixXf& distM,
			                    const std::vector<int>& clusterNode, const float& rangeValue,
								const MetricPreparation& object, const int& normOption, const MatrixXf& array,
								int& index);

	// compute the Sc by input range value for PCA case only
	const float get_Sc_by_range(const Eigen::MatrixXf& distM, const std::vector<int>& clusterNode,
													 const float& rangeValue, const MatrixXf& array, int& index);
};

#endif /* SRC_COMMON_VALIDITYMEASUREMENT_H_ */
