/*
 * ValidityMeasurement.h
 *
 *  Created on: Jun 24, 2018
 *      Author: lieyu
 */

#ifndef SRC_COMMON_VALIDITYMEASUREMENT_H_
#define SRC_COMMON_VALIDITYMEASUREMENT_H_

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <assert.h>
#include <tuple>
#include "Distance.h"

using namespace std;


/*
 * @brief The C++ implementation for the paper https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=4761242
 */
class ValidityMeasurement {
public:

	/*
	 * @brief f_c=h(DDc)*g(Sc)
	 */
	float f_c;


	/*
	 * @brief Default constructor
	 */
	ValidityMeasurement();


	/*
	 * @brief Destructor
	 */
	virtual ~ValidityMeasurement();


	/*
	 * @brief Compute the validity measurement for the clustering result of general norm option
	 *
	 * @param[in] normOption The norm option
	 * @param[in] array The matrix coordinates of streamlines
	 * @param[in] group The labels of different streamlines
	 * @param[in] object The MetricPreparation class object
	 * @param[in] isPBF The bool tag to tell whether it is a PBF data set or not
	 */
	void computeValue(const int& normOption, const MatrixXf& array, const std::vector<int>& group,
			const MetricPreparation& object, const bool& isPBF);


	/*
	 * @brief Compute the validity measurement for the clustering result of PCA case only
	 *
	 * @param[in] array The matrix coordinates of streamlines
	 * @param[in] group The labels of different streamlines
	 */
	void computeValue(const MatrixXf& array, const std::vector<int>& group);

private:

	/*
	 * @brief min and max of S_c
	 */
	float min_Sc, max_Sc;


	/*
	 * @brief Get the MST (minimal spanning tree) results for the clustering with general norm option
	 *
	 * @param[out] values The output value vector to be updated
	 * @param[in] clusterNode The candidate index in the vector
	 * @param[in] object MetricPreparation class object
	 * @param[in] normOption The norm option
	 * @param[in] array The matrix coordinates
	 * @param[in] isPBF whether it is a PBF data set or not
	 */
	void getMST_Parent_Node(std::tuple<float, float, float>& values,
							const std::vector<int>& clusterNode,
							const MetricPreparation& object,
							const int& normOption,
							const MatrixXf& array,
							const bool& isPBF);


	/*
	 * @brief Get the MST (minimal spanning tree) results for the clustering with PCA
	 *
	 * @param[out] values The output value vector to be updated
	 * @param[in] clusterNode The candidate index in the vector
	 * @param[in] array The matrix coordinates
	 */
	void getMST_Parent_Node(std::tuple<float, float, float>& values,
				const std::vector<int>& clusterNode, const MatrixXf& array);


	/*
	 * @brief Get the Sc value for the clustering with general norm option
	 *
	 * @param[in] isPBF Whether the data set is PBF or not
	 * @param[in] distM The distance matrix as input
	 * @param[in] clusterNode The candidate index in the vector
	 * @param[in] rangeValue The ranged value as input
	 * @param[in] object MetricPreparation class object
	 * @param[in] normOption The norm option
	 * @param[in] array The matrix coordinates
	 * @param[out] index The cluster index to be updated
	 */
	const float get_Sc_by_range(const bool& isPBF, const Eigen::MatrixXf& distM,
			                    const std::vector<int>& clusterNode, const float& rangeValue,
								const MetricPreparation& object, const int& normOption, const MatrixXf& array,
								int& index);

	/*
	 * @brief Get the Sc value for the clustering with PCA only
	 *
	 * @param[in] distM The distance matrix as input
	 * @param[in] clusterNode The candidate index in the vector
	 * @param[in] rangeValue The ranged value as input
	 * @param[in] array The matrix coordinates
	 * @param[out] index The cluster index to be updated
	 */
	const float get_Sc_by_range(const Eigen::MatrixXf& distM, const std::vector<int>& clusterNode,
						        const float& rangeValue, const MatrixXf& array, int& index);
};

#endif /* SRC_COMMON_VALIDITYMEASUREMENT_H_ */
