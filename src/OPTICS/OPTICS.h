/*
 * @brief The class for OPTICS clustering
 * @details
 * 	The algorithm is implemented based on https://en.wikipedia.org/wiki/OPTICS_algorithm
 * @author Lieyu Shi
 */

#ifndef _OPTICS_H_
#define _OPTICS_H_

#include "Predefined.h"
#include "ValidityMeasurement.h"
#include <string>
#include <queue>


/*
 * @brief The DensityClustering class to perform OPTICS clustering on integral curves
 */
class DensityClustering
{
public:

	/*
	 * @brief The constructor with parameters of arguments
	 *
	 * @param[in] argc The count of arguments
	 * @param[in] argv The strings of arguments
	 */
	DensityClustering(const int& argc,
					  char **argv);


	/*
	 * @brief The destructor
	 */
	~DensityClustering();


	/*
	 * @brief Perform the OPTICS clustering on the input data set
	 */
	void performClustering();

private:

	/*
	 * @brief final output vector or ordered list
	 */
	vector<int> orderedList;

	/*
	 * @brief The vectors of nodes for the point for OPTICS clustering
	 */
	vector<PointNode> nodeVec;

	/*
	 * @brief The MetricPreparation object for computing the distance
	 */
	MetricPreparation object;

	/*
	 * @brief The norm option
	 */
	int normOption;

	/*
	 * @brief The data set
	 */
	DataSet ds;

	/*
	 * @brief Whether it is a PBF or not
	 */
	bool isPBF;

	/*
	 * @brief Whether it is pathlines or not
	 */
	bool isPathlines;


	/*
	 * @brief Set the data set from the argument and set necessary parameters for the sampling
	 *
	 * @param[in] argc The count of arguments
	 * @param[in] argv The strings of arguments
	 */
	void setDataset(const int& argc,
				    char **argv);


	/*
	 * @brief Set the norm option
	 */
	void setNormOption();


	/*
	 * @brief Perform the OPTICS clustering with input parameters
	 *
	 * @param[in] radius_eps The radius of neighborhood search for the OPTICS
	 * @param[in] minPts The minimal points for the neighborhood search
	 */
	void OPTICS(const float& radius_eps,
				const int& minPts);


	/*
	 * @brief Update the priority queue with the epsilon-neighborhood of two points
	 *
	 * @param[in] index The index of point
	 * @param[in] neighbor The neighborhood points
	 * @param[out] seeds The LinkedList object
	 * @param[in] radius_eps The radius of neighborhood search for the OPTICS
	 * @param[in] minPts The minimal points for the neighborhood search
	 */
	void update(const int& index,
				const vector<int>& neighbor,
				LinkedList& seeds,
				const float& radius_eps,
				const int& minPts);


	/*
	 * @brief It is to query about the neighborhood candidates for a given point
	 *
	 * @param[in] index The index of the target point
	 * @param[in] radius_eps The radius of neighborhood search for the OPTICS
	 * @return vector<int> object that contains the neighboring candidates
	 */
	const vector<int> regionQuery(const int& index,
								  const float& radius_eps);

	/*
	 * @brief Get the distance range
	 *
	 * @param[out] minDist The minimal distance
	 * @param[out] maxDist The maximal distance
	 */
	void getDistRange(float& minDist, 
					  float& maxDist);

	/*
	 * @brief Set the minimal points as parameter
	 */
	const int setMinPts();


	/*
	 * @brief Set and return the epsilong with ratio to the maxDist
	 *
	 * @param[in] minDist The minimal distance
	 * @param[in] maxDist The maximal distance
	 * @param The selected radius as parameter for the OPTICS clustering
	 */
	const float setTimesMin(const float& minDist, 
					  		const float& maxDist);


	/*
	 * @brief Get the reachability distance for the target
	 *
	 * @param[in] first The index of first point
	 * @param[in] target The index of target point
	 * @param[in] minPts The min points of the clustering
	 * @return A float value for the distance
	 */
	const float getReachability(const int& first,
								const int& target,
								const int& minPts);


	/*
	 * @brief Extract closest and furthest representatives and perform clustering evaluation for the clustering results
	 *
	 * @param[in] radius_eps The radius of neighborhood search for the OPTICS
	 * @param[in] minPts The minimal points for the neighborhood search
	 */
	void extractFeatures(const float& radius_eps,
							   const int& minPts);

	/*
	 * @brief Compute the cored distance for the candidates
	 *
	 * @param[in] radius_eps The radius of neighborhood search for the OPTICS
	 * @param[in] minPts The minimal points for the neighborhood search
	 */
	void computeCoredDistance(const float& radius_eps,
							  const int& minPts);


	/*
	 * @brief how to get group information based on reachability-plot
	 *
	 * @param[in] radius_eps The radius of neighborhood search for the OPTICS
	 */
	void getGroup(const float& radius_eps);

	/*
	 * @brief Write reachability plot data in the txt for further evaluation
	 */
	void writeReachability();


	/*
	 * @brief set the eps as averaged minPt-th dist
	 *
	 * @param[in] radius_eps The radius of neighborhood search for the OPTICS
	 */
	const float getMinPt_thDist(const int& minPts);
};

#endif
