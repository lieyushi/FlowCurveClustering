/*
 * @brief The class to perform the k-medoids clustering
 * @author Lieyu Shi
 */


#ifndef _KMEDOIDS_H
#define _KMEDOIDS_H


#include "IOHandler.h"
#include "Initialization.h"
#include "Silhouette.h"
#include "ValidityMeasurement.h"


/*
 * @brief The struct to record the parameters for k-medoids
 */
struct Parameter
{
	int initialization;
	bool isSample;
	Parameter(const int& initialization, const bool& isSample)
	: initialization(initialization), isSample(isSample)
	{}
	Parameter()
	{}
	~Parameter()
	{}
};


/*
 * @brief The class object to record the event and time
 */
struct TimeRecorder
{
	std::vector<string> eventList;
	std::vector<string> timeList;
};


/*
 * @brief The class to perform the k-medoids clustering technique on the integral curves
 */
class KMedoids
{
public:

	/*
	 * @brief The k-medoids class constructor with parameter and data set input
	 *
	 * @param[in] pm The Parameter class object for the user input control
	 * @param[in] data The matrix coordinates of the streamlines
	 * @param[in] numOfClusters The number of clusters as input
	 */
	KMedoids(const Parameter& pm,
			 const Eigen::MatrixXf& data,
			 const int& numOfClusters);


	/*
	 * @brief The destructor
	 */
	~KMedoids();


	/*
	 * @brief Perform the k-medoids clustering technique and calculate the clustering evaluation metrics
	 *
	 * @param[out] fline The feature line extraction result (centroid, closest, furthest candidates)
	 * @param[in] normOption The norm option
	 * @param[out] sil Silhouette The Silhouette class object to record clustering evaluation metrics
	 * @param[out] tr The TimeRecorder class object
	 */
	void getMedoids(FeatureLine& fline,
					const int& normOption,
					Silhouette& sil,
					TimeRecorder& tr) const;

	/*
	 * @brief The number of clusters as result
	 */
	int numOfClusters;

private:

	/*
	 * @brief The initialization of k-medoids
	 * @details
	 * 	It has three options: 1. random coordinate initialization, 2. randomly chosen from samples, 3. chosen with
	 * 	k-means++ (sampled based on uniform distribution probability)
	 */
	int initialStates;

	/*
	 * @brief Whether the medoid of the cluster is calculated from the sample or by iteration.
	 * @details
	 * 	True means medoid should be from samples, false means medoids are calculated by iterative calculation
	 */
	bool isSample;

	/*
	 * @brief The matrix coordinates
	 */
	Eigen::MatrixXf data;


	/*
	 * @brief Get the initial centroid by k-medoids initialization
	 *
	 * @param[out] initialCenter The initialized center coordinates
	 * @param[in] object The MetricPreparation object for calculating the distance values
	 * @param[in] normOption The norm option
	 */
	void getInitCenter(MatrixXf& initialCenter,
					   const MetricPreparation& object,
					   const int& normOption) const;


	/*
	 * @brief Compute the medoids by either calculating the median or iteration
	 *
	 * @param[out] centerTemp The medoid coordinates to be updated
	 * @param[in] neighborVec The neighboring candidates belonging to a cluster
	 * @param[in] normOption The norm option
	 * @param[in] object The MetricPreparation object for distance calculation
	 */
	void computeMedoids(MatrixXf& centerTemp, 
						const vector<vector<int> >& neighborVec, 
						const int& normOption, 
						const MetricPreparation& object) const;
};


#endif
