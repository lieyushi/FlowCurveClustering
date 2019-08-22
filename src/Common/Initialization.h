/*
 * @brief This is the source cpp for the class Initialization.h. It is for the k-means initialization
 * @author Lieyu Shi
 */


#ifndef _INITIALIZATION_H_
#define _INITIALIZATION_H_

#include <eigen3/Eigen/Dense>
#include <vector>
#include <ctime>
#include <cassert>
#include "Distance.h"


using namespace std;
using namespace Eigen;


/*
 * @brief The class to initialize the samples for the k-means clustering technique
 */
class Initialization
{
public:


	/*
	 * @brief To generate the random coordinates for the k-means initialization
	 *
	 * @param[out] clusterCenter The random initialization to be updated
	 * @param[in] column The column size
	 * @param[in] cArray The input matrix coordinates
	 * @param[in] Cluster The number of centroids
	 */
	static void generateRandomPos(MatrixXf& clusterCenter,
								  const int& column,
								  const MatrixXf& cArray,
								  const int& Cluster);


	/*
	 * @brief To generate the initialization from the samples
	 *
	 * @param[out] clusterCenter The initialized centroid coordinates
	 * @param[in] column The size of column
	 * @param[in] cArray The original matrix coordinates as input
	 * @paramp[in] Cluster The count of clusters
	 */
	static void generateFromSamples(MatrixXf& clusterCenter,
								    const int& column,
								    const MatrixXf& cArray,
								    const int& Cluster);


	/*
	 * @brief This is the k-means++ initialization
	 *
	 * @param[out] clusterCenter The cluster centroid to be initialized
	 * @param[in] column The size of column
	 * @param[in] cArray The matrix coordinates of the lines
	 * @param[in] normOption The norm option
	 * @param[in] object The MetricPreparation
	 */
	static void generateFarSamples(MatrixXf& clusterCenter,
								   const int& column,
								   const MatrixXf& cArray,
								   const int& Cluster,
								   const int& normOption,
								   const MetricPreparation& object);

};


#endif
