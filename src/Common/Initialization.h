/*
 * This class is the centroid initialization for k-means clustering technqiue
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


class Initialization
{
public:

	// generate the centroids from random coordinates
	static void generateRandomPos(MatrixXf& clusterCenter,
								  const int& column,
								  const MatrixXf& cArray,
								  const int& Cluster);

	// generate the centroids from the random samples
	static void generateFromSamples(MatrixXf& clusterCenter,
								    const int& column,
								    const MatrixXf& cArray,
								    const int& Cluster);

	// generate the centroids by k-means++
	static void generateFarSamples(MatrixXf& clusterCenter,
								   const int& column,
								   const MatrixXf& cArray,
								   const int& Cluster,
								   const int& normOption,
								   const MetricPreparation& object);

};


#endif
