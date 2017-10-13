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
	static void generateRandomPos(MatrixXf& clusterCenter,
								  const int& column,
								  const MatrixXf& cArray,
								  const int& Cluster);

	static void generateFromSamples(MatrixXf& clusterCenter,
								    const int& column,
								    const MatrixXf& cArray,
								    const int& Cluster);

	static void generateFarSamples(MatrixXf& clusterCenter,
								   const int& column,
								   const MatrixXf& cArray,
								   const int& Cluster,
								   const int& normOption,
								   const MetricPreparation& object);

};


#endif