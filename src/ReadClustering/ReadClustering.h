/*
 * ReadClustering.h
 *
 *  Created on: Mar 13, 2018
 *      Author: lieyu
 */

#ifndef SRC_READCLUSTERING_READCLUSTERING_H_
#define SRC_READCLUSTERING_READCLUSTERING_H_

#include "IOHandler.h"
#include "Initialization.h"
#include "Silhouette.h"
#include <unordered_map>
#include "ValidityMeasurement.h"
#include <sstream>

using namespace std;


/*
 * @brief The data set class object that contains the information
 * @details
 *  It includes original coordinates, the matrix coordinates, the map of string and candidates, the max group hash map
 *  with index, the number of elements and the neighbor vector including candidates of the clusters
 */
struct Dataset
{
	/* original coordinates */
	std::vector<std::vector<float> > dataVec;

	/* sampled array */
	Eigen::MatrixXf array;

	/* label information */
	unordered_map<string, std::vector<int> > groupAggregate;

	/* cluster number */
	unordered_map<string, int> maxGroup;

	/* number of elements inside */
	int numOfElements;

	std::vector<std::vector<int> > neighborVec;
};


/*
 * @brief The read clustering class object to read the clustering data and calculate the clustering evaluation metrics
 * from the input vtk file
 */
class ReadClustering {
public:

	/*
	 * @brief The default constructor
	 */
	ReadClustering();


	/*
	 * @brief The destructor
	 */
	virtual ~ReadClustering();


	/*
	 * @brief The public function called by main.cpp to calculate the clustering evaluation
	 *
	 * @param[in] fileName The name of the vtk file
	 */
	void getEvaluation(const char* fileName);

private:

	/*
	 * @brief The activityList vector to store event
	 */
	std::vector<string> activityList;

	/*
	 * @brief The timeList vector to store time information
	 */
	std::vector<string> timeList;

	/*
	 * @brief Dataset object
	 */
	Dataset ds;

	/*
	 * @brief the max element
	 */
	int maxElements;

	/*
	 * @brief judge whether it is a PBF or not
	 */
	bool isPBF;

	/*
	 * @brief read data from file and store it into Dataset ds
	 *
	 * @param[in] fileName The name of the vtk file as input
	 */
	void readData(const char* fileName);


	/*
	 * @brief compute four analysis evaluation measures
	 */
	void computeEvaluation();


	/*
	 * @breif write those analysis framework into the local file
	 */
	void writeAnalysis();


	/*
	 * @brief Compute evaluation based on norm option
	 *
	 * @param[out] The iterator of the hash map for recording the norm option and the input number
	 */
	void computeEvaluation(std::unordered_map<string, std::vector<int> >::const_iterator& iter);

	/*
	 * @brief Perform SVD decomposition for equal-sized streamlines
	 *
	 * @param[out] cArray The matrix coordinates after the dimensionality reduction process
	 * @param[in] data The matrix coordinates of the input integral curves
	 * @param[in] Row The row size
	 * @param[in] Column The column size
	 * @param[out] PC_Number The number of PCs
	 */
	void performSVD(MatrixXf& cArray, const Eigen::MatrixXf& data,
		 const int& Row, const int& Column, int& PC_Number);

};

#endif /* SRC_READCLUSTERING_READCLUSTERING_H_ */
