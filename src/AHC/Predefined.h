/*
 * @brief Preliminary class definition before the AHC clustering algorithms
 * @author Lieyu Shi
 */

#ifndef _PREDEFINED_H
#define _PREDEFINED_H

#include "IOHandler.h"
#include "Initialization.h"
#include "Silhouette.h"


/*
 * @brief The class to store data set related information
 * @details
 * 	It has coordinates, matrix representation, max elements, vertex count, dimension, string name, full vtk name
 * 	and data set name
 */
struct DataSet
{
	vector<vector<float> > dataVec;	//original dataset
	Eigen::MatrixXf dataMatrix;	//sampled dataset
	int maxElements = -1;
	int vertexCount = -1;
	int dimension = -1;

	string strName;
	string fullName;
	string dataName;

};


/*
 * @brief The class to store the cluster and its inclusive candidates inside
 */
struct Ensemble
{
	int index = -1;

	/* to alleviate the computational cost to traverse all node elements */
	std::vector<int> element;

	Ensemble(const int& index): index(index)
	{}

	Ensemble()
	{}
};


/*
 * @brief Delete the elements in the vector
 *
 * @param[out] Original The vector to be operated on
 * @param[in] first The first index to be delete
 * @param[in] second The second index to be deleted
 */
template <class T>
void deleteVecElements(std::vector<T>& origine, const T& first, const T& second);


/*
 * @brief The class that contains the distance between two nodes i and j
 * @details
 * 	we will use a min-heap to perserve sorted distance for hirarchical clustering
 */
struct DistNode
{
	int first = -1, second = -1;
	float distance = -1.0;

	DistNode(const int& first, const int& second, const float& dist):first(first), second(second), distance(dist)
	{}

	DistNode()
	{}
};

#endif
