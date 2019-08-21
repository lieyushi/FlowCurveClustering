/*
 * @brief The preliminary class before the AHC clustering starts
 * @author Lieyu Shi
 */


#ifndef _PREDEFINED_H
#define _PREDEFINED_H

#include "IOHandler.h"
#include "Initialization.h"
#include "Silhouette.h"


/*
 * @brief The class to store the information of the data set
 * @details
 * 	Record the coordinates, matrix, max elements, vertex count, dimension, string name, full vtk name and data name of
 * 	the data set
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
 * @brief The node that contains the candidates included in the cluster
 */
struct Ensemble
{
	int index = -1;
	bool merged = false;
	std::vector<int> element;

	Ensemble(const int& index): index(index)
	{}

	Ensemble()
	{}
};


/*
 * @brief Delete the two elements in the vector
 *
 * @param[out] original The input vector to be operated on
 * @param[in] first The first element to be deleted
 * @param[in] second The second element to be deleted
 */
template <class T>
void deleteVecElements(std::vector<T>& origine, const T& first, const T& second);

#endif
