/*
 * @brief The file for spectral clustering needed to be decided which should be contained
 * @author Lieyu Shi
 */


#ifndef _PREDEFINED_H
#define _PREDEFINED_H

#include "IOHandler.h"
#include "Initialization.h"
#include "Silhouette.h"


/*
 * @brief The data set class object to contain basic information
 * @details
 * 	It includes original coordinates, sampled matrix coordinates, max elements, vertex count, dimension and names
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
 * @brief The struct that has the candidates and size of clusters
 */
struct Ensemble
{
	int size;
	std::vector<int> element;
};


/*
 * @brief Remove two elements in template vector
 *
 * @param[out] original The vector to be operated on
 * @param[in] first The first index
 * @param[in] second The second index
 */
template <class T>
void deleteVecElements(std::vector<T>& original, const T& first, const T& second);


#endif
