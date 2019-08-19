#ifndef _PREDEFINED_H
#define _PREDEFINED_H

#include "IOHandler.h"
#include "Initialization.h"
#include "Silhouette.h"


// the dataset class to record relevant information
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


// define a treeNode structure to store AHC clustering tree
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


// remove two elements in template vector
template <class T>
void deleteVecElements(std::vector<T>& origine, const T& first, const T& second);

#endif
