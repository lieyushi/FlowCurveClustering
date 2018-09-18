/*
 * readDistRange.cpp
 *
 *  Created on: Sep 17, 2018
 *      Author: lieyu
 */
#include <iostream>
#include <fstream>
#include <string>
#include <float.h>
#include <limits>
#include "IOHandler.h"
#include "Initialization.h"

using namespace std;
using namespace Eigen;

struct Dataset
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

void getDistRange(const Dataset& ds);
/* set dataset from user command */
void setDataset(Dataset& ds, const int& argc, char **argv);


int main(int argc, char* argv[])
{
	if(argc!=3)
	{
		std::cout << "parameter option is not right!" << std::endl;
		exit(1);
	}

	Dataset ds;
	setDataset(ds, argc, argv);

	getDistRange(ds);

	return 0;
}

void getDistRange(const Dataset& ds)
{
	for(int i=0; i<16; ++i)
	{
		if(i!=0 && i!=1 && i!=2 && i!=4 && i!=12 && i!=13 && i!=14 && i!=15)
			continue;
		/* very hard to decide whether needed to perform such pre-processing */
		MetricPreparation object = MetricPreparation(ds.dataMatrix.rows(), ds.dataMatrix.cols());
		object.preprocessing(ds.dataMatrix, ds.dataMatrix.rows(), ds.dataMatrix.cols(), i);

		deleteDistanceMatrix(ds.dataMatrix.rows());
		getDistanceMatrix(ds.dataMatrix, i, object);

		const int& Row = ds.dataMatrix.rows();
		float min_dist = numeric_limits<float>::max(), max_dist = numeric_limits<float>::min();

	#pragma omp parallel for reduction(min:min_dist) num_threads(8)
		for (int i = 0; i < Row; ++i)
		{
			for (int j = 0; j < Row; ++j)
			{
				if(i==j)
					continue;
				min_dist = std::min(min_dist, distanceMatrix[i][j]);
			}
		}

	#pragma omp parallel for reduction(max:max_dist) num_threads(8)
		for (int i = 0; i < Row; ++i)
		{
			for (int j = 0; j < Row; ++j)
			{
				if(i==j)
					continue;
				max_dist = std::max(max_dist, distanceMatrix[i][j]);
			}
		}

		std::cout << "norm " << i << " has min " << min_dist << " and max " << max_dist << std::endl;

		std::ofstream readme("../dataset/dist_range", ios::app | ios::out);
		if(readme.fail())
		{
			std::cout << "Error for opening readme!" << std::endl;
			exit(1);
		}

		readme << "For norm " << i << ", min is " << min_dist << ", max is " << max_dist << ", and (max - min) is " <<
				(max_dist-min_dist) << std::endl;
		readme << std::endl;

		readme.close();
	}
}

/* set dataset from user command */
void setDataset(Dataset& ds, const int& argc, char **argv)
{
	if(argc!=3)
	{
		std::cout << "Input argument should have 3!" << endl
		          << "./cluster inputFile_name(in dataset folder) "
		          << "data_dimension(3)" << endl;
		exit(1);
	}
	ds.strName = string("../dataset/")+string(argv[1]);
	ds.dataName = string(argv[1]);
	ds.dimension = atoi(argv[2]);

	int sampleOption;
    std::cout << "choose a sampling method for the dataset?" << std::endl
	    	  << "1.directly filling with last vertex; 2. uniform sampling." << std::endl;
	std::cin >> sampleOption;
	assert(sampleOption==1||sampleOption==2);

	IOHandler::readFile(ds.strName,ds.dataVec,ds.vertexCount,ds.dimension,ds.maxElements);

	ds.fullName = ds.strName+"_full.vtk";
	IOHandler::printVTK(ds.fullName, ds.dataVec, ds.vertexCount, ds.dimension);

	if(sampleOption==1)
		IOHandler::expandArray(ds.dataMatrix,ds.dataVec,ds.dimension,ds.maxElements);
	else if(sampleOption==2)
		IOHandler::sampleArray(ds.dataMatrix,ds.dataVec,ds.dimension,ds.maxElements);
}





