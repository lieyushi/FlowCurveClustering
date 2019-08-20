#ifndef _IOHANDLER_H_
#define _IOHANDLER_H_

#include <fstream>
#include <vector>
#include <iostream>
#include <cstring>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <climits>
#include <cassert>
#include <float.h>
#include <unordered_map>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/SVD>


#include "Silhouette.h"


using namespace std;
using namespace Eigen;


// the class of extractedLine to record the linenumber, cluster index and coordinates
struct ExtractedLine
{
	int lineNum;
	int cluster;
	ExtractedLine(const int& pointIndex,
				  const int& cluster)
				 :lineNum(pointIndex),cluster(cluster)
	{}
};


// the class of the centroid lines
struct MeanLine
{
	std::vector<float> minCenter;
	int cluster;
	MeanLine(const std::vector<float>& minCenter,
			 const int& cluster)
	        :minCenter(minCenter), cluster(cluster)
	{}
};


// the class of streamline query results
struct StringQuery
{
	int index;
	std::vector<int> neighbor;
	StringQuery()
	{   }
	StringQuery(const int& index,
				const std::vector<int>& neighbor): 
				index(index), neighbor(neighbor)
	{	}
};


// the class of FeatureLine for clustering result
struct FeatureLine
{
	std::vector<MeanLine> centerMass;
	std::vector<ExtractedLine> closest;
	std::vector<ExtractedLine> furthest;
	std::vector<int> group;
	std::vector<int> totalNum;

	FeatureLine()
	{

	}

	FeatureLine(const std::vector< std::vector<float> >& dataVec)
	{
		group = std::vector<int>(dataVec.size());
		totalNum = std::vector<int>(dataVec.size());
	}
};


// most of the functions related to I/O operation are contained inside
class IOHandler
{
	
public:
	
	/*
	 * @brief Read the data from the file into the vector<vector<>>
	 * @param fileName: The file name of the data set
	 * @param dataVec: The streamline coordinates to be updated
	 * @param vertexCount: The total vertex count
	 * @param dimension: 2 or 3 indicates it's 2D or 3D point
	 * @param maxElement: The max dimension of streamlines
	 */
	static void readFile(const string& fileName, 
						 std::vector< std::vector<float > >& dataVec, 
						 int& vertexCount, 
						 const int& dimension,
						 int& maxElement);

	/*
	 * @brief Read the particle-based data from the given frame number
	 * @param fileName: The file of the particle based data
	 * @param dataVec: The streamline coordinates
	 * @param vertexCount: The total point count
	 * @param dimension: 2 (2d) or 3 (3d) points
	 * @param trajectoryNum: Number of trajectories
	 * @param Frame: Number of frames
	 */
	static void readFile(const string& fileName, 
						 std::vector< std::vector<float > >& dataVec, 
						 int& vertexCount, 
						 const int& dimension,
						 const int& trajectoryNum, 
						 const int& Frame);

	/*
	 * @brief Print the vtk for streamlines
	 * @param fileName: The given vtk file name and position
	 * @param dataVec: The streamline coordinates
	 * @param vertexCount: the count of vertices
	 * @param dimension: 2 or 3 indicates 2D or 3D
	 * @param clusterNumber: the labels for each streamline
	 * @param sCluster: The float scalar value for each streamline
	 */
	static void printVTK(const string& fileName, 
						 const std::vector< std::vector<float > >& dataVec, 
						 const int& vertexCount, 
						 const int& dimension,
						 const std::vector<int>& clusterNumber,
						 const std::vector<float>& sCluster);


	/*
	 * @brief Print the vtk of the streamlines
	 * @param fileName: The name of the data set file
	 * @param dataVec: The streamline coordinates
	 * @param vertexCount: total count of vertices
	 * @param dimension: 2 or 3
	 */
	static void printVTK(const string& fileName, 
						 const std::vector< std::vector<float > >& dataVec, 
						 const int& vertexCount, 
						 const int& dimension);


	/*
	 * @brief Print the vtk file of the streamlines
	 * @param fileName: the name of the data set
	 * @param dataVec: The streamline coordinates
	 * @param vertexCount: the count of vertices
	 * @param dimension: 2 or 3
	 * @param sCluster: the scalar value for streamlines
	 */
	static void printVTK(const string& fileName, 
						 const std::vector<MeanLine>& dataVec, 
						 const int& vertexCount, 
						 const int& dimension,
						 const std::vector<float>& sCluster);

	static void printToFull(const std::vector< std::vector<float> >& dataVec, 
							const std::vector<int>& group, 
				 			const std::vector<int>& totalNum, 
				 			const string& groupName, 
				 			const string& fullName, 
				 			const int& dimension);

	static void printToFull(const std::vector< std::vector<float> >& dataVec, 
							const std::vector<float>& sData, 
				 			const string& groupName, 
				 			const string& fullName, 
				 			const int& dimension);

	static void printToFull(const std::vector< std::vector<float> >& origin,
				 			const std::vector<int>& group, 
				 			const string& fullName,
				 			const string& groupName,
				 			const int& dimension);

	static void printToFull(const std::vector< std::vector<float> >& dataVec, 
							const std::vector<int>& group, 
							const std::vector<float>& sCluster, 
				 			const string& groupName, 
				 			const string& fullName, 
				 			const int& dimension);

	static void writeReadme(const double& PCA_KMeans_delta, 
							const double& KMeans_delta);

	static void writeReadme(const string& comment,
							const std::vector<float>& sAverage);

	static void writeReadme(const std::vector<string>& timeName, 
							const std::vector<double>& timeDiff,
							const int& cluster);

	static void writeReadme(const std::vector<string>& timeName,
					 const std::vector<string>& timeDiff,
					 const int& cluster);

	static void writeReadme(const std::vector<ExtractedLine>& closest, 
							const std::vector<ExtractedLine>& furthest, 
							const int& normOption);

	static void writeReadme(const std::vector<ExtractedLine>& closest, 
							const std::vector<ExtractedLine>& furthest);

/* write the average rotation of closest and furthest extraction */
	static void writeReadme(const float& closestAverage, const float& furthestAverage);

	static void writeReadme(const string& comments);

/* write value of the silhouette class */
	static void writeReadme(const float& entropy, const Silhouette& sil, const string& norm_str);

/* print information into README */
	static void writeReadMe(const float& value, const string& dataSet, const string& clustering,
							const string& value_name);

	static void writeGroupSize(const std::vector<int>& storage);


/* expand array to the longest size so that we can perform entrywise comparison */
	static void expandArray(MatrixXf& data, 
							const std::vector< std::vector<float> >& dataVec, 
							const int& dimension, 
							const int& maxElements);

/* directly align vector object to Eigen::MatrixXf. Map vector to vector as well */
	static void expandArray(std::vector<std::vector<float> >& equalArray,
							const std::vector<std::vector<float> >& trajectories, 
						 	const int& dimension,
						 	const int& maxRowNum);

/* Uniformly sample array along streamlines instead of filling by the last index */
	static void sampleArray(MatrixXf& data, 
							const std::vector< std::vector<float> >& dataVec, 
							const int& dimension, 
							const int& maxElements);	

/* form array directly copied by vector for distribution-based comparison */
	static void formArray(float ***data, 
						  const std::vector< std::vector<float> >& dataVec, 
						  const int& dimension);

/* sample equal-sized array by equal arcs given a numOfVertices count */
	static void uniformArcSampling(MatrixXf& data,
								   const std::vector< std::vector<float> >& dataVec,
								   const int& dimension,
								   const int& maxElements);


	static void deleteArray(float **data, 
							const int& row);

	static void assignVec(std::vector<std::vector<float> >& closestStreamline,
						  std::vector<int>& cluster, 
						  const std::vector<ExtractedLine>& closest, 
						  int& pointNumber,
						  const std::vector< std::vector<float> >& dataVec);

	static void assignVec(std::vector<int>& cluster,
						  const std::vector<MeanLine>& centerMass);

	static void writeGroup(const std::vector<int>& group, 
						   const std::vector< std::vector<float> >& dataVec);

	static void printQuery(const int& normOption,
						   const int& order,
						   const StringQuery& queryResult, 
						   const std::vector<std::vector<float> >& dataVec);

	static void printTXT(float **data,
						 const int& Row,
						 const int& Column);


/* Birch code working header files for file reading */

	static void printFeature(const string& fileName,
				 			 const std::vector<std::vector<float> >& array,
				 			 const std::vector<float>& sCluster,
				 			 const int& dimension);

	static void printFeature(const string& fileName,
							 const std::vector<std::vector<float> >& array,
							 const std::vector<float>& sCluster,
							 const std::vector<float>& rotation,
							 const int& dimension);


	static void printClusters(const std::vector< std::vector<float> >& dataVec, 
							  const std::vector<int>& group, 
				 			  const std::vector<int>& storage, 
				 			  const string& groupName, 
				 			  const string& fullName, 
				 			  const int& dimension);

	static void printClustersNoise(const std::vector< std::vector<float> >& dataVec, 
							  const std::vector<int>& group, 
				 			  const std::vector<int>& storage, 
				 			  const string& groupName, 
				 			  const string& fullName, 
				 			  const int& dimension);

	static void generateReadme(const std::vector<string>& activityList,
							   const std::vector<double>& timeList,
							   const int& normOption,
							   const int& numClusters,
							   const float& sValue,
							   const float& threshold);

	static void generateReadme(const std::vector<string>& activityList,
							   const std::vector<string>& timeList);

	/* need to store each label for elements for NID computation */
	static void generateGroups(const std::vector<std::vector<int> >& storage);

	/* need to store each label for elements for NID computation */
	static void generateGroups(const std::vector<std::vector<int> >& storage, const string& fileName);

	/* read clustering number as a dictionary */
	static void readClusteringNumber(std::unordered_map<int,int>& clusMap, const string& fileName);
};


#endif
