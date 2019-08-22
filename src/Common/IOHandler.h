/*
 * @brief This class contains the I/O functions for the program
 * @author Lieyu Shi
 */

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


/*
 * @brief the class of extractedLine to record the linenumber, cluster index and coordinates
 */
struct ExtractedLine
{
	int lineNum;
	int cluster;
	ExtractedLine(const int& pointIndex,
				  const int& cluster)
				 :lineNum(pointIndex),cluster(cluster)
	{}
};


/*
 * @brief the class of the centroid lines
 */
//
struct MeanLine
{
	std::vector<float> minCenter;
	int cluster;
	MeanLine(const std::vector<float>& minCenter,
			 const int& cluster)
	        :minCenter(minCenter), cluster(cluster)
	{}
};


/*
 * @brief the class of streamline query results
 */
//
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


/*
 * @brief the class of FeatureLine for clustering result
 */
//
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


/*
 * @brief most of the functions related to I/O operation are contained inside
 */
class IOHandler
{
	
public:
	
	/*
	 * @brief Read the data from the file into the vector<vector<>>
	 *
	 * @param[in] fileName The file name of the data set
	 * @param[out] dataVec The streamline coordinates to be updated
	 * @param[out] vertexCount The total vertex count
	 * @param[in] dimension 2 or 3 indicates it's 2D or 3D point
	 * @param[out] maxElement The max dimension of streamlines
	 */
	static void readFile(const string& fileName, 
						 std::vector< std::vector<float > >& dataVec, 
						 int& vertexCount, 
						 const int& dimension,
						 int& maxElement);

	/*
	 * @brief Read the particle-based data from the given frame number
	 *
	 * @param[in] fileName The file of the particle based data
	 * @param[out] dataVec The streamline coordinates
	 * @param[in] vertexCount The total point count
	 * @param[in] dimension 2 (2d) or 3 (3d) points
	 * @param[in] trajectoryNum Number of trajectories
	 * @param[in] Frame Number of frames
	 */
	static void readFile(const string& fileName, 
						 std::vector< std::vector<float > >& dataVec, 
						 int& vertexCount, 
						 const int& dimension,
						 const int& trajectoryNum, 
						 const int& Frame);

	/*
	 * @brief Print the vtk for streamlines
	 *
	 * @param[in] fileName The given vtk file name and position
	 * @param[in] dataVec The streamline coordinates
	 * @param[in] vertexCount the count of vertices
	 * @param[in] dimension 2 or 3 indicates 2D or 3D
	 * @param[in] clusterNumber the labels for each streamline
	 * @param[in] sCluster The float scalar value for each streamline
	 */
	static void printVTK(const string& fileName, 
						 const std::vector< std::vector<float > >& dataVec, 
						 const int& vertexCount, 
						 const int& dimension,
						 const std::vector<int>& clusterNumber,
						 const std::vector<float>& sCluster);


	/*
	 * @brief Print the vtk of the streamlines
	 *
	 * @param[in] fileName The name of the data set file
	 * @param[in] dataVec The streamline coordinates
	 * @param[in] vertexCount total count of vertices
	 * @param[in] dimension 2 or 3
	 */
	static void printVTK(const string& fileName, 
						 const std::vector< std::vector<float > >& dataVec, 
						 const int& vertexCount, 
						 const int& dimension);


	/*
	 * @brief Print the vtk file of the streamlines
	 *
	 * @param[in] fileName the name of the data set
	 * @param[in] dataVec The streamline coordinates
	 * @param[in] vertexCount the count of vertices
	 * @param[in] dimension 2 or 3
	 * @param[in] sCluster the scalar value for streamlines
	 */
	static void printVTK(const string& fileName, 
						 const std::vector<MeanLine>& dataVec, 
						 const int& vertexCount, 
						 const int& dimension,
						 const std::vector<float>& sCluster);


	/*
	 * @brief Print the new group labels to the existing vtk file
	 *
	 * @param[in] dataVec The streamline coordinates as input
	 * @param[in] group The label vector for each streamline
	 * @param[in] totalNum The size of each cluster just in case to show dominant clusters
	 * @param[in] groupName The string text for the labels
	 * @param[in] fullName The full name of the vtk to be written in
	 * @param[in] dimension 2 or 3
	 */
	static void printToFull(const std::vector< std::vector<float> >& dataVec, 
							const std::vector<int>& group, 
				 			const std::vector<int>& totalNum, 
				 			const string& groupName, 
				 			const string& fullName, 
				 			const int& dimension);


	/*
	 * @brief Print scalar values to the existing vtk file
	 *
	 * @param[in] dataVec The streamline data set
	 * @param[in] sData The scalar values on the streamlines
	 * @param[in] groupName The string text for the scalar values
	 * @param[in] fullName The full name of the .vtk file
	 * @param[in] dimension 2 or 3
	 */
	static void printToFull(const std::vector< std::vector<float> >& dataVec, 
							const std::vector<float>& sData, 
				 			const string& groupName, 
				 			const string& fullName, 
				 			const int& dimension);


	/*
	 * @brief print the group information to the existing full vtk file with label names
	 *
	 * @param[in] dataVec The streamline coordinates
	 * @param[in] group The labels of the streamlines
	 * @param[in] fullName The full name of the primary vtk file
	 * @param[in] groupName The group name of the labels
	 * @param[in] dimension The dimension of the points
	 */
	static void printToFull(const std::vector< std::vector<float> >& origin,
				 			const std::vector<int>& group, 
				 			const string& fullName,
				 			const string& groupName,
				 			const int& dimension);


	/*
	 * @brief Print scalar and label values to the existing vtk file
	 *
	 * @param[in] dataVec The streamline data set
	 * @param[in] group The group labels for the streamlines
	 * @param[in] sCluster The scalar values on the streamlines
	 * @param[in] groupName The string text for the scalar values
	 * @param[in] fullName The full name of the .vtk file
	 * @param[in] dimension 2 or 3
	 */
	static void printToFull(const std::vector< std::vector<float> >& dataVec, 
							const std::vector<int>& group, 
							const std::vector<float>& sCluster, 
				 			const string& groupName, 
				 			const string& fullName, 
				 			const int& dimension);


	/*
	 * @brief Write two double values in the readme
	 *
	 * @param[in] PCA_KMeans_delta The time for PCA k-means
	 * @param[in] KMeans_delta The time for k-means
	 */
	static void writeReadme(const double& PCA_KMeans_delta, 
							const double& KMeans_delta);


	/*
	 * @brief Write string and float vector in the readme
	 *
	 * @param[in] comment The string to be put in the readme
	 * @param[in] sAverage The float vector to be written
	 */
	static void writeReadme(const string& comment,
							const std::vector<float>& sAverage);


	/*
	 * @brief Write the double and string vector and cluster in the readme
	 *
	 * @param[in] timeName The string vector
	 * @param[in] timeDiff The double vector
	 * @param[in] cluster The number of clusters
	 */
	static void writeReadme(const std::vector<string>& timeName, 
							const std::vector<double>& timeDiff,
							const int& cluster);


	/*
	 * @brief Write two string vectors in the readme for printing information
	 *
	 * @param[in] timeName The string vector
	 * @param[in] timeDiff The string vector
	 * @param[in] cluster The size of cluster
	 */
	static void writeReadme(const std::vector<string>& timeName,
					 const std::vector<string>& timeDiff,
					 const int& cluster);


	/*
	 * @brief Write the extracted line information w.r.t. the norm in the readme
	 *
	 * @param[in] closest The closest extracted lines in the readme
	 * @param[in] furthest The furthest extracted lines in the readme
	 * @param[in] normOption The norm option
	 */
	static void writeReadme(const std::vector<ExtractedLine>& closest, 
							const std::vector<ExtractedLine>& furthest, 
							const int& normOption);


	/*
	 * @brief Write closest and furthest extracted lines in the readme
	 *
	 * @param[in] closest The closest extracted lines
	 * @param[in] furthest The furthest extracted lines
	 */
	static void writeReadme(const std::vector<ExtractedLine>& closest, 
							const std::vector<ExtractedLine>& furthest);


	/*
	 * @brief Write two float values in the readme
	 *
	 * @param[in] closestAverage The float to be written
	 * @param[in] furthestAverage The float to be written
	 */
	static void writeReadme(const float& closestAverage, const float& furthestAverage);


	/*
	 * @brief write string arrays in the readme
	 *
	 * @param[in] comments The string to be put in the readme
	 */
	static void writeReadme(const string& comments);


	/*
	 * @brief Write the entropy, silhoutte and string in the readme
	 *
	 * @param[in] entropy The entropy value
	 * @param[in] sil The Silhouette class object
	 * @param[in] norm_str The string object
	 */
	static void writeReadme(const float& entropy, const Silhouette& sil, const string& norm_str);


	/*
	 * @brief Print float and string information into the readme
	 *
	 * @param[in] value The float value
	 * @param[in] dataSet The string type
	 * @param[in] clustering The string of the clustering technique
	 * @param[in] value_name The string for the value name
	 */
	static void writeReadMe(const float& value, const string& dataSet, const string& clustering,
							const string& value_name);


	/*
	 * @brief write the finalized group size in the storage
	 *
	 * @param[in] storage The size of clusters from the clustering results
	 */
	static void writeGroupSize(const std::vector<int>& storage);


	/*
	 * @brief Direct repeating the last point of the streamlines
	 *
	 * @param[out] data The matrix coordinates to be updated
	 * @param[in] dataVec The line coordinates of streamlines as input
	 * @param[in] dimension 2 or 3
	 * @param[in] maxElement The max dimension of streamlines
	 */
	static void expandArray(MatrixXf& data, 
							const std::vector< std::vector<float> >& dataVec, 
							const int& dimension, 
							const int& maxElements);


	/*
	 * @brief Directly align vector object to Eigen::MatrixXf. Map vector to vector as well, or called direct repeating
	 *
	 * @param[out] equalArray The matrix to be assigned the coordinates
	 * @param[in] trajectories The trajectory coordinates
	 * @param[in] dimension The dimension
	 * @param[in] maxElement The max element
	 */
	static void expandArray(std::vector<std::vector<float> >& equalArray,
							const std::vector<std::vector<float> >& trajectories, 
						 	const int& dimension,
						 	const int& maxRowNum);


	/*
	 * @brief Sample the streamlines on the intervals to preservce thier geometry shapes
	 *
	 * @param[out] data The matrix coordinates to be updated
	 * @param[in] dataVec The streamline coordinates as input
	 * @param[in] dimension 2 or 3
	 * @param[in] maxElements The maximal length of streamlines
	 */
	static void sampleArray(MatrixXf& data, 
							const std::vector< std::vector<float> >& dataVec, 
							const int& dimension, 
							const int& maxElements);	


	/*
	 * @brief Assign vector values to the pointer type
	 *
	 * @param[out] data The pointer to be updated
	 * @param[in] dataVec The streamline coordinates
	 * @param[in] dimension 2 or 3
	 */
	static void formArray(float ***data, 
						  const std::vector< std::vector<float> >& dataVec, 
						  const int& dimension);


	/*
	 * @brief Sample the streamlines with equal arc length
	 *
	 * @param[out] data The matrix coordinates to be updated
	 * @param[in] dataVec The input streamline coordinates
	 * @param[in] dimension 2 or 3
	 * @param[in] maxElements The max dimension of streamlines
	 */
	static void uniformArcSampling(MatrixXf& data,
								   const std::vector< std::vector<float> >& dataVec,
								   const int& dimension,
								   const int& maxElements);


	/*
	 * @brief delete the array
	 *
	 * @param[out] data The pointer type
	 * @param[in] row The number of rows
	 */
	static void deleteArray(float **data, 
							const int& row);


	/*
	 * @brief Assign the ExtractedLine information to the closest streamline coordinates and cluster
	 *
	 * @param[out] closestStreamline The closest extracted streamlines
	 * @param[out] cluster The number of clusters
	 * @param[in] closest The closest ExtractedLine object as input
	 * @param[out] pointNumber The number of points
	 * @param[in] dataVec The input coordinates of the streamlines
	 */
	static void assignVec(std::vector<std::vector<float> >& closestStreamline,
						  std::vector<int>& cluster, 
						  const std::vector<ExtractedLine>& closest, 
						  int& pointNumber,
						  const std::vector< std::vector<float> >& dataVec);


	/*
	 * @brief Assign the size of clustes to the center
	 *
	 * @param[out] cluster The size of clusters
	 * @param[in] centerMass The extracted mean line
	 */
	static void assignVec(std::vector<int>& cluster,
						  const std::vector<MeanLine>& centerMass);


	/*
	 * @brief Write group information into the local file as backup
	 *
	 * @param[in] The labels of group for the streamlines
	 * @param[in] The streamline coordinates
	 */
	static void writeGroup(const std::vector<int>& group, 
						   const std::vector< std::vector<float> >& dataVec);


	/*
	 * @brief Print the query information from the streamline query
	 *
	 * @param[in] normOption The norm option
	 * @param[in] order The int of order
	 * @param[in] queryResult The string query result
	 * @param[in] dataVec The streamline coordinates
	 */
	static void printQuery(const int& normOption,
						   const int& order,
						   const StringQuery& queryResult, 
						   const std::vector<std::vector<float> >& dataVec);


	/*
	 * @brief print the coordinates of the streamlines into the local txt file
	 *
	 * @param[in] The float pointer
	 * @param[in] Row The size of rows
	 * @param[in] Column The size of columns
	 */
	static void printTXT(float **data,
						 const int& Row,
						 const int& Column);


	/*
	 * @brief print the streamlines with attached scalar values in the vtk format
	 *
	 * @param[in] fileName The file name of the vtk
	 * @param[in] array The streamline coordinates
	 * @param[in] sCluster The float vector as input
	 * @param[in] dimension The dimension of the points
	 */
	static void printFeature(const string& fileName,
				 			 const std::vector<std::vector<float> >& array,
				 			 const std::vector<float>& sCluster,
				 			 const int& dimension);


	/*
	 * @brief print the streamlines with several scalar values in the vtk file
	 *
	 * @param[in] fileName The name of the vtk file
	 * @param[in] array The streamline coordinates
	 * @param[in] sCluster The float vector
	 * @param[in] rotation The float vector
	 * @param[in] dimension The dimension of the points
	 */
	static void printFeature(const string& fileName,
							 const std::vector<std::vector<float> >& array,
							 const std::vector<float>& sCluster,
							 const std::vector<float>& rotation,
							 const int& dimension);


	/*
	 * @brief print the label and size of clusters for the streamlines
	 *
	 * @param[in] dataVec The streamline coordinates
	 * @param[in] group The group labels of the streamlines
	 * @param[in] storage The size of clusters
	 * @param[in] groupName The string of the label name
	 * @param[in] fullName The full name of the primary vtk file
	 * @param[in] dimension The dimension of points
	 */
	static void printClusters(const std::vector< std::vector<float> >& dataVec, 
							  const std::vector<int>& group, 
				 			  const std::vector<int>& storage, 
				 			  const string& groupName, 
				 			  const string& fullName, 
				 			  const int& dimension);


	/*
	 * @brief print the clustered streamlines with noise information generated from DBSCAN
	 *
	 * @param[in] dataVec The streamline coordinates
	 * @param[in] group The labels of streamlines
	 * @param[in] storage The size of clusters
	 * @param[in] groupName The string of the labels
	 * @param[in] fullName The name of the primary vtk
	 * @param[in] dimension The dimension of the points
	 */
	static void printClustersNoise(const std::vector< std::vector<float> >& dataVec, 
							  const std::vector<int>& group, 
				 			  const std::vector<int>& storage, 
				 			  const string& groupName, 
				 			  const string& fullName, 
				 			  const int& dimension);


	/*
	 * @brief generate readme txt with number of clusters, input threshold and activity list, mostly for BIRCH
	 *
	 * @param[in] activityList The vector to store the activities
	 * @param[in] timeList The vector to store the float of time
	 * @param[in] normOption The norm option
	 * @param[in] numClusters The number of clusters generated
	 * @param[in] sValue The silhouette value
	 * @param[in] threshold The distance threshold as input for the BIRCH clustering
	 */
	static void generateReadme(const std::vector<string>& activityList,
							   const std::vector<double>& timeList,
							   const int& normOption,
							   const int& numClusters,
							   const float& sValue,
							   const float& threshold);


	/*
	 * @brief print the readme information with activity list and time list that records some status during the computation
	 *
	 * @param[in] activityList The vector that has the relative activities
	 * @param[in] timeList The vector of time for respective status
	 */
	static void generateReadme(const std::vector<string>& activityList,
							   const std::vector<string>& timeList);


	/*
	 * @brief Write the candidates of each cluster into the local file
	 *
	 * @param[in] storage The candidates belonging to each cluster
	 */
	static void generateGroups(const std::vector<std::vector<int> >& storage);


	/*
	 * @brief Generate and store the candidates of all the clusters in the local file
	 *
	 * @param[in] storage The candidates belonging to each cluster
	 * @param[in] fileName The name of the file for storing the candidiate information
	 */
	static void generateGroups(const std::vector<std::vector<int> >& storage, const string& fileName);


	/*
	 * @brief Read the number of clusters as input for different similarity measures from the local file
	 *
	 * @param[out] clusMap The hash map that records the number of cluster and its norm option
	 * @param[in] fileName The name of the file to be read in
	 */
	static void readClusteringNumber(std::unordered_map<int,int>& clusMap, const string& fileName);
};


#endif
