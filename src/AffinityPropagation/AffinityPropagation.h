/* @brief Affinity propagation is a newly emerging clustering techniques based on "message passing" between data points
 *
 * @details
 * 	The wikipage could be referenced at https://en.wikipedia.org/wiki/Affinity_propagation
 * A sample C++ code could be seen at https://github.com/jincheng9/AffinityPropagation/blob/master/affinity_propagation.cpp,
 * or https://github.com/nojima/affinity-propagation-sparse/blob/master/ap.cpp
 *
 * @author Lieyu Shi
 */

#ifndef _AFFINITY_PROPAGATION_H_
#define _AFFINITY_PROPAGATION_H_

#include "Predefined.h"
#include "ValidityMeasurement.h"
#include <unordered_set>
#include <map>
#include <string>


#define LAMBDA 0.5


/*
 * @brief Struct of storing the parameters for the AP clustering techniques
 *
 * @details
 * 	It stores the sample strategy (1. directly filling with last vertex, 2. uniform sampling, 3. equal-arc sampling),
 * 	extraction strategy (1. centroid, closest and furthest, 2. median), and the maxIteration.
 */
struct Para
{

	/* 1.directly filling with last vertex; 2. uniform sampling, 3. equal-arc sampling */
	int sampled;

	/* extraction option, 1. centroid, closest and furthest, 2. median */
	int extractOption;

	/* max iteration for AP clustering */
	int maxIteration;
};


/*
 * @brief The class for performing the affinity propagation clustering technique
 */
class AffinityPropagation
{

public:

	/*
	 * @brief This is a default constructor for AffinityPropagation class object
	 */
	AffinityPropagation();


	/*
	 * @brief Create an AffinityPropagation object with command line parameters
	 *
	 * @details
	 * 	Create an AffinityPropagation class object from argument strings and pre-set parameters. This constructor will be
	 * 	ready for the following calculation of the AP clustering algorithm
	 *
	 * @param[in] argc The count of argument as input
	 * @param[in] argv The argument line with string type of data set names and dimension count
	 * @param[in] p The Para object with some pre-defined parameter values
	 * @param[in] automatic The bool object for automatic parameter setting or not
	 */
	AffinityPropagation(const int& argc, char **argv, const Para& p, bool& automatic);


	/*
	 * @brief A destructor of the class
	 *
	 * @details
	 * 	Delete the global pointer variable distanceMatrix
	 */
	~AffinityPropagation();


	/*
	 * @brief The member function to perform the affinity propagation clustering on similarity measures
	 *
	 * @details
	 * 	Perform the AP clustering techniques on selected similarity measures and evaluation metric calculation. The selected
	 * 	similarity measures can be revised for any choice
	 */
	void performClustering();

private:


	/*
	 * @brief metric preparation object to be stored ahead of time
	 */
	MetricPreparation object;

	/*
	 * @brief input norm option
	 */
	int normOption = -1;

	/*
	 * @brief group information
	 */
	std::vector<int> group;

	/*
	 * @brief activityList vector to store event
	 */
	std::vector<string> activityList;

	/*
	 * @brief timeList vector to store time information
	 */
	std::vector<string> timeList;

	/*
	 * @brief store dataset information
	 */
	DataSet ds;

	/*
	 * @brief how many clusters to be needed
	 */
	int numberOfClusters = -1;

	/*
	 * @brief extraction option, 1. centroid, closest and furthest, 2. median, 3. statistical representation
	 */
	int extractOption = -1;

	/*
	 * @brief max iteration
	 */
	int maxIteration = -1;

	/*
	 * @brief tell whether it is a tag for PBF dataset
	 */
	bool isPBF;

	/*
	 * @brief tell whether it is a pathline data set or not
	 */
	bool isPathlines;

	/*
	 * @brief S[i,i] initialization format
	 */
	int initialOption;

	/*
	 * @brief whether use two-staged AP or not
	 */
	bool useTwoStage;


	/*
	 * @brief Extract the cluster representatives of clusters and calculate the evaluation metrics for the clustering results
	 *
	 * @details
	 * 	Extract the features and calculate the clustering evaluation metrics for the clustering result. Check the cluster size
	 * 	output information to guarantee the ascending order of cluster sizes. Print the extracted vtk files for visual
	 * 	analysis.
	 *
	 * @param[in] storage The number of candidates included in each cluster
	 * @param[in] neighborVec The candidates (only streamline index included) for each streamline cluster
	 * @param[in] centroid The centroids of streamline clusters
	 */
	void extractFeatures(const std::vector<int>& storage, const std::vector<std::vector<int> >& neighborVec,
            			 const Eigen::MatrixXf& centroid);


	/*
	 * @brief A member function to read the geometric coordinates from given file name and position
	 *
	 * @param[in] argc The count of arguments
	 * @param[in] argv The argument char array type that includes data set name and dimension count
	 */
	void setDataset(const int& argc, char **argv);


	/*
	 * @brief It provides console for user parameter input
	 */
	void getParameterUserInput();


	/*
	 * @brief It performs some necessary operations given the input parameters
	 *
	 * @param[in] p A Para object that contains parameters for pre-processing and activation for two-staged AP
	 */
	void setParameterAutomatic(const Para& p);


	/*
	 * @brief It performs the AP clustering algorithm on selected similarity measure norm (number tag)
	 *
	 * @details
	 * 	Perform the AP clustering on a given similarity measure norm. It includes first-stage (and second-staged AP for
	 * 	option) AP clustering, the group labeling, clustering evaluation calculation and representative extraction procedures.
	 *
	 * @param[in] norm The number representing the similarity measure
	 */
	void clusterByNorm(const int& norm);


	/*
	 * @brief Set the labels of clustering technique for all the individual lines
	 *
	 * @details
	 * 	This function will set the group labels in ascending order of cluster size and furthermore calculate the candidate
	 * 	vectors, centroid vectors and cluster size vector
	 *
	 * @param[out] neighborVec The assemble of lines that belong to the same cluster
	 * @param[out] storage The number of assemble candidates for each cluster
	 * @param[out] centroid The center coordinates of each cluster
	 * @param[out] groupTag The cluster label for each candididate individual line
	 */
	void setLabel(vector<vector<int> >& neighborVec, vector<int>& storage, Eigen::MatrixXf& centroid,
			std::vector<int>& groupTag);


	/*
	 * @brief This function is to calculate the normalized entropy for the clustering result
	 *
	 * @param[in] storage The size of each streamline cluster
	 * @pagram[out] EntropyRatio The normalized entropy value to be assigned from calculation
	 */
	void getEntropyRatio(const std::vector<int>& storage, float& EntropyRatio);


	/*
	 * @brief This function is to perform the AP clustering
	 *
	 * @param[out] matrixS The matrix S to be updated
	 * @param[out] matrixR The matrix R to be updated
	 * @param[out] matrixA The matrix A to be updated
	 * @param[out] distMatrix The distance matrix as input
	 * @param[in] coordinates The coordinate of streamlines/pathlines as input
	 */
	void performAPClustering(Eigen::MatrixXf& matrixS, Eigen::MatrixXf& matrixR,
		Eigen::MatrixXf& matrixA, float** distMatrix, const Eigen::MatrixXf& coordinates);


	/*
	 * @brief Calculate the matrix S given a distance matrix and the coordinates of lines
	 *
	 * @param[out] matrixS The matrix S to be assigned values
	 * @param[in] distMatrix The distance matrix
	 * @param[in] coordinates The coordinate matrix for the streamlines/pathlines
	 */
	void getMatrixS(Eigen::MatrixXf& matrixS, float** distMatrix, const Eigen::MatrixXf& coordinates);


	/*
	 * @brief Initialize the matrix S, R and A for AP clustering
	 *
	 * @param[out] matrixS The matrix S to be initialized
	 * @param[out] matrixR The matrix R to be initialized
	 * @param[out] matrixA The matrix A to be initialized
	 * @param[in] rows The row size of the matrix
	 */
	void initializeMatrices(Eigen::MatrixXf& matrixS, Eigen::MatrixXf& matrixR, Eigen::MatrixXf& matrixA,
			const int& rows);


	/*
	 * @brief Update the responsibility matrix R by input matrix A and S
	 *
	 * @param[out] matrixR The matrix R to be updated
	 * @param[in] matrixA The matrix A as input
	 * @param[in] matrixR The matrix S as input
	 */
	void updateResponsibility(Eigen::MatrixXf& matrixR, const Eigen::MatrixXf& matrixA,
							  const Eigen::MatrixXf& matrixS);


	/*
	 * @brief Update the availability matrix A
	 *
	 * @param[out] matrixA the availability matrix A to be updated
	 * @param[in] matrixR the responsibility matrix R as input
	 */
	void updateAvailability(Eigen::MatrixXf& matrixA, const Eigen::MatrixXf& matrixR);


	/*
	 * @brief Get the group label for all the streamlines/pathlines
	 *
	 * @param[in] matrixR The responsibility matrix R as input
	 * @param[in] matrixA The availability matrix A as input
	 * @param[in] matrixS The matrix S as input
	 * @param[out] neighborVec The neighborhood vector for each cluster
	 * @param[out] storage The int vector for recording the size of each cluster
	 * @param[out] groupTag The labels for each streamlines/pathlines
	 */
	void getGroupAssignment(const Eigen::MatrixXf& matrixR, const Eigen::MatrixXf& matrixA,
			  	  	  	  	const Eigen::MatrixXf& matrixS, std::vector<std::vector<int> >& neighborVec,
							std::vector<int>& storage, std::vector<int>& groupTag);

	/*
	 * @brief The function is to get the distance matrix for centroid streamlines/pathlines
	 *
	 * @param[out] centroidDistMatrix The distance matrix for centroid streamlines/pathlines
	 * @param[in] norm The input similarity label
	 * @param[in] centroid The coordinate matrix for the centroids
	 */
	void getDistMatrixForCentroids(float*** centroidDistMatrix, const int& norm, const Eigen::MatrixXf& centroid);


	/*
	 * @brief This function is to read distance matrix from the local file if it exists
	 *
	 * @param[in] norm It is the index for similarity measure
	 */
	void getDistanceMatrixFromFile(const int& norm);


	/*
	 * @brief This function is to calculate the centroids by AP clustering labels
	 *
	 * @param[out] storage The vector of int for the size of each cluster
	 * @param[out] neighborVec The candidate of each cluster
	 * @param[out] centroid The centroid streamline to be updated
	 * @param[out] groupTag The labels for each streamline
	 * @param[in] centroidGroup The number of candidates for each cluster
	 * @param[in] groupSize The size of groups generated
	 */
	void getHierarchicalClusters(std::vector<int>& storage, std::vector<std::vector<int> >& neighborVec,
		Eigen::MatrixXf& centroid, std::vector<int>& group, const std::vector<int>& centroidGroup,
		const int& groupSize);

};


#endif
