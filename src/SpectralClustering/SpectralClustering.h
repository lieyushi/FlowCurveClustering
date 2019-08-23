/*
 * @brief The SpectralClustering class to perform the spectral clustering on input data set
 *
 * @detais
 * 	Spectral clustering is a graph-based technique to map original streamlines to a spectral embedding space.
 * The problem itself is a NP-hard and we used instead a relaxed versions with Graph Laplacians.
 *
 *  Detailed procedures can be referred at https://tarekmamdouh.wordpress.com/2014/09/28/spectral-clustering/
 * and TVCG paper http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=6702500
 * local scaling for Gaussian kernel size might be defined as 0.05*totalCount as in
 * Blood Flow Clustering and Applications in Virtual Stenting of Intracranial Aneurysms
 */


#ifndef _SPECTRAL_CLUSTERING_H_
#define _SPECTRAL_CLUSTERING_H_

#include "Predefined.h"
#include "Evrot.h"
#include "ValidityMeasurement.h"
#include <unordered_set>
#include <map>
#include <queue>
#include <string>


/*
 * @brief update date size for gradient descent
 */
#ifndef GradientStep
	#define GradientStep 0.3
#endif


/*
 * @brief The Parameter struct to enable parameter tuning for spectral clustering algorithm
 */
struct Para
{

	/*
	 * @brief 1.directly filling with last vertex; 2. uniform sampling, 3. equal-arc sampling
	 */
	int sampled;

	/*
	 * @brief Laplacian option: 1.Normalized Laplacian, 2.Unsymmetric Laplacian
	 */
	int LaplacianOption;

	/*
	 * @brief local scaling by sorted distance: true, false
	 */
	bool isDistSorted;

	/*
	 * @brief post-processing method: 1.k-means, 2.eigenvector rotation
	 */
	int postProcessing;

	/*
	 * @brief derivative method for eigen rotation: 1.numerical derivative, 2.true derivative
	 */
	int mMethod;

	/*
	 * @brief extraction option, 1. centroid, closest and furthest, 2. median, 3. statistical representation
	 */
	int extractOption;
};


/*
 * @brief The spectral clustering class for streamline/pathline clustering
 */
class SpectralClustering
{

public:

	/*
	 * @brief default constructor
	 */
	SpectralClustering();

	/*
	 * @brief The argument constructor with argc and argv
	 *
	 * @param[in] argc The count of arguments
	 * @param[in] argv The char* array of arguments
	 * @param[in] p The Para object
	 * @param[out] automatic The bool flag
	 */
	SpectralClustering(const int& argc, char **argv, const Para& p, bool& automatic);


	/*
	 * @brief destructor
	 */
	~SpectralClustering();


	/*
	 * @brief perform spectral clustering function
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
	 * @brief k-means initialization option
	 */
	int initializationOption = -1;

	/*
	 * @brief distance range vector
	 */
	std::vector<float> distRange;

	/*
	 * @brief Gaussian kernel radius for generating adjacency matrix
	 */
	std::vector<float> sigmaVec;

	/*
	 * @brief Laplacian option, 1: Unnormalized Laplacian, 2: normalized Laplacian, 3: Random Walk Laplacian
	 */
	int LaplacianOption = -1;

	/*
	 * @brief what kind of 5-th neighbor point would be obtained?
	 */
	bool isDistSorted = -1;

	/*
	 * @brief what kind of post-processing is to be chosen
	 */
	int postProcessing = -1;

	/*
	 * @brief extraction option, 1. centroid, closest and furthest (this is what I implemented), 2. median
	 */
	int extractOption = -1;

	/*
	 * @brief scaling factor for spectral clustering to decide Gaussian kernel size
	 */
	int SCALING;

	/*
	 * @brief whether used to find optimal number of clustering
	 */
	bool isOptimal;

	/*
	 * @brief preset number
	 */
	int presetNumber;

	/*
	 * @brief whether read cluster
	 */
	bool readCluster;

	/*
	 * @brief whether it is a pathlines
	 */
	bool isPathlines;


	/*
	 * @brief extract features from datasets as representative curves and calculate the clustering evaluation
	 *
	 * @param[in] storage The size of different clusters
	 * @param[in] neighborVec The candidates included in each cluster
	 * @param[in] centroid The centroid coordinates of cluster
	 */
	void extractFeatures(const std::vector<int>& storage, const std::vector<std::vector<int> >& neighborVec,
            			 const Eigen::MatrixXf& centroid);

	/*
	 * @brief set data set from user command
	 *
	 * @param[in] argc The count of arguments
	 * @param[in] argv The char* array of argument string
	 */
	void setDataset(const int& argc, char **argv);


	/*
	 * @brief set necessary parameter
	 */
	void getParameterUserInput();


	/*
	 * @brief set automatic parameter from the Para object input
	 *
	 * @param[in] p The Para object for the parameters
	 */
	void setParameterAutomatic(const Para& p);


	/*
	 * @brief run spectral clustering based on different norm input
	 *
	 * @param[in] norm The norm option
	 */
	void clusterByNorm(const int& norm);


	/*
	 * @brief perform group-labeling information for all the streamlines
	 *
	 * @param[in] neighborVec The neighboring vector of candidates belonging to the clusters
	 * @param[out] storage The individual size of clusters
	 * @param[out] centroid The centroid coordinates of the clusters
	 */
	void setLabel(vector<vector<int> >& neighborVec, vector<int>& storage, Eigen::MatrixXf& centroid);


	/*
	 * @brief get weighted adjacency matrix by Gaussian kernel
	 *
	 * @param[out] adjacencyMatrix The weighted adjacency matrix computed from the Gaussian graph algorithm
	 */
	void getAdjacencyMatrix(Eigen::MatrixXf& adjacencyMatrix);


	/*
	 * @brief get degree matrix by the adjacency matrix
	 *
	 * @param[in] adjacencyMatrix The adjacency matrix as input
	 * @param[ou] degreeMatrix The degree matrix (diagonal matrix) as output
	 */
	void getDegreeMatrix(const Eigen::MatrixXf& adjacencyMatrix, Eigen::DiagonalMatrix<float,Dynamic>& degreeMatrix);


	/*
	 * @brief get Laplacian matrix
	 *
	 * @param[in] adjacencyMatrix The adjacency matrix as input
	 * @param[out] degreeMatrix The degree matrix to be updated
	 * @param[out] laplacianMatrix The Laplacian matrix to be calculated
	 */
	void getLaplacianMatrix(const Eigen::MatrixXf& adjacencyMatrix, Eigen::DiagonalMatrix<float,Dynamic>& degreeMatrix,
							Eigen::MatrixXf& laplacianMatrix);


	/*
	 * @brief decide optimal cluster number by eigenvectors of Laplacian matrix
	 *
	 * @param[in] laplacianMatrix The input Laplacian matrix for the eigen-rotation minimization
	 * @param[in] norm The norm option
	 */
	void getEigenClustering(const Eigen::MatrixXf& laplacianMatrix, const int& norm);


	/*
	 * @brief get local scaling from NIPS 2002 paper
	 */
	void getSigmaList();


	/*
	 * @brief get entropy ratio from size of clusters as input
	 *
	 * @param[in] storage The size of clusters as input
	 * @param[out] EntropyRatio The entropy ratio to calculate
	 */
	void getEntropyRatio(const std::vector<int>& storage, float& EntropyRatio);


	/*
	 * @brief record the preset number of clusters
	 *
	 * @param[in] number The number of clusters as input
	 */
	void recordPreset(const int& number);


	/*
	 * @brief record the result of optimal number of clusters
	 *
	 * @param[in] normOption The norm option
	 * @param[in] clusNum The number of clusters
	 */
	void recordOptimalResult(const int& normOption, const int& clusNum);


/* url: https://www.cs.cmu.edu/~aarti/Class/10701/readings/Luxburg06_TR.pdf */
/********************************** Perform k-means clustering *********************************************/
	/*
	 * @brief normalize each row first
	 *
	 * @param[out] eigenVec The matrix with eigen-vectors to be updated
	 */
	void normalizeEigenvec(Eigen::MatrixXf& eigenVec);


	/*
	 * @brief perform k-means clustering for the normalized eigen vector matrix
	 *
	 * @param[in] eigenVec The input eigen vector as matrix
	 * @param[out] storage The size of clusters to be updated
	 * @param[out] neighborVec The vector of candidates belonging to clusters
	 */
	void performKMeans(const Eigen::MatrixXf& eigenVec,
					   std::vector<int>& storage,
					   std::vector<std::vector<int> >& neighborVec);


/********************************** Vector Rotation from library *********************************************
 ********************************** from library https://github.com/pthimon/clustering ***********************/

	/*
	 * @brief The max quality initialized value
	 */
	float mMaxQuality = 0;

	/*
	 * @brief The initialized value of mMethod
	 */
	int mMethod = -1;

	/*
	 * @brief get cluster information based on eigenvector rotation
	 *
	 * @param[out] storage The size of clusters
	 * @param[out] neighborVec The candidates belonging to clusters
	 * @param[out] clusterCenter The centroid coordinates of clusters
	 * @param[in] X The matrix X
	 */
	void getEigvecRotation(std::vector<int>& storage, std::vector<std::vector<int> >& neighborVec,
			               Eigen::MatrixXf& clusterCenter, const Eigen::MatrixXf& X);

};


/*
 * @brief Calculate the matrix^powNumber, and it would be inv if powNumber is -1
 *
 * @param[out] matrix The original matrix
 * @param[in] powNumber The exponential index
 */
void getMatrixPow(Eigen::DiagonalMatrix<float,Dynamic>& matrix, const float& powNumber);


#endif
