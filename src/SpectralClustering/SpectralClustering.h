/* Spectral clustering is a graph-based technique to map original streamlines to a spectral embedding space.
 * The problem itself is a NP-hard and we used instead a relaxed versions with Graph Laplacians.
 * Detailed procedures can be referred at https://tarekmamdouh.wordpress.com/2014/09/28/spectral-clustering/
 * and TVCG paper http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=6702500
 */


#ifndef _SPECTRAL_CLUSTERING_H_
#define _SPECTRAL_CLUSTERING_H_

#include "Predefined.h"
#include "Evrot.h"
#include <unordered_set>
#include <map>
#include <string>

/* local scaling for Gaussian kernel size */
#define SCALING 5


/* update date size for gradient descent */
#ifndef GradientStep
#define GradientStep 0.3
#endif

class SpectralClustering
{

public:

/* default constructor */
	SpectralClustering();

/* argument constructor with argc and argv */
	SpectralClustering(const int& argc, char **argv);

/* destructor */
	~SpectralClustering();

/* perform clustering function */
	void performClustering();

private:

/**********************************************************************************************************
 **************************************   Private member variables   **************************************
 **********************************************************************************************************/

/* metric preparation object to be stored ahead of time */
	MetricPreparation object;

/* input norm option */
	int normOption;

/* group information */
	std::vector<int> group;

/* activityList vector to store event */
	std::vector<string> activityList;

/* timeList vector to store time information */
	std::vector<string> timeList;

/* store dataset information */
	DataSet ds;

/* how many clusters to be needed */
	int numberOfClusters;

/* k-means initialization option */
	int initializationOption;

/* distance range vector */
	std::vector<float> distRange;

/* Gaussian kernel radius for generating adjacency matrix */
	std::vector<float> sigmaVec;

/* Laplacian option, 1: Unnormalized Laplacian, 2: normalized Laplacian, 3: Random Walk Laplacian */
	int LaplacianOption;

/* what kind of 5-th neighbor point would be obtained? */
	bool isDistSorted;

/* what kind of post-processing is to be chosen */
	int postProcessing;

/**********************************************************************************************************
 **************************************   Private member functions   **************************************
 **********************************************************************************************************/

/* extract features from datasets as representative curves */
	void extractFeatures(const std::vector<int>& storage, const std::vector<std::vector<int> >& neighborVec,
            			 const Eigen::MatrixXf& centroid);

/* set dataset from user command */
	void setDataset(const int& argc, char **argv);

/* set norm option, must be within 0-12 */
	void setNormOption();

/* perform group-labeling information */
	void setLabel(vector<vector<int> >& neighborVec, vector<int>& storage, Eigen::MatrixXf& centroid);

/* get weighted adjacency matrix by Gaussian kernel */
	void getAdjacencyMatrix(Eigen::MatrixXf& adjacencyMatrix);

/* get degree matrix */
	void getDegreeMatrix(const Eigen::MatrixXf& adjacencyMatrix, Eigen::DiagonalMatrix<float,Dynamic>& degreeMatrix);

/* get Laplacian matrix */
	void getLaplacianMatrix(const Eigen::MatrixXf& adjacencyMatrix, Eigen::DiagonalMatrix<float,Dynamic>& degreeMatrix,
							Eigen::MatrixXf& laplacianMatrix);

/* decide optimal cluster number by eigenvectors of Laplacian matrix */
	void getEigenClustering(const Eigen::MatrixXf& laplacianMatrix);

/* get local scaling from NIPS 2002 paper */
	void getSigmaList();

/* get entropy ratio */
	void getEntropyRatio(const std::vector<int>& storage, float& EntropyRatio);

/* url: https://www.cs.cmu.edu/~aarti/Class/10701/readings/Luxburg06_TR.pdf */
/********************************** Perform k-means clustering *********************************************/
	/* normalize each row first */
	void normalizeEigenvec(Eigen::MatrixXf& eigenVec);

	/* perform k-means clustering */
	void performKMeans(const Eigen::MatrixXf& eigenVec,
					   std::vector<int>& storage,
					   std::vector<std::vector<int> >& neighborVec);


/* http://lihi.eew.technion.ac.il/files/Publications/SelfTuningClustering.pdf */
/********************************** Vector Rotation *********************************************/
	/* lexicographical list of (i,j) */
	vector<std::pair<int,int> > lexicogList;

	/* theta list for gradient descent method */
	vector<float> thetaList;

	/* generate lexicogList */
	void setLexicogList();

	/* generate theta list */
	void setThetaList();

	/* get Givens rotation matrix given k and theta */
	Eigen::MatrixXf getGivensRotation(const int& k, const float& theta);

	/* get Matrix U */
	Eigen::MatrixXf getMatrixU(const int& a, const int& b);

	/* get Matrix V */
	Eigen::MatrixXf getMatrixV(const int& k);

	/* get Matrix A */
	Eigen::MatrixXf getMatrixA(const int& k, const Eigen::MatrixXf& X);

	/* get dJ/d_theta the gradient */
	const float getGradientToTheta(const int& k,  const Eigen::MatrixXf& X);


/********************************** Vector Rotation from library *********************************************
 ********************************** from library https://github.com/pthimon/clustering ***********************/
	float mMaxQuality = 0;
	int mMethod;

	/* get cluster information based on eigenvector rotation */
	void getEigvecRotation(std::vector<int>& storage, std::vector<std::vector<int> >& neighborVec,
			               Eigen::MatrixXf& clusterCenter, const Eigen::MatrixXf& X);

	/* get derivate method as input, 1. numerical derivative, 2. true derivative */
	void getDerivateMethod();
};

void getMatrixPow(Eigen::DiagonalMatrix<float,Dynamic>& matrix, const float& powNumber);


#endif
