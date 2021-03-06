/*
 * @brief The PCA-based clustering and k-means clustering technique on the input integral curves
 * @author Lieyu Shi
 */


#include "PCA_Cluster.h"
#include <sys/time.h>

using namespace std;


/*
 * @brief The PCA and k-means clustering technique to perform the clustering on data sets from the argument
 *
 * @param[in] argc The count of arguments
 * @param[in] argv The argv* array of arguments
 */
void featureExtraction(const int& argc,
					   char **argv);


/*
 * @brief Perform the PCA clustering technique on the integral curves
 *
 * @param[in] fileName The file name of the data set
 * @param[in] dataVec The coordinates of the streamlines
 * @param[in] cluster The number of clusters as input
 * @param[in] dimension The dimension of 2 or 3
 * @param[in] fullName The full name of the primary vtk file
 * @param[in] maxElements The max elements of the streamlines
 * @param[in] data The matrix coordinates after re-sampling
 */
void performPCA_Cluster(const string& fileName, 
						const std::vector< std::vector<float> >& dataVec, 
						const int& cluster, 
						const int& dimension, 
						const string& fullName, 
						const int& maxElements, 
						const Eigen::MatrixXf& data,
						TimeRecorder& tr,
						Silhouette& sil);


/*
 * @brief Perform direct k-means on the input data set
 *
 * @param[in] fileName The file name of the data set
 * @param[in] dataVec The coordinates of the streamlines
 * @param[in] cluster The number of clusters as input
 * @param[in] dimension The dimension of 2 or 3
 * @param[in] fullName The full name of the primary vtk file
 * @param[in] maxElements The max elements of the streamlines
 * @param[in] data The matrix coordinates after re-sampling
 * @param[in] normOption The norm option
 * @param[out] tr The TimeRecorder object
 * @param[out] sil The Silhouette class object
 */
void performK_Means(const string& fileName, 
					const std::vector< std::vector<float> >& dataVec, 
					const int& cluster, 
					const int& dimension, 
					const string& fullName, 
					const int& maxElements, 
					const Eigen::MatrixXf& data,
					const int& normOption,
					TimeRecorder& tr,
					Silhouette& sil);


/*
 * @brief The initialization option for k-means
 */
int initializationOption;

/*
 * @brief The tag whether it is PBF or not
 */
bool isPBF;

/*
 * @brief The post processing strategy
 */
int post_processing;

/*
 * @breif Whether read the number of clusters from local file
 */
bool readCluster;


int main(int argc, char* argv[])
{
	featureExtraction(argc, argv);
	return 0;
}


/*
 * @brief The PCA and k-means clustering technique to perform the clustering on data sets from the argument
 *
 * @param[in] argc The count of arguments
 * @param[in] argv The argv* array of arguments
 */
void featureExtraction(const int& number,
					   char **argv)
{
	while(number!=3)
	{
		std::cout << "Input argument should have 3!" << endl
		          << "./cluster inputFile_name(in dataset folder) "
		          << "data_dimension(3)" << endl;
		exit(1);
	}
	const string& strName = string("../dataset/")+string(argv[1]);
	//const string& strName = "../dataset/pbfDataset";
	const int& dimension = atoi(argv[2]);
	//const string& pbfPath = "/media/lieyu/Seagate Backup Plus Drive/PBF_2013Macklin/pbf_velocitySeparate/source_data/Frame ";

	//std::cout << strName << std::endl;
	//fullName = "../dataset/streamlines_cylinder_9216_full.vtk";
	//const string& strName = "../dataset/streamlines_cylinder_9216";
	//const string& strName = "../dataset/pbf_data";
	
	std::cout << "It is a PBF dataset? 1.Yes, 0.No" << std::endl;
	int PBFjudgement;
	std::cin >> PBFjudgement;
	assert(PBFjudgement==1||PBFjudgement==0);
	isPBF = (PBFjudgement==1);

	/* input for judge whether it is a pathline data set so that MCP can be called */
	bool isPathlines;
	std::cout << "It is pathlines? 1.Yes, 0.No" << std::endl;
	std::cin >> PBFjudgement;
	assert(PBFjudgement==1||PBFjudgement==0);
	isPathlines = (PBFjudgement==1);

	/* set how many clusters and max vertex count of the data set */
	int cluster, vertexCount;

	/* choose k-means initialization method, 2 is often adopted providing better visualization effect */
	std::cout << "Please choose initialization option for seeds:" << std::endl
			  << "1.chose random positions, 2.Chose from samples, 3.k-means++ sampling" << endl;
	std::cin >> initializationOption;
	assert(initializationOption==1 || initializationOption==2 
		   || initializationOption==3);

	int samplingMethod;
	if(isPathlines)
		samplingMethod = 1;

	else
	{
		/* select sampling strategy, and 2 is often for geometric clustering */
		std::cout << "Please choose sampling strategy: " << std::endl
				  << "1.directly filling, 2.uniformly sampling" << std::endl;
		std::cin >> samplingMethod;
	}

    assert(samplingMethod==1 || samplingMethod==2);

    /* whether number of clusters is read from user input or from ../dataset/cluster_number */
    std::cout << "Please choose cluster number method, 0.user input, 1.read clustering: " << std::endl;
    int clusterInput;
    std::cin >> clusterInput;
    assert(clusterInput==0 || clusterInput==1);
    readCluster = (clusterInput==1);

	std::unordered_map<int,int> clusterMap;
	if(readCluster)
	{
		IOHandler::readClusteringNumber(clusterMap, "cluster_number");
	}


    TimeRecorder tr;

	Silhouette sil;

	struct timeval start, end;
	double timeTemp;
	int maxElements;

	gettimeofday(&start, NULL);
	std::vector< std::vector<float> > dataVec;
	IOHandler::readFile(strName, dataVec, vertexCount, dimension, maxElements);
	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec  - start.tv_sec) * 1000000u 
			   + end.tv_usec - start.tv_usec) / 1.e6;
	tr.eventList.push_back("I-O file reader takes: ");
	tr.timeList.push_back(to_string(timeTemp)+"s");

	if(!readCluster)
	{
		std::cout << "Please input a cluster number (>=2) among [2, " << dataVec.size() << "]: " << std::endl;
		std::cin >> cluster;
	}
	else
	{
		cluster = clusterMap[0];
	}

	stringstream ss;
	ss << strName << "_differentNorm_full.vtk";
	const string& fullName = ss.str();
	IOHandler::printVTK(ss.str(), dataVec, vertexCount, dimension);
	ss.str("");

	Eigen::MatrixXf data;

	/* PCA computation is always using brute-force filling arrays by last point */
	IOHandler::expandArray(data, dataVec, dimension, maxElements);

	std::cout << "PCA-based clustering starts..." << std::endl;
	ss << strName << "_PCAClustering";
	gettimeofday(&start, NULL);
	performPCA_Cluster(ss.str(), dataVec, cluster, dimension, fullName, maxElements, data, tr, sil);
	std::cout << "Max element is " << maxElements << std::endl;
	ss.str("");
	ss.clear();
	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec  - start.tv_sec) * 1000000u 
			   + end.tv_usec - start.tv_usec) / 1.e6;

	tr.eventList.push_back("PCA+KMeans takes: ");
	tr.timeList.push_back(std::to_string(timeTemp)+"s");

	IOHandler::writeReadme(tr.eventList, tr.timeList, cluster);

	tr.eventList.clear();
	tr.timeList.clear();


	sil.reset();

	/*  0: Euclidean Norm
		1: Fraction Distance Metric
		2: piece-wise angle average
		3: Bhattacharyya metric for rotation
		4: average rotation
		5: signed-angle intersection
		6: normal-direction multivariate distribution
		7: Bhattacharyya metric with angle to a fixed direction
		8: Piece-wise angle average \times standard deviation
		9: normal-direction multivariate un-normalized distribution
		10: x*y/|x||y| borrowed from machine learning
		11: cosine similarity
		12: Mean-of-closest point distance (MCP)
		13: Hausdorff distance min_max(x_i,y_i)
		14: Signature-based measure taken from http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6231627
		15: Procrustes distance taken from http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6787131
		16: entropy-based distance metric taken from http://vis.cs.ucdavis.edu/papers/pg2011paper.pdf
		17: time-series MCP distance from https://www.sciencedirect.com/science/article/pii/S0097849318300128
			for pathlines only
	*/
	if(samplingMethod==2)
		IOHandler::sampleArray(data, dataVec, dimension, maxElements);

	for(int i = 0;i<=17;i++)
	{
		if(isPathlines)
		{
			if(i!=0 && i!=1 && i!=2 && i!=4 && i!=12 && i!=13 && i!=14 && i!=15 && i!=17)
				continue;
		}

		else
		{
			if(i!=0 && i!=1 && i!=2 && i!=4 && i!=12 && i!=13 && i!=14 && i!=15)
				continue;
		}

		if(readCluster)
			cluster = clusterMap[i];
		else
		{
			std::cout << "Please input a cluster number (>=2) for norm " << i << " in [2, "
					<< dataVec.size() << "]: " << std::endl;
			std::cin >> cluster;
		}

		std::cout << "Kmeans on norm " << i << " starts..." << std::endl;
		gettimeofday(&start, NULL);
		ss << strName << "_KMeans";
		performK_Means(ss.str(), dataVec, cluster, dimension, fullName, maxElements, data,i, tr, sil);

		ss.str("");
		gettimeofday(&end, NULL);
		timeTemp = ((end.tv_sec  - start.tv_sec) * 1000000u 
					+ end.tv_usec - start.tv_usec) / 1.e6;
		tr.eventList.push_back("K-means on norm "+to_string(i)+" takes: ");
		tr.timeList.push_back(to_string(timeTemp)+"s");

		IOHandler::writeReadme(tr.eventList, tr.timeList, cluster);

		tr.eventList.clear();
		tr.timeList.clear();
		sil.reset();
	}
}


/*
 * @brief Perform the PCA clustering technique on the integral curves
 *
 * @param[in] fileName The file name of the data set
 * @param[in] dataVec The coordinates of the streamlines
 * @param[in] cluster The number of clusters as input
 * @param[in] dimension The dimension of 2 or 3
 * @param[in] fullName The full name of the primary vtk file
 * @param[in] maxElements The max elements of the streamlines
 * @param[in] data The matrix coordinates after re-sampling
 */
void performPCA_Cluster(const string& fileName, 
					    const std::vector< std::vector<float> >& dataVec, 
					    const int& cluster, 
						const int& dimension, 
						const string& fullName, 
						const int& maxElements, 
						const Eigen::MatrixXf& data,
						TimeRecorder& tr,
						Silhouette& sil)
{

	std::vector<MeanLine> centerMass;
	std::vector<int> group(dataVec.size());
	std::vector<ExtractedLine> closest;
	std::vector<ExtractedLine> furthest;
	std::vector<int> totalNum(dataVec.size());

	// choose an appropriate post processing technique for PCA rank space
	std::cout << "Please select a post-processing: 1. k-means, 2. AHC-average." << std::endl;
	std::cin >> post_processing;
	assert(post_processing==1 || post_processing==2);

	PCA_Cluster::performPCA_Clustering(data, dataVec.size(), maxElements, centerMass,
			           group, totalNum, closest, furthest, cluster, tr, sil);

	std::vector<std::vector<float> > closestStreamline;
	std::vector<std::vector<float> > furthestStreamline;
	std::vector<int> closestCluster, furthestCluster, meanCluster;
	int closestPoint, furthestPoint;

	IOHandler::assignVec(closestStreamline, closestCluster, closest, closestPoint, dataVec);

	IOHandler::assignVec(furthestStreamline, furthestCluster, furthest, 
						 furthestPoint, dataVec);

/* get the average rotation of the extraction */
	std::vector<float> closestRotation, furthestRotation;
	const float& closestAverage = getRotation(closestStreamline, closestRotation);
	const float& furthestAverage = getRotation(furthestStreamline, furthestRotation);

	tr.eventList.push_back("Average rotation of closest for PCA clustering is: ");
	tr.timeList.push_back(to_string(closestAverage));

	tr.eventList.push_back("Average rotation of furthest for PCA clustering is: ");
	tr.timeList.push_back(to_string(furthestAverage));
/* finish the rotation computation */

	IOHandler::assignVec(meanCluster, centerMass);

	IOHandler::printVTK(fileName+string("_PCA_closest.vtk"), closestStreamline, 
						closestPoint/dimension, dimension, closestCluster,
						sil.sCluster);

	IOHandler::printVTK(fileName+string("_PCA_furthest.vtk"), furthestStreamline, 
						furthestPoint/dimension, dimension, furthestCluster,
						sil.sCluster);

	IOHandler::printVTK(fileName+string("_PCA_mean.vtk"), centerMass, 
						centerMass.size()*centerMass[0].minCenter.size()/dimension, 
						dimension, sil.sCluster);

	std::cout << "Finish printing vtk for pca-clustering result!" << std::endl;

	if(post_processing==1)
		IOHandler::printToFull(dataVec, group, totalNum, string("PCA_KMeans"), fullName, dimension);
	else if(post_processing==2)
		IOHandler::printToFull(dataVec, group, totalNum, string("PCA_AHC"), fullName, dimension);

	//IOHandler::writeReadme(closest, furthest);

	IOHandler::printToFull(dataVec, sil.sData, "PCA_SValueLine", fullName, 3);

	IOHandler::printToFull(dataVec, group, sil.sCluster, "PCA_SValueCluster", 
						   fullName, 3);
}


/*
 * @brief Perform direct k-means on the input data set
 *
 * @param[in] fileName The file name of the data set
 * @param[in] dataVec The coordinates of the streamlines
 * @param[in] cluster The number of clusters as input
 * @param[in] dimension The dimension of 2 or 3
 * @param[in] fullName The full name of the primary vtk file
 * @param[in] maxElements The max elements of the streamlines
 * @param[in] data The matrix coordinates after re-sampling
 * @param[in] normOption The norm option
 * @param[out] tr The TimeRecorder object
 * @param[out] sil The Silhouette class object
 */
void performK_Means(const string& fileName, 
					const std::vector< std::vector<float> >& dataVec, 
					const int& cluster, 
					const int& dimension, 
					const string& fullName, 
					const int& maxElements, 
					const Eigen::MatrixXf& data, 
					const int& normOption,
					TimeRecorder& tr,
					Silhouette& sil)
{
	std::vector<MeanLine> centerMass;
	std::vector<ExtractedLine> closest;
	std::vector<ExtractedLine> furthest;
	std::vector<int> group(dataVec.size());
	std::vector<int> totalNum(dataVec.size());
	PCA_Cluster::performDirectK_Means(data, dataVec.size(), maxElements, 
									  centerMass, group, totalNum, 
									  closest, furthest, cluster, normOption, tr, sil);

	std::vector<std::vector<float> > closestStreamline, furthestStreamline;
	std::vector<int> closestCluster, furthestCluster, meanCluster;
	int closestPoint, furthestPoint;
	IOHandler::assignVec(closestStreamline, closestCluster, closest, closestPoint, dataVec);
	IOHandler::assignVec(furthestStreamline, furthestCluster, furthest, 
						 furthestPoint, dataVec);

/* get the average rotation of the extraction */
	std::vector<float> closestRotation, furthestRotation;
	const float& closestAverage = getRotation(closestStreamline, closestRotation);
	const float& furthestAverage = getRotation(furthestStreamline, furthestRotation);

	tr.eventList.push_back("Average rotation of closest for K-means clustering on norm "
			               + to_string(normOption) + " is: ");
	tr.timeList.push_back(to_string(closestAverage));

	tr.eventList.push_back("Average rotation of furthest for K-means clustering on norm "
			               + to_string(normOption) + " is: ");
	tr.timeList.push_back(to_string(furthestAverage));
/* finish the rotation computation */


	IOHandler::assignVec(meanCluster, centerMass);
	IOHandler::printVTK(fileName+string("_norm")+to_string(normOption)+string("_mean.vtk"), 
						centerMass, 
						centerMass.size()*centerMass[0].minCenter.size()/dimension, 
						dimension, sil.sCluster);
	IOHandler::printVTK(fileName+"_norm"+to_string(normOption)+"_closest.vtk", 
						closestStreamline, closestPoint/dimension, dimension, 
						closestCluster, sil.sCluster);
	IOHandler::printVTK(fileName+"_norm"+to_string(normOption)+"_furthest.vtk", 
						furthestStreamline, furthestPoint/dimension, 
						dimension, furthestCluster, sil.sCluster);
	std::cout << "Finish printing vtk for k-means clustering result!" << std::endl;

	IOHandler::printToFull(dataVec, group, totalNum, string("norm")+to_string(normOption)
						   +string("_KMeans"), fullName, dimension);

	//IOHandler::writeReadme(closest, furthest, normOption);

	IOHandler::printToFull(dataVec, sil.sData, "norm"+to_string(normOption)+"_SValueLine", 
						   fullName, 3);
	IOHandler::printToFull(dataVec, group, sil.sCluster, "norm"+to_string(normOption)+"_SValueCluster", fullName, 3);

	centerMass.clear();
	closest.clear();
	furthest.clear();
	group.clear();
	totalNum.clear();
}


