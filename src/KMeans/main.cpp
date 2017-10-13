#include "PCA_Cluster.h"
#include <sys/time.h>

using namespace std;

void featureExtraction(const int& argc,
					   char **argv);

void performPCA_Cluster(const string& fileName, 
						const std::vector< std::vector<float> >& dataVec, 
						const int& cluster, 
						const int& dimension, 
						const string& fullName, 
						const int& maxElements, 
						const Eigen::MatrixXf& data,
						float& entropy,
						Silhouette& sil);

void performK_Means(const string& fileName, 
					const std::vector< std::vector<float> >& dataVec, 
					const int& cluster, 
					const int& dimension, 
					const string& fullName, 
					const int& maxElements, 
					const Eigen::MatrixXf& data,
					const int& normOption,
					float& entropy,
					Silhouette& sil);

int initializationOption;

int main(int argc, char* argv[])
{
	featureExtraction(argc, argv);
	return 0;
}

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
	
	int cluster, vertexCount;
	std::cout << "Please input a cluster number (>=2):" << std::endl;
	std::cin >> cluster;

	std::cout << "Please choose initialization option for seeds:" << std::endl
			  << "1.chose random positions, 2.Chose from samples, 3.k-means++ sampling" << endl;
	std::cin >> initializationOption;
	assert(initializationOption==1 || initializationOption==2 
		   || initializationOption==3);

	std::cout << "Please choose sampling strategy: " << std::endl
			  << "1.directly filling, 2.uniformly sampling" << std::endl;
    int samplingMethod;
    std::cin >> samplingMethod;
    assert(samplingMethod==1 || samplingMethod==2);

	std::vector<string> timeName;
	std::vector<double> timeDiff;
	std::vector<float> entropyVec;
	float entropy;

	/* a Silhouette method to estimate the clustering effect */
	Silhouette silhou;

	struct timeval start, end;
	double timeTemp;
	int maxElements;

	gettimeofday(&start, NULL);
	std::vector< std::vector<float> > dataVec;
	IOHandler::readFile(strName, dataVec, vertexCount, dimension, maxElements);
	//IOHandler::readFile(pbfPath, dataVec, vertexCount, dimension, 128000, 1500);
	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec  - start.tv_sec) * 1000000u 
			   + end.tv_usec - start.tv_usec) / 1.e6;
	timeName.push_back("I-O file reader time");
	timeDiff.push_back(timeTemp);

	stringstream ss;
	ss << strName << "_differentNorm_full.vtk";
	const string& fullName = ss.str();
	IOHandler::printVTK(ss.str(), dataVec, vertexCount, dimension);
	ss.str("");

	Eigen::MatrixXf data;

	/* PCA computation is always using brute-force filling arrays by last point */
	IOHandler::expandArray(data, dataVec, dimension, maxElements);

	
	ss << strName << "_PCAClustering";
	gettimeofday(&start, NULL);
	//ss << "TACs_PCAClustering";
	performPCA_Cluster(ss.str(), dataVec, cluster, dimension, 
					   fullName, maxElements, data, entropy, silhou);
	entropyVec.push_back(entropy);
	ss.str("");
	ss.clear();
	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec  - start.tv_sec) * 1000000u 
			   + end.tv_usec - start.tv_usec) / 1.e6;
	timeName.push_back("PCA+K_Means operation time");
	timeDiff.push_back(timeTemp);

	std::vector<float> averageS;
	averageS.push_back(silhou.sAverage);

	silhou.reset();

	//IOHandler::writeReadme("PCA Silhouette of each cluster is ", silhou.sCluster);

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
	*/
	if(samplingMethod==2)
		IOHandler::sampleArray(data, dataVec, dimension, maxElements);

	for(int i = 0;i<12;i++)
	{
		//if(i!=7)
		//	continue;
		gettimeofday(&start, NULL);
		ss << strName << "_KMeans";
		performK_Means(ss.str(), dataVec, cluster, dimension, 
					   fullName, maxElements, data,i, entropy, silhou);
		entropyVec.push_back(entropy);
		ss.str("");
		gettimeofday(&end, NULL);
		timeTemp = ((end.tv_sec  - start.tv_sec) * 1000000u 
					+ end.tv_usec - start.tv_usec) / 1.e6;
		timeName.push_back("Direct K_Means operation time for norm "+to_string(i));
		timeDiff.push_back(timeTemp);
		if(silhou.sData.empty())
			silhou.sAverage = 0;
		averageS.push_back(silhou.sAverage);

		silhou.reset();
	}

	IOHandler::writeReadme(timeName, timeDiff, cluster, entropyVec);

	IOHandler::writeReadme("Average Silhouette value ", averageS);
}


void performPCA_Cluster(const string& fileName, 
					    const std::vector< std::vector<float> >& dataVec, 
					    const int& cluster, 
						const int& dimension, 
						const string& fullName, 
						const int& maxElements, 
						const Eigen::MatrixXf& data,
						float& entropy,
						Silhouette& sil)
{

	std::vector<MeanLine> centerMass;
	std::vector<int> group(dataVec.size());
	std::vector<ExtractedLine> closest;
	std::vector<ExtractedLine> furthest;
	std::vector<int> totalNum(dataVec.size());

	PCA_Cluster::performPCA_Clustering(data, dataVec.size(), maxElements, 
						centerMass, group, totalNum, closest, furthest, cluster, entropy,
						sil);

	std::vector<std::vector<float> > closestStreamline;
	std::vector<std::vector<float> > furthestStreamline;
	std::vector<int> closestCluster, furthestCluster, meanCluster;
	int closestPoint, furthestPoint;

	IOHandler::assignVec(closestStreamline, closestCluster, closest, closestPoint, dataVec);

	IOHandler::assignVec(furthestStreamline, furthestCluster, furthest, 
						 furthestPoint, dataVec);

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

	IOHandler::printToFull(dataVec, group, totalNum, string("PCA_KMeans"), 
						   fullName, dimension);

	IOHandler::writeReadme(closest, furthest);

	IOHandler::printToFull(dataVec, sil.sData, "PCA_SValueLine", fullName, 3);

	IOHandler::printToFull(dataVec, group, sil.sCluster, "PCA_SValueCluster", 
						   fullName, 3);
}


void performK_Means(const string& fileName, 
					const std::vector< std::vector<float> >& dataVec, 
					const int& cluster, 
					const int& dimension, 
					const string& fullName, 
					const int& maxElements, 
					const Eigen::MatrixXf& data, 
					const int& normOption,
					float& entropy,
					Silhouette& sil)
{
	std::vector<MeanLine> centerMass;
	std::vector<ExtractedLine> closest;
	std::vector<ExtractedLine> furthest;
	std::vector<int> group(dataVec.size());
	std::vector<int> totalNum(dataVec.size());
	PCA_Cluster::performDirectK_Means(data, dataVec.size(), maxElements, 
									  centerMass, group, totalNum, 
									  closest, furthest, cluster, normOption, entropy,
									  sil);
	std::vector<std::vector<float> > closestStreamline, furthestStreamline;
	std::vector<int> closestCluster, furthestCluster, meanCluster;
	int closestPoint, furthestPoint;
	IOHandler::assignVec(closestStreamline, closestCluster, closest, closestPoint, dataVec);
	IOHandler::assignVec(furthestStreamline, furthestCluster, furthest, 
						 furthestPoint, dataVec);
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
	IOHandler::writeReadme(closest, furthest, normOption);

	IOHandler::printToFull(dataVec, sil.sData, "norm"+to_string(normOption)+"_SValueLine", 
						   fullName, 3);
	IOHandler::printToFull(dataVec, group, sil.sCluster, 
		      "norm"+to_string(normOption)+"_SValueCluster", fullName, 3);

	centerMass.clear();
	closest.clear();
	furthest.clear();
	group.clear();
	totalNum.clear();
}


