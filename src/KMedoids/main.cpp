#include "KMedoids.h"
#include <sys/time.h>

using namespace std;

void featureExtraction(const int& argc,
					   char **argv);

void performKMedoids(const string& fileName, 
					 const std::vector< std::vector<float> >& dataVec,  
					 const int& dimension, 
					 const string& fullName, 
					 const KMedoids& kmedoid,
					 const int& normOption,
					 float& entropy,
					 Silhouette& sil);

void recordInitilization(const Parameter& pm,
						 const int& sampleOption);

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
	const int& dimension = atoi(argv[2]);

	int numOfClusters, vertexCount;
	std::cout << "Please input a cluster number (>=2):" << std::endl;
	std::cin >> numOfClusters;


/*-------------------------------------Input parameter choice-------------------------*/
	Parameter pm;

	std::cout << "Please choose initialization option for seeds:" << std::endl
			  << "1.chose random positions, 2.Chose from samples, 3.k-means++ sampling" << endl;
	std::cin >> pm.initialization;
	assert(pm.initialization==1||pm.initialization==2||pm.initialization==3);

	int sampleOption;
	std::cout << "Please choose sample strategy option for centroids: " << std::endl
			  << "1.from sample, 2.by iterations" << std::endl;
    std::cin >> sampleOption;
    assert(sampleOption==1||sampleOption==2);
    if(sampleOption==1)
    	pm.isSample = true;
    else if(sampleOption==2)
    	pm.isSample = false;

    std::cout << "choose a sampling method for the dataset?" << std::endl
	    	  << "1.directly filling with last vertex; 2. uniform sampling." << std::endl;
	std::cin >> sampleOption;
	assert(sampleOption==1||sampleOption==2);
/*-------------------------------------Finish parameter choice-------------------------*/

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
	std::vector<float> averageS;

	if(sampleOption==1)
		IOHandler::expandArray(data, dataVec, dimension, maxElements);
	else if(sampleOption==2)
		IOHandler::sampleArray(data, dataVec, dimension, maxElements);

	/*
	    0: Euclidean Norm
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
		11: cosine similarity cos(x*y/|x||y|)/pi
	*/

	KMedoids kmedoid(pm, data, numOfClusters);

	for(int i = 0;i<12;i++)
	{
		gettimeofday(&start, NULL);
		ss << strName << "_KMeans";
		performKMedoids(ss.str(), dataVec, dimension, 
					   fullName, kmedoid, i, entropy, silhou);
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

	IOHandler::writeReadme(timeName, timeDiff, kmedoid.getNumOfClusters(), entropyVec);

	IOHandler::writeReadme("Average Silhouette value ", averageS);

	recordInitilization(pm, sampleOption);
}


void performKMedoids(const string& fileName, 
					 const std::vector< std::vector<float> >& dataVec,  
					 const int& dimension, 
					 const string& fullName, 
					 const KMedoids& kmedoid, 
					 const int& normOption,
					 float& entropy,
					 Silhouette& sil)
{
	FeatureLine fl(dataVec);
	kmedoid.getMedoids(fl, normOption, entropy, sil);

	std::vector<std::vector<float> > closestStreamline, furthestStreamline;
	std::vector<int> closestCluster, furthestCluster, meanCluster;
	int closestPoint, furthestPoint;
	IOHandler::assignVec(closestStreamline, closestCluster, fl.closest, 
						 closestPoint, dataVec);
	IOHandler::assignVec(furthestStreamline, furthestCluster, fl.furthest, 
						 furthestPoint, dataVec);
	IOHandler::assignVec(meanCluster, fl.centerMass);
	IOHandler::printVTK(fileName+string("_norm")+to_string(normOption)+string("_mean.vtk"), 
						fl.centerMass, 
						fl.centerMass.size()*fl.centerMass[0].minCenter.size()/dimension, 
						dimension, sil.sCluster);
	IOHandler::printVTK(fileName+"_norm"+to_string(normOption)+"_closest.vtk", 
						closestStreamline, closestPoint/dimension, dimension, 
						closestCluster, sil.sCluster);
	IOHandler::printVTK(fileName+"_norm"+to_string(normOption)+"_furthest.vtk", 
						furthestStreamline, furthestPoint/dimension, 
						dimension, furthestCluster, sil.sCluster);
	std::cout << "Finish printing vtk for k-means clustering result!" << std::endl;

	IOHandler::printToFull(dataVec, fl.group, fl.totalNum, string("norm")+to_string(normOption)
						   +string("_KMeans"), fullName, dimension);
	IOHandler::writeReadme(fl.closest, fl.furthest, normOption);

	IOHandler::printToFull(dataVec, sil.sData, "norm"+to_string(normOption)+"_SValueLine", 
						   fullName, 3);
	IOHandler::printToFull(dataVec, fl.group, sil.sCluster, 
		      "norm"+to_string(normOption)+"_SValueCluster", fullName, 3);
}


void recordInitilization(const Parameter& pm,
						 const int& sampleOption)
{
	std::ofstream readme("../dataset/README", ios::out | ios::app);
	if(readme.fail())
	{
		std::cout << "cannot create README file!" << std::endl;
		exit(1);
	}

	readme << std::endl;
	readme << "Initial centroid is: ";
	if(pm.initialization==1)
		readme << pm.initialization << ".random initialization" 
			   << std::endl;
	else if(pm.initialization==2)
		readme << pm.initialization << ".sample initialization" 
			   << std::endl;
	else if(pm.initialization==3)
		readme << pm.initialization << ".kmeans++ initialization" 
			   << std::endl;
    
    readme << "Medoid is: ";
	if(pm.isSample)
		readme << pm.isSample << ".inside samples" << std::endl;
	else
		readme << pm.isSample << ".from iterations" << std::endl;

	readme << "Sampling is: ";
	if(sampleOption==1)
		readme << sampleOption << ".directly filling" << std::endl;
	else if(sampleOption==2)
		readme << sampleOption << ".uniformly sampling" << std::endl;

	readme.close();

}