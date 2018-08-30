#include "KMedoids.h"
#include <sys/time.h>

using namespace std;

bool isPBF;


void featureExtraction(const int& argc,
					   char **argv);

void performKMedoids(const string& fileName, 
					 const std::vector< std::vector<float> >& dataVec,  
					 const int& dimension, 
					 const string& fullName, 
					 const KMedoids& kmedoid,
					 const int& normOption,
					 Silhouette& sil,
					 EvaluationMeasure& measure,
					 TimeRecorder& tr);

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

	std::cout << "It is a PBF dataset? 1. Yes, 0. No." << std::endl;
	int isPBFInput;
	std::cin >> isPBFInput;
	assert(isPBFInput==1||isPBFInput==0);
	isPBF = (isPBFInput==1);

	int numOfClusters, vertexCount;

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

	EvaluationMeasure measure;

	TimeRecorder tr;


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
	tr.eventList.push_back("I-O file reader takes: ");
	tr.timeList.push_back(to_string(timeTemp)+"s");

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
		14: Signature-based measure from http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6231627
		15: Procrustes distance take from http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6787131
	*/

	KMedoids kmedoid(pm, data, numOfClusters);

	for(int i = 0;i<16;i++)
	{
		if(i!=0 && i!=1 && i!=2 && i!=4 && i!=12 && i!=13 && i!=14 && i!=15)
			continue;

		std::cout << "Please input a cluster number (>=2) for norm " << i << " in [2, "
				<< dataVec.size() << "]: " << std::endl;
		std::cin >> numOfClusters;

		gettimeofday(&start, NULL);
		ss << strName << "_KMeans";
		performKMedoids(ss.str(), dataVec, dimension, fullName, kmedoid, i, silhou, measure, tr);
		ss.str("");
		gettimeofday(&end, NULL);
		timeTemp = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
		tr.eventList.push_back("Direct K_Means operation time for norm "+to_string(i)+" takes: ");
		tr.timeList.push_back(to_string(timeTemp)+"s");
		if(silhou.sData.empty())
			silhou.sAverage = 0;

		silhou.reset();
	}

	IOHandler::writeReadme(tr.eventList, tr.timeList, kmedoid.getNumOfClusters());

	recordInitilization(pm, sampleOption);

	/* print silhouette values in the readme file */
	IOHandler::writeReadme("Average Silhouette value is ", measure.silVec);

	/* print gamma statistics vector in the readme */
	IOHandler::writeReadme("Average Gamma statistics value is ", measure.gammaVec);

	/* print entropy value in the readme file */
	IOHandler::writeReadme("Average Entropy value is ", measure.entropyVec);

	/* print DB index values in the readme file */
	IOHandler::writeReadme("Average DB index is ", measure.dbIndexVec);
}


void performKMedoids(const string& fileName, 
					 const std::vector< std::vector<float> >& dataVec,  
					 const int& dimension, 
					 const string& fullName, 
					 const KMedoids& kmedoid, 
					 const int& normOption,
					 Silhouette& sil,
					 EvaluationMeasure& measure,
					 TimeRecorder& tr)
{
	FeatureLine fl(dataVec);
	kmedoid.getMedoids(fl, normOption, sil, measure, tr);

	std::vector<std::vector<float> > closestStreamline, furthestStreamline;
	std::vector<int> closestCluster, furthestCluster, meanCluster;
	int closestPoint, furthestPoint;
	IOHandler::assignVec(closestStreamline, closestCluster, fl.closest, 
						 closestPoint, dataVec);
	IOHandler::assignVec(furthestStreamline, furthestCluster, fl.furthest, 
						 furthestPoint, dataVec);


/* get the average rotation of the extraction */
	std::vector<float> closestRotation, furthestRotation;
	const float& closestAverage = getRotation(closestStreamline, closestRotation);
	const float& furthestAverage = getRotation(furthestStreamline, furthestRotation);

	tr.eventList.push_back("Average rotation of closest for K-medoids clustering on norm "
						   + to_string(normOption) + " is: ");
	tr.timeList.push_back(to_string(closestAverage));

	tr.eventList.push_back("Average rotation of furthest for K-medoids clustering on norm "
						   + to_string(normOption) + " is: ");
	tr.timeList.push_back(to_string(furthestAverage));
/* finish the rotation computation */



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
