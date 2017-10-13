#include "DensityClustering.h"

std::vector<string> activityList;
std::vector<string> timeList;

float multiTimes;
int minPts;

DensityClustering::DensityClustering(const int& argc,
									 char **argv)
{
	setDataset(argc, argv);
	setNormOption();

	object = MetricPreparation(ds.dataMatrix.rows(), ds.dataMatrix.cols());
	object.preprocessing(ds.dataMatrix, ds.dataMatrix.rows(), ds.dataMatrix.cols(), normOption);
	
	nodeVec = vector<PointNode>(ds.dataMatrix.rows(),PointNode());
}


DensityClustering::~DensityClustering()
{

}


void DensityClustering::performClustering()
{
	float minDist, maxDist;
	getDistRange(minDist, maxDist);
	std::cout << "Distance range is [" << minDist << ", "
			  << maxDist << "]." << std::endl;
	minPts = setMinPts();
	multiTimes = setTimesMin(minDist, maxDist);

	struct timeval start, end;
	double timeTemp;
	gettimeofday(&start, NULL);

	DBSCAN(maxDist*multiTimes, minPts);

	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec  - start.tv_sec) * 1000000u 
			   + end.tv_usec - start.tv_usec) / 1.e6;
	activityList.push_back("DBSCAN clustering takes: ");
	timeList.push_back(to_string(timeTemp)+" s");

	extractFeatures(minDist*multiTimes, minPts);
}


void DensityClustering::DBSCAN(const float& radius_eps,
							   const int& minPts)
{
	int C = 0;
	for(int i=0;i<ds.dataMatrix.rows();++i)
	{
		if(nodeVec[i].visited)
			continue;
		nodeVec[i].visited = true;
		vector<int> neighbor = regionQuery(i,radius_eps);
		if(neighbor.size()<minPts)
			nodeVec[i].type = NOISE;
		else
		{
			expandCluster(i, neighbor, C, radius_eps, minPts);
			++C;
		}
	}
}


void DensityClustering::expandCluster(const int& index,
									  vector<int>& neighbor,
									  const int& cluster_id,
									  const float& radius_eps,
									  const int& minPts)
{
	nodeVec[index].group = cluster_id;
	int insideElement;
	for (int i = 0; i < neighbor.size(); ++i)
	{
		insideElement = neighbor[i];
		if(!nodeVec[insideElement].visited)
		{
			nodeVec[insideElement].visited = true;
			vector<int> newNeighbor = regionQuery(insideElement,radius_eps);
			if(newNeighbor.size()>=minPts)
			{
				neighbor.insert(neighbor.end(), newNeighbor.begin(), newNeighbor.end());
			}
		}
		if(nodeVec[insideElement].group==-1)
			nodeVec[insideElement].group = cluster_id;
	}
}


const vector<int> DensityClustering::regionQuery(const int& index,
												 const float& radius_eps)
{
	vector<int> neighborArray;
	neighborArray.push_back(index);
	float tempDist;
	for (int i = 0; i < ds.dataMatrix.rows(); ++i)
	{
		if(i==index)
			continue;
		tempDist=getDisimilarity(ds.dataMatrix.row(index),
				 ds.dataMatrix.row(i),index,i,normOption, object);
		if(tempDist<=radius_eps)
			neighborArray.push_back(i);
	}
	return neighborArray;
}


void DensityClustering::setDataset(const int& argc,
								   char **argv)
{
	if(argc!=3)
	{
		std::cout << "Input argument should have 3!" << endl
		          << "./cluster inputFile_name(in dataset folder) "
		          << "data_dimension(3)" << endl;
		exit(1);
	}
	ds.strName = string("../dataset/")+string(argv[1]);
	ds.dimension = atoi(argv[2]);

	int sampleOption;
    std::cout << "choose a sampling method for the dataset?" << std::endl
	    	  << "1.directly filling with last vertex; 2. uniform sampling." << std::endl;
	std::cin >> sampleOption;
	assert(sampleOption==1||sampleOption==2);

	IOHandler::readFile(ds.strName,ds.dataVec,ds.vertexCount,ds.dimension,ds.maxElements);

	ds.fullName = ds.strName+"_differentNorm_full.vtk";
	IOHandler::printVTK(ds.fullName, ds.dataVec, ds.vertexCount, ds.dimension);

	if(sampleOption==1)
		IOHandler::expandArray(ds.dataMatrix,ds.dataVec,ds.dimension,ds.maxElements);
	else if(sampleOption==2)
		IOHandler::sampleArray(ds.dataMatrix,ds.dataVec,ds.dimension,ds.maxElements);
}


void DensityClustering::setNormOption()
{
	std::cout << "Choose a norm from 0-11!" << std::endl;
	std::cin >> normOption;
	std::cout << std::endl;
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
		11: cosine similarity
	*/
	bool found = false;
	for (int i = 0; i < 12&&!found; ++i)
	{
		if(normOption==i)
		{
			found = true;
			break;
		}
	}
	if(!found)
	{
		std::cout << "Cannot find the norm!" << std::endl;
		exit(1);
	}
}


void DensityClustering::getDistRange(float& minDist, 
					                 float& maxDist)
{
	const float& Percentage = 0.1;
	const int& Rows = ds.dataMatrix.rows();
	const int& chosen = int(Percentage*Rows);
	minDist = FLT_MAX;
	maxDist = -1.0;
#pragma omp parallel num_threads(8)
	{
	#pragma omp for nowait
		for (int i = 0; i < chosen; ++i)
		{
			for (int j = 0; j < Rows; ++j)
			{
				if(i==j)
					continue;
				float dist = getDisimilarity(ds.dataMatrix.row(i),
					  ds.dataMatrix.row(j),i,j,normOption,object);
			#pragma omp critical 
				{
					if(dist<minDist)
						minDist=dist;
					if(dist>maxDist)
						maxDist=dist;
				}
			}
		}	
	}
	std::cout << minDist << " " << maxDist << std::endl;
}


const int DensityClustering::setMinPts()
{
	std::cout << std::endl;
	std::cout << "Input the minPts for DBSCAN in [0" << ", "
			  << ds.dataMatrix.rows() << "]:" << std::endl;
	int minPts;
	std::cin >> minPts;
	if(minPts<=0 || minPts>=ds.dataMatrix.rows())
	{
		std::cout << "Error for out-of-range minPts!" << std::endl;
		exit(1);
	}
	return minPts;
}


const float DensityClustering::setTimesMin(const float& minDist, 
					  					   const float& maxDist)
{
	std::cout << std::endl;
	float lowerBound = minDist/maxDist;
	std::cout << "Input the multiplication for DBSCAN radius in ["
	          << lowerBound << ",1.0]:" << std::endl;
	float multiTimes;
	std::cin >> multiTimes;
	if(multiTimes>=1.0 || multiTimes<=lowerBound)
	{
		std::cout << "Error for out-of-range minPts!" << std::endl;
		exit(1);
	}
	return multiTimes;
}


void DensityClustering::extractFeatures(const float& radius_eps,
							   			const int& minPts)
{
	int maxGroup = -INT_MAX+1;
#pragma omp parallel num_threads(8)
	{
	#pragma omp for nowait
		for (int i = 0; i < nodeVec.size(); ++i)
		{			
			int groupID = nodeVec[i].group;
		#pragma omp critical
			{ 
				if(groupID!=-1 && groupID>maxGroup)
					maxGroup = groupID;
			}
		}	
	}	
	std::cout << "Max group is: " << maxGroup << std::endl;

/* re-index the group id by increasing number */
	int numClusters = maxGroup+1;
	std::vector<int> container(numClusters,0);
	for (int i = 0; i < nodeVec.size(); ++i)
	{
		if(nodeVec[i].group!=-1)
			++container[nodeVec[i].group];
	}

	int increasingOrder[numClusters];
	std::multimap<int,int> groupMap;

	for (int i = 0; i < numClusters; ++i)
		groupMap.insert(std::pair<int,int>(container[i],i));

	std::fill(container.begin(), container.end(), 0);
	int groupNo = 0;
	for (std::multimap<int,int>::iterator it=groupMap.begin();it!=groupMap.end();++it)
	{
		if(it->first>0)
		{
			increasingOrder[it->second] = groupNo;
			container[groupNo] = it->first;
			++groupNo;
		}
	}

	numClusters = groupNo+1;	/* plus -1 as group */

#pragma omp parallel for schedule(dynamic) num_threads(8)
	for (int i = 0; i < nodeVec.size(); ++i)
	{
		if(nodeVec[i].group!=-1)
			nodeVec[i].group=increasingOrder[nodeVec[i].group];
	}

	/* in case -1, we use 0 to record number of -1 as noise */

	std::vector<int> item_cids(nodeVec.size());
	std::vector<std::vector<int> > storage(numClusters);
	/* -1 group as group[0] */
	for (int i = 0; i < nodeVec.size(); ++i)
	{
		item_cids[i] = nodeVec[i].group;
		storage[nodeVec[i].group+1].push_back(i);
	}

	container.insert(container.begin(),storage[0].size());

	IOHandler::printClustersNoise(ds.dataVec,item_cids,container, 
		 "norm"+to_string(normOption),ds.fullName,ds.dimension);

	struct timeval start, end;
	double timeTemp;

	numClusters-=1;
	gettimeofday(&start, NULL);
	Silhouette sil;
	sil.computeValue(normOption,ds.dataMatrix,ds.dataMatrix.rows(),
		ds.dataMatrix.cols(),item_cids,object,numClusters);
	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec  - start.tv_sec) * 1000000u 
			   + end.tv_usec - start.tv_usec) / 1.e6;
	activityList.push_back("Silhouette calculation takes: ");
	timeList.push_back(to_string(timeTemp)+" s");

	const int& numNoise = storage[0].size();
	storage.erase(storage.begin());
	/* compute the centroid coordinates of each clustered group */

	gettimeofday(&start, NULL);

	Eigen::MatrixXf centroid = MatrixXf::Zero(numClusters,ds.dataMatrix.cols());
	vector<vector<float> > cenVec(numClusters);
#pragma omp parallel for schedule(dynamic) num_threads(8)
	for (int i=0;i<numClusters;++i)
	{
		const std::vector<int>& groupRow = storage[i];
		for (int j = 0; j < groupRow.size(); ++j)
		{
			centroid.row(i)+=ds.dataMatrix.row(groupRow[j]);
		}		
		centroid.row(i)/=groupRow.size();
		const Eigen::VectorXf& vec = centroid.row(i);
		cenVec[i] = vector<float>(vec.data(), vec.data()+ds.dataMatrix.cols());
	}

	vector<vector<float> > closest(numClusters);
	vector<vector<float> > furthest(numClusters);

#pragma omp parallel for schedule(dynamic) num_threads(8)
	for (int i=0;i<numClusters;++i)
	{
		float minDist = FLT_MAX;
		float maxDist = -10;
		int minIndex = -1, maxIndex = -1;
		const std::vector<int>& groupRow = storage[i];
		const Eigen::VectorXf& eachCentroid = centroid.row(i);
		for (int j = 0; j < groupRow.size(); ++j)
		{
			float distance = getDisimilarity(eachCentroid,ds.dataMatrix,groupRow[j],normOption,object);
			if(minDist>distance)
			{
				minDist = distance;
				minIndex = groupRow[j];
			}
			if(maxDist<distance)
			{				
				maxDist = distance;
				maxIndex = groupRow[j];
			}
		}
		closest[i] = ds.dataVec[minIndex];
		furthest[i] = ds.dataVec[maxIndex];
	}

	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec  - start.tv_sec) * 1000000u 
			   + end.tv_usec - start.tv_usec) / 1.e6;
	activityList.push_back("Feature extraction takes: ");
	timeList.push_back(to_string(timeTemp)+" s");


	std::cout << "Finishing extracting features!" << std::endl;	
	IOHandler::printFeature("norm"+to_string(normOption)+"_closest.vtk", 
			closest, sil.sCluster, ds.dimension);
	IOHandler::printFeature("norm"+to_string(normOption)+"_furthest.vtk",
			furthest, sil.sCluster, ds.dimension);
	IOHandler::printFeature("norm"+to_string(normOption)+"_centroid.vtk", 
			cenVec, sil.sCluster,ds.dimension);

	IOHandler::printToFull(ds.dataVec, sil.sData, 
			"norm"+to_string(normOption)+"_SValueLine", ds.fullName, ds.dimension);
	IOHandler::printToFull(ds.dataVec, item_cids, sil.sCluster, 
		      "norm"+to_string(normOption)+"_SValueCluster", ds.fullName, ds.dimension);

	activityList.push_back("Norm option is: ");
	timeList.push_back(to_string(normOption));

	activityList.push_back("numCluster is: ");
	timeList.push_back(to_string(numClusters));

	activityList.push_back("Average Silhouette is: ");
	timeList.push_back(to_string(sil.sAverage));

	activityList.push_back("Noise number is: ");
	timeList.push_back(to_string(numNoise));

	activityList.push_back("radius eps is: ");
	timeList.push_back(to_string(multiTimes));

	activityList.push_back("MinPts is: ");
	timeList.push_back(to_string(minPts));

	IOHandler::generateReadme(activityList,timeList);
}