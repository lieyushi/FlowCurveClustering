#include "AHC.h"

/* default constructor */
AHC::AHC()
{

}

/* argument constructor with argc and argv */
AHC::AHC(const int& argc, char **argv)
{
	setDataset(argc, argv);
	setNormOption();

	/* very hard to decide whether needed to perform such pre-processing */
	object = MetricPreparation(ds.dataMatrix.rows(), ds.dataMatrix.cols());
	object.preprocessing(ds.dataMatrix, ds.dataMatrix.rows(), ds.dataMatrix.cols(), normOption);

	/* would store distance matrix instead because it would save massive time */
	struct timeval start, end;
	double timeTemp;
	gettimeofday(&start, NULL);

	if(!getDistanceMatrix(ds.dataMatrix, normOption, object))
	{
		std::cout << "Failure to compute distance matrix!" << std::endl;
	}

	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec  - start.tv_sec) * 1000000u 
			   + end.tv_usec - start.tv_usec) / 1.e6;
	activityList.push_back("Distance matrix computing takes: ");
	timeList.push_back(to_string(timeTemp)+" s");

	getDistRange();
}

/* destructor */
AHC::~AHC()
{
	deleteDistanceMatrix(ds.dataMatrix.rows());
}

/* perform clustering function */
void AHC::performClustering()
{
	std::vector<Ensemble> nodeVec;

	/* perform hierarchical clustering */
	std::cout << "---------------------" << std::endl;
	std::cout << "1. clustering by a fixed group, 2. clustering by a distance threshold." << std::endl;
	int clusteringOption;
	std::cin >> clusteringOption;
	assert(clusteringOption==1 || clusteringOption==2);

	if(clusteringOption==1)
		bottomUp_byGroup(nodeVec);
	else if(clusteringOption==2)
		bottomUp_byThreshold(nodeVec);

	vector<vector<int>> neighborVec(numberOfClusters);
	// element size for all groups
	vector<int> storage(numberOfClusters);

	// geometric center
	Eigen::MatrixXf centroid = Eigen::MatrixXf::Zero(numberOfClusters,ds.dataMatrix.cols());


	// set label information
	setLabel(nodeVec, neighborVec, storage, centroid);

	nodeVec.clear();

	extractFeatures(storage, neighborVec, centroid);
}


/* perform hierarchical clustering by given a group */
void AHC::bottomUp_byGroup(std::vector<Ensemble>& nodeVec)
{
	const int& Row = ds.dataMatrix.rows();
	std::cout << "-------------------------------------------------------------------------------" << std::endl;
	std::cout << "Expected number of clusters from [0, " << Row << "]:";
	std::cin >> expectedClusters;
	assert(expectedClusters>0 && expectedClusters<Row/10);

	/* would store distance matrix instead because it would save massive time */
	struct timeval start, end;
	double timeTemp;
	gettimeofday(&start, NULL);

	int clusterCount = 0;
	const int& minExpected = 0.8*expectedClusters;
	const int& maxExpected = 1.2*expectedClusters;
	float minDist = distRange[0], maxDist = distRange[1]/4.0;
	int iteration = 0;
	std::cout << ".." << std::endl;
	std::cout << ".." << std::endl;
	std::cout << "Binary search starts!" << std::endl;
	while(true)
	{
		distanceThreshold = (minDist+maxDist)/2.0;
		hierarchicalMerging(nodeVec);
		clusterCount = nodeVec.size();
		std::cout << "Iteration " << (++iteration) << " finds " << clusterCount << " groups!" << std::endl;
		if(clusterCount>=minExpected && clusterCount<=maxExpected)
			break;
		else if(clusterCount<minExpected)
			maxDist = distanceThreshold;
		else 
			minDist = distanceThreshold;
	}

	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec  - start.tv_sec) * 1000000u 
			   + end.tv_usec - start.tv_usec) / 1.e6;
	stringstream ss;
	ss << expectedClusters;
	string cluster_str = ss.str();

	ss.str("");
	ss << iteration;

	activityList.push_back("To achieve "+cluster_str+" groups will take " + ss.str() + " binary search and take: ");
	timeList.push_back(to_string(timeTemp)+" s");

}


/* perform hierarchical clustering by given a threshold */
void AHC::bottomUp_byThreshold(std::vector<Ensemble>& nodeVec)
{
	std::cout << "-------------------------------------------------------------------------------" << std::endl;
	std::cout << "Input threshold: [" << distRange[0] << "," << distRange[1] <<"]: ";
	std::cin >> distanceThreshold;
	assert(distanceThreshold>distRange[0] && distanceThreshold<distRange[1]);

	/* would store distance matrix instead because it would save massive time */
	struct timeval start, end;
	double timeTemp;
	gettimeofday(&start, NULL);

	hierarchicalMerging(nodeVec);

	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec  - start.tv_sec) * 1000000u 
			   + end.tv_usec - start.tv_usec) / 1.e6;
	stringstream ss;
	ss << distanceThreshold;
	activityList.push_back("To cluster by distance "+ss.str()+" will take: ");
	timeList.push_back(to_string(timeTemp)+" s");

}


/* perform AHC merging by given a distance threshold */
void AHC::hierarchicalMerging(std::vector<Ensemble>& nodeVec)
{
	const int Row = ds.dataMatrix.rows();

	nodeVec.clear();
//could have used vector, but since there're too many operations inside so should use set
	nodeVec = std::vector<Ensemble>(Row);
//create node in forest structure
#pragma omp parallel for schedule(dynamic) num_threads(8)
	for(int i=0;i<nodeVec.size();++i)
	{
		nodeVec[i].index = i;
		nodeVec[i].element.push_back(i);
	}

	//two iterators to record positions of set
	std::vector<Ensemble>::iterator iter_i, iter_j;

	//vector to store new node
	std::vector<Ensemble> newNodeList;
	do
	{
		//insert new node obtained from previous step
		if(!newNodeList.empty())
		{
			nodeVec.insert(nodeVec.end(), newNodeList.begin(), newNodeList.end());
			newNodeList.clear();
		}
		//iter_i prior, iter_j consecutive
		iter_i = nodeVec.begin();
		iter_j = iter_i;
		++iter_j;
		int mergedCount = 0;
		while(iter_i!=nodeVec.end())
		{
			// j reaches end or i is merged, should move i forward
			if(iter_j==nodeVec.end() || (*iter_i).merged)
			{
				++iter_i;
				iter_j = iter_i;
				++iter_j;
			}
			// j node already merged, so no longer consideration
			else if((*iter_j).merged)
				++iter_j;
			// move j and calculate distance for mutual pairs
			else
			{
				//compute distance between two clusters by single/complete/average linkages
				const float& linkageDist = getDistAtNodes((*iter_i).element, (*iter_j).element, linkageOption);

				//larger distance than threshold, then move forward
				if(linkageDist>distanceThreshold)
					++iter_j;
				//merge two clusters into one cluster if smaller than threshold
				else
				{
					//add merged node whose index is total element size
					vector<int> first = (*iter_i).element, second = (*iter_j).element;
					Ensemble newNode = Ensemble(first.size()+second.size());
					newNode.element = first;
					newNode.element.insert(newNode.element.begin(), second.begin(), second.end());
					newNodeList.push_back(newNode);

					(*iter_i).merged = true;
					(*iter_j).merged = true;

					++iter_i;
					iter_j = iter_i;
					++iter_j;

					mergedCount+=2;
				}
			}
		}

		/* erase would cost so much time so we'd better directly use copy
		for (auto iter=nodeVec.begin(); iter!=nodeVec.end();)
		{
			if((*iter).merged)
				iter=nodeVec.erase(iter);
			else
				++iter;
		}*/


		/* use copy and backup to delete those merged elements */
		assert(nodeVec.size()>=mergedCount);
		std::vector<Ensemble> copyNode(nodeVec.size()-mergedCount);
		int c_i = 0;
		for(int i=0;i<nodeVec.size();++i)
		{
			if(!nodeVec[i].merged)
				copyNode[c_i++] = nodeVec[i];
		}
		nodeVec.clear();
		nodeVec = copyNode;
		copyNode.clear();

		mergedCount = 0;

	}while(!newNodeList.empty());	//merging happens constantly

	newNodeList.clear();

	numberOfClusters = nodeVec.size();

	/* use alpha function to sort the group by its size */
	std::sort(nodeVec.begin(), nodeVec.end(), [](const Ensemble& e1, const Ensemble& e2)
	{return e1.element.size()<e2.element.size() ||(e1.element.size()==e2.element.size()&&e1.index<e2.index);});
}


/* perform group-labeling information */
void AHC::setLabel(const std::vector<Ensemble>& nodeVec, vector<vector<int> >& neighborVec,
			      vector<int>& storage, Eigen::MatrixXf& centroid)
{
// group tag by increasing order
	int groupID = 0;

	// element list for each group
	vector<int> eachContainment;

	// find group id and neighboring vec
	for(auto iter = nodeVec.begin(); iter!=nodeVec.end();++iter)
	{
		eachContainment = (*iter).element;
		neighborVec[groupID] = eachContainment;
	#pragma omp parallel num_threads(8)
		{
		#pragma omp for nowait
			for(int i=0;i<eachContainment.size();++i)
			{
				group[eachContainment[i]] = groupID;
			#pragma omp critical
				centroid.row(groupID) += ds.dataMatrix.row(eachContainment[i]);
			}
		}
		storage[groupID] = (*iter).element.size();
		centroid.row(groupID)/=eachContainment.size();
		++groupID;
		eachContainment.clear();
	}
}



/* extract features from datasets as representative curves */
void AHC::extractFeatures(const std::vector<int>& storage, const std::vector<std::vector<int> >& neighborVec,
		                  const Eigen::MatrixXf& centroid)
{
	const int& Row = ds.dataMatrix.rows();
	const int& Column = ds.dataMatrix.cols();

	std::cout << "Final group number information: " << std::endl;
	for (int i = 0; i < storage.size(); ++i)
	{
		std::cout << storage[i] << " ";
	}
	std::cout << std::endl;

	IOHandler::printClustersNoise(ds.dataVec,group,storage,
			 "norm"+to_string(normOption),ds.fullName,ds.dimension);

	struct timeval start, end;
	double timeTemp;

	gettimeofday(&start, NULL);
	Silhouette sil;
	sil.computeValue(normOption,ds.dataMatrix,ds.dataMatrix.rows(),ds.dataMatrix.cols(),group,object,numberOfClusters);
	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec  - start.tv_sec) * 1000000u 
			   + end.tv_usec - start.tv_usec) / 1.e6;
	activityList.push_back("Silhouette calculation takes: ");
	timeList.push_back(to_string(timeTemp)+" s");

	/* compute the centroid coordinates of each clustered group */

	gettimeofday(&start, NULL);

	vector<vector<float> > closest(numberOfClusters);
	vector<vector<float> > furthest(numberOfClusters);

	/* extract the closest and furthest streamlines to centroid */

#pragma omp parallel for schedule(dynamic) num_threads(8)
	for (int i=0;i<numberOfClusters;++i)
	{
		float minDist = FLT_MAX;
		float maxDist = -10;
		int minIndex = -1, maxIndex = -1;
		const std::vector<int>& groupRow = neighborVec[i];
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

	std::vector<std::vector<float> > center_vec(numberOfClusters, vector<float>(Column));
#pragma omp parallel for schedule(dynamic) num_threads(8)
	for (int i = 0; i < center_vec.size(); ++i)
	{
		for (int j = 0; j < Column; ++j)
		{
			center_vec[i][j] = centroid(i,j);
		}
	}


	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec  - start.tv_sec) * 1000000u 
			   + end.tv_usec - start.tv_usec) / 1.e6;
	activityList.push_back("Feature extraction takes: ");
	timeList.push_back(to_string(timeTemp)+" s");

	std::cout << "Finishing extracting features!" << std::endl;	

	stringstream ss;
	ss << "norm_" << normOption;

	string linkage = getLinkageStr();

	IOHandler::printFeature(ds.dataName+"_AHC_"+linkage+"_closest_"+ss.str()+".vtk", closest, sil.sCluster, ds.dimension);
	IOHandler::printFeature(ds.dataName+"_AHC_"+linkage+"_furthest_"+ss.str()+".vtk", furthest, sil.sCluster, ds.dimension);
	IOHandler::printFeature(ds.dataName+"_AHC_"+linkage+"_centroid_"+ss.str()+".vtk", center_vec, sil.sCluster,ds.dimension);

	IOHandler::printToFull(ds.dataVec, sil.sData, "AHC_SValueLine_"+ss.str(), ds.fullName, ds.dimension);
	IOHandler::printToFull(ds.dataVec, group, sil.sCluster, "AHC_SValueCluster_"+ss.str(), ds.fullName, ds.dimension);

	activityList.push_back("numCluster is: ");
	timeList.push_back(to_string(numberOfClusters));

	activityList.push_back("Average Silhouette is: ");
	timeList.push_back(to_string(sil.sAverage));

	IOHandler::generateReadme(activityList,timeList);

	IOHandler::writeReadme("Linkage: "+linkage+", ");
}

/* set dataset from user command */
void AHC::setDataset(const int& argc, char **argv)
{
	if(argc!=3)
	{
		std::cout << "Input argument should have 3!" << endl
		          << "./cluster inputFile_name(in dataset folder) "
		          << "data_dimension(3)" << endl;
		exit(1);
	}
	ds.strName = string("../dataset/")+string(argv[1]);
	ds.dataName = string(argv[1]);
	ds.dimension = atoi(argv[2]);

	int sampleOption;
    std::cout << "choose a sampling method for the dataset?" << std::endl
	    	  << "1.directly filling with last vertex; 2. uniform sampling." << std::endl;
	std::cin >> sampleOption;
	assert(sampleOption==1||sampleOption==2);

	IOHandler::readFile(ds.strName,ds.dataVec,ds.vertexCount,ds.dimension,ds.maxElements);

	ds.fullName = ds.strName+"_full.vtk";
	IOHandler::printVTK(ds.fullName, ds.dataVec, ds.vertexCount, ds.dimension);

	if(sampleOption==1)
		IOHandler::expandArray(ds.dataMatrix,ds.dataVec,ds.dimension,ds.maxElements);
	else if(sampleOption==2)
		IOHandler::sampleArray(ds.dataMatrix,ds.dataVec,ds.dimension,ds.maxElements);

	group = std::vector<int>(ds.dataMatrix.rows());

	std::cout << "---------------------------" << std::endl;
	std::cout << "Input linkage option: 0.single linkage, 1.complete linkage, 2.average linkage" << std::endl;
	std::cin >> linkageOption;
	assert(linkageOption==0||linkageOption==1||linkageOption==2);
}

/* set norm option, must be within 0-12 */
void AHC::setNormOption()
{
	std::cout << "Input a norm option 0-12!" << std::endl;
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
		12: mean of closest point distance
	*/
	bool found = false;
	for (int i = 0; i < 13&&!found; ++i)
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


/* set threshold for AHC function */
void AHC::getDistRange()
{
	const float& Percentage = 0.05;
	const int& Row = ds.dataMatrix.rows();

	distRange = vector<float>(2);
	distRange[0] = FLT_MAX, distRange[1] = FLT_MIN;
	const int& totalSize = int(Percentage*Row);
#pragma omp parallel num_threads(8)
	{
	#pragma omp for nowait
		for (int i = 0; i < totalSize; ++i)
		{
			float tempDist, i_min = FLT_MAX, i_max = FLT_MIN;
			for (int j = 0; j < Row; ++j)
			{
				if(distanceMatrix)
					tempDist = distanceMatrix[i][j];
				else
					tempDist = getDisimilarity(ds.dataMatrix, i, j, normOption, object);
				if(tempDist<i_min)
				{
					i_min = tempDist;
				}
				if(tempDist>i_max)
				{
					i_max = tempDist;
				}
			}

		#pragma omp critical
			{
				distRange[0] = std::min(distRange[0], i_min);
				distRange[1] = std::max(distRange[1], i_max);
			}
		}
	}

	std::cout << "Distance threshold is: [" << distRange[0] << ", " << distRange[1] << "]." << std::endl;
}


const float AHC::getDistAtNodes(const vector<int>& firstList, const vector<int>& secondList, const int& Linkage)
{
	const int& m = firstList.size();
	const int& n = secondList.size();
	assert(m!=0);
	assert(n!=0);
	/* 0: single linkage, min(x_i,y_j)
	 * 1: complete linkdage, max(x_i,y_j)
	 * 2: average linkage, sum/x_i/y_j
	 */
	float result, value;
	switch(Linkage)
	{
	case 0:	//single linkage
		{
			result = FLT_MAX;
		#pragma omp parallel for reduction(min:result) num_threads(8)
			for(int i=0;i<m;++i)
			{
				for(int j=0;j<n;++j)
				{
					if(distanceMatrix)
						value = distanceMatrix[firstList[i]][secondList[j]];
					else
						value = getDisimilarity(ds.dataMatrix, firstList[i], secondList[j], normOption, object);
					result = std::min(result, value);
				}
			}
		}
		break;

	case 1:	//complete linkage
		{
			result = FLT_MIN;
		#pragma omp parallel for reduction(max:result) num_threads(8)
			for(int i=0;i<m;++i)
			{
				for(int j=0;j<n;++j)
				{
					if(distanceMatrix)
						value = distanceMatrix[firstList[i]][secondList[j]];
					else
						value = getDisimilarity(ds.dataMatrix, firstList[i], secondList[j], normOption, object);
					result = std::max(result, value);
				}
			}
		}
		break;

	case 2:
		{
			result = 0;
		#pragma omp parallel for reduction(+:result) num_threads(8)
			for(int i=0;i<m;++i)
			{
				for(int j=0;j<n;++j)
				{
					if(distanceMatrix)
						value = distanceMatrix[firstList[i]][secondList[j]];
					else
						value = getDisimilarity(ds.dataMatrix, firstList[i], secondList[j], normOption, object);
					result+=value;
				}
			}
			result/=m*n;
		}
		break;

	default:
		std::cout << "error!" << std::endl;
		exit(1);
	}
	return result;
}


/* get string for linkage type */
string AHC::getLinkageStr()
{
	string result;
	switch(linkageOption)
	{
	case 0:
		result = "single";
		break;

	case 1:
		result = "complete";
		break;

	case 2:
		result = "average";
		break;
	}
	return result;
}