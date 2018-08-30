#include "OPTICS.h"

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
	
	/* if the dataset is not PBF, then should record distance matrix for Gamma matrix compution */
	if(!isPBF)
	{
		deleteDistanceMatrix(ds.dataMatrix.rows());

		if(!getDistanceMatrix(ds.dataMatrix, normOption, object))
		{
			std::cout << "Failure to compute distance matrix!" << std::endl;
		}
	}

	nodeVec = vector<PointNode>(ds.dataMatrix.rows(),PointNode());
}


DensityClustering::~DensityClustering()
{

}


void DensityClustering::performClustering()
{
	float radius_eps;

	int minPts = setMinPts();

	int epsOption = 2;
	/*std::cout << "Choose eps selection method. 1. user input of multiplication, 2. minPt-th dist." << std::endl;
	std::cin >> epsOption;
	assert(epsOption==1||epsOption==2);*/

	if(epsOption==1)
	{
		float minDist, maxDist;
		getDistRange(minDist, maxDist);
		std::cout << "Distance range is [" << minDist << ", "
				  << maxDist << "]." << std::endl;
		multiTimes = setTimesMin(minDist, maxDist);

		radius_eps = maxDist*multiTimes;
	}
	else if(epsOption==2)
	{
		radius_eps = getMinPt_thDist(minPts);
	}

	struct timeval start, end;
	double timeTemp;
	gettimeofday(&start, NULL);

	OPTICS(radius_eps, minPts);

	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec  - start.tv_sec) * 1000000u 
			   + end.tv_usec - start.tv_usec) / 1.e6;
	activityList.push_back("OPTICS clustering takes: ");
	timeList.push_back(to_string(timeTemp)+" s");

	getGroup(radius_eps);
	extractFeatures(radius_eps, minPts);
}


void DensityClustering::OPTICS(const float& radius_eps,
							   const int& minPts)
{
	computeCoredDistance(radius_eps, minPts);

	for(int i=0;i<ds.dataMatrix.rows();++i)
	{
		if(nodeVec[i].visited)
			continue;
		const vector<int>& neighbor = nodeVec[i].neighbor;
			  /*regionQuery(i,radius_eps)*/;
		nodeVec[i].visited = true;
		orderedList.push_back(i);
		if(nodeVec[i].core_distance!=-1.0)
		{
			/* linkedList is good, but would be O(n), so delete it */
			 LinkedList seeds;

			/* should use priority queue */
			update(i, neighbor, seeds, radius_eps, minPts);
			pointNode *temp = seeds.start;
			while(temp)
			{
				if(nodeVec[temp->value.index].visited)
				{
					temp = temp->next;
					continue;
				}
				const vector<int>& neighborChild = nodeVec[temp->value.index].neighbor;
					  /*regionQuery(temp->value.index, radius_eps)*/;
				nodeVec[temp->value.index].visited = true;
				orderedList.push_back(temp->value.index);
				if(nodeVec[temp->value.index].core_distance!=-1.0)
					update(temp->value.index, neighborChild, seeds, radius_eps, minPts);
				temp = seeds.start;
			}
		}
	}
}


void DensityClustering::update(const int& index, const vector<int>& neighbor,
	 LinkedList& seeds, //a min heap for orderedPoint
	 const float& radius_eps, const int& minPts)
{
	const float& coredist = nodeVec[index].core_distance;
	for(int i=0;i<neighbor.size();++i)
	{
		if(!nodeVec[neighbor[i]].visited)
		{
			float dist_toCenter;
			if(distanceMatrix)
				dist_toCenter = distanceMatrix[neighbor[i]][index];
			else
				dist_toCenter = getDisimilarity(ds.dataMatrix.row(neighbor[i]),
				  ds.dataMatrix.row(index), neighbor[i], index, normOption, object);
			const float& biggerDist = std::max(coredist, dist_toCenter);
			if(nodeVec[neighbor[i]].reachabilityDist==-1.0)
			{
				nodeVec[neighbor[i]].reachabilityDist=biggerDist;
				seeds.insertNode(new pointNode(OrderedPoint(neighbor[i], biggerDist)));
			}
			else
			{
				float& reachDist = nodeVec[neighbor[i]].reachabilityDist;
				if(biggerDist<reachDist)
				{
					reachDist = biggerDist;
					seeds.updateNode(neighbor[i], reachDist);
				}
			}
		}
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
		if(distanceMatrix)
			tempDist = distanceMatrix[index][i];
		else
			tempDist=getDisimilarity(ds.dataMatrix.row(index),ds.dataMatrix.row(i),index,i,normOption, object);
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

	/* need to judge whether it is a PBF dataset or not */
	std::cout << "It is a PBF dataset? 1.Yes, 0.No." << std::endl;
	int PBFInput;
	std::cin >> PBFInput;
	assert(PBFInput==1||PBFInput==0);
	isPBF = (PBFInput==1);

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
	std::cout << "Choose a norm from 0-15!" << std::endl;
	std::cin >> normOption;
	std::cout << std::endl;
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
	bool found = false;
	for (int i = 0; i < 16&&!found; ++i)
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
				float dist;
				if(distanceMatrix)
					dist = distanceMatrix[i][j];
				else
					dist = getDisimilarity(ds.dataMatrix.row(i),ds.dataMatrix.row(j),i,j,normOption,object);
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
}


const int DensityClustering::setMinPts()
{
	/*std::cout << std::endl;
	std::cout << "Input the minPts for OPTICS in [0" << ", "
			  << ds.dataMatrix.rows() << "], 6 is recommended:" << std::endl;
			  */
	int minPts = 6;
	//std::cin >> minPts;

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
	std::cout << "Input the multiplication for OPTICS radius in ["
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

#pragma omp parallel for schedule(static) num_threads(8)
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

	const int& Row = ds.dataMatrix.rows();
	float entropy = 0.0, probability;
	for(int i=0;i<container.size();++i)
	{
		probability = float(container[i])/float(Row);
		entropy+=probability*log2f(probability);
	}
	entropy = -entropy/log2f(numClusters);


	IOHandler::printClustersNoise(ds.dataVec,item_cids,container, 
		 "norm"+to_string(normOption),ds.fullName,ds.dimension);

	struct timeval start, end;
	double timeTemp;

	numClusters-=1;

	gettimeofday(&start, NULL);
	Silhouette sil;
	sil.computeValue(normOption,ds.dataMatrix,ds.dataMatrix.rows(),ds.dataMatrix.cols(),
					 item_cids,object,numClusters, isPBF);

	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec  - start.tv_sec) * 1000000u 
			   + end.tv_usec - start.tv_usec) / 1.e6;
	activityList.push_back("Silhouette calculation takes: ");
	timeList.push_back(to_string(timeTemp)+" s");


	IOHandler::writeReadme(entropy, sil);

	const int& numNoise = storage[0].size();
	storage.erase(storage.begin());


	/* record labeling information */
	IOHandler::generateGroups(storage);

	/* compute the centroid coordinates of each clustered group */

	gettimeofday(&start, NULL);

	Eigen::MatrixXf centroid = MatrixXf::Zero(numClusters,ds.dataMatrix.cols());
	vector<vector<float> > cenVec(numClusters);
#pragma omp parallel for schedule(static) num_threads(8)
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

#pragma omp parallel for schedule(static) num_threads(8)
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

/* measure closest and furthest rotation */
	std::vector<float> closestRot, furthestRot;
	const float& closestAverage = getRotation(closest, closestRot);
	const float& furthestAverage = getRotation(furthest, furthestRot);

	IOHandler::writeReadme(closestAverage, furthestAverage);

	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
	activityList.push_back("Feature extraction takes: ");
	timeList.push_back(to_string(timeTemp)+" s");

	ValidityMeasurement vm;
	vm.computeValue(normOption, ds.dataMatrix, item_cids, object, isPBF);
	activityList.push_back("Optics Validity measure is: ");
	timeList.push_back(to_string(vm.f_c));

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

	activityList.push_back("Noise number is: ");
	timeList.push_back(to_string(numNoise));

	activityList.push_back("radius eps is: ");
	timeList.push_back(to_string(multiTimes));

	activityList.push_back("MinPts is: ");
	timeList.push_back(to_string(minPts));

	IOHandler::generateReadme(activityList,timeList);
}

/* compute neighbor information, update core-distance and store neighbor index information*/
void DensityClustering::computeCoredDistance(const float& radius_eps,
							  				 const int& minPts)
{
#pragma omp parallel for schedule(static) num_threads(8)
	for (int i = 0; i < ds.dataMatrix.rows(); ++i)
	{
		vector<float> distRecord;	//record distance value
		for (int j = 0; j < ds.dataMatrix.rows(); ++j)
		{
			if(j==i)
				continue;
			float tempDist;
			if(distanceMatrix)
				tempDist = distanceMatrix[i][j];
			else
				tempDist = getDisimilarity(ds.dataMatrix.row(i),ds.dataMatrix.row(j),i, j, normOption, object);
			if(tempDist<=radius_eps)
			{
				nodeVec[i].neighbor.push_back(j);
				distRecord.push_back(tempDist);
			}
		}
		/* find minPts-th smallest element in vector by linear traversal */
		if(distRecord.size()>=minPts)
		{
			/* A k*n complex algorithm */
			/*
			vector<float> smallestRange(minPts,FLT_MAX);	//update to get minPts-th smallest
			for(int k=0;k<distRecord.size();++k)
			{
				if(distRecord[k]<smallestRange[minPts-1])
					smallestRange[minPts-1]=distRecord[k];
				for(int l=minPts-1;l>=1;--l)
				{
					if(smallestRange[l]>smallestRange[l-1])
						std::swap(smallestRange[l], smallestRange[l-1]);
				}
			}
			nodeVec[i].core_distance = smallestRange[minPts-1];
			*/

			/* instead we shall apply a n*logk algorithm */
			std::priority_queue<float> smallestRange;
			for(int k=0;k<distRecord.size();++k)
			{
				smallestRange.push(distRecord[k]);
				if(smallestRange.size()>minPts)
					smallestRange.pop();
			}
			nodeVec[i].core_distance = smallestRange.top();
		}
	}	

	std::cout << "Precomputing for cored-distance is done!" << std::endl;
}


const float DensityClustering::getReachability(const int& first,
											   const int& target,
											   const int& minPts)
{
	const int& size = nodeVec[target].neighbor.size();
	if(size<minPts)
		return -1.0;
	float dist;

	if(distanceMatrix)
		dist = distanceMatrix[first][target];
	else
		dist = getDisimilarity(ds.dataMatrix.row(first),ds.dataMatrix.row(target),first, target, normOption, object);
	return std::max(nodeVec[target].core_distance, dist);
}


/* how to get group information based on reachability-plot */
void DensityClustering::getGroup(const float& radius_eps)
{
	std::cout << "----Parameter regime----" << std::endl;
	int continueOption;
	do
	{
		/* group tag information */
		int tag = 0;
		std::cout << "Input threshold for OPTICS reachability-plot: ";
		float threshold;
		std::cin >> threshold;
		threshold*=radius_eps;
		std::cout << threshold << std::endl;
		bool findSummit = false;
		for(int i=0;i<orderedList.size();++i)
		{
			if(nodeVec[orderedList[i]].reachabilityDist>=threshold)
			{
				findSummit = true;
			}
			else
			{
				if(findSummit)
				{
					findSummit = false;
					++tag;
				}
				nodeVec[orderedList[i]].group = tag;
			}
		}
		std::cout << "Finally it forms " << (tag+1) << " clusters!" << std::endl;
		std::cout << "Want to continue with parameter? 1. Yes, 0. No." << std::endl;
		std::cin >> continueOption;
		assert(continueOption==1||continueOption==0);
	}while(continueOption==1);
	writeReachability();
}


void DensityClustering::writeReachability()
{
	ofstream ofile("../dataset/reachability.txt", ios::out);
	if(!ofile)
	{
		std::cout << "Error creating the file!" << std::endl;
		exit(1);
	}
	for (int i = 0; i < orderedList.size(); ++i)
	{
		ofile << orderedList[i] << " " << nodeVec[orderedList[i]].reachabilityDist 
		<< " " << nodeVec[orderedList[i]].group << std::endl;
	}
	ofile << std::endl;
	ofile.close();
}


/* set the eps as averaged minPt-th dist */
const float DensityClustering::getMinPt_thDist(const int& minPts)
{
	float result = 0.0;
	const int& rowSize = ds.dataMatrix.rows();
	const int& seletedRow = 0.2*rowSize;
#pragma omp parallel num_threads(8)
	{
	#pragma omp for nowait
		for (int i = 0; i < seletedRow; ++i)
		{
			/* a linear k*n implementation by directly using a vector for linear mapping */
			/*
			std::vector<float> minDistVec(minPts, FLT_MAX);
			float tempDist;
			for (int j=0;j<rowSize;++j)
			{
				if(i==j)
					continue;
				if(distanceMatrix)
					tempDist = distanceMatrix[i][j];
				else
					tempDist=getDisimilarity(ds.dataMatrix.row(i), ds.dataMatrix.row(j),i,j,normOption, object);

				if(tempDist<minDistVec[minPts-1])
					minDistVec[minPts-1]=tempDist;
				for(int l=minPts-1;l>=1;--l)
				{
					if(minDistVec[l]>minDistVec[l-1])
						std::swap(minDistVec[l], minDistVec[l-1]);
				}
			}*/

			/* use a priority_queue<float> with n*logk time complexity */
			std::priority_queue<float> minDistArray;
			float tempDist;
			for (int j=0;j<rowSize;++j)
			{
				if(i==j)
					continue;
				if(distanceMatrix)
					tempDist = distanceMatrix[i][j];
				else
					tempDist=getDisimilarity(ds.dataMatrix.row(i), ds.dataMatrix.row(j),i,j,normOption, object);

				minDistArray.push(tempDist);
				if(minDistArray.size()>minPts)
					minDistArray.pop();

			}

		#pragma omp critical
			result += minDistArray.top();
		}
	}
	return result/seletedRow;
}
