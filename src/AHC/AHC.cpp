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

	object = MetricPreparation(ds.dataMatrix.rows(), ds.dataMatrix.cols());
	object.preprocessing(ds.dataMatrix, ds.dataMatrix.rows(), ds.dataMatrix.cols(), normOption);

	getDistanceMatrix(ds.dataMatrix, normOption, object);
	setThreshold();
}

/* destructor */
AHC::~AHC()
{
	deleteDistanceMatrix(ds.dataMatrix.rows());
}

/* perform clustering function */
void AHC::performClustering()
{
	const int& Row = ds.dataMatrix.rows();

//could have used vector, but since there're too many operations inside so should use set
	std::unordered_set<TreeNode*> nodeSet;
//create node in forest structure
	for(int i=0;i<ds.dataMatrix.rows();++i)
		nodeSet.insert(new TreeNode(i));

	//two iterators to record positions of set 
	std::unordered_set<TreeNode*>::iterator iter_i, iter_j;

	//vector to store new node 
	std::vector<TreeNode*> newNodeList;
	do
	{
		//insert new node obtained from previous step
		if(!newNodeList.empty())
		{
			nodeSet.insert(newNodeList.begin(), newNodeList.end());
			newNodeList.clear();
		}
		//iter_i prior, iter_j consecutive
		iter_i = nodeSet.begin();
		iter_j = iter_i;
		++iter_j;
		while(iter_i!=nodeSet.end())
		{
			// j reaches end, should move i forward
			if(iter_j==nodeSet.end())
			{
				++iter_i;
				iter_j = iter_i;
				++iter_j;
			}

			// move j and calculate distance for mutual pairs
			else
			{
				vector<int> firstList, secondList;

				//dfs traversal to get element vector for two nodes
				dfsTraversal(firstList, *iter_i);
				dfsTraversal(secondList, *iter_j);

				//compute distance between two clusters by single/complete/average linkages
				const float& linkageDist = getDistAtNodes(firstList, secondList, linkageOption);

				//larger distance than threshold, then move forward
				if(linkageDist>distanceThreshold)
					++iter_j;
				//merge two clusters into one cluster if smaller than threshold
				else
				{
					TreeNode* left = *iter_i, *right = *iter_j;
					if(!left||!right)
					{
						std::cout << "Error found NULL in agglomerative hierarchical clustering!" << std::endl;\
						exit(1);
					}

					//add merged node whose index is total element size
					TreeNode* newNode = new TreeNode(firstList.size()+secondList.size());
					newNode->left=left;
					newNode->right=right;
					newNodeList.push_back(newNode);

					//remove two nodes from orginal vector
					nodeSet.erase(iter_i);
					nodeSet.erase(iter_j);

					//calculate j-i
					int diff = std::distance(iter_i, iter_j);
					//more than one, then move i to i+1
					if(diff>1)
					{
						++iter_i;
						iter_j=iter_i;
						++iter_j;
					}
					//i and j are adjacent, then i = j+1, j = i+1
					else
					{
						++iter_j;
						iter_i=iter_j;
						++iter_j;
					}
				}
			}
		}

	}while(!newNodeList.empty());	//merging happens constantly

	//create an ordered_map to sort the node by size of contained elements
	std::map<int, TreeNode*> increasingOrder;
	for(auto iter = nodeSet.begin();iter!=nodeSet.end();++iter)
	{
		increasingOrder.insert(make_pair((*iter)->index, *iter));
	}

	nodeSet.clear();

	// group tag by increasing order
	int groupID = 0;

	// element lists for all groups
	vector<vector<int>> neighborVec(increasingOrder.size());

	// element list for each group
	vector<int> eachContainment;

	// element size for all groups
	vector<int> storage(increasingOrder.size());

	// find centroid matrix
	Eigen::MatrixXf centroid = Eigen::MatrixXf::Zero(increasingOrder.size(),ds.dataMatrix.cols());

	// find group id and neighboring vec
	for(auto iter = increasingOrder.begin(); iter!=increasingOrder.end(); ++iter)
	{
		dfsTraversal(eachContainment, (*iter).second);
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
			eachContainment.clear();
			storage[groupID] = (*iter).first;
			centroid.row(groupID)/=eachContainment.size();
			++groupID;
		}
	}
	numberOfClusters = groupID;
	increasingOrder.clear();

	extractFeatures(storage, neighborVec, centroid);
}

/* extract features from datasets as representative curves */
void AHC::extractFeatures(const std::vector<int>& storage, const std::vector<std::vector<int> >& neighborVec,
		                  const Eigen::MatrixXf& centroid)
{
	const int& Row = ds.dataMatrix.rows();
	const int& Column = ds.dataMatrix.cols();

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
	IOHandler::printFeature("optimization_closest.vtk", closest, sil.sCluster, ds.dimension);
	IOHandler::printFeature("optimization_furthest.vtk", furthest, sil.sCluster, ds.dimension);
	IOHandler::printFeature("optimization_centroid.vtk", center_vec, sil.sCluster,ds.dimension);

	IOHandler::printToFull(ds.dataVec, sil.sData, "optimization_SValueLine", ds.fullName, ds.dimension);
	IOHandler::printToFull(ds.dataVec, group, sil.sCluster, "optimization_SValueCluster", ds.fullName, ds.dimension);

	activityList.push_back("numCluster is: ");
	timeList.push_back(to_string(numberOfClusters));

	activityList.push_back("Average Silhouette is: ");
	timeList.push_back(to_string(sil.sAverage));

	IOHandler::generateReadme(activityList,timeList);
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
void AHC::setThreshold()
{
	const float& Percentage = 0.05;
	const int& Row = ds.dataMatrix.rows();
	float minDist = FLT_MAX, maxDist = FLT_MIN;
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
				minDist = std::min(minDist, i_min);
				maxDist = std::max(maxDist, i_max);
			}
		}
	}

	std::cout << "Distance threshold is: [" << minDist << ", " << maxDist << "]." << std::endl;
	std::cout << "Input threshold: ";
	std::cin >> distanceThreshold;
	assert(distanceThreshold>minDist && distanceThreshold<maxDist);
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

