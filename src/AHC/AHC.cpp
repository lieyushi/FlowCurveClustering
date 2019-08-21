/*
 * @brief This is the source cpp for implementing the member functions for the class AHC
 * @author Lieyu Shi
 */

#include "AHC.h"


/*
 * @brief The default constructor
 */
AHC::AHC()
{

}


/*
 * @brief A constructor with parameters
 * @details
 * 	To set the data set and perform some reading operation into the member variables from the file
 *
 * @param[in] argc The count of argv
 * @param[in] argv The argument string line
 */
AHC::AHC(const int& argc, char **argv)
{
	setDataset(argc, argv);
}


/*
 * @brief The destructor of the class AHC
 * @details
 * 	Delete the global pointer distance matrix
 *
 */
AHC::~AHC()
{
	deleteDistanceMatrix(ds.dataMatrix.rows());
}


/*
 * @brief The function to perform AHC clustering by a given norm
 *
 * @details
 * 	The function will first judge whether the local storage of distance matrix exists or not. If it exists, the program
 * 	will read in the distance matrix from the local file, otherwise it will calculate the distance matrix. In this format
 * 	the time for calculating distance matrix can be saved for different clustering techniques.
 * 	Then it will read in number of clusters as input for the clustering operation.
 */
void AHC::performClustering_by_norm()
{
	/* very hard to decide whether needed to perform such pre-processing, but still create a
	 * MetricPreparation object in case some pre-calculation for similarity meaures can be ready
	 * before the pairwise distance matrix calculation
	 */
	object = MetricPreparation(ds.dataMatrix.rows(), ds.dataMatrix.cols());
	object.preprocessing(ds.dataMatrix, ds.dataMatrix.rows(), ds.dataMatrix.cols(), normOption);

	/* would store distance matrix instead because it would save massive time */
	struct timeval start, end;
	double timeTemp;
	gettimeofday(&start, NULL);

	// check whether the file for distance matrix exists or not
	std::ifstream distFile(("../dataset/"+to_string(normOption)).c_str(), ios::in);
	if(distFile.fail())	// not exist, will calculate the distance matrix and store them in local files
	{
		distFile.close();
		// calculate the distance matrix
		getDistanceMatrix(ds.dataMatrix, normOption, object);
		// store the distance matrix values in the local files
		std::ofstream distFileOut(("../dataset/"+to_string(normOption)).c_str(), ios::out);
		for(int i=0;i<ds.dataMatrix.rows();++i)
		{
			for(int j=0;j<ds.dataMatrix.rows();++j)
			{
				distFileOut << distanceMatrix[i][j] << " ";
			}
			distFileOut << std::endl;
		}
		distFileOut.close();
	}
	else	// the file for distance matrix exists, then directly reads in the pair-wise values
	{
		std::cout << "read distance matrix..." << std::endl;

		distanceMatrix = new float*[ds.dataMatrix.rows()];
	#pragma omp parallel for schedule(static) num_threads(8)
		for (int i = 0; i < ds.dataMatrix.rows(); ++i)
		{
			distanceMatrix[i] = new float[ds.dataMatrix.rows()];
		}

		// read the distance values from the .txt file
		int i=0, j;
		string line;
		stringstream ss;
		while(getline(distFile, line))
		{
			j=0;
			ss.str(line);
			while(ss>>line)
			{
				if(i==j)
					distanceMatrix[i][j]=0;
				else
					distanceMatrix[i][j] = std::atof(line.c_str());
				++j;
			}
			++i;
			ss.str("");
			ss.clear();
		}
		distFile.close();
	}

	// record the time for distance matrix computation time
	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec  - start.tv_sec) * 1000000u
			   + end.tv_usec - start.tv_usec) / 1.e6;
	activityList.push_back("Distance matrix computing for norm "+to_string(normOption)+" takes: ");
	timeList.push_back(to_string(timeTemp)+" s");

	// create node-related parameters for AHC clustering
	std::unordered_map<int, Ensemble> node_map;
	std::vector<DistNode> dNodeVec;
	std::vector<Ensemble> nodeVec;

	/* set the ditNode vector */
	setValue_merge(dNodeVec, node_map);

	/* perform hiarchical clustering where within each step would merge two nodes */
	hierarchicalMerging(node_map, dNodeVec, nodeVec);

	if(!lMethod)	// perform the AHC clustering with lMethod not activated
	{
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
}


/*
 * @brief Perform the clustering for selected similarity measure labels w.r.t. user input and data set type
 * @details
 *	If the number of clusters are pre-stored in the "cluster_number" file, the code will read the numbers first.
 *	Then for different similarity measures, it will decide whether the L-method is activated or not. If L-method
 *	is activated, the number of clusters is set to be 1, otherwise it will be set as user input. Hierarchical merging
 *	operation for the tree is called after parameter setting is finished.
 */
void AHC::performClustering()
{
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
		16: entropy-based distance metric taken from http://vis.cs.ucdavis.edu/papers/pg2011paper.pdf
		17: time-series MCP distance from https://www.sciencedirect.com/science/article/pii/S0097849318300128
			for pathlines only
	*/

	// read the input number for different similarity measures from the file "cluster_number"
	std::unordered_map<int,int> clusterMap;
	if(readCluster)
	{
		IOHandler::readClusteringNumber(clusterMap, "cluster_number");
	}

	// std::vector<int> cluster_array;
	//for(int i=2; i<=100; ++i)
	//{
		//cluster_array.push_back(i);
		for(normOption=0;normOption<=17;++normOption)
		{
			if(isPathlines)	// for pathlines, will consider d_T (17)
			{
				if(normOption!=0 && normOption!=1 && normOption!=2 && normOption!=4 && normOption!=12
				   && normOption!=13 && normOption!=14 && normOption!=15 && normOption!=17)
					continue;
			}
			else	// for streamlines, will not consider d_T (17)
			{
				if(normOption!=0 && normOption!=1 && normOption!=2 && normOption!=4 && normOption!=12
				&& normOption!=13 && normOption!=14 && normOption!=15)
					continue;
			}
			std::cout << "----------------------------------------------------------" << std::endl;
			std::cout << "norm " << normOption << " starts......" << std::endl;
			timeList.clear();
			activityList.clear();

			/* L-method is not performed. It's a normal AHC procedure */
			if(!lMethod)
			{
				const int& Row = ds.dataMatrix.rows();
				if(readCluster)
				{
					numberOfClusters = clusterMap[normOption];
				}
				else
				{
					std::cout << "---------------------------------------" << std::endl;
					std::cout << "Input cluster number among [0, " << Row << "] for norm " << normOption << ": ";
					std::cin >> numberOfClusters;
					assert(numberOfClusters>0 && numberOfClusters<Row);
				}
				//numberOfClusters = i;
				assert(numberOfClusters>0 && numberOfClusters<Row);
			}
			/* perform L-method for detecting optimal num of clusters */
			else if(lMethod)
			{
				numberOfClusters = 1;
			}
			// perform clustering by given input of norm option
			performClustering_by_norm();
		}
	//}

	/*std::ofstream curve("../dataset/curveValue.txt", ios::out);
	for(int i=0; i<4; ++i)
	{
		for(int j=0; j<curveValue[0].size(); ++j)
		{
			curve << curveValue[i][j] << " ";
		}
		curve << std::endl;
	}
	curve.close();*/
}


/*
 * @brief Perform hierarchical merge for the trees by a given required cluster number
 * @details
 * 	Hiarachically merge the nodes until the number of cluster is reached. Then based on whether L-method is activated
 * 	or not, the posterior operation will be called on either finding the clustering information or finding the optimal
 * 	number of clusters
 *
 * @param[out] node_map The initial node with each node representing one streamline/pathline
 * @param[out] dNodeVec The DistNode vector which has the indices of two nodes and their distance
 * @param[out] nodeVec The vector of Ensemble which has candidate index of the cluster
 */
void AHC::hierarchicalMerging(std::unordered_map<int, Ensemble>& node_map, std::vector<DistNode>& dNodeVec,
							  std::vector<Ensemble>& nodeVec)
{
	std::map<int, float> dist_map;

	/* would store distance matrix instead because it would save massive time */
	struct timeval start, end;
	double timeTemp;
	gettimeofday(&start, NULL);

	const int Row = ds.dataMatrix.rows();

	DistNode poped;

	/* find node-pair with minimal distance */
	float minDist = FLT_MAX;
	int target = -1;
	for (int i = 0; i < dNodeVec.size(); ++i)
	{
		if(dNodeVec[i].distance<minDist)
		{
			target = i;
			minDist = dNodeVec[i].distance;
		}
	}

	// find which distNode is to be popped
	poped = dNodeVec[target];

	int index = Row, currentNumber;
	do	// perform iterative hierarchical merging until the final cluster reaches the given input number
	{
		if(lMethod)	// if the l-method is enabled, record the number of clusters and merged distance
		{
			dist_map.insert(std::make_pair(node_map.size(), poped.distance));
		}
		//create new node merged and input it into hash unordered_map
		vector<int> first = (node_map[poped.first]).element;
		vector<int> second = (node_map[poped.second]).element;

		/* index would be starting from Row */
		Ensemble newNode(index);
		newNode.element = first;
		newNode.element.insert(newNode.element.end(), second.begin(), second.end());
		node_map.insert(make_pair(index, newNode));

		//delete two original nodes
		node_map.erase(poped.first);
		node_map.erase(poped.second);

		/* the difficulty lies how to update the min-heap with linkage
		 * This would take 2NlogN.
		 * Copy all node-pairs that are not relevant to merged nodes to new vec.
		 * For relevant, would update the mutual distance by linkage
		 */

		/* how many clusters exist */
		currentNumber = node_map.size();

		target = -1, minDist = FLT_MAX;

		// create new distNode vector
		std::vector<DistNode> tempVec(currentNumber*(currentNumber-1)/2);

		// update and find the minimal distance for next merging
		int current = 0, i_first, i_second;
		for(int i=0;i<dNodeVec.size();++i)
		{
			i_first=dNodeVec[i].first, i_second=dNodeVec[i].second;
			/* not relevant, directly copied to new vec */
			if(i_first!=poped.first&&i_first!=poped.second&&i_second!=poped.first&&i_second!=poped.second)
			{
				tempVec[current]=dNodeVec[i];
				if(tempVec[current].distance<minDist)
				{
					target = current;
					minDist = tempVec[current].distance;
				}
				++current;
			}
		}

		// merge two nodes and update the node-distance relative to these two nodes
		for (auto iter=node_map.begin();iter!=node_map.end();++iter)
		{
			if((*iter).first!=newNode.index)
			{
				tempVec[current].first = (*iter).first;
				tempVec[current].second = newNode.index;
				tempVec[current].distance=getDistAtNodes(newNode.element,(*iter).second.element, linkageOption);
				if(tempVec[current].distance<minDist)
				{
					target = current;
					minDist = tempVec[current].distance;
				}
				++current;
			}
		}

		if(target>=0 && tempVec.size()>=1)
		{
			poped = tempVec[target];

			/* judge whether current is assigned to right value */
			assert(current==tempVec.size());
			dNodeVec.clear();
			dNodeVec = tempVec;
			tempVec.clear();
			++index;
		}

	}while(node_map.size()!=numberOfClusters);	//merging happens whenever requested cluster is not met

	if(lMethod)	// invoke the l-method to find the optimal number of clusters
	{
		/* perform L-method computation to detect optimal number of AHC */
		DetermClusterNum dcn;
		dcn.iterativeRefinement(dist_map);
		std::cout << "Otimal number of clusters by L-Method is " << dcn.getFinalNumOfClusters() << std::endl;
		dcn.recordLMethodResult(normOption);
	}

	else	// otherwise, just perform the AHC clustering
	{
		nodeVec=std::vector<Ensemble>(node_map.size());
		int tag = 0;
		for(auto iter=node_map.begin();iter!=node_map.end();++iter)
			nodeVec[tag++]=(*iter).second;

		gettimeofday(&end, NULL);
		timeTemp = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;

		activityList.push_back("Hirarchical clustering of norm "+to_string(normOption)+" for "+
							   to_string(numberOfClusters)+" groups takes: ");
		timeList.push_back(to_string(timeTemp)+" s");
		/* task completed, would delete memory contents */
		dNodeVec.clear();
		node_map.clear();
		/* use alpha function to sort the group by its size in ascending order */
		std::sort(nodeVec.begin(), nodeVec.end(), [](const Ensemble& e1, const Ensemble& e2)
		{return e1.element.size()<e2.element.size()||(e1.element.size()==e2.element.size()&&e1.index<e2.index);});
	}
}


/*
 * @brief Set the labels and compute the centroid and cluster related information
 * @details
 *	With generated node information, the cluster size, cluster centroids and candidates belonging to the same cluster
 *	will be determined for further clustering evaluation metric calculation.
 *
 * @param[in] nodeVec The remained node vector
 * @param[out] neighborVec The candidate vector of each cluster to be updated
 * @param[out] storage The size of each cluster to be updated
 * @param[out] centroid The centroid coordinates to be updated
 */
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
				// update the label index for each streamline candidates
				group[eachContainment[i]] = groupID;
			#pragma omp critical
				// update the centroid coordinates of the cluster
				centroid.row(groupID) += ds.dataMatrix.row(eachContainment[i]);
			}
		}
		storage[groupID] = (*iter).element.size();
		centroid.row(groupID)/=eachContainment.size();
		++groupID;
		eachContainment.clear();
	}
}


/*
 * @brief Extract the features and calculate the evaluation metrics for clustering results
 * @details
 * 	Based on the clustering result, the cluster representatives will be extracted first for each cluster based on the
 * 	closest/furthest candidate to the cluster centroid.
 * 	The clustering evaluation metrics will be computed for the quantitative analysis of the clustering result.
 * 	All the information (cluster representatives stored in .vtk file, clustering evaluation metrics stored in readme)
 * 	will be recorded and stored in the designated folders for further batch processing.
 *
 * @param[in] storage The size of each cluster as input as input
 * @param[in] neighborVec The candidate vector for each cluster as input
 * @param[in] centroid The centroid for each cluster as input
 */
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

	/* record labeling information */
	// IOHandler::generateGroups(neighborVec);

	IOHandler::printClusters(ds.dataVec,group,storage,"norm"+to_string(normOption),ds.fullName,ds.dimension);

	struct timeval start, end;
	double timeTemp;

	/* compute the centroid coordinates of each clustered group */

	gettimeofday(&start, NULL);

	vector<vector<float> > closest(numberOfClusters);
	vector<vector<float> > furthest(numberOfClusters);

	/* extract the closest and furthest streamlines to centroid */
#pragma omp parallel for schedule(static) num_threads(8)
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

	// re-assign centroid coordinates to the vector<vector<float>>
	std::vector<std::vector<float> > center_vec(numberOfClusters, vector<float>(Column));
#pragma omp parallel for schedule(static) num_threads(8)
	for (int i = 0; i < center_vec.size(); ++i)
	{
		for (int j = 0; j < Column; ++j)
		{
			center_vec[i][j] = centroid(i,j);
		}
	}

	// calculate the normalized entropy
	float EntropyRatio;
	getEntropyRatio(storage, EntropyRatio);

	std::cout << "Entropy ratio is: " << EntropyRatio << std::endl;

	// record the time for feature extraction
	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec  - start.tv_sec) * 1000000u 
			   + end.tv_usec - start.tv_usec) / 1.e6;
	activityList.push_back("Feature extraction for norm "+to_string(normOption)+ " takes: ");
	timeList.push_back(to_string(timeTemp)+" s");

	// calculate the normalized validity measurement
	ValidityMeasurement vm;
	vm.computeValue(normOption, ds.dataMatrix, group, object, isPBF);
	activityList.push_back("AHC Validity measure is: ");
	stringstream fc_ss;
	fc_ss << vm.f_c;
	timeList.push_back(fc_ss.str());

	// calculate the silhouette, the Gamma statistics and DB index
	gettimeofday(&start, NULL);
	Silhouette sil;
	sil.computeValue(normOption,ds.dataMatrix,ds.dataMatrix.rows(),ds.dataMatrix.cols(),
			         group,object,numberOfClusters, isPBF);
	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec  - start.tv_sec) * 1000000u
			   + end.tv_usec - start.tv_usec) / 1.e6;
	activityList.push_back("Silhouette calculation for norm " +to_string(normOption)+" takes: ");
	timeList.push_back(to_string(timeTemp)+" s");

	std::cout << "Finishing extracting features!" << std::endl;	

	stringstream ss;
	ss << "norm_" << normOption;

	/* measure closest and furthest rotation */
	std::vector<float> closestRotation, furthestRotation;
	const float& closestAverage = getRotation(closest, closestRotation);
	const float& furthestAverage = getRotation(furthest, furthestRotation);

	// record the linkage type, norm option and normalized entropy
	string linkage = getLinkageStr();
	string normStr = getNormStr();
	string entropyStr = getEntropyStr(EntropyRatio);

	// create the .vtk for streamline labels and cluster representatives
	IOHandler::printFeature(ds.dataName+"_AHC_"+linkage+"_closest_"+ss.str()+".vtk", closest, sil.sCluster,
			                closestRotation, ds.dimension);
	IOHandler::printFeature(ds.dataName+"_AHC_"+linkage+"_furthest_"+ss.str()+".vtk", furthest, sil.sCluster,
							furthestRotation, ds.dimension);
	IOHandler::printFeature(ds.dataName+"_AHC_"+linkage+"_centroid_"+ss.str()+".vtk", center_vec, sil.sCluster,ds.dimension);

	IOHandler::printToFull(ds.dataVec, sil.sData, "AHC_SValueLine_"+ss.str(), ds.fullName, ds.dimension);
	IOHandler::printToFull(ds.dataVec, group, sil.sCluster, "AHC_SValueCluster_"+ss.str(), ds.fullName, ds.dimension);

	// generate README for evaluation metrics
	activityList.push_back("numCluster is: ");
	timeList.push_back(to_string(numberOfClusters));

	activityList.push_back("Average Silhouette is: ");
	timeList.push_back(to_string(sil.sAverage));

	activityList.push_back("Average rotation of closest is: ");
	timeList.push_back(to_string(closestAverage));

	activityList.push_back("Average rotation of furthest is: ");
	timeList.push_back(to_string(furthestAverage));

	IOHandler::generateReadme(activityList,timeList);
	IOHandler::writeReadme("Linkage: "+linkage+", "+"norm option is "+normStr+", ");
	IOHandler::writeGroupSize(storage);

	/* print entropy value for the clustering algorithm */
	IOHandler::writeReadme(EntropyRatio, sil, "For norm "+to_string(normOption));
	IOHandler::writeReadme(closestAverage, furthestAverage);

	//curveValue[0].push_back(sil.sAverage);
	//curveValue[1].push_back(sil.gammaStatistic);
	//curveValue[2].push_back(sil.dbIndex);
	//curveValue[3].push_back(vm.f_c);
}


/*
 * @brief Set the data set and perform necessary operations with user parameter input
 * @details
 * 	The function will read in the coordinates of the streamlines/pathlines from the given argument.
 * 	Then it will provide necessary sampling strategy based on user input and data set type
 * 	Furthmore, parameter input will be enforced from the console.
 *
 * @param[in] argc The count of argument string
 * @param[in] argv The char* array for argument
 */
void AHC::setDataset(const int& argc, char **argv)
{
	if(argc!=3)
	{
		std::cout << "Input argument should have 3!" << endl
		          << "./cluster inputFile_name(in dataset folder) "
		          << "data_dimension(3)" << endl;
		exit(1);
	}
	// get the attribute for data set
	ds.strName = string("../dataset/")+string(argv[1]);
	ds.dataName = string(argv[1]);
	ds.dimension = atoi(argv[2]);

	/* get the bool tag for isPBF */
	std::cout << "It is a PBF dataset? 1.Yes, 0.No" << std::endl;
	int PBFjudgement;
	std::cin >> PBFjudgement;
	assert(PBFjudgement==1||PBFjudgement==0);
	isPBF = (PBFjudgement==1);

	/* get the bool tag for isPBF */
	std::cout << "It is a pathline dataset? 1.Yes, 0.No" << std::endl;
	std::cin >> PBFjudgement;
	assert(PBFjudgement==1||PBFjudgement==0);
	isPathlines = (PBFjudgement==1);

	// set up the sample option by user input and data set type (pathlines or streamlines)
	int sampleOption;
	if(isPathlines)	// default direct-repeating for pathlines to match the time steps
		sampleOption = 1;
	else	// streamline sample option can be versatile
	{
		std::cout << "choose a sampling method for the dataset?" << std::endl
				  << "1.directly filling with last vertex; 2. uniform sampling." << std::endl;
		std::cin >> sampleOption;
		assert(sampleOption==1||sampleOption==2);
	}

	// read the coordinates into the member variales
	IOHandler::readFile(ds.strName,ds.dataVec,ds.vertexCount,ds.dimension,ds.maxElements);

	// get the path of full name and print the streamlines vtk
	ds.fullName = ds.strName+"_full.vtk";
	IOHandler::printVTK(ds.fullName, ds.dataVec, ds.vertexCount, ds.dimension);

	// perform sampling operation with user parameters
	if(sampleOption==1)
		IOHandler::expandArray(ds.dataMatrix,ds.dataVec,ds.dimension,ds.maxElements);
	else if(sampleOption==2)
		IOHandler::sampleArray(ds.dataMatrix,ds.dataVec,ds.dimension,ds.maxElements);

	// create the label index for each individual streamline/pathline
	group = std::vector<int>(ds.dataMatrix.rows());

	// whether to activate the L-method or not
	std::cout << "Perform L-method to detect optimal num of clusters? 0: No, 1: Yes! " << std::endl;
	std::cin >> lMethod;
	assert(lMethod==0 || lMethod==1);
	lMethod = (lMethod==1);

	// which linkage type to be selected
	std::cout << "---------------------------" << std::endl;
	std::cout << "Input linkage option: 0.single linkage, 1.complete linkage, 2.average linkage" << std::endl;
	std::cin >> linkageOption;
	assert(linkageOption==0||linkageOption==1||linkageOption==2);

	// lMethod is not activated, so will ask for number of cluster as input
	if(!lMethod)
	{
		std::cout << "---------------------------" << std::endl;
		std::cout << "Choose cluster number input method: 0.user input, 1.read from file: " << std::endl;
		int clusterInput;
		std::cin >> clusterInput;
		assert(clusterInput==0||clusterInput==1);
		readCluster = (clusterInput==1);
	}
}


/*
 * @brief Get the distance between two nodes with a given linkage type
 *
 * @param[in] firstList The first node with candidates
 * @param[in] secondList The second node with candidates
 * @param[in] Linkage The linkage type, 0 for single, 1 for complete and 2 for average
 * @return The distance value between two nodes in AHC clustering
 */
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
			result = -1.0;
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

	case 2:	// average linkage
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

	default:	// no linkage option, should print "Error" information and exit the program
		std::cout << "error!" << std::endl;
		exit(1);
	}
	return result;
}


/*
 * @brief Get the string for linkage type
 * @return A string type for linkage
 */
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


/*
 * @brief Calculate the normalized entropy
 *
 * @param[in] storage The size of different clusters
 * @param[out] EntropyRatio The normalized entropy to be updated
 */
void AHC::getEntropyRatio(const std::vector<int>& storage, float& EntropyRatio)
{
	EntropyRatio = 0;
	const int& Row = ds.dataMatrix.rows();
	for (int i = 0; i < storage.size(); ++i)
	{
		float ratio = float(storage[i])/float(Row);
		EntropyRatio-=ratio*log2f(ratio);
	}
	EntropyRatio/=log2f(storage.size());
}


/*
 * @brief Get the string type of input similarity measure
 * @return A string type
 */
string AHC::getNormStr()
{
	stringstream ss;
	ss << normOption;
	return ss.str();
}


/*
 * @brief Get the string type for entropy value
 *
 * @param[out] EntropyRatio The normalized entropy value
 * @return The string of the float value
 */
string AHC::getEntropyStr(const float& EntropyRatio)
{
	stringstream ss;
	ss << EntropyRatio;
	return ss.str();
}	


/*
 * @brief Set the merged nodes and perform necessary merge operations before the starting of AHC
 *
 * @param[out] dNodeVec The node vector to be updated
 * @param[out] node_map The map to record the index and node
 */
void AHC::setValue_merge(std::vector<DistNode>& dNodeVec, std::unordered_map<int, Ensemble>& node_map)
{
	const int& Row = ds.dataMatrix.rows();

	/* find the node of closest distance */
	std::vector<int> miniNode(Row);
#pragma omp parallel for schedule(static) num_threads(8)
	for(int i=0;i<Row;++i)
	{
		float miniDist = FLT_MAX, dist;
		int index = -1;
		for(int j=0;j<Row;++j)
		{
			if(i==j)
				continue;
			if(distanceMatrix)
				dist = distanceMatrix[i][j];
			else
				dist = getDisimilarity(ds.dataMatrix, i, j, normOption, object);

			if(miniDist>dist)
			{
				miniDist=dist;
				index=j;
			}
		}
		miniNode[i]=index;
	}

	std::vector<bool> isIn(Row, false);

	// update the map for node
	int tag = 0;
	for(int i=0;i<Row;++i)
	{
		if(!isIn[i])
		{
			Ensemble en;
			if(miniNode[miniNode[i]]==i)
			{
				en.element.push_back(i);
				en.element.push_back(miniNode[i]);
				isIn[i]=true;
				isIn[miniNode[i]]=true;
				node_map[tag] = en;
			}
			else
			{
				en.element.push_back(i);
				isIn[i]=true;
				node_map[tag] = en;
			}
			++tag;
		}
	}

	// update the dNodeVec from the newly create nodes
	const int& mapSize = node_map.size();
	dNodeVec = std::vector<DistNode>(mapSize*(mapSize-1)/2);

	tag = 0;
	for(auto start = node_map.begin(); start!=node_map.end(); ++start)
	{
		for(auto end = node_map.begin(); end!=node_map.end() && end!=start; ++end)
		{
			dNodeVec[tag].first = start->first;
			dNodeVec[tag].second = end->first;
			dNodeVec[tag].distance = getDistAtNodes(start->second.element,end->second.element, linkageOption);
			++tag;
		}
	}
	assert(tag==dNodeVec.size());
}


/*
 * @brief Set value for the dNodeVec and node_map as initialization of the AHC procedure
 *
 * @param[out] dNodeVec The vector of nodes to be updated
 * @param[out] node_map The map for recording index and candidate streamlines in the cluster
 */
void AHC::setValue(std::vector<DistNode>& dNodeVec, std::unordered_map<int, Ensemble>& node_map)
{
	const int& Row = ds.dataMatrix.rows();
	dNodeVec = std::vector<DistNode>(Row*(Row-1)/2);
	int tag = 0;
	// record the node i, node j and their distance into the vector
	for(int i=0;i<Row-1;++i)
	{
		for(int j=i+1;j<Row;++j)
		{
			dNodeVec[tag].first = i;
			dNodeVec[tag].second = j;
			if(distanceMatrix)
				dNodeVec[tag].distance = distanceMatrix[i][j];
			else
				dNodeVec[tag].distance = getDisimilarity(ds.dataMatrix, i, j, normOption, object);
			++tag;
		}
	}
	assert(tag==dNodeVec.size());
	// record the index of the node
	for(int i=0;i<Row;++i)
	{
		node_map[i].element.push_back(i);
	}
}


