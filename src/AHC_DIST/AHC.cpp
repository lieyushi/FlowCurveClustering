/*
 * @brief This is the source cpp for the AHC clustering with distance threshold as input
 * @details
 * 	Different from the standard AHC which requires a number of cluster as input, the AHC_DIDS accepts a distance
 * 	threshold as input and merges any two nodes that have the distance value within this range. It is very similar
 * 	to BIRCH clustering algorithm and it is not used in flow visualization. However, we still place it here with
 * 	documentation just in case anyone feels interest in it.
 * @author Lieyu Shi
 */

#include "AHC.h"


/*
 * @brief The default constructor
 */
AHC::AHC() {

}


/*
 * @brief The constructor with parameters
 * @details
 *	It firstly sets up the data set from the argument string and reads in from the local file
 *	then create the distance matrix for AHC clustering
 *
 * @param[in] argc Count of argument
 * @param[in] argv Char* array of argument
 */
AHC::AHC(const int& argc, char **argv) {

	// set the data set and norm option
	setDataset(argc, argv);
	setNormOption();

	/* very hard to decide whether needed to perform such pre-processing but still create a class
	 * object as cached before the pairwise distance matrix computation
	 */
	object = MetricPreparation(ds.dataMatrix.rows(), ds.dataMatrix.cols());
	object.preprocessing(ds.dataMatrix, ds.dataMatrix.rows(), ds.dataMatrix.cols(), normOption);

	/* would store distance matrix instead because it would save massive time */
	struct timeval start, end;
	double timeTemp;
	gettimeofday(&start, NULL);

	// calculate the distance matrix
	getDistanceMatrix(ds.dataMatrix, normOption, object);

	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec
			- start.tv_usec) / 1.e6;
	activityList.push_back("Distance matrix computing takes: ");
	timeList.push_back(to_string(timeTemp) + " s");

	getDistRange();
}


/*
 * @brief The destructor for the class
 */
AHC::~AHC() {
	// delete the distance matrix
	deleteDistanceMatrix(ds.dataMatrix.rows());
}


/*
 * @brief Perform AHC clustering with distance input
 * @details
 * 	It will select the type of AHC clustering, either by input of a cluster number, or a distance threshold. Then
 * 	perform the hiararchical clustering by bottom-up merge of the tree. Then posterior calculation on feature
 * 	extraction and clustering evaluation is performed for quantitative and visual analysis
 *
 */
void AHC::performClustering() {

	std::vector<Ensemble> nodeVec;

	/* perform hierarchical clustering */
	std::cout << "---------------------" << std::endl;
	std::cout
			<< "1. clustering by a fixed group, 2. clustering by a distance threshold."
			<< std::endl;
	int clusteringOption;
	std::cin >> clusteringOption;
	assert(clusteringOption == 1 || clusteringOption == 2);

	// choose AHC by distance threshold or fixed group
	if (clusteringOption == 1)
		bottomUp_byGroup(nodeVec);
	else if (clusteringOption == 2)
		bottomUp_byThreshold(nodeVec);

	vector<vector<int>> neighborVec(numberOfClusters);
	// element size for all groups
	vector<int> storage(numberOfClusters);

	// geometric center
	Eigen::MatrixXf centroid = Eigen::MatrixXf::Zero(numberOfClusters,
			ds.dataMatrix.cols());

	// set label information
	setLabel(nodeVec, neighborVec, storage, centroid);

	nodeVec.clear();

	extractFeatures(storage, neighborVec, centroid);
}


/*
 * @breif Perform bottom-up clustering with cluster number as input
 * @details
 * 	Iteratively merge the hierarchical tree until the tree size is approximately close to the required number
 *
 * @param[out] nodeVec A vector for Ensemble objects
 */
void AHC::bottomUp_byGroup(std::vector<Ensemble>& nodeVec) {
	const int& Row = ds.dataMatrix.rows();
	std::cout << "-------------------------------------------------------------------------------" << std::endl;
	std::cout << "Expected number of clusters from [0, " << Row << "]:";
	std::cin >> expectedClusters;
	assert(expectedClusters > 0 && expectedClusters < Row / 10);

	/* would store distance matrix instead because it would save massive time */
	struct timeval start, end;
	double timeTemp;
	gettimeofday(&start, NULL);

	// perform the clustering algorithm until the approximate clusters are reached
	int clusterCount = 0;
	const int& minExpected = 0.8 * expectedClusters;
	const int& maxExpected = 1.2 * expectedClusters;
	float minDist = distRange[0], maxDist = distRange[1] / 4.0;
	int iteration = 0;
	std::cout << ".." << std::endl;
	std::cout << ".." << std::endl;
	std::cout << "Binary search starts!" << std::endl;
	while (true && iteration <= 20) {
		distanceThreshold = (minDist + maxDist) / 2.0;
		hierarchicalMerging(nodeVec);
		clusterCount = nodeVec.size();
		std::cout << "Iteration " << (++iteration) << " finds " << clusterCount
				<< " groups!" << std::endl;
		if (clusterCount >= minExpected && clusterCount <= maxExpected)
			break;
		else if (clusterCount < minExpected)
			maxDist = distanceThreshold;
		else
			minDist = distanceThreshold;
	}

	// record relevant information
	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec
			- start.tv_usec) / 1.e6;
	stringstream ss;
	ss << expectedClusters;
	string cluster_str = ss.str();

	ss.str("");
	ss << iteration;

	activityList.push_back("To achieve " + cluster_str + " groups will take " + ss.str()
					+ " binary search and take: ");
	timeList.push_back(to_string(timeTemp) + " s");

}


/*
 * @brief Perform the bottom-up clustering by distance threshold of the linkage type
 *
 * @param[out] nodeVec The Ensemble vector to be updated
 */
void AHC::bottomUp_byThreshold(std::vector<Ensemble>& nodeVec) {
	std::cout
			<< "-------------------------------------------------------------------------------"
			<< std::endl;
	std::cout << "Input threshold: [" << distRange[0] << "," << distRange[1]
			<< "]: ";
	std::cin >> distanceThreshold;
	assert(distanceThreshold > distRange[0] && distanceThreshold < distRange[1]);

	/* would store distance matrix instead because it would save massive time */
	struct timeval start, end;
	double timeTemp;
	gettimeofday(&start, NULL);

	// merge the nodes
	hierarchicalMerging(nodeVec);

	// record relevant information
	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec
			- start.tv_usec) / 1.e6;
	stringstream ss;
	ss << distanceThreshold;
	activityList.push_back("To cluster by distance " + ss.str() + " will take: ");
	timeList.push_back(to_string(timeTemp) + " s");
}


/*
 * @brief Perform hierarchical merging for AHC if the merged distance is lower than the distance threshold
 *
 * @param[out] nodeVec The Ensemble vector to be updated
 */
void AHC::hierarchicalMerging(std::vector<Ensemble>& nodeVec) {
	const int Row = ds.dataMatrix.rows();

	nodeVec.clear();
	//could have used vector, but since there're too many operations inside so should use set
	nodeVec = std::vector<Ensemble>(Row);
	//create node in forest structure
#pragma omp parallel for schedule(static) num_threads(8)
	for (int i = 0; i < nodeVec.size(); ++i) {
		nodeVec[i].index = i;
		nodeVec[i].element.push_back(i);
	}

	//two iterators to record positions of set
	std::vector<Ensemble>::iterator iter_i, iter_j;

	//vector to store new node
	std::vector<Ensemble> newNodeList;
	do {
		//insert new node obtained from previous step
		if (!newNodeList.empty()) {
			nodeVec.insert(nodeVec.end(), newNodeList.begin(),
					newNodeList.end());
			newNodeList.clear();
		}
		//iter_i prior, iter_j consecutive
		iter_i = nodeVec.begin();
		iter_j = iter_i;
		++iter_j;
		int mergedCount = 0;
		while (iter_i != nodeVec.end()) {
			// j reaches end or i is merged, should move i forward
			if (iter_j == nodeVec.end() || (*iter_i).merged) {
				++iter_i;
				iter_j = iter_i;
				++iter_j;
			}
			// j node already merged, so no longer consideration
			else if ((*iter_j).merged)
				++iter_j;
			// move j and calculate distance for mutual pairs
			else {
				//compute distance between two clusters by single/complete/average linkages
				const float& linkageDist = getDistAtNodes((*iter_i).element,
						(*iter_j).element, linkageOption);

				//larger distance than threshold, then move forward
				if (linkageDist > distanceThreshold)
					++iter_j;
				//merge two clusters into one cluster if smaller than threshold
				else {
					//add merged node whose index is total element size
					vector<int> first = (*iter_i).element, second =
							(*iter_j).element;
					Ensemble newNode = Ensemble(first.size() + second.size());
					newNode.element = first;
					newNode.element.insert(newNode.element.begin(),
							second.begin(), second.end());
					newNodeList.push_back(newNode);

					(*iter_i).merged = true;
					(*iter_j).merged = true;

					++iter_i;
					iter_j = iter_i;
					++iter_j;

					mergedCount += 2;
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
		assert(nodeVec.size() >= mergedCount);
		std::vector<Ensemble> copyNode(nodeVec.size() - mergedCount);
		int c_i = 0;
		for (int i = 0; i < nodeVec.size(); ++i) {
			if (!nodeVec[i].merged)
				copyNode[c_i++] = nodeVec[i];
		}
		nodeVec.clear();
		nodeVec = copyNode;
		copyNode.clear();

		mergedCount = 0;

	} while (!newNodeList.empty());	//merging happens constantly

	newNodeList.clear();

	numberOfClusters = nodeVec.size();

	/* use alpha function to sort the group by its size */
	std::sort(nodeVec.begin(), nodeVec.end(), [](const Ensemble& e1, const Ensemble& e2)
			{	return e1.element.size()<e2.element.size() ||(e1.element.size()==e2.element.size()&&e1.index<e2.index);});
}


/*
 * @brief Set label for streamlines from the hierarchical tree
 *
 * @param[in] nodeVec The Ensemble vector as input
 * @param[out] neighborVec The candidate vector for each cluster
 * @param[out] storage The size of clusters to be updated
 * @param[out] centroid The centroids of streamlines to be updated
 */
void AHC::setLabel(const std::vector<Ensemble>& nodeVec, vector<vector<int> >& neighborVec,
		vector<int>& storage, Eigen::MatrixXf& centroid)
{
// group tag by increasing order
	int groupID = 0;

	// element list for each group
	vector<int> eachContainment;

	// find group id and neighboring vec
	for (auto iter = nodeVec.begin(); iter != nodeVec.end(); ++iter)
	{
		eachContainment = (*iter).element;
		neighborVec[groupID] = eachContainment;
	#pragma omp parallel num_threads(8)
		{
		#pragma omp for nowait
			for (int i = 0; i < eachContainment.size(); ++i) {
				group[eachContainment[i]] = groupID;
			#pragma omp critical
				centroid.row(groupID) += ds.dataMatrix.row(eachContainment[i]);
			}
		}
		storage[groupID] = (*iter).element.size();
		centroid.row(groupID) /= eachContainment.size();
		++groupID;
		eachContainment.clear();
	}
}


/*
 * @brief Extract the features and compute the evaluation metrics
 *
 * @param[in] storage size of clusters as input
 * @param[in] neighborVec candidate vectors for each cluster
 * @param[in] centroid The centroids for each cluster
 */
void AHC::extractFeatures(const std::vector<int>& storage, const std::vector<std::vector<int> >& neighborVec,
		const Eigen::MatrixXf& centroid)
{
	const int& Row = ds.dataMatrix.rows();
	const int& Column = ds.dataMatrix.cols();

	std::cout << "Final group number information: " << std::endl;
	for (int i = 0; i < storage.size(); ++i) {
		std::cout << storage[i] << " ";
	}
	std::cout << std::endl;

	/* record labeling information */
	// IOHandler::generateGroups(neighborVec);

	IOHandler::printClusters(ds.dataVec, group, storage, "norm" + to_string(normOption),
			ds.fullName, ds.dimension);

	struct timeval start, end;
	double timeTemp;

	// calculate the evaluation metrics
	gettimeofday(&start, NULL);
	Silhouette sil;
	sil.computeValue(normOption, ds.dataMatrix, ds.dataMatrix.rows(), ds.dataMatrix.cols(), group,
			object, numberOfClusters, isPBF, neighborVec);
	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec
			- start.tv_usec) / 1.e6;
	activityList.push_back("Silhouette calculation takes: ");
	timeList.push_back(to_string(timeTemp) + " s");

	/* compute the centroid coordinates of each clustered group */
	gettimeofday(&start, NULL);

	vector<vector<float> > closest(numberOfClusters);
	vector<vector<float> > furthest(numberOfClusters);

	/* extract the closest and furthest streamlines to centroid */

#pragma omp parallel for schedule(static) num_threads(8)
	for (int i = 0; i < numberOfClusters; ++i) {
		float minDist = FLT_MAX;
		float maxDist = -10;
		int minIndex = -1, maxIndex = -1;
		const std::vector<int>& groupRow = neighborVec[i];
		const Eigen::VectorXf& eachCentroid = centroid.row(i);
		for (int j = 0; j < groupRow.size(); ++j) {
			float distance = getDisimilarity(eachCentroid, ds.dataMatrix,
					groupRow[j], normOption, object);
			if (minDist > distance) {
				minDist = distance;
				minIndex = groupRow[j];
			}
			if (maxDist < distance) {
				maxDist = distance;
				maxIndex = groupRow[j];
			}
		}
		closest[i] = ds.dataVec[minIndex];
		furthest[i] = ds.dataVec[maxIndex];
	}

	std::vector<std::vector<float> > center_vec(numberOfClusters,
			vector<float>(Column));
#pragma omp parallel for schedule(static) num_threads(8)
	for (int i = 0; i < center_vec.size(); ++i) {
		for (int j = 0; j < Column; ++j) {
			center_vec[i][j] = centroid(i, j);
		}
	}

	// calculate the normalized entropy ratio
	float EntropyRatio;
	getEntropyRatio(storage, EntropyRatio);

	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec
			- start.tv_usec) / 1.e6;
	activityList.push_back("Feature extraction takes: ");
	timeList.push_back(to_string(timeTemp) + " s");

	// calculate the normalized validity measurement
	ValidityMeasurement vm;
	vm.computeValue(normOption, ds.dataMatrix, group, object, isPBF);
	activityList.push_back("AHC Validity measure is: ");
	stringstream fc_ss;
	fc_ss << vm.f_c;
	timeList.push_back(fc_ss.str());

	std::cout << "Finishing extracting features!" << std::endl;

	// record relevant information
	stringstream ss;
	ss << "norm_" << normOption;

	string linkage = getLinkageStr();
	string normStr = getNormStr();

	// print the featured information as result
	IOHandler::printFeature(
			ds.dataName + "_AHC_Dist_" + linkage + "_closest_" + ss.str()
					+ ".vtk", closest, sil.sCluster, ds.dimension);
	IOHandler::printFeature(
			ds.dataName + "_AHC_Dist_" + linkage + "_furthest_" + ss.str()
					+ ".vtk", furthest, sil.sCluster, ds.dimension);
	IOHandler::printFeature(
			ds.dataName + "_AHC_Dist_" + linkage + "_centroid_" + ss.str()
					+ ".vtk", center_vec, sil.sCluster, ds.dimension);

	IOHandler::printToFull(ds.dataVec, sil.sData,
			"AHC_Dist_SValueLine_" + ss.str(), ds.fullName, ds.dimension);
	IOHandler::printToFull(ds.dataVec, group, sil.sCluster,
			"AHC_Dist_SValueCluster_" + ss.str(), ds.fullName, ds.dimension);

	// generate necessary readme file
	activityList.push_back("numCluster is: ");
	timeList.push_back(std::to_string(numberOfClusters));

	IOHandler::generateReadme(activityList, timeList);

	IOHandler::writeReadme(
			"Linkage: " + linkage + ", " + "norm option is " + normStr);

	IOHandler::writeGroupSize(storage);

	/* print entropy value for the clustering algorithm */
	IOHandler::writeReadme(EntropyRatio, sil, "For norm "+std::to_string(normOption));

	/* measure closest and furthest rotation */
	std::vector<float> closestRot, furthestRot;
	const float& closestAverage = getRotation(closest, closestRot);
	const float& furthestAverage = getRotation(furthest, furthestRot);

	IOHandler::writeReadme(closestAverage, furthestAverage);
}


/*
 * @brief Set data set and perform necessary operations with user parameters
 *
 * @param[in] argc Count of arguments
 * @param[in] argv Char* array of arguments
 */
void AHC::setDataset(const int& argc, char **argv) {
	if (argc != 3) {
		std::cout << "Input argument should have 3!" << endl
				<< "./cluster inputFile_name(in dataset folder) "
				<< "data_dimension(3)" << endl;
		exit(1);
	}
	ds.strName = string("../dataset/") + string(argv[1]);
	ds.dataName = string(argv[1]);
	ds.dimension = atoi(argv[2]);

	/* get the bool tag for isPBF */
	std::cout << "It is a PBF dataset? 1.Yes, 0.No" << std::endl;
	int PBFjudgement;
	std::cin >> PBFjudgement;
	assert(PBFjudgement == 1 || PBFjudgement == 0);
	isPBF = (PBFjudgement == 1);

	// set the sampling option
	int sampleOption;
	std::cout << "choose a sampling method for the dataset?" << std::endl
			<< "1.directly filling with last vertex; 2. uniform sampling."
			<< std::endl;
	std::cin >> sampleOption;
	assert(sampleOption == 1 || sampleOption == 2);

	// read from the file
	IOHandler::readFile(ds.strName, ds.dataVec, ds.vertexCount, ds.dimension,
			ds.maxElements);

	ds.fullName = ds.strName + "_full.vtk";
	IOHandler::printVTK(ds.fullName, ds.dataVec, ds.vertexCount, ds.dimension);

	// perform sampling
	if (sampleOption == 1)
		IOHandler::expandArray(ds.dataMatrix, ds.dataVec, ds.dimension,
				ds.maxElements);
	else if (sampleOption == 2)
		IOHandler::sampleArray(ds.dataMatrix, ds.dataVec, ds.dimension,
				ds.maxElements);

	group = std::vector<int>(ds.dataMatrix.rows());

	// choose linkage type
	std::cout << "---------------------------" << std::endl;
	std::cout
			<< "Input linkage option: 0.single linkage, 1.complete linkage, 2.average linkage"
			<< std::endl;
	std::cin >> linkageOption;
	assert(linkageOption == 0 || linkageOption == 1 || linkageOption == 2);
}


/*
 * @brief Set norm option
 */
void AHC::setNormOption() {
	std::cout << "Input a norm option 0-12!" << std::endl;
	std::cin >> normOption;
	std::cout << std::endl;

	// choose distance metrics according to number
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
	bool found = false;
	for (int i = 0; i < 16 && !found; ++i) {
		if (normOption == i) {
			found = true;
			break;
		}
	}
	if (!found) {
		std::cout << "Cannot find the norm!" << std::endl;
		exit(1);
	}
}


/*
 * @brief Calculate the range of distance matrix
 */
void AHC::getDistRange() {
	const float& Percentage = 0.05;
	const int& Row = ds.dataMatrix.rows();

	distRange = vector<float>(2);
	distRange[0] = FLT_MAX, distRange[1] = FLT_MIN;
	const int& totalSize = int(Percentage * Row);
#pragma omp parallel num_threads(8)
	{
	#pragma omp for nowait
		for (int i = 0; i < totalSize; ++i) {
			float tempDist, i_min = FLT_MAX, i_max = FLT_MIN;
			for (int j = 0; j < Row; ++j) {
				if (distanceMatrix)
					tempDist = distanceMatrix[i][j];
				else
					tempDist = getDisimilarity(ds.dataMatrix, i, j, normOption,
							object);
				if (tempDist < i_min) {
					i_min = tempDist;
				}
				if (tempDist > i_max) {
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

	std::cout << "Distance threshold is: [" << distRange[0] << ", " << distRange[1] << "]."
			<< std::endl;
}


/*
 * @brief Get distance between nodes by linkage type
 *
 * @param[in] firstList The first node
 * @param[in] secondList The second node
 * @param[in] Linkage The linkage type
 * @return A float value for the distance
 */
const float AHC::getDistAtNodes(const vector<int>& firstList, const vector<int>& secondList,
		const int& Linkage)
{
	const int& m = firstList.size();
	const int& n = secondList.size();
	assert(m != 0);
	assert(n != 0);
	/* 0: single linkage, min(x_i,y_j)
	 * 1: complete linkdage, max(x_i,y_j)
	 * 2: average linkage, sum/x_i/y_j
	 */
	float result, value;
	switch (Linkage)
	{
	case 0:	//single linkage
		{
			result = FLT_MAX;
		#pragma omp parallel for reduction(min:result) num_threads(8)
			for (int i = 0; i < m; ++i) {
				for (int j = 0; j < n; ++j) {
					if (distanceMatrix)
						value = distanceMatrix[firstList[i]][secondList[j]];
					else
						value = getDisimilarity(ds.dataMatrix, firstList[i],
								secondList[j], normOption, object);
					result = std::min(result, value);
				}
			}
		}
		break;

	case 1:	//complete linkage
		{
			result = FLT_MIN;
		#pragma omp parallel for reduction(max:result) num_threads(8)
			for (int i = 0; i < m; ++i) {
				for (int j = 0; j < n; ++j) {
					if (distanceMatrix)
						value = distanceMatrix[firstList[i]][secondList[j]];
					else
						value = getDisimilarity(ds.dataMatrix, firstList[i],
								secondList[j], normOption, object);
					result = std::max(result, value);
				}
			}
		}
		break;

	case 2: 	// average linkage
		{
			result = 0;
		#pragma omp parallel for reduction(+:result) num_threads(8)
			for (int i = 0; i < m; ++i) {
				for (int j = 0; j < n; ++j) {
					if (distanceMatrix)
						value = distanceMatrix[firstList[i]][secondList[j]];
					else
						value = getDisimilarity(ds.dataMatrix, firstList[i],
								secondList[j], normOption, object);
					result += value;
				}
			}
			result /= m * n;
		}
		break;

	default:	// error
		std::cout << "error!" << std::endl;
		exit(1);
	}
	return result;
}


/*
 * @brief Get string for the linkage type
 * @return String of linkage type
 */
string AHC::getLinkageStr()
{
	string result;
	switch (linkageOption)
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
 * @brief Get entropy ratio
 * @param[in] storage Size of clusters
 * @param[out] EntropyRatio The normalized entropy to be updated
 */
void AHC::getEntropyRatio(const std::vector<int>& storage, float& EntropyRatio)
{
	EntropyRatio = 0;
	const int& Row = ds.dataMatrix.rows();
	for (int i = 0; i < storage.size(); ++i) {
		float ratio = float(storage[i]) / float(Row);
		EntropyRatio -= ratio * log2f(ratio);
	}
	EntropyRatio /= log2f(storage.size());
}


/*
 * @brief Get string for norm
 * @return String of norm
 */
string AHC::getNormStr() {
	stringstream ss;
	ss << normOption;
	return ss.str();
}


/*
 * @brief Get the string for entropy value
 * @param[in] EntropyRatio The entropy value
 * @return String type
 */
string AHC::getEntropyStr(const float& EntropyRatio)
{
	stringstream ss;
	ss << EntropyRatio;
	return ss.str();
}
