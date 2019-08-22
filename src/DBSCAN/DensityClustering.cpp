/*
 * @brief The source cpp for implementing the clustering algorithm DBSCAN for flow visualization
 * @author Lieyu Shi
 */

#include "DensityClustering.h"

std::vector<string> activityList;
std::vector<string> timeList;


/*
 * @brief The ratio w.r.t the maximal distance in (0,1)
 */
float multiTimes;

/*
 * @brief The minPts to decide the neighbor number
 */
int minPts;


/*
 * @brief DBSCAN constructor with parameters
 *
 * @param[in] argc The count of arguments
 * @param[in] argv The char* array with data sets
 */
DensityClustering::DensityClustering(const int& argc, char **argv)
{
	struct timeval start, end;
	double timeTemp;
	gettimeofday(&start, NULL);

	// set data set and norm option
	setDataset(argc, argv);
	setNormOption();

	// create the object for distance matrix computation
	object = MetricPreparation(ds.dataMatrix.rows(), ds.dataMatrix.cols());
	object.preprocessing(ds.dataMatrix, ds.dataMatrix.rows(),
			ds.dataMatrix.cols(), normOption);

	/* if the dataset is not PBF, then should record distance matrix for Gamma matrix compution */
	if (!isPBF)
	{
		deleteDistanceMatrix(ds.dataMatrix.rows());

		std::ifstream distFile(("../dataset/"+to_string(normOption)).c_str(), ios::in);
		if(distFile.fail())	// the distance matrix file not exists, will create new one
		{
			distFile.close();
			getDistanceMatrix(ds.dataMatrix, normOption, object);
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
		else	// distance matrix file exists, directly read from the file
		{
			std::cout << "read distance matrix..." << std::endl;

			distanceMatrix = new float*[ds.dataMatrix.rows()];
		#pragma omp parallel for schedule(static) num_threads(8)
			for (int i = 0; i < ds.dataMatrix.rows(); ++i)
			{
				distanceMatrix[i] = new float[ds.dataMatrix.rows()];
			}

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
						distanceMatrix[i][j] = 0;
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
	}

	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec
			- start.tv_usec) / 1.e6;
	activityList.push_back("Computing distance matrix for norm "+to_string(normOption)+" takes: ");
	timeList.push_back(to_string(timeTemp) + " s");

	nodeVec = vector<PointNode>(ds.dataMatrix.rows(), PointNode());
}


/*
 * @brief Destructor
 */
DensityClustering::~DensityClustering() {

}


/*
 * @brief Perform DBSCAN clustering
 */
void DensityClustering::performClustering() {

	struct timeval start, end;
	double timeTemp;
	gettimeofday(&start, NULL);

	// read in the minPts as a parameter
	minPts = setMinPts();
	float distThreshold = getDistThreshold(minPts);

	// perform DBSCAN clustering
	DBSCAN(distThreshold, minPts);

	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec
			- start.tv_usec) / 1.e6;
	activityList.push_back("DBSCAN clustering for norm "+to_string(normOption)+" takes: ");
	timeList.push_back(to_string(timeTemp) + " s");

	// extract features and calculate the clustering evaluation metrics
	extractFeatures(distThreshold, minPts);
}


/*
 * @brief Get the distance threshold
 *
 * @param[in] minPts The min points for the DBSCAN
 */
const float DensityClustering::getDistThreshold(const int& minPts)
{
	int distOption = 2;	// set the default parameter type
	/*
	 std::cout << "Choose distThreshold setup option: 1.user input, 2.minPts-th dist." << std::endl;
	 std::cin >> distOption;*/

	assert(distOption == 1 || distOption == 2);

	if (distOption == 1) 	// if input for the radius, should let the user know the distance range
	{
		float minDist, maxDist;
		getDistRange(minDist, maxDist);
		std::cout << "Distance range is [" << minDist << ", " << maxDist << "]."
				<< std::endl;
		multiTimes = setTimesMin(minDist, maxDist);
		return multiTimes * maxDist;
	} else if (distOption == 2) {	// otherwise, directly get the average distance
		/* should be pointed as average distance of minPts-th dist */
		return getAverageDist(minPts);
	}
}


/*
 * @brief Compute the minPts-th dist for all candidates
 *
 * @param[in] minPts The min point numbers for the DBSCAN clustering
 */
const float DensityClustering::getAverageDist(const int& minPts)
{
	float result = 0.0;
	const int& rowSize = ds.dataMatrix.rows();
#pragma omp parallel num_threads(8)
	{
	#pragma omp for nowait
		for (int i = 0; i < rowSize; ++i) {
			/* use a priority_queue<float> with n*logk time complexity */
			std::priority_queue<float> minDistArray;
			float tempDist;
			for (int j = 0; j < rowSize; ++j) {
				if (i == j)
					continue;
				if (distanceMatrix)
					tempDist = distanceMatrix[i][j];
				else
					tempDist = getDisimilarity(ds.dataMatrix.row(i),
							ds.dataMatrix.row(j), i, j, normOption, object);

				minDistArray.push(tempDist);
				if (minDistArray.size() > minPts)
					minDistArray.pop();
			}

		#pragma omp critical
			result += minDistArray.top();
		}
	}
	return result / rowSize;
}


/*
 * @brief Perform the DBSCAN clustering with given parameters
 *
 * @param[in] radius_eps The radius to check the neighboring information
 * @param[in] minPts The minPts for the DBSCAN clustering
 */
void DensityClustering::DBSCAN(const float& radius_eps, const int& minPts) {
	int C = 0;
	for (int i = 0; i < ds.dataMatrix.rows(); ++i)
	{
		if (nodeVec[i].visited)
			continue;
		nodeVec[i].visited = true;
		vector<int> neighbor = regionQuery(i, radius_eps);
		if (neighbor.size() < minPts)
			nodeVec[i].type = NOISE;
		else
		{
			expandCluster(i, neighbor, C, radius_eps, minPts);
			++C;
		}
	}
}


/*
 * @brief Expand the cluster with candidates that lie within range of the target
 *
 * @param[in] index The index of candidate streamlines
 * @param[in] neighbor The neighborhood candidates found
 * @param[in] cluster_id The label for the clusters
 * @param[in] radius_eps The radius for searching around the neighborhood
 * @param[in] minPts The min number of points for the DBSCAN clustering
 */
void DensityClustering::expandCluster(const int& index, vector<int>& neighbor,
	const int& cluster_id, const float& radius_eps, const int& minPts)
{
	nodeVec[index].group = cluster_id;
	int insideElement;
	for (int i = 0; i < neighbor.size(); ++i) {
		insideElement = neighbor[i];
		if (!nodeVec[insideElement].visited) {
			nodeVec[insideElement].visited = true;
			vector<int> newNeighbor = regionQuery(insideElement, radius_eps);
			if (newNeighbor.size() >= minPts) {
				neighbor.insert(neighbor.end(), newNeighbor.begin(),
						newNeighbor.end());
			}
		}
		if (nodeVec[insideElement].group == -1)
			nodeVec[insideElement].group = cluster_id;
	}
}


/*
 * @brief Perform the region-based query for the target streamline within a given radius
 *
 * @param[in] index The index of the target streamlines
 * @param[in] radius_eps The radius of neighborhood checking
 * @return A vector<int> object that contains the region candidates
 */
const vector<int> DensityClustering::regionQuery(const int& index,
		const float& radius_eps) {
	vector<int> neighborArray;
	neighborArray.push_back(index);
	float tempDist;
	for (int i = 0; i < ds.dataMatrix.rows(); ++i) {
		if (i == index)
			continue;

		/* in case somebody uses distance matrix */
		if (distanceMatrix)
			tempDist = distanceMatrix[index][i];
		else
			tempDist = getDisimilarity(ds.dataMatrix.row(index),
					ds.dataMatrix.row(i), index, i, normOption, object);
		if (tempDist <= radius_eps)
			neighborArray.push_back(i);
	}
	return neighborArray;
}


/*
 * @brief Set the data set from the arguments
 *
 * @param[in] argc Count of arguments
 * @param[in] argv The argument string type
 */
void DensityClustering::setDataset(const int& argc, char **argv) {
	if (argc != 3) {
		std::cout << "Input argument should have 3!" << endl
				<< "./cluster inputFile_name(in dataset folder) "
				<< "data_dimension(3)" << endl;
		exit(1);
	}
	ds.strName = string("../dataset/") + string(argv[1]);
	ds.dimension = atoi(argv[2]);

	/* get the bool tag for isPBF */
	std::cout << "It is a PBF dataset? 1.Yes, 0.No" << std::endl;
	int PBFjudgement;
	std::cin >> PBFjudgement;
	assert(PBFjudgement == 1 || PBFjudgement == 0);
	isPBF = (PBFjudgement == 1);

	// check whether it is pathlines or not
	std::cout << "It is a pathlines dataset? 1.Yes, 0.No" << std::endl;
	std::cin >> PBFjudgement;
	assert(PBFjudgement == 1 || PBFjudgement == 0);
	isPathlines = (PBFjudgement == 1);

	// decide the sampling strategy and operation for the given data sets
	int sampleOption;

	if(isPathlines)
		sampleOption = 1;
	else
	{
		std::cout << "choose a sampling method for the dataset?" << std::endl
				<< "1.directly filling with last vertex; 2. uniform sampling (recommended!); "
			    << "3. equal-arc sampling. " << std::endl;
		std::cin >> sampleOption;
	}
	assert(sampleOption == 1 || sampleOption == 2 || sampleOption == 3);

	// read the coordinates from the file
	IOHandler::readFile(ds.strName, ds.dataVec, ds.vertexCount, ds.dimension, ds.maxElements);

	ds.fullName = ds.strName + "_differentNorm_full.vtk";
	IOHandler::printVTK(ds.fullName, ds.dataVec, ds.vertexCount, ds.dimension);

	if (sampleOption == 1)
		IOHandler::expandArray(ds.dataMatrix, ds.dataVec, ds.dimension, ds.maxElements);
	else if (sampleOption == 2)
		IOHandler::sampleArray(ds.dataMatrix, ds.dataVec, ds.dimension, ds.maxElements);
	else if (sampleOption == 3)
		IOHandler::uniformArcSampling(ds.dataMatrix, ds.dataVec, ds.dimension, ds.maxElements);
}


/*
 * @brief Read the norm option from user input
 */
void DensityClustering::setNormOption() {

	if(isPathlines)
	{
		std::cout << "Choose a norm from 0-17!" << std::endl;
		std::cin >> normOption;
		assert(normOption>=0 && normOption<=17);
	}
	else
	{
		std::cout << "Choose a norm from 0-15!" << std::endl;
		std::cin >> normOption;
		assert(normOption>=0 && normOption<=15);
	}

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
}


/*
 * @brief Calculate the minimal and maximal distance range for the user input of radius
 *
 * @param[in] minDist The minimal distance to be updated
 * @param[in] maxDist The maximal distance to be updated
 */
void DensityClustering::getDistRange(float& minDist, float& maxDist)
{
	const float& Percentage = 0.1;
	const int& Rows = ds.dataMatrix.rows();
	const int& chosen = int(Percentage * Rows);
	minDist = FLT_MAX;
	maxDist = -1.0;
#pragma omp parallel num_threads(8)
	{
	#pragma omp for nowait
		for (int i = 0; i < chosen; ++i) {
			float tempDist;
			for (int j = 0; j < Rows; ++j) {
				if (i == j)
					continue;
				if (distanceMatrix)
					tempDist = distanceMatrix[i][j];
				else
					tempDist = getDisimilarity(ds.dataMatrix.row(i),
							ds.dataMatrix.row(j), i, j, normOption, object);
				#pragma omp critical
				{
					if (tempDist < minDist)
						minDist = tempDist;
					if (tempDist > maxDist)
						maxDist = tempDist;
				}
			}
		}
	}
	std::cout << minDist << " " << maxDist << std::endl;
}


/*
 * @brief Set the minPts parameter, default value 6 is already enough for our paper
 */
const int DensityClustering::setMinPts() {
	/*std::cout << std::endl;
	std::cout << "Input the minPts for DBSCAN in [0" << ", "
			<< ds.dataMatrix.rows() << "], 6 preferred: " << std::endl;*/
	int minPts = 6;
	//std::cin >> minPts;
	if (minPts <= 0 || minPts >= ds.dataMatrix.rows()) {
		std::cout << "Error for out-of-range minPts!" << std::endl;
		exit(1);
	}
	return minPts;
}


/*
 * @brief Select a ratio between [0,1] for the radius input
 *
 * @param[in] minDist The minimal distance value as input
 * @param[in] maxDist The maximal distance value as input
 * @return A radius for the DBSCAN clustering after the user input
 */
const float DensityClustering::setTimesMin(const float& minDist, const float& maxDist)
{
	std::cout << std::endl;
	float lowerBound = minDist / maxDist;
	std::cout << "Input the multiplication for DBSCAN radius in [" << lowerBound
			<< ",1.0]:" << std::endl;
	float multiTimes;
	std::cin >> multiTimes;
	if (multiTimes >= 1.0 || multiTimes <= lowerBound) {
		std::cout << "Error for out-of-range minPts!" << std::endl;
		exit(1);
	}
	return multiTimes;
}


/*
 * @brief Extract features for the clustering results and calculate the clustering evaluation metrics
 *
 * @param[in] radius_eps The radius for DBSCAN clustering
 * @param[in] minPts The minPts parameter for DBSCAN
 */
void DensityClustering::extractFeatures(const float& radius_eps, const int& minPts)
{
	// find the maximal cluster labels with openmp critical operation, could be disabled
	int maxGroup = -INT_MAX + 1;
#pragma omp parallel num_threads(8)
	{
	#pragma omp for nowait
		for (int i = 0; i < nodeVec.size(); ++i) {
			int groupID = nodeVec[i].group;
		#pragma omp critical
			{
				if (groupID != -1 && groupID > maxGroup)
					maxGroup = groupID;
			}
		}
	}
	std::cout << "Max group is: " << maxGroup << std::endl;

	/* re-index the group id by increasing number */
	int numClusters = maxGroup + 1;
	std::vector<int> container(numClusters, 0);
	for (int i = 0; i < nodeVec.size(); ++i) {
		if (nodeVec[i].group != -1)
			++container[nodeVec[i].group];
	}

	int increasingOrder[numClusters];
	std::multimap<int, int> groupMap;

	for (int i = 0; i < numClusters; ++i)
		groupMap.insert(std::pair<int, int>(container[i], i));

	std::fill(container.begin(), container.end(), 0);
	int groupNo = 0;
	for (std::multimap<int, int>::iterator it = groupMap.begin();
			it != groupMap.end(); ++it) {
		if (it->first > 0) {
			increasingOrder[it->second] = groupNo;
			container[groupNo] = it->first;
			++groupNo;
		}
	}

	numClusters = groupNo + 1; /* plus -1 as group */

#pragma omp parallel for schedule(static) num_threads(8)
	for (int i = 0; i < nodeVec.size(); ++i) {
		if (nodeVec[i].group != -1)
			nodeVec[i].group = increasingOrder[nodeVec[i].group];
	}

	/* in case -1, we use 0 to record number of -1 as noise */
	std::vector<int> item_cids(nodeVec.size());
	std::vector<std::vector<int> > storage(numClusters);
	/* -1 group as group[0] */
	for (int i = 0; i < nodeVec.size(); ++i) {
		item_cids[i] = nodeVec[i].group;
		storage[nodeVec[i].group + 1].push_back(i);
	}

	container.insert(container.begin(), storage[0].size());

	/* compute balanced Entropy value for the clustering algorithm */
	const int& Row = ds.dataMatrix.rows();
	float entropy = 0.0, probability;
	for (int i = 0; i < container.size(); ++i) {
		probability = float(container[i]) / float(Row);
		entropy += probability * log2f(probability);
	}
	entropy = -entropy / log2f(numClusters);

	IOHandler::printClustersNoise(ds.dataVec, item_cids, container, "norm" + to_string(normOption),
			ds.fullName, ds.dimension);

	struct timeval start, end;
	double timeTemp;

	numClusters -= 1;

	const int& numNoise = storage[0].size();
	storage.erase(storage.begin());

	/* record labeling information */
	// IOHandler::generateGroups(storage);

	/* compute the centroid coordinates of each clustered group */

	gettimeofday(&start, NULL);

	// compute the centroid coordinates for the clusters
	Eigen::MatrixXf centroid = MatrixXf::Zero(numClusters, ds.dataMatrix.cols());
	vector<vector<float> > cenVec(numClusters);
#pragma omp parallel for schedule(static) num_threads(8)
	for (int i = 0; i < numClusters; ++i)
	{
		const std::vector<int>& groupRow = storage[i];
		for (int j = 0; j < groupRow.size(); ++j)
		{
			centroid.row(i) += ds.dataMatrix.row(groupRow[j]);
		}
		centroid.row(i) /= groupRow.size();
		const Eigen::VectorXf& vec = centroid.row(i);
		cenVec[i] = vector<float>(vec.data(), vec.data() + ds.dataMatrix.cols());
	}

	// extract the streamlines closest and furthest to the centroids for each cluster
	vector<vector<float> > closest(numClusters);
	vector<vector<float> > furthest(numClusters);

#pragma omp parallel for schedule(static) num_threads(8)
	for (int i = 0; i < numClusters; ++i) {
		float minDist = FLT_MAX;
		float maxDist = -10;
		int minIndex = -1, maxIndex = -1;
		const std::vector<int>& groupRow = storage[i];
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

	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec
			- start.tv_usec) / 1.e6;
	activityList.push_back("Feature extraction takes: ");
	timeList.push_back(to_string(timeTemp) + " s");

	// calculate the normalized validity measurement
	ValidityMeasurement vm;
	vm.computeValue(normOption, ds.dataMatrix, item_cids, object, isPBF);
	activityList.push_back("Validity measure is: ");
	stringstream fc_ss;
	fc_ss << vm.f_c;
	timeList.push_back(fc_ss.str());

	// calculate the silhouette, db index and gamma statistics for the evaluation
	gettimeofday(&start, NULL);
	Silhouette sil;
	sil.computeValue(normOption, ds.dataMatrix, ds.dataMatrix.rows(), ds.dataMatrix.cols(),
			item_cids, object, numClusters, isPBF);
	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec
			- start.tv_usec) / 1.e6;
	activityList.push_back("Silhouette calculation takes: ");
	timeList.push_back(to_string(timeTemp) + " s");

	// print the cluster representatives
	std::cout << "Finishing extracting features!" << std::endl;
	IOHandler::printFeature("norm" + to_string(normOption) + "_closest.vtk",
			closest, sil.sCluster, ds.dimension);
	IOHandler::printFeature("norm" + to_string(normOption) + "_furthest.vtk",
			furthest, sil.sCluster, ds.dimension);
	IOHandler::printFeature("norm" + to_string(normOption) + "_centroid.vtk",
			cenVec, sil.sCluster, ds.dimension);

	IOHandler::printToFull(ds.dataVec, sil.sData,
			"norm" + to_string(normOption) + "_SValueLine", ds.fullName,
			ds.dimension);
	IOHandler::printToFull(ds.dataVec, item_cids, sil.sCluster,
			"norm" + to_string(normOption) + "_SValueCluster", ds.fullName,
			ds.dimension);

	// record some time for readme
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

	IOHandler::generateReadme(activityList, timeList);

	IOHandler::writeReadme(entropy, sil, "For norm "+to_string(normOption));

	/* measure closest and furthest rotation */
	std::vector<float> closestRot, furthestRot;
	const float& closestAverage = getRotation(closest, closestRot);
	const float& furthestAverage = getRotation(furthest, furthestRot);

	IOHandler::writeReadme(closestAverage, furthestAverage);

}
