/*
 * This is the source file for the implementation of member functions for AffinityPropagation class
 */

#include "AffinityPropagation.h"


/*
 * @brief This is a default constructor for AffinityPropagation class object
 * @param No parameter input
 * @return No return type
 */

AffinityPropagation::AffinityPropagation()
{

}


/*
 * @brief Create an AffinityPropagation object with command line parameters
 * @param argc: the count of argument
 * @param argv: the argument line with string type of data set names and dimension count
 * @param p: a Para object with some pre-defined parameter values
 * @param automatic: a bool object for automatic parameter setting or not
 * @return Nothing but to create an AffinityPropagation object
 */
AffinityPropagation::AffinityPropagation(const int& argc, char **argv, const Para& p, bool& automatic)
{
	// set the data set information from the provided data set string name
	setDataset(argc, argv);

	if(automatic)	// automate the parameter setting
		setParameterAutomatic(p);

	else	// manually input the parameter
		getParameterUserInput();

	/* select how to initialize the matrixS elements with preference value */
	std::cout << "Please select a MatrixS initialization? 1.median value, 2.minimal value (recommended!)." << std::endl;
	std::cin >> initialOption;
	assert(initialOption==1||initialOption==2);
}


/*
 * @brief A destructor of the class
 * @param No parameter
 * @return No return values
 */
AffinityPropagation::~AffinityPropagation()
{
	// clear the cache memory for distance matrix
	deleteDistanceMatrix(ds.dataMatrix.rows());
}


/*
 * @brief The member function to perform the affinity propagation clustering on similarity measures
 * @param No parameter provided
 * @return Nothing for return
 */
void AffinityPropagation::performClustering()
{
	//distance metric type
	/*  0: Euclidean Norm, d(a,b) = (\sum_(a-b)^2)^(1/2).
		1: Fraction Distance Metric, d(a,b) = (\sum_(a-b)^p)^(1/p), we choose p==0.5
		2: piece-wise angle average, from http://www2.cs.uh.edu/~chengu/Publications/3DFlowVis/curveClustering.pdf
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

	for(int i=0;i<=17;++i)
	{
		if(isPathlines)	// for pathlines, it will call similarity measure d_T (17)
		{
			if(i!=0 && i!=1 && i!=2 && i!=4 && i!=12 && i!=13 && i!=14 && i!=15 && i!=17)
				continue;
		}
		else	// for streamlines, d_T (17) will not be involved
		{
			/* don't want to deal with many too naive metrics */
			if(i!=0 && i!=1 && i!=2 && i!=4 && i!=12 && i!=13 && i!=14 && i!=15)
				continue;
		}

		std::cout << "----------------------------------------------------" << std::endl;
		std::cout << "Experiment on norm " << i << " starts!--------------" << std::endl;

		// clear out the recorded string information
		activityList.clear();
		timeList.clear();

		// perform clustering on the selected similarity measure i
		clusterByNorm(i);

		std::cout << std::endl;
	}
}


/*
 * @brief It performs the AP clustering algorithm on selected similarity measure norm (number tag)
 * @param norm: The number representing the similarity measure
 * @return Nothing to return. Generate some necessary evaluation metrics and data files
 */
void AffinityPropagation::clusterByNorm(const int& norm)
{
	// The parameters to record time needed for calculation
	struct timeval start, end;
	double timeTemp;

	// calculate the distance matrix given the similarity measure type
	getDistanceMatrixFromFile(norm);

	Eigen::MatrixXf matrixR, matrixA, matrixS;

	gettimeofday(&start, NULL);

	/*-------------------------First-level Affinity Propagation----------------------------*/

	// perform the AP clustering based on given distance matrix and matrix S, R and A
	performAPClustering(matrixS, matrixR, matrixA, distanceMatrix, ds.dataMatrix);

	// calculate and record the time for first-level AP clustering
	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec  - start.tv_sec) * 1000000u+ end.tv_usec - start.tv_usec) / 1.e6;

	activityList.push_back("First-level affinity propagation takes: ");
	timeList.push_back(to_string(timeTemp)+" s");

	// some parameters for two-level AP clustering algorithm
	std::vector<std::vector<int> > neighborVec;
	std::vector<int> storage;
	Eigen::MatrixXf centroid;

	// get exemplary examples from the first-level AP
	getGroupAssignment(matrixR, matrixA, matrixS, neighborVec, storage, group);

	// set the labels of initial samples by first-level AP
	setLabel(neighborVec, storage, centroid, group);

	activityList.push_back("First-level affinity propagation generates: ");
	timeList.push_back(to_string(storage.size())+" groups");

	if(useTwoStage)	// two-staged AP is activated
	{

	/*----------------------Second-level Affinity Propagation ---------------------------------------
	 * Use the centroid of the first level and then apply affinipty propagation once again ----------
	 */
		gettimeofday(&start, NULL);

		/* get distance matrix for the centroids */
		float** centroidDistMatrix = NULL;
		getDistMatrixForCentroids(&centroidDistMatrix, normOption, centroid);

		// perform second-level Affinity Propagation on centroids of the streamlines/pathlines
		performAPClustering(matrixS, matrixR, matrixA, centroidDistMatrix, centroid);

		// release the memory of centroidDistMatrix
	#pragma omp parallel for schedule(static) num_threads(8)
		for(int i=0; i<centroid.rows(); ++i)
		{
			delete[] centroidDistMatrix[i];
			centroidDistMatrix[i] = NULL;
		}
		delete[] centroidDistMatrix;
		centroidDistMatrix = NULL;

		// record the time into the README
		gettimeofday(&end, NULL);
		timeTemp = ((end.tv_sec  - start.tv_sec) * 1000000u+ end.tv_usec - start.tv_usec) / 1.e6;
		activityList.push_back("Second-level affinity propagation takes: ");
		timeList.push_back(to_string(timeTemp)+" s");

		/* extract the group information */
		std::vector<std::vector<int> > secondNeighborVec;
		std::vector<int> secondStorage;
		Eigen::MatrixXf secondCentroid;

		std::vector<int> centroidGroup(centroid.rows());

		/* get exemplary examples */
		getGroupAssignment(matrixR, matrixA, matrixS, secondNeighborVec, secondStorage, centroidGroup);

		// get the label of each candidate lines by two-level AP clustering
		setLabel(secondNeighborVec, secondStorage, secondCentroid, centroidGroup);

		secondNeighborVec.clear();

		// record the consumed time
		activityList.push_back("Second-level affinity propagation generates: ");
		timeList.push_back(to_string(secondStorage.size())+" groups");

	/*------------------------ Get the true group id by hierarchical affinity propagation -----------------*/
		// should re-calculate the centroid, storage and neighborVec for new clusters
		getHierarchicalClusters(storage, neighborVec, centroid, group, centroidGroup, secondStorage.size());
	}

	// begin to calculate the evaluation metrics and cluster representatives
	extractFeatures(storage, neighborVec, centroid);

}


/*
 * @brief Set the labels of clustering technique for all the individual lines
 * @param neighborVec: The assemble of lines that belong to the same cluster
 * @param storage: The number of assemble candidates for each cluster
 * @param centroid: The center coordinates of each cluster
 * @param groupTag: The cluster label for each candididate individual line
 */
void AffinityPropagation::setLabel(vector<vector<int> >& neighborVec, vector<int>& storage, Eigen::MatrixXf& centroid,
		std::vector<int>& groupTag)
{
	// record the pair {cluster size, cluster candidate index}
	std::vector<Ensemble> nodeVec;

	for(int i=0;i<storage.size();++i)
	{
		if(storage[i]==0)
			continue;
		nodeVec.push_back({storage[i], neighborVec[i]});
	}

	numberOfClusters = nodeVec.size();

	std::cout << "Cluster label setting begins with " << nodeVec.size() << " clusters..." << std::endl;

	/* sort group index by size of elements containd inside to make sure that, 0 cluster has the
	 * smallest size of candidates
	 */
	std::sort(nodeVec.begin(), nodeVec.end(), [](const Ensemble& first, const Ensemble& second)
	{return first.size<second.size|| (first.size==second.size&&first.element[0]<second.element[0]);});

	// re-define the neighborVec, storage and centroid coordinates given the new cluster index
	neighborVec = std::vector<std::vector<int> >(nodeVec.size());
	storage = std::vector<int>(nodeVec.size());
	centroid = Eigen::MatrixXf(nodeVec.size(), ds.dataMatrix.cols());

	// re-calculate the coordinates of the cluster centroids
#pragma omp parallel for schedule(static) num_threads(8)
	for(int i=0;i<nodeVec.size();++i)
	{
		neighborVec[i] = nodeVec[i].element;
		storage[i] = nodeVec[i].size;
		Eigen::VectorXf tempVec = Eigen::VectorXf::Zero(ds.dataMatrix.cols());
		for(int j=0;j<storage[i];++j)
		{
			tempVec+=ds.dataMatrix.row(i).transpose();
			/* don't forget to re-compute the group tag */
			groupTag[neighborVec[i][j]]=i;
		}
		centroid.row(i) = tempVec/storage[i];
	}

	std::cout << "Cluster label setting ends..." << std::endl;
}


/*
 * @brief Extract the cluster representatives of clusters and calculate the evaluation metrics for the clustering results
 * @param storage: The number of candidates included in each cluster
 * @param neighborVec: The candidates (only streamline index included) for each streamline cluster
 * @param centroid: The centroids of streamline clusters
 * @return Nothing. It will output the cluster representatives and evaluation metric values
 */
void AffinityPropagation::extractFeatures(const std::vector<int>& storage, const std::vector<std::vector<int> >& neighborVec,
		                  const Eigen::MatrixXf& centroid)
{
	const int& Row = ds.dataMatrix.rows();
	const int& Column = ds.dataMatrix.cols();

	/* record labeling information */
	// IOHandler::generateGroups(neighborVec);

	// Output the number of candidates inside each streamline cluster
	std::cout << "Final group number information: " << std::endl;
	for (int i = 0; i < storage.size(); ++i)
	{
		std::cout << storage[i] << " ";
	}
	std::cout << std::endl;

	// calculate the normalized entropy to check the balance of cluster size
	float EntropyRatio;
	getEntropyRatio(storage, EntropyRatio);

	// print the cluster labels in the primary .vtk file
	IOHandler::printClusters(ds.dataVec,group,storage,"AP_norm"+to_string(normOption),ds.fullName,ds.dimension);

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

	// convert the centroid matrix into vector<vector<float>> type. It is not necessary actually
	std::vector<std::vector<float> > center_vec(numberOfClusters, vector<float>(Column));
#pragma omp parallel for schedule(static) num_threads(8)
	for (int i = 0; i < center_vec.size(); ++i)
	{
		for (int j = 0; j < Column; ++j)
		{
			center_vec[i][j] = centroid(i,j);
		}
	}

	// Record the time for extracting the cluster representative lines
	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec  - start.tv_sec) * 1000000u 
			   + end.tv_usec - start.tv_usec) / 1.e6;
	activityList.push_back("Feature extraction takes: ");
	timeList.push_back(to_string(timeTemp)+" s");

	// calculate the normalized validity measurement metric for clustering evaluation
	ValidityMeasurement vm;
	vm.computeValue(normOption, ds.dataMatrix, group, object, isPBF);
	activityList.push_back("Validity measure is: ");
	stringstream fc_ss;
	fc_ss << vm.f_c;
	timeList.push_back(fc_ss.str());

	std::cout << "Finishing extracting features!" << std::endl;	

	// calculate silhouette, the Gamma statistics and DB index for clustering evaluation
	gettimeofday(&start, NULL);
	Silhouette sil;
	sil.computeValue(normOption,ds.dataMatrix,ds.dataMatrix.rows(),ds.dataMatrix.cols(),group,object,
					 numberOfClusters, isPBF, neighborVec);
	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec  - start.tv_sec) * 1000000u
			   + end.tv_usec - start.tv_usec) / 1.e6;
	activityList.push_back("Silhouette calculation takes: ");
	timeList.push_back(to_string(timeTemp)+" s");

	stringstream ss;
	ss << "norm_" << normOption;

	/* measure closest and furthest rotation */
	std::vector<float> closestRotation, furthestRotation;
	const float& closestAverage = getRotation(closest, closestRotation);
	const float& furthestAverage = getRotation(furthest, furthestRotation);

	/* save closest, furthest and centroid representative streamlines */
	IOHandler::printFeature(ds.dataName+"_AP_closest_"+ss.str()+".vtk", closest, sil.sCluster,
			closestRotation, ds.dimension);
	IOHandler::printFeature(ds.dataName+"_AP_furthest_"+ss.str()+".vtk", furthest, sil.sCluster,
			furthestRotation, ds.dimension);
	IOHandler::printFeature(ds.dataName+"_AP_centroid_"+ss.str()+".vtk", center_vec, sil.sCluster,ds.dimension);

	IOHandler::printToFull(ds.dataVec, sil.sData, "AP_SValueLine_"+ss.str(), ds.fullName, ds.dimension);
	IOHandler::printToFull(ds.dataVec, group, sil.sCluster, "AP_SValueCluster_"+ss.str(), ds.fullName, ds.dimension);

	// record the clustering evaluation metric values in the txt file
	activityList.push_back("numCluster is: ");
	timeList.push_back(to_string(numberOfClusters));

	activityList.push_back("Norm option is: ");
	timeList.push_back(to_string(normOption));

	IOHandler::generateReadme(activityList,timeList);

	/* print entropy value for the clustering algorithm */
	IOHandler::writeReadme(EntropyRatio, sil, "For norm "+to_string(normOption));

	IOHandler::writeReadme(closestAverage, furthestAverage);
}


/*
 * @brief A member function to read the geometric coordinates from given file name and position
 * @param argc: The count of arguments
 * @param argv: The argument char array type that includes data set name and dimension count
 * @return Nothing. Just assign the coordinate matrix to the member variables of the class
 */
void AffinityPropagation::setDataset(const int& argc, char **argv)
{
	// the argc should be 3, e.g., ./ap cylinder 3
	if(argc!=3)
	{
		std::cout << "Input argument should have 3!" << endl
		          << "./cluster inputFile_name(in dataset folder) "
		          << "data_dimension(3)" << endl;
		exit(1);
	}

	// extract the required information from argument string
	ds.strName = string("../dataset/")+string(argv[1]);
	ds.dataName = string(argv[1]);
	ds.dimension = atoi(argv[2]);

	/* get the bool tag for variable isPBF */
	std::cout << "It is a PBF dataset? 1.Yes, 0.No" << std::endl;
	int PBFjudgement;
	std::cin >> PBFjudgement;
	assert(PBFjudgement==1||PBFjudgement==0);
	isPBF = (PBFjudgement==1);

	/* check whether it is a Pathline data set or not */
	std::cout << "It is a Pathline? 1.Yes, 0. No" << std::endl;
	std::cin >> PBFjudgement;
	assert(PBFjudgement==1||PBFjudgement==0);
	isPathlines = (PBFjudgement==1);

	// read from the file into the member variables
	IOHandler::readFile(ds.strName,ds.dataVec,ds.vertexCount,ds.dimension,ds.maxElements);

	// print the streamline/pathline vtk file
	ds.fullName = ds.strName+"_full.vtk";
	IOHandler::printVTK(ds.fullName, ds.dataVec, ds.vertexCount, ds.dimension);
}


/*
 * @brief This function is to calculate the normalized entropy for the clustering result
 * @param storage: The size of each streamline cluster
 * @pagram EntropyRatio: The normalized entropy value to be assigned from calculation
 * @return Nothing, just perform the clustering evaluation calculation
 */
void AffinityPropagation::getEntropyRatio(const std::vector<int>& storage, float& EntropyRatio)
{
	// the formula is -s[i]/S * log(s[i]/S), and then normalized by log(numOfClusters)
	EntropyRatio = 0;
	const int& Row = ds.dataMatrix.rows();
	for (int i = 0; i < storage.size(); ++i)
	{
		float ratio = float(storage[i])/float(Row);
		EntropyRatio-=ratio*log2f(ratio);
	}
	/* the higher value shows that the final clusters are balanced and almost equal sized, while the
		low value shows the contrary
	*/
	EntropyRatio/=log2f(storage.size());
}


/*
 * @brief It performs some necessary operations given the input parameters
 * @param p: A Para object that contains parameters for pre-processing and activation for two-staged AP
 * @return Nothing, just perform the sampling and assignment of member variables
 */
void AffinityPropagation::setParameterAutomatic(const Para& p)
{
	// if the data set is pathline, will direct expand the array on the back
	if(isPathlines)
		IOHandler::expandArray(ds.dataMatrix,ds.dataVec,ds.dimension,ds.maxElements);
	else	// it is streamline
	{
		if(p.sampled==1)	// sampling is to directly expand the array from the back
			IOHandler::expandArray(ds.dataMatrix,ds.dataVec,ds.dimension,ds.maxElements);
		else if(p.sampled==2)	// sample the array on the intervals without change of geometric shape
			IOHandler::sampleArray(ds.dataMatrix,ds.dataVec,ds.dimension,ds.maxElements);
		else if(p.sampled==3)	// sample the array with equal arcs such that
			IOHandler::uniformArcSampling(ds.dataMatrix,ds.dataVec,ds.dimension,ds.maxElements);
	}

	// ceate a label vector for each candidate line
	group = std::vector<int>(ds.dataMatrix.rows());

	// assign the parameters for AP clustering
	extractOption = p.extractOption;
	maxIteration = p.maxIteration;

	/* whether to activate two-staged AP or not, see Jun Tao FlowString TVCG 2016 paper for details */
	std::cout << "Whether to activate two-staged AP or not? 1.Yes, 2.No," << std::endl;
	int twoStageOption;
	std::cin >> twoStageOption;
	assert(twoStageOption==1 || twoStageOption==2);
	useTwoStage = (twoStageOption==1);
}


/*
 * @brief It provides console for user parameter input
 * @param No parameters for the function
 * @return Nothing, just set up parameters and provide relevant computation
 */
void AffinityPropagation::getParameterUserInput()
{
	// User input for streamline/pathline sampleOption
	int sampleOption;
	std::cout << "choose a sampling method for the dataset?" << std::endl
			  << "1.directly filling with last vertex; 2. uniform sampling." << std::endl;
	std::cin >> sampleOption;
	assert(sampleOption==1||sampleOption==2);

	if(isPathlines)	// if is pathlines, directly repeat the last vertex of pathlines
		IOHandler::expandArray(ds.dataMatrix,ds.dataVec,ds.dimension,ds.maxElements);
	else	// for streamlines, there are multiple options for that
	{
		if(sampleOption==1)	// direct repeat the last vertex
			IOHandler::expandArray(ds.dataMatrix,ds.dataVec,ds.dimension,ds.maxElements);
		else if(sampleOption==2)	// sample the array on the intervals
			IOHandler::sampleArray(ds.dataMatrix,ds.dataVec,ds.dimension,ds.maxElements);
		else if(sampleOption==3)	// sample the array with equal arc
			IOHandler::uniformArcSampling(ds.dataMatrix,ds.dataVec,ds.dimension,ds.maxElements);
	}

	group = std::vector<int>(ds.dataMatrix.rows());

	// select cluster represnetative strategy, and 1 is recommended
	std::cout << "Select extraction method: 1.centroid, closest and furthest (recommended!), 2.median."
			<< std::endl;
	std::cin >> extractOption;
	assert(extractOption==1||extractOption==2);

	// Input the maximal iteration for AP clustering algorithm
	std::cout << "Input max iteration for affinity propagation: " << std::endl;
	std::cin >> maxIteration;
	assert(maxIteration>0);
}


/*
 * @brief Calculate the matrix S given a distance matrix and the coordinates of lines
 * @param matrixS: The matrix S to be assigned values
 * @param distMatrix: The distance matrix
 * @param coordinates: The coordinate matrix for the streamlines/pathlines
 * @return Nothing, just initialize the matrix S for further AP clustering operation
 */
void AffinityPropagation::getMatrixS(Eigen::MatrixXf& matrixS, float** distMatrix, const Eigen::MatrixXf& coordinates)
{
	std::cout << "Start initializing matrix S..." << std::endl;

	const int& rows = matrixS.rows();

	/* define a vector to store pair-wise distance vector and get the median */
	const int& distVecSize = rows*(rows-1)/2;
	std::vector<float> distVec(distVecSize);
	int count = 0;

	/* find the minimal dissimilarity value from the distance matrix */
	float minV = (float)FLT_MAX;
	float tempDist;
	for(int i=0;i<rows-1;++i)
	{
		for(int j=i+1;j<rows;++j)
		{
			if(distMatrix)	// if distance matrix exists, direct fetch the cached value
				tempDist = distMatrix[i][j];
			else	// otherwise, has to calculate the distance matrix
				tempDist = getDisimilarity(coordinates, i, j, normOption, object);

			/* conventionally we assign -d*d as non-diagonal entries for matrix S */
			matrixS(i,j) = -tempDist;
			matrixS(j,i) = matrixS(i,j);

			minV = std::min(minV, matrixS(i,j));
			distVec[count++] = matrixS(i,j);
		}
	}

	std::cout << "min Value is " << minV << std::endl;
	assert(count==distVecSize);

	float initialValue;
	if(initialOption==1)	// the initialization is by median of distance matrix values
	{
		/* get median value to be assigned for S(i,i) */
		float medianValue, leftMedian, rightMedian;

		/* odd size, just pick mid index */
		if(distVecSize%2==1)
			medianValue = select(distVec, 0, distVecSize-1, distVecSize/2);
		/* even size, choose average of left and right */
		else if(distVecSize%2==0)
		{
			leftMedian = select(distVec, 0, distVecSize-1, (distVecSize-1)/2);
			rightMedian = select(distVec, 0, distVecSize-1, distVecSize/2);
			medianValue = (leftMedian+rightMedian)/2.0;
		}
		// assign the preference value as median of the distance matrix values
		initialValue = medianValue;
	}
	else if(initialOption==2)	// the initialization is by minimal dissimilarity value
	{
		initialValue = minV;
	}
	std::cout << "Initial value is " << initialValue << std::endl;

	/* assign the initialValue to diagonal matrix element */
#pragma omp parallel for schedule(static) num_threads(8)
	for(int i=0;i<rows;++i)
		matrixS(i,i) = initialValue;

	std::cout << "Finish initializing matrix S..." << std::endl;
}


/*
 * @brief Initialize the matrix S, R and A for AP clustering
 * @param matrixS: The matrix S to be initialized
 * @param matrixR: The matrix R to be initialized
 * @param matrixA: The matrix A to be initialized
 * @param rows: The row size of the matrix
 * @return Nothing, just resize the memory for three matrices
 */
void AffinityPropagation::initializeMatrices(Eigen::MatrixXf& matrixS, Eigen::MatrixXf& matrixR,
											Eigen::MatrixXf& matrixA, const int& rows)
{
	/* initialize all three matrices as zero entry */
	matrixS = Eigen::MatrixXf::Zero(rows, rows);
	matrixR = Eigen::MatrixXf::Zero(rows, rows);
	matrixA = Eigen::MatrixXf::Zero(rows, rows);
}


/*
 * @brief Update the responsibility matrix R by input matrix A and S
 * @param matrixR: The matrix R to be updated
 * @param matrixA: The matrix A as input
 * @param matrixR: The matrix S as input
 * @return Nothing, just update the entry-wise value of matrix R
 */
void AffinityPropagation::updateResponsibility(Eigen::MatrixXf& matrixR, const Eigen::MatrixXf& matrixA,
											  const Eigen::MatrixXf& matrixS)
{
	const int& rows = matrixR.rows();
	// update the R with relaxed value of S and R
#pragma omp parallel for schedule(static) num_threads(8)
	for(int i=0;i<rows;++i)
	{
		for(int k=0;k<rows;++k)
		{
			/* don't use FLT_MIN because FLT_MIN == 0.0 */
			float maxValue = -FLT_MAX;
			for(int kk=0;kk<rows;++kk)
			{
				if(kk==k)
					continue;
				maxValue = std::max(maxValue, matrixS(i,kk)+matrixA(i,kk));
			}
			/* in wikipage it's update by R[i,k] = S[i][k]-maxValue, but here use a Laplace smoothor for convergence */
			matrixR(i,k) = (1-LAMBDA)*(matrixS(i,k)-maxValue)+LAMBDA*matrixR(i,k);
		}
	}
}


/*
 * @brief Update the availability matrix A
 * @param matrixA: the availability matrix A to be updated
 * @param matrixR: the responsibility matrix R as input
 * @return Nothing but just update the matrix A
 */
void AffinityPropagation::updateAvailability(Eigen::MatrixXf& matrixA, const Eigen::MatrixXf& matrixR)
{
	const int& rows = matrixR.rows();
#pragma omp parallel for schedule(static) num_threads(8)
	for(int i=0;i<rows;++i)
	{
		for(int k=0;k<rows;++k)
		{
			/* for diagonal matrix, update by summation of non-diagonal entries in the row */
			if(i==k)
			{
				float summation = 0.0;
				for(int ii=0;ii<rows;++ii)
				{
					if(ii==i)
						continue;
					summation+=std::max((float)0.0, matrixR(ii,k));
				}

				/* smoothing update instead of direct assignment */
				matrixA(i,k)=(1-LAMBDA)*summation+LAMBDA*matrixA(i,k);
			}
			else
			{
				float summation = 0.0;
				for(int ii=0;ii<rows;++ii)
				{
					if(ii==i||ii==k)
						continue;
					summation+=std::max((float)0.0, matrixR(ii,k));
				}
				matrixA(i,k)=(1-LAMBDA)*std::min((float)0.0, matrixR(k,k)+summation)+LAMBDA*matrixA(i,k);
			}
		}
	}
}


/*
 * @brief Get the group label for all the streamlines/pathlines
 * @param matrixR: The responsibility matrix R as input
 * @param matrixA: The availability matrix A as input
 * @param matrixS: The matrix S as input
 * @param neighborVec: The neighborhood vector for each cluster
 * @param storage: The int vector for recording the size of each cluster
 * @param groupTag: The labels for each streamlines/pathlines
 * @return Nothing, will update the cluster labels for all the streamlines/pathlines
 */
void AffinityPropagation::getGroupAssignment(const Eigen::MatrixXf& matrixR, const Eigen::MatrixXf& matrixA,
			  	  	  	  	const Eigen::MatrixXf& matrixS, std::vector<std::vector<int> >& neighborVec,
							std::vector<int>& storage, std::vector<int>& groupTag)
{
	std::vector<int> centerVec;
	const int& rows = matrixR.rows();

	/* store the candidate whose diagonal summation is positive */
	float diagonalSum;
	for(int i=0;i<rows;++i)
	{
		diagonalSum=matrixR(i,i)+matrixA(i,i);
		if(diagonalSum>0)
		{
			centerVec.push_back(i);
		}
	}

	const int& centerSize = centerVec.size();
	/* get group tag information for each candidate streamline */
#pragma omp parallel for schedule(static) num_threads(8)
	for(int i=0;i<rows;++i)
	{
		int index, element;
		float maxSim = -FLT_MAX;
		for(int j=0;j<centerSize;++j)
		{
			element = centerVec[j];
			if(matrixS(i,element)>maxSim)
			{
				maxSim = matrixS(i,element);
				index = element;
			}
		}
		groupTag[i]=index;
	}

	/* output group information and cluster size */
	std::map<int,int> groupMap;
	for(int i=0;i<rows;++i)
	{
		/* group tag not int the hash map */
		if(groupMap.find(groupTag[i])==groupMap.end())
		{
			groupMap.insert(make_pair(groupTag[i],0));
		}
	}

	/* give them new index starting from 0 */
	int count = 0;
	for(auto iter = groupMap.begin();iter!=groupMap.end();++iter)
	{
		iter->second = count++;
	}

	numberOfClusters = groupMap.size();

	/* assign contained element and size */
	neighborVec = std::vector<std::vector<int> >(numberOfClusters);
	storage = std::vector<int>(numberOfClusters);
	for(int i=0;i<rows;++i)
	{
		count = groupMap[group[i]];
		neighborVec[count].push_back(i);
	}

	/* assign the storage vector */
	for(int i=0;i<storage.size();++i)
	{
		storage[i] = neighborVec[i].size();
	}
}


/*
 * @brief The function is to get the distance matrix for centroid streamlines/pathlines
 * @param centroidDistMatrix: The distance matrix for centroid streamlines/pathlines
 * @param norm: The input similarity label
 * @param centroid: The coordinate matrix for the centroids
 * @return Nothing, update the distance matrix
 */
void AffinityPropagation::getDistMatrixForCentroids(float*** centroidDistMatrix, const int& norm,
		const Eigen::MatrixXf& centroid)
{
	const int& rows = centroid.rows();
	*centroidDistMatrix = new float*[rows];

	/* in order to calculate the distance matrix given norm, we need to calculate the object first. This object
	 * is to pre-calculate some preliminary stuff for distance matrix computation. I know it is redundant but in
	 * practice it can help to accelerate the performance a little bit
	 */

	MetricPreparation centroidObj = MetricPreparation(centroid.rows(), centroid.cols());
	centroidObj.preprocessing(centroid, centroid.rows(), centroid.cols(), norm);

	// calculate the distance matrix among centroid matrix coordinates
#pragma omp parallel for schedule(static) num_threads(8)
	for(int i=0; i<rows; ++i)
	{
		(*centroidDistMatrix)[i] = new float[rows];
		for(int j=0; j<rows; ++j)
		{
			(*centroidDistMatrix)[i][j] = getDisimilarity(centroid, i, j, norm, centroidObj);
		}
	}
}


/*
 * @brief This function is to read distance matrix from the local file if it exists
 * @param norm: It is the index for similarity measure
 * @return Nothing, just update the global variable distanceMatrix
 */
void AffinityPropagation::getDistanceMatrixFromFile(const int& norm)
{
	normOption = norm;

	/* very hard to decide whether needed to perform such pre-processing, but recommended
	 * to create a cached object for further pair-wise distance matrix calculation
	 */
	object = MetricPreparation(ds.dataMatrix.rows(), ds.dataMatrix.cols());
	object.preprocessing(ds.dataMatrix, ds.dataMatrix.rows(), ds.dataMatrix.cols(), normOption);

	/* would store distance matrix instead because it would save massive time */
	struct timeval start, end;
	double timeTemp;
	gettimeofday(&start, NULL);

	// in case the distance matrix already exists for other similarity, will clean it first
	deleteDistanceMatrix(ds.dataMatrix.rows());

	// read distance matrix from the local file in ../dataset/
	std::ifstream distFile(("../dataset/"+to_string(normOption)).c_str(), ios::in);

	// the local file of distance matrix does not exist, then will create the file
	if(distFile.fail())
	{
		distFile.close();
		// calculate the distance matrix from norm option
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
	else // the local file for distance matrix computation exists, then directly read in
	{
		std::cout << "read distance matrix..." << std::endl;

		// create the distance matrix and read in the content
		distanceMatrix = new float*[ds.dataMatrix.rows()];
	#pragma omp parallel for schedule(static) num_threads(8)
		for (int i = 0; i < ds.dataMatrix.rows(); ++i)
		{
			distanceMatrix[i] = new float[ds.dataMatrix.rows()];
		}
		int i=0, j;
		string line;
		stringstream ss;
		// extract the distance values from the file
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

	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec  - start.tv_sec) * 1000000u
			   + end.tv_usec - start.tv_usec) / 1.e6;
	activityList.push_back("Distance matrix computing for norm "+to_string(normOption)+" takes: ");
	timeList.push_back(to_string(timeTemp)+" s");
}


/*
 * @brief This function is to perform the AP clustering
 * @param matrixS: The matrix S to be updated
 * @param matrixR: The matrix R to be updated
 * @param matrixA: The matrix A to be updated
 * @param distMatrix: The distance matrix as input
 * @param coordinates: The coordinate of streamlines/pathlines as input
 * @return Nothing, just perform the AP clustering
 */
void AffinityPropagation::performAPClustering(Eigen::MatrixXf& matrixS, Eigen::MatrixXf& matrixR,
		Eigen::MatrixXf& matrixA, float** distMatrix, const Eigen::MatrixXf& coordinates)
{
	/* initialize S, R, A */
	initializeMatrices(matrixS, matrixR, matrixA, coordinates.rows());

	/* get S */
	getMatrixS(matrixS, distMatrix, coordinates);

	int current = 0;
	while(current++<maxIteration)
	{
		std::cout << "Iteration " << current << std::endl;

		/* update responsibility */
		updateResponsibility(matrixR, matrixA, matrixS);

		/* update availability */
		updateAvailability(matrixA, matrixR);

	}
}


/*
 * @brief This function is to calculate the centroids by AP clustering labels
 * @param storage: The vector of int for the size of each cluster
 * @param neighborVec: The candidate of each cluster
 * @param centroid: The centroid streamline to be updated
 * @param groupTag: The labels for each streamline
 * @param centroidGroup: The number of candidates for each cluster
 * @param groupSize: The size of groups generated
 * @return Nothing, just to update the neighborVec, storage and centroid
 */
void AffinityPropagation::getHierarchicalClusters(std::vector<int>& storage, std::vector<std::vector<int> >& neighborVec,
		Eigen::MatrixXf& centroid, std::vector<int>& groupTag, const std::vector<int>& centroidGroup,
		const int& groupSize)
{
	neighborVec.clear();
	neighborVec.resize(groupSize);
	storage.resize(groupSize);
	centroid = Eigen::MatrixXf::Zero(groupSize, centroid.cols());

	int groupID;
	for(int i=0; i<groupTag.size(); ++i)
	{
		groupID = centroidGroup[groupTag[i]];
		groupTag[i] = groupID;
		neighborVec[groupID].push_back(i);
		centroid.row(groupID)+=ds.dataMatrix.row(i);
	}

#pragma omp parallel for schedule(static) num_threads(8)
	for(int i=0; i<groupSize; ++i)
	{
		centroid.row(i)/=neighborVec[i].size();
		storage[i] = neighborVec[i].size();
	}
}
