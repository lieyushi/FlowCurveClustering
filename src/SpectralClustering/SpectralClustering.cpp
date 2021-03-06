/*
 * @brief The SpectralClustering class to perform the spectral clustering on input data set
 *
 * @detais
 * 	Spectral clustering is a graph-based technique to map original streamlines to a spectral embedding space.
 * The problem itself is a NP-hard and we used instead a relaxed versions with Graph Laplacians.
 *
 *  Detailed procedures can be referred at https://tarekmamdouh.wordpress.com/2014/09/28/spectral-clustering/
 * and TVCG paper http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=6702500
 * local scaling for Gaussian kernel size might be defined as 0.05*totalCount as in
 * Blood Flow Clustering and Applications in Virtual Stenting of Intracranial Aneurysms
 */


#include "SpectralClustering.h"


/*
 * @brief default constructor
 */
SpectralClustering::SpectralClustering()
{

}


/*
 * @brief The argument constructor with argc and argv
 *
 * @param[in] argc The count of arguments
 * @param[in] argv The char* array of arguments
 * @param[in] p The Para object
 * @param[out] automatic The bool flag
 */
SpectralClustering::SpectralClustering(const int& argc, char **argv, const Para& p, bool& automatic)
{
	setDataset(argc, argv);

	if(automatic)
		setParameterAutomatic(p);
	else
		getParameterUserInput();

}


/*
 * @brief destructor
 */
SpectralClustering::~SpectralClustering()
{
	deleteDistanceMatrix(ds.dataMatrix.rows());
}


/*
 * @brief perform spectral clustering function
 */
void SpectralClustering::performClustering()
{
	//distance metric type
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
	if(postProcessing==2)
	{
		std::cout << "Find optimal activated? 0. No, 1. Yes: " << std::endl;
		int optimalOption;
		std::cin >> optimalOption;
		assert(optimalOption==0 || optimalOption==1);
		isOptimal = (optimalOption==1);

		std::cout << "Please input the preset number of clusters in [2, " << ds.dataVec.size() << "]: " << std::endl;
		std::cin >> presetNumber;
		assert(presetNumber>=2 && presetNumber<=ds.dataVec.size());

		/* record initial number of clusters of user input */
		recordPreset(presetNumber);

		readCluster = false;
	}
	else if(postProcessing==1)
	{
		std::cout << "Please choose cluster number method, 0.user input, 1.read clustering: " << std::endl;
		int clusterInput;
		std::cin >> clusterInput;
		assert(clusterInput==0 || clusterInput==1);
		readCluster = (clusterInput==1);
	}

	std::unordered_map<int,int> clusterMap;
	if(readCluster)
	{
		IOHandler::readClusteringNumber(clusterMap, "cluster_number");
	}


	for(int i=0;i<=17;++i)
	{
		if(isPathlines)
		{
			/* don't want to deal with many too naive metrics */
			if(i!=0 && i!=1 && i!=2 && i!=4 && i!=12 && i!=13 && i!=14 && i!=15 && i!=17)
				continue;
		}
		else
		{
			if(i!=0 && i!=1 && i!=2 && i!=4 && i!=12 && i!=13 && i!=14 && i!=15)
				continue;
		}

		if(postProcessing==1)
		{
			if(readCluster)
				presetNumber = clusterMap[i];
			else
			{
				std::cout << "Please input the preset number of clusters for norm " << i << " among [2, "
						<< ds.dataVec.size() << "]: " << std::endl;
				std::cin >> presetNumber;
			}
			assert(presetNumber>=2 && presetNumber<=ds.dataVec.size());
		}

		std::cout << "----------------------------------------------------" << std::endl;
		std::cout << "Experiment on norm " << i << " starts!--------------" << std::endl;

		activityList.clear();
		timeList.clear();

		activityList.push_back("Preset numOfClusters for norm "+to_string(i) +" is: ");
		timeList.push_back(to_string(presetNumber));

		struct timeval start, end;
		double timeTemp;
		gettimeofday(&start, NULL);

		clusterByNorm(i);

		gettimeofday(&end, NULL);
		timeTemp = ((end.tv_sec  - start.tv_sec) * 1000000u
						   + end.tv_usec - start.tv_usec) / 1.e6;
		activityList.push_back("SC for "+to_string(i) +" takes: ");
		timeList.push_back(to_string(timeTemp)+"s");

		std::cout << std::endl;
	}
}


/*
 * @brief run spectral clustering based on different norm input
 *
 * @param[in] norm The norm option
 */
void SpectralClustering::clusterByNorm(const int& norm)
{
	normOption = norm;

	/* very hard to decide whether needed to perform such pre-processing */
	object = MetricPreparation(ds.dataMatrix.rows(), ds.dataMatrix.cols());
	object.preprocessing(ds.dataMatrix, ds.dataMatrix.rows(), ds.dataMatrix.cols(), normOption);

	/* would store distance matrix instead because it would save massive time */
	struct timeval start, end;
	double timeTemp;
	gettimeofday(&start, NULL);

	deleteDistanceMatrix(ds.dataMatrix.rows());

	std::ifstream distFile(("../dataset/"+to_string(normOption)).c_str(), ios::in);
	if(distFile.fail())
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
	else
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
	activityList.push_back("Distance matrix computing for norm "+to_string(norm)+" takes: ");
	timeList.push_back(to_string(timeTemp)+" s");

	getSigmaList();


	Eigen::MatrixXf adjacencyMatrix, laplacianMatrix;
	Eigen::DiagonalMatrix<float,Dynamic> degreeMatrix;

	/* get weighted adjacency matrix by Gaussian kernel */
	getAdjacencyMatrix(adjacencyMatrix);

	/* get degree matrix */
	getDegreeMatrix(adjacencyMatrix, degreeMatrix);

	/* get Laplacian matrix */
	getLaplacianMatrix(adjacencyMatrix, degreeMatrix, laplacianMatrix);

	getEigenClustering(laplacianMatrix, norm);
}


/*
 * @brief perform group-labeling information for all the streamlines
 *
 * @param[in] neighborVec The neighboring vector of candidates belonging to the clusters
 * @param[out] storage The individual size of clusters
 * @param[out] centroid The centroid coordinates of the clusters
 */
void SpectralClustering::setLabel(vector<vector<int> >& neighborVec, vector<int>& storage, Eigen::MatrixXf& centroid)
{
	std::vector<Ensemble> nodeVec(storage.size());
	std::cout << "Cluster label setting begins with " << numberOfClusters << " clusters..." << std::endl;
#pragma omp parallel for schedule(static) num_threads(8)
	for(int i=0;i<nodeVec.size();++i)
	{
		nodeVec[i].size = storage[i];
		nodeVec[i].element = neighborVec[i];
	}

	/* sort group index by size of elements containd inside */
	std::sort(nodeVec.begin(), nodeVec.end(), [](const Ensemble& first, const Ensemble& second)
	{return first.size<second.size|| (first.size==second.size&&first.element[0]<second.element[0]);});

	neighborVec = std::vector<std::vector<int> >(nodeVec.size());
	storage = std::vector<int>(nodeVec.size());
	centroid = Eigen::MatrixXf(nodeVec.size(), ds.dataMatrix.cols());

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
			group[neighborVec[i][j]]=i;
		}
		centroid.row(i) = tempVec/storage[i];
	}

	std::cout << "Cluster label setting ends..." << std::endl;
}


/*
 * @brief extract features from datasets as representative curves and calculate the clustering evaluation
 *
 * @param[in] storage The size of different clusters
 * @param[in] neighborVec The candidates included in each cluster
 * @param[in] centroid The centroid coordinates of cluster
 */
void SpectralClustering::extractFeatures(const std::vector<int>& storage, const std::vector<std::vector<int> >& neighborVec,
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

	string pprocessing;
	switch(postProcessing)
	{
	case 1:
		pprocessing="Kmeans";
		break;

	case 2:
		pprocessing="EigenRot";
		break;
	}

	float EntropyRatio;
	getEntropyRatio(storage, EntropyRatio);

	/* record labeling information */
	// IOHandler::generateGroups(neighborVec);


	IOHandler::printClusters(ds.dataVec,group,storage,"SC_"+pprocessing+"_norm"+to_string(normOption),ds.fullName,ds.dimension);

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

	std::vector<std::vector<float> > center_vec(numberOfClusters, vector<float>(Column));
#pragma omp parallel for schedule(static) num_threads(8)
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

	ValidityMeasurement vm;
	vm.computeValue(normOption, ds.dataMatrix, group, object, false);
	activityList.push_back("SC Validity measure is: ");
	stringstream fc_ss;
	fc_ss << vm.f_c;
	timeList.push_back(fc_ss.str());

	gettimeofday(&start, NULL);
	Silhouette sil;
	sil.computeValue(normOption,ds.dataMatrix,ds.dataMatrix.rows(),ds.dataMatrix.cols(),group,object,
			         numberOfClusters, false, neighborVec);
	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
	activityList.push_back("Silhouette calculation takes: ");
	timeList.push_back(to_string(timeTemp)+" s");

	std::cout << "Finishing extracting features!" << std::endl;	

	stringstream ss;
	ss << "norm_" << normOption;

	std::vector<float> closestRotation, furthestRotation;
	const float& closestAverage = getRotation(closest, closestRotation);
	const float& furthestAverage = getRotation(furthest, furthestRotation);

	/* save closest, furthest and centroid representative streamlines */
	IOHandler::printFeature(ds.dataName+"_SC"+pprocessing+"_closest_"+ss.str()+".vtk", closest, sil.sCluster,
			closestRotation, ds.dimension);
	IOHandler::printFeature(ds.dataName+"_SC"+pprocessing+"_furthest_"+ss.str()+".vtk", furthest, sil.sCluster,
			furthestRotation, ds.dimension);
	IOHandler::printFeature(ds.dataName+"_SC"+pprocessing+"_centroid_"+ss.str()+".vtk", center_vec, sil.sCluster,ds.dimension);

	IOHandler::printToFull(ds.dataVec, sil.sData, "SC"+pprocessing+"_SValueLine_"+ss.str(), ds.fullName, ds.dimension);
	IOHandler::printToFull(ds.dataVec, group, sil.sCluster, "SC"+pprocessing+"_SValueCluster_"+ss.str(), ds.fullName, ds.dimension);

	activityList.push_back("numCluster is: ");
	timeList.push_back(to_string(numberOfClusters));

	activityList.push_back("Norm option is: ");
	timeList.push_back(to_string(normOption));

	activityList.push_back("SC post-processing is: ");
	switch(postProcessing)
	{
	case 1:
		timeList.push_back("k-means");
		break;

	case 2:
		timeList.push_back("vector rotation");
		break;
	}

	activityList.push_back("Average rotation of closest is: ");
	timeList.push_back(to_string(closestAverage));

	activityList.push_back("Average rotation of furthest is: ");
	timeList.push_back(to_string(furthestAverage));

	IOHandler::generateReadme(activityList,timeList);

	IOHandler::writeReadme(EntropyRatio, sil, "For norm "+to_string(normOption));
}


/*
 * @brief set data set from user command
 *
 * @param[in] argc The count of arguments
 * @param[in] argv The char* array of argument string
 */
void SpectralClustering::setDataset(const int& argc, char **argv)
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

	IOHandler::readFile(ds.strName,ds.dataVec,ds.vertexCount,ds.dimension,ds.maxElements);

	ds.fullName = ds.strName+"_full.vtk";
	IOHandler::printVTK(ds.fullName, ds.dataVec, ds.vertexCount, ds.dimension);

	SCALING = 0.05*ds.dataVec.size();
}


/*
 * @brief get local scaling from NIPS 2002 paper
 */
void SpectralClustering::getSigmaList()
{
	const int& Row = ds.dataMatrix.rows();
	sigmaVec = std::vector<float>(Row);

	if(isDistSorted)
	{
		/* get SCALING-th smallest dist */
	#pragma omp parallel for schedule(static) num_threads(8)
		for(int i=0;i<Row;++i)
		{
			/* instead we implement a n*logk priority_queue method for finding k-th smallest element */
			std::priority_queue<float> limitQueue;
			float tempDist;
			for(int j=0;j<Row;++j)
			{
				if(i==j)
					continue;
				if(distanceMatrix)
					tempDist = distanceMatrix[i][j];
				else
					tempDist = getDisimilarity(ds.dataMatrix, i, j, normOption, object);
				// element is even larger than the biggest
				limitQueue.push(tempDist);
				if(limitQueue.size()>SCALING)
					limitQueue.pop();
			}

			sigmaVec[i] = limitQueue.top();
		}
	}
	else
	{
		/* directly by index since in both papers only mention i-th neighboring point */
	#pragma omp parallel for schedule(static) num_threads(8)
		for(int i=0;i<Row;++i)
		{
			if(i<SCALING)
			{
				if(distanceMatrix)
					sigmaVec[i]=distanceMatrix[i][SCALING];
				else
					sigmaVec[i]=getDisimilarity(ds.dataMatrix, i, SCALING, normOption, object);
			}
			else
			{
				if(distanceMatrix)
					sigmaVec[i]=distanceMatrix[i][SCALING-1];
				else
					sigmaVec[i]=getDisimilarity(ds.dataMatrix, i, SCALING-1, normOption, object);
			}
		}
	}
	std::cout << "Finish local scaling..." << std::endl;
}


/*
 * @brief get entropy ratio from size of clusters as input
 *
 * @param[in] storage The size of clusters as input
 * @param[out] EntropyRatio The entropy ratio to calculate
 */
void SpectralClustering::getEntropyRatio(const std::vector<int>& storage, float& EntropyRatio)
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
 * @brief get weighted adjacency matrix by Gaussian kernel
 *
 * @param[out] adjacencyMatrix The weighted adjacency matrix computed from the Gaussian graph algorithm
 */
void SpectralClustering::getAdjacencyMatrix(Eigen::MatrixXf& adjacencyMatrix)
{
	//in case of diagonal matrix element is not assigned
	adjacencyMatrix = Eigen::MatrixXf::Zero(ds.dataMatrix.rows(), ds.dataMatrix.rows());
#pragma omp parallel for schedule(static) num_threads(8)
	for(int i=0;i<adjacencyMatrix.rows();++i)
	{
		for(int j=0;j<adjacencyMatrix.cols();++j)
		{
			float dist_ij;
			if(i==j)
				continue;
			else if(distanceMatrix)
			{
				dist_ij = distanceMatrix[i][j];
			}
			else
				dist_ij = getDisimilarity(ds.dataMatrix, i, j, normOption, object);
			adjacencyMatrix(i,j)=exp(-dist_ij*dist_ij/sigmaVec[i]/sigmaVec[j]);
		}
	}

	std::cout << "Finish computing adjacency matrix!" << std::endl;
}


/*
 * @brief get degree matrix by the adjacency matrix
 *
 * @param[in] adjacencyMatrix The adjacency matrix as input
 * @param[ou] degreeMatrix The degree matrix (diagonal matrix) as output
 */
void SpectralClustering::getDegreeMatrix(const Eigen::MatrixXf& adjacencyMatrix, Eigen::DiagonalMatrix<float,Dynamic>& degreeMatrix)
{
	degreeMatrix = Eigen::DiagonalMatrix<float,Dynamic>(ds.dataMatrix.rows());
	Eigen::VectorXf v = VectorXf::Zero(ds.dataMatrix.rows());
#pragma omp parallel for schedule(static) num_threads(8)
	for(int i=0;i<v.size();++i)
	{
		float summation = 0;
		for(int j=0;j<adjacencyMatrix.cols();++j)
		{
			summation+=adjacencyMatrix(i,j);
		}
		v(i) = summation;
	}

	degreeMatrix.diagonal() = v;

	std::cout << "Fnish computing degree matrix!" << std::endl;
}


/*
 * @brief get Laplacian matrix
 *
 * @param[in] adjacencyMatrix The adjacency matrix as input
 * @param[out] degreeMatrix The degree matrix to be updated
 * @param[out] laplacianMatrix The Laplacian matrix to be calculated
 */
void SpectralClustering::getLaplacianMatrix(const Eigen::MatrixXf& adjacencyMatrix,
		                                    Eigen::DiagonalMatrix<float,Dynamic>& degreeMatrix,
											Eigen::MatrixXf& laplacianMatrix)
{
	switch(LaplacianOption)
	{
	default:
	case 1:
	/* L = D^(-1)A */
		getMatrixPow(degreeMatrix, -1.0);
		laplacianMatrix=degreeMatrix*adjacencyMatrix;
		break;

	case 2:
		Eigen::MatrixXf dMatrix = Eigen::MatrixXf(adjacencyMatrix.rows(),adjacencyMatrix.cols());
		const Eigen::VectorXf& m_v = degreeMatrix.diagonal();
		for(int i=0;i<dMatrix.rows();++i)
			dMatrix(i,i) = m_v(i);
		laplacianMatrix = dMatrix-adjacencyMatrix;
		break;
	}
}


/*
 * @brief decide optimal cluster number by eigenvectors of Laplacian matrix
 *
 * @param[in] laplacianMatrix The input Laplacian matrix for the eigen-rotation minimization
 * @param[in] norm The norm option
 */
void SpectralClustering::getEigenClustering(const Eigen::MatrixXf& laplacianMatrix, const int& norm)
{
	struct timeval start, end;
	gettimeofday(&start, NULL);

	/* eigen decomposition for Hermite matrix (real and symmetric matrix) */
	std::cout << "Eigen decomposition starts!..." << std::endl;
	SelfAdjointEigenSolver<MatrixXf> eigensolver(laplacianMatrix);
	std::cout << "Eigen decomposition ends!..." << std::endl;

	gettimeofday(&end, NULL);
	float timeTemp = ((end.tv_sec-start.tv_sec)*1000000u+end.tv_usec-start.tv_usec)/1.e6;
	activityList.push_back("Eigen decomposition takes: ");
	timeList.push_back(to_string(timeTemp)+" s");

	const int& eigenRows = presetNumber;
	std::cout << "Eigen rows are: " << eigenRows << std::endl;
	//const int& eigenRows = 5;

	Eigen::MatrixXf eigenVec(eigenRows, ds.dataMatrix.rows());

	const int& Row = laplacianMatrix.rows();

	/* from paper we know it should get largest eigenvalues, and from eigen library we know it's latter */
	for(int i=Row-1;i>Row-eigenRows-1;--i)
		eigenVec.row(Row-1-i) = eigensolver.eigenvectors().col(i).transpose();
	eigenVec.transposeInPlace();

	/* how many elements in each cluster */
	std::vector<int> storage;

	/* which elements stored in each cluster */
	std::vector<std::vector<int> > neighborVec;

	/* centroid cluster */
	Eigen::MatrixXf clusterCenter;

	/* k-means as a post-processing */
	if(postProcessing==1)
	{
		normalizeEigenvec(eigenVec);

		performKMeans(eigenVec,storage,neighborVec);

		setLabel(neighborVec, storage, clusterCenter);

		extractFeatures(storage,neighborVec,clusterCenter);
	}
	/* eigenvector rotation */
	else if(postProcessing==2)
	{
		getEigvecRotation(storage,neighborVec,clusterCenter,eigenVec);

		if(neighborVec.empty())
			return;

		setLabel(neighborVec, storage, clusterCenter);

		if(isOptimal)
			recordOptimalResult(norm, neighborVec.size());
		else
			extractFeatures(storage,neighborVec,clusterCenter);
	}
}


/*
 * @brief Calculate the matrix^powNumber, and it would be inv if powNumber is -1
 *
 * @param[out] matrix The original matrix
 * @param[in] powNumber The exponential index
 */
void getMatrixPow(Eigen::DiagonalMatrix<float,Dynamic>& matrix, const float& powNumber)
{
	Eigen::VectorXf& m_v = matrix.diagonal();
#pragma omp parallel for schedule(static) num_threads(8)
	for(int i=0;i<m_v.size();++i)
		m_v(i) = pow(m_v(i), powNumber);
}


/*
 * @brief normalize each row first
 *
 * @param[out] eigenVec The matrix with eigen-vectors to be updated
 */
void SpectralClustering::normalizeEigenvec(Eigen::MatrixXf& eigenVec)
{
	const int& rows = eigenVec.rows();
#pragma omp parallel for schedule(static) num_threads(8)
	for(int i=0;i<rows;++i)
	{
		eigenVec.row(i)/=eigenVec.row(i).norm();
	}
}


/*
 * @brief perform k-means clustering for the normalized eigen vector matrix
 *
 * @param[in] eigenVec The input eigen vector as matrix
 * @param[out] storage The size of clusters to be updated
 * @param[out] neighborVec The vector of candidates belonging to clusters
 */
void SpectralClustering::performKMeans(const Eigen::MatrixXf& eigenVec,
									   std::vector<int>& storage,
									   std::vector<std::vector<int> >& neighborVec)
{

	const int& Row = eigenVec.rows();
	const int& Column = eigenVec.cols();

	float moving=1000, tempMoving, before;

	numberOfClusters = presetNumber;

	storage = std::vector<int>(numberOfClusters);

	/* centerTemp is temporary term for storing centroid position, clusterCenter is permanent */
	MatrixXf centerTemp, clusterCenter;

	/* chosen from sample for initialization of k-means */
	Initialization::generateFromSamples(clusterCenter,Column,eigenVec,numberOfClusters);

	int tag = 0;

	neighborVec=std::vector< std::vector<int> >(numberOfClusters);

	float PCA_KMeans_delta, KMeans_delta;

	std::cout << "...k-means started!" << std::endl;

	struct timeval start, end;
	gettimeofday(&start, NULL);

	do
	{
		before = moving;
		/* preset cluster number recorder */
		std::fill(storage.begin(), storage.end(), 0);

		centerTemp = MatrixXf::Zero(numberOfClusters, Column);

	#pragma omp parallel for schedule(static) num_threads(8)
		for (int i = 0; i < numberOfClusters; ++i)
		{
			neighborVec[i].clear();
		}

	#pragma omp parallel num_threads(8)
		{
		#pragma omp for nowait
			for (int i = 0; i < Row; ++i)
			{
				float dist = FLT_MAX;
				float temp;
				int clusTemp;
				for (int j = 0; j < numberOfClusters; ++j)
				{
					temp = (eigenVec.row(i)-clusterCenter.row(j)).norm();
					if(temp<dist)
					{
						dist = temp;
						clusTemp = j;
					}
				}

			#pragma omp critical
				{
					storage[clusTemp]++;
					neighborVec[clusTemp].push_back(i);
					group[i] = clusTemp;
					centerTemp.row(clusTemp)+=eigenVec.row(i);
				}
			}
		}

		moving = FLT_MIN;

	#pragma omp parallel for reduction(max:moving) num_threads(8)
		for (int i = 0; i < numberOfClusters; ++i)
		{
			if(storage[i]>0)
			{
				centerTemp.row(i)/=storage[i];
				tempMoving = (centerTemp.row(i)-clusterCenter.row(i)).norm();
				clusterCenter.row(i) = centerTemp.row(i);
				if(moving<tempMoving)
					moving = tempMoving;
			}
		}
		std::cout << "K-means iteration " << ++tag << " completed, and moving is "
		<< moving << "!" << std::endl;
	}while(abs(moving-before)/before >= 1.0e-3 && tag < 50 && moving>0.01);

	gettimeofday(&end, NULL);
	float timeTemp = ((end.tv_sec-start.tv_sec)*1000000u+end.tv_usec-start.tv_usec)/1.e6;
	activityList.push_back("K-means takes: ");
	timeList.push_back(to_string(timeTemp)+" s");

	for(auto iter=storage.begin(); iter!=storage.end();)
	{
		if(*iter==0)
			storage.erase(iter);
		else
			++iter;
	}
	numberOfClusters = storage.size();
}


/*
 * @brief get cluster information based on eigenvector rotation
 *
 * @param[out] storage The size of clusters
 * @param[out] neighborVec The candidates belonging to clusters
 * @param[out] clusterCenter The centroid coordinates of clusters
 * @param[in] X The matrix X
 */
void SpectralClustering::getEigvecRotation(std::vector<int>& storage, std::vector<std::vector<int> >& neighborVec,
        								   Eigen::MatrixXf& clusterCenter, const Eigen::MatrixXf& X)
{
	mMaxQuality = 0;
	Eigen::MatrixXf vecRot;
	Eigen::MatrixXf vecIn = X.block(0,0,X.rows(),2);
	Evrot *e = NULL;

	struct timeval start, end;
	gettimeofday(&start, NULL);

	const int& xCols = X.cols();

	std::cout << "Eigenvector rotation starts within " << xCols << " columns..." << std::endl;
	for (int g=2; g <= xCols; g++)
	{
		// make it incremental (used already aligned vectors)
		std::cout << "column " << g << ":";
		if( g > 2 )
		{
			vecIn.resize(X.rows(),g);
			vecIn.block(0,0,vecIn.rows(),g-1) = e->getRotatedEigenVectors();
			vecIn.block(0,g-1,X.rows(),1) = X.block(0,g-1,X.rows(),1);
			delete e;
		}
		//perform the rotation for the current number of dimensions
		e = new Evrot(vecIn, mMethod);

		//save max quality
		if (e->getQuality() > mMaxQuality)
		{
			mMaxQuality = e->getQuality();
		}

		if(isnan(e->getQuality())||isinf(e->getQuality()))
		{
			std::cout << "Meet with nan or inf! Stop! " << std::endl;
			return;
		}

		std::cout << " max quality is " << mMaxQuality << ", Evrot has quality " << e->getQuality() << std::endl;
		//save cluster data for max cluster or if we're near the max cluster (so prefer more clusters)
		if ((e->getQuality() > mMaxQuality) || (mMaxQuality - e->getQuality() <= 0.001))
		{
			neighborVec = e->getClusters();
			vecRot = e->getRotatedEigenVectors();
		}
	}

	gettimeofday(&end, NULL);
	float timeTemp = ((end.tv_sec-start.tv_sec)*1000000u+end.tv_usec-start.tv_usec)/1.e6;
	activityList.push_back("Eigenvector rotation takes: ");
	timeList.push_back(to_string(timeTemp)+" s");

	if(neighborVec.empty())
		return;

	clusterCenter = Eigen::MatrixXf::Zero(neighborVec.size(),vecRot.cols());
	storage = std::vector<int>(neighborVec.size());

#pragma omp parallel for schedule(static) num_threads(8)
	for (unsigned int i=0; i < neighborVec.size(); i++)
	{
		storage[i] = neighborVec[i].size();
		for (unsigned int j=0; j < neighborVec[i].size(); j++)
		{
			//sum points within cluster
			clusterCenter.row(i) += vecRot.row(neighborVec[i][j]);
		}
	}

#pragma omp parallel for schedule(static) num_threads(8)
	for (unsigned int i=0; i < neighborVec.size(); i++) {
		//find average point within cluster
		clusterCenter.row(i) = clusterCenter.row(i) / neighborVec[i].size();
	}

	numberOfClusters = neighborVec.size();
}


/*
 * @brief set automatic parameter from the Para object input
 *
 * @param[in] p The Para object for the parameters
 */
void SpectralClustering::setParameterAutomatic(const Para& p)
{
	std::cout << "It is a pathline data set? 1.Yes, 0.No." << std::endl;
	int pathlineOption;
	std::cin >> pathlineOption;
	assert(pathlineOption==1||pathlineOption==0);
	isPathlines = (pathlineOption==1);

	if(isPathlines)
	{
		IOHandler::expandArray(ds.dataMatrix,ds.dataVec,ds.dimension,ds.maxElements);
	}
	else
	{
		if(p.sampled==1)
			IOHandler::expandArray(ds.dataMatrix,ds.dataVec,ds.dimension,ds.maxElements);
		else if(p.sampled==2)
			IOHandler::sampleArray(ds.dataMatrix,ds.dataVec,ds.dimension,ds.maxElements);
		else if(p.sampled==3)
			IOHandler::uniformArcSampling(ds.dataMatrix,ds.dataVec,ds.dimension,ds.maxElements);
	}

	group = std::vector<int>(ds.dataMatrix.rows());

	/* the default value for streamline clustering is 2 normalized Laplacian */
	LaplacianOption = p.LaplacianOption;

	isDistSorted = p.isDistSorted;

	postProcessing = p.postProcessing;

	mMethod = p.mMethod;

	extractOption = p.extractOption;

}


/*
 * @brief set necessary parameter
 */
void SpectralClustering::getParameterUserInput()
{
	std::cout << "It is a pathline data set? 1.Yes, 0.No." << std::endl;
	int pathlineOption;
	std::cin >> pathlineOption;
	assert(pathlineOption==1||pathlineOption==0);
	isPathlines = (pathlineOption==1);

	int sampleOption;

	if(isPathlines)
		sampleOption = 1;
	else
	{
		std::cout << "choose a sampling method for the dataset?" << std::endl
				  << "1.directly filling with last vertex; 2. uniform sampling." << std::endl;
		std::cin >> sampleOption;
	}
	assert(sampleOption==1||sampleOption==2);

	if(sampleOption==1)
		IOHandler::expandArray(ds.dataMatrix,ds.dataVec,ds.dimension,ds.maxElements);
	else if(sampleOption==2)
		IOHandler::sampleArray(ds.dataMatrix,ds.dataVec,ds.dimension,ds.maxElements);
	else if(sampleOption==3)
		IOHandler::uniformArcSampling(ds.dataMatrix,ds.dataVec,ds.dimension,ds.maxElements);

	group = std::vector<int>(ds.dataMatrix.rows());

	/* the default value for streamline clustering is 2 normalized Laplacian */
	std::cout << "---------------------------" << std::endl;
	std::cout << "Laplacian option: 1.Normalized Laplacian, 2.Unsymmetric Laplacian" << std::endl;
	std::cout << "..And in streamline clustering people tend to choose 1.Normalized Laplacian!-----------" << std::endl;
	std::cin >> LaplacianOption;
	assert(LaplacianOption==1||LaplacianOption==2);


	int sortedOption;
	std::cout << "Please choose whether local scaling by sorted distance: 1. yes, 2. no: " << std::endl;
	std::cin >> sortedOption;
	assert(sortedOption==1||sortedOption==2);
	if(sortedOption==1)
		sortedOption = true;
	else if(sortedOption==2)
		sortedOption = false;

	std::cout << "-----------------------------------------------------------------------" << std::endl;
	std::cout << "Input a desired cluster number among [1, " << ds.dataMatrix.rows() << "]: ";
	std::cin >> presetNumber;
	assert(presetNumber>1 && presetNumber<ds.dataMatrix.rows()/10);

	std::cout << "-----------------------------------------------------------------------" << std::endl;
	std::cout << "Input a post-processing method: 1.k-means, 2.eigenvector rotation: " << std::endl;
	std::cin >> postProcessing;
	assert(postProcessing==1||postProcessing==2);

	if(postProcessing==2)
	{
		std::cout << "------------------------------------------------" << std::endl;
		std::cout << "Please input derivative method: 1.numerical derivative, 2.true derivative." << std::endl;
		std::cin >> mMethod;
		assert(mMethod==1 || mMethod==2);
	}

    std::cout << "Please choose cluster number method, 0.user input, 1.read clustering: " << std::endl;
    int clusterInput;
    std::cin >> clusterInput;
    assert(clusterInput==0 || clusterInput==1);
    readCluster = (clusterInput==1);

}


/*
 * @brief record the preset number of clusters
 *
 * @param[in] number The number of clusters as input
 */
void SpectralClustering::recordPreset(const int& number)
{
	std::ofstream readme("../dataset/optimal.txt",ios::out | ios::app);
	if(!readme)
	{
		std::cout << "Error creating readme!" << std::endl;
		exit(1);
	}
	readme << "Preset cluster number is: " << number << std::endl;
	readme << std::endl;
	readme.close();
}


/*
 * @brief record the result of optimal number of clusters
 *
 * @param[in] normOption The norm option
 * @param[in] clusNum The number of clusters
 */
void SpectralClustering::recordOptimalResult(const int& normOption, const int& clusNum)
{
	std::ofstream readme("../dataset/optimal.txt",ios::out | ios::app);
	if(!readme)
	{
		std::cout << "Error creating readme!" << std::endl;
		exit(1);
	}
	readme << "Optimal number of cluster for norm " << normOption << " with sc eigen-rotation minimization is "
		   << clusNum << std::endl;
	readme << std::endl;
	readme.close();
}
