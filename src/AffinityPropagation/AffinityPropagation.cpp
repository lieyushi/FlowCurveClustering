#include "AffinityPropagation.h"

/* default constructor */
AffinityPropagation::AffinityPropagation()
{

}

/* argument constructor with argc and argv */
AffinityPropagation::AffinityPropagation(const int& argc, char **argv, const Para& p, bool& automatic)
{
	setDataset(argc, argv);

	if(automatic)
		setParameterAutomatic(p);

	else
		getParameterUserInput();

}

/* destructor */
AffinityPropagation::~AffinityPropagation()
{
	deleteDistanceMatrix(ds.dataMatrix.rows());
}

/* perform clustering function */
void AffinityPropagation::performClustering()
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
	*/
	for(int i=0;i<=15;++i)
	{
		/* don't want to deal with many too naive metrics */
		if(i!=0 && i!=1 && i!=2 && i!=4 && i!=12 && i!=14 && i!=15)
			continue;


		std::cout << "----------------------------------------------------" << std::endl;
		std::cout << "Experiment on norm " << i << " starts!--------------" << std::endl;

		activityList.clear();
		timeList.clear();

		clusterByNorm(i);

		std::cout << std::endl;
	}
}


/* run clustering based on different norm */
void AffinityPropagation::clusterByNorm(const int& norm)
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

	if(!getDistanceMatrix(ds.dataMatrix, normOption, object))
	{
		std::cout << "Failure to compute distance matrix!" << std::endl;
	}

	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec  - start.tv_sec) * 1000000u
			   + end.tv_usec - start.tv_usec) / 1.e6;
	activityList.push_back("Distance matrix computing takes: ");
	timeList.push_back(to_string(timeTemp)+" s");

	Eigen::MatrixXf matrixR, matrixA, matrixS;

	gettimeofday(&start, NULL);

	/* initialize S, R, A */
	initializeMatrices(matrixS, matrixR, matrixA);

	/* get S */
	getMatrixS(matrixS);

	int current = 0;
	while(current++<maxIteration)
	{
		std::cout << "Iteration " << current << std::endl;

		/* update responsibility */
		updateResponsibility(matrixR, matrixA, matrixS);

		/* update availability */
		updateAvailability(matrixA, matrixR);

	}

	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec  - start.tv_sec) * 1000000u+ end.tv_usec - start.tv_usec) / 1.e6;

	activityList.push_back("Affinity propagation takes: ");
	timeList.push_back(to_string(timeTemp)+" s");

	std::vector<std::vector<int> > neighborVec;
	std::vector<int> storage;
	Eigen::MatrixXf centroid;

	/* get exemplary examples */
	getGroupAssignment(matrixR, matrixA, matrixS, neighborVec, storage);

	setLabel(neighborVec, storage, centroid);

	activityList.push_back("Affinity propagation generates: ");
	timeList.push_back(to_string(storage.size())+" groups");

	extractFeatures(storage, neighborVec, centroid);

}


/* perform group-labeling information */
void AffinityPropagation::setLabel(vector<vector<int> >& neighborVec, vector<int>& storage, Eigen::MatrixXf& centroid)
{
	std::vector<Ensemble> nodeVec(storage.size());

	std::cout << "Cluster label setting begins with " << nodeVec.size() << " clusters..." << std::endl;
#pragma omp parallel for schedule(dynamic) num_threads(8)
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

#pragma omp parallel for schedule(dynamic) num_threads(8)
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



/* extract features from datasets as representative curves */
void AffinityPropagation::extractFeatures(const std::vector<int>& storage, const std::vector<std::vector<int> >& neighborVec,
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

	float EntropyRatio;
	getEntropyRatio(storage, EntropyRatio);

	IOHandler::printClusters(ds.dataVec,group,storage,"AP_norm"+to_string(normOption),ds.fullName,ds.dimension);

	/* if the dataset is not PBF, then should record distance matrix for Gamma matrix compution */
	if(!isPBF)
	{
		deleteDistanceMatrix(ds.dataMatrix.rows());

		if(!getDistanceMatrix(ds.dataMatrix, normOption, object))
		{
			std::cout << "Failure to compute distance matrix!" << std::endl;
		}
	}


	struct timeval start, end;
	double timeTemp;

	gettimeofday(&start, NULL);
	Silhouette sil;
	sil.computeValue(normOption,ds.dataMatrix,ds.dataMatrix.rows(),ds.dataMatrix.cols(),group,object,
			         numberOfClusters, isPBF, neighborVec);
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

	activityList.push_back("numCluster is: ");
	timeList.push_back(to_string(numberOfClusters));

	activityList.push_back("Norm option is: ");
	timeList.push_back(to_string(normOption));


	IOHandler::generateReadme(activityList,timeList);

/* print entropy value for the clustering algorithm */
	IOHandler::writeReadme(EntropyRatio, sil);

	IOHandler::writeReadme(closestAverage, furthestAverage);

}


/* set dataset from user command */
void AffinityPropagation::setDataset(const int& argc, char **argv)
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

/* get the bool tag for isPBF */
	std::cout << "It is a PBF dataset? 1.Yes, 0.No" << std::endl;
	int PBFjudgement;
	std::cin >> PBFjudgement;
	assert(PBFjudgement==1||PBFjudgement==0);
	isPBF = (PBFjudgement==1);

	IOHandler::readFile(ds.strName,ds.dataVec,ds.vertexCount,ds.dimension,ds.maxElements);

	ds.fullName = ds.strName+"_full.vtk";
	IOHandler::printVTK(ds.fullName, ds.dataVec, ds.vertexCount, ds.dimension);
}


/* get entropy ratio, lower value tells dinstinguishable cluster while higher value tells a more uniformality. */
void AffinityPropagation::getEntropyRatio(const std::vector<int>& storage, float& EntropyRatio)
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



/* set automatic parameter */
void AffinityPropagation::setParameterAutomatic(const Para& p)
{

	if(p.sampled==1)
		IOHandler::expandArray(ds.dataMatrix,ds.dataVec,ds.dimension,ds.maxElements);
	else if(p.sampled==2)
		IOHandler::sampleArray(ds.dataMatrix,ds.dataVec,ds.dimension,ds.maxElements);
	else if(p.sampled==3)
		IOHandler::uniformArcSampling(ds.dataMatrix,ds.dataVec,ds.dimension,ds.maxElements);

	group = std::vector<int>(ds.dataMatrix.rows());

	extractOption = p.extractOption;

	maxIteration = p.maxIteration;

}



/* set parameter */
void AffinityPropagation::getParameterUserInput()
{
	int sampleOption;
	std::cout << "choose a sampling method for the dataset?" << std::endl
			  << "1.directly filling with last vertex; 2. uniform sampling." << std::endl;
	std::cin >> sampleOption;
	assert(sampleOption==1||sampleOption==2);

	if(sampleOption==1)
		IOHandler::expandArray(ds.dataMatrix,ds.dataVec,ds.dimension,ds.maxElements);
	else if(sampleOption==2)
		IOHandler::sampleArray(ds.dataMatrix,ds.dataVec,ds.dimension,ds.maxElements);
	else if(sampleOption==3)
		IOHandler::uniformArcSampling(ds.dataMatrix,ds.dataVec,ds.dimension,ds.maxElements);

	group = std::vector<int>(ds.dataMatrix.rows());

	std::cout << "Select extraction method: 1.centroid, closest and furthest, 2.median, 3.statiscal representation." << std::endl;
	std::cin >> extractOption;
	assert(extractOption==1||extractOption==2||extractOption==3);

	std::cout << "Input max iteration for affinity propagation: " << std::endl;
	std::cin >> maxIteration;
	assert(maxIteration>0);


}


/* get matrix S from distance matrix */
void AffinityPropagation::getMatrixS(Eigen::MatrixXf& matrixS)
{
	std::cout << "Start initializing matrix S..." << std::endl;

	const int& rows = matrixS.rows();

	/* define a vector to store pair-wise distance vector and get the median */
	const int& distVecSize = rows*(rows-1)/2;
	std::vector<float> distVec(distVecSize);
	int count = 0;

	/* fill the matrix S */
	float tempDist;
	for(int i=0;i<rows-1;++i)
	{
		for(int j=i+1;j<rows;++j)
		{
			if(distanceMatrix)
				tempDist = distanceMatrix[i][j];
			else
				tempDist = getDisimilarity(ds.dataMatrix, i, j, normOption, object);

			/* conventionally we assign -d*d as non-diagonal entries for matrix S */
			matrixS(i,j) = -tempDist*tempDist;
			matrixS(i,j) = matrixS(j,i);

			distVec[count++] = matrixS(i,j);
		}
	}

	assert(count==distVecSize);

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

	/* assign medianValue to diagonal matrix element */
#pragma omp parallel for schedule(dynamic) num_threads(8)
	for(int i=0;i<rows;++i)
		matrixS(i,i) = medianValue;

	std::cout << "Finish initializing matrix S..." << std::endl;

}


/* initialize matrix S, R, A */
void AffinityPropagation::initializeMatrices(Eigen::MatrixXf& matrixS, Eigen::MatrixXf& matrixR,
											Eigen::MatrixXf& matrixA)
{
	const int& rows = ds.dataMatrix.rows();

	/* initialize all three matrices as zero entry */
	matrixS = Eigen::MatrixXf::Zero(rows, rows);
	matrixR = Eigen::MatrixXf::Zero(rows, rows);
	matrixA = Eigen::MatrixXf::Zero(rows, rows);

}


/* update responsibility matrix R */
void AffinityPropagation::updateResponsibility(Eigen::MatrixXf& matrixR, const Eigen::MatrixXf& matrixA,
											  const Eigen::MatrixXf& matrixS)
{
	const int& rows = matrixR.rows();
#pragma omp parallel for schedule(dynamic) num_threads(8)
	for(int i=0;i<rows;++i)
	{
		for(int k=0;k<rows;++k)
		{
			/* don't use FLT_MIN because FLT_MIN == 0.0 */
			float maxValue = float(-INT_MIN);
			for(int kk=0;kk<rows;++kk)
			{
				if(kk==k)
					continue;
				maxValue = std::max(maxValue, matrixR(i,kk)+matrixA(i,kk));
			}
			/* in wikipage it's update by R[i,k] = S[i][k]-maxValue, but here use a Laplace smoothor for convergence */
			matrixR(i,k) = (1-LAMBDA)*(matrixS(i,k)-maxValue)+LAMBDA*matrixR(i,k);
		}
	}

}


/* update availability matrix */
void AffinityPropagation::updateAvailability(Eigen::MatrixXf& matrixA, const Eigen::MatrixXf& matrixR)
{
	const int& rows = matrixR.rows();
#pragma omp parallel for schedule(dynamic) num_threads(8)
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


/* get assignment by three matrices */
void AffinityPropagation::getGroupAssignment(const Eigen::MatrixXf& matrixR, const Eigen::MatrixXf& matrixA,
			  	  	  	  	const Eigen::MatrixXf& matrixS, std::vector<std::vector<int> >& neighborVec,
							std::vector<int>& storage)
{
	std::vector<int> centerVec;
	const int& rows = matrixR.rows();

	/* store the candidate whose diagonal summation is positive */
	float diagonalSum;
	for(int i=0;i<rows;++i)
	{
		diagonalSum=matrixR(i,i)+matrixA(i,i);
		if(diagonalSum>0)
			centerVec.push_back(i);
	}

	const int& centerSize = centerVec.size();
	/* get group tag information for each candidate streamline */
#pragma omp parallel for schedule(dynamic) num_threads(8)
	for(int i=0;i<rows;++i)
	{
		int index, element;
		float maxSim = (float)(-INT_MIN);
		for(int j=0;j<centerSize;++j)
		{
			element = centerVec[j];
			if(matrixS(i,element)>maxSim)
			{
				maxSim = matrixS(i,element);
				index = element;
			}
		}
		group[i]=element;
	}

	/* output group information and cluster size */
	std::map<int,int> groupMap;
	for(int i=0;i<rows;++i)
	{
		/* group tag not int the hash map */
		if(groupMap.find(group[i])==groupMap.end())
			groupMap.insert(make_pair(group[i],0));
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

	for(int i=0;i<storage.size();++i)
		storage[i] = neighborVec[i].size();
}
