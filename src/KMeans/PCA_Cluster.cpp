#include "PCA_Cluster.h"

const float& TOR_1 = 0.999;
const int& CLUSTER = 8;

extern int initializationOption;

/* an external value to judge whether it is a PBF or not */
extern bool isPBF;

void PCA_Cluster::performPCA_Clustering(const Eigen::MatrixXf& data, 
										const int& Row, 
										const int& Column, 
										std::vector<MeanLine>& massCenter,
									    std::vector<int>& group, 
									    std::vector<int>& totalNum, 
									    std::vector<ExtractedLine>& closest,
									    std::vector<ExtractedLine>& furthest,
										TimeRecorder& tr,
										Silhouette& sil)
{
	MatrixXf cArray, SingVec;
	VectorXf meanTrajectory(Column);
	int PC_Number;



	performSVD(cArray, data, Row, Column, PC_Number, SingVec, meanTrajectory, tr);
	performPC_KMeans(cArray, Row, Column, PC_Number, SingVec, meanTrajectory, 
					 massCenter, CLUSTER, group, totalNum, closest, furthest, data, tr, sil);
}



void PCA_Cluster::performSVD(MatrixXf& cArray, 
							 const Eigen::MatrixXf& data, 
							 const int& Row, 
							 const int& Column,
							 int& PC_Number, 
							 MatrixXf& SingVec, 
							 VectorXf& meanTrajectory,
							 TimeRecorder& tr)
{
	Eigen::MatrixXf temp = data;

#pragma omp parallel for schedule(dynamic) num_threads(8)
	for (int i = 0; i < Column; ++i)
	{
		meanTrajectory(i) = temp.transpose().row(i).mean();
	}
#pragma omp parallel for schedule(dynamic) num_threads(8)
	for (int i = 0; i < Row; ++i)
	{
		temp.row(i) = temp.row(i)-meanTrajectory.transpose();
	}

	struct timeval start, end;
	gettimeofday(&start, NULL);
	/* perform SVD decomposition for temp */
	JacobiSVD<MatrixXf> svd(temp, ComputeThinU | ComputeThinV);
	//const VectorXf& singValue = svd.singularValues();
	SingVec = svd.matrixV();
	gettimeofday(&end, NULL);
	const double& delta = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;

	tr.eventList.push_back("SVD takes ");
	tr.timeList.push_back(to_string(delta)+"s");

	/* compute new attribute space based on principal component */
	MatrixXf coefficient = temp*SingVec;
	/*  decide first r dorminant PCs with a threshold */
	const float& varianceSummation = coefficient.squaredNorm();
	float tempSum = 0.0;
	const float& threshold = TOR_1*varianceSummation;
	
	for (int i = 0; i < Column; ++i)
	{
		tempSum+=(coefficient.transpose().row(i)).squaredNorm();
		if(tempSum>threshold)
		{
			PC_Number = i;
			break;
		}
	}

	cArray = MatrixXf(Row, PC_Number);
#pragma omp parallel for schedule(dynamic) num_threads(8)
	for (int i = 0; i < PC_Number; ++i)
	{
		cArray.transpose().row(i) = coefficient.transpose().row(i);
	}

	std::cout << "SVD completed!" << std::endl;

	SingVec.transposeInPlace();
}



void PCA_Cluster::performPC_KMeans(const MatrixXf& cArray, 
								   const int& Row, 
								   const int& Column, 
								   const int& PC_Number, 
				 				   const MatrixXf& SingVec, 
				 				   const VectorXf& meanTrajectory, 
				 				   std::vector<MeanLine>& massCenter, 
				 				   const int& Cluster, 
				 				   std::vector<int>& group, 
				 				   std::vector<int>& totalNum, 
				 				   std::vector<ExtractedLine>& closest,
				 				   std::vector<ExtractedLine>& furthest, 
				 				   const Eigen::MatrixXf& data,
								   TimeRecorder& tr,
								   Silhouette& sil)
{
	MetricPreparation object(Row, Column);
	object.preprocessing(data, Row, Column, 0);
/* perform K-means clustering */
	MatrixXf clusterCenter;

	switch(initializationOption)
	{
	case 1:
		Initialization::generateRandomPos(clusterCenter, PC_Number, cArray, Cluster);
		break;

	case 2:
		Initialization::generateFromSamples(clusterCenter, PC_Number, cArray, Cluster);
		break;

	case 3:
		Initialization::generateFarSamples(clusterCenter, PC_Number, cArray, 
										   Cluster, 0, object);
		break;
	}

	std::cout << "Initialization completed!" << std::endl;

	float moving=1000, tempMoving, before;
	int storage[Cluster];

	MatrixXf centerTemp;  //store provisional center coordinate

	int tag = 0;

	std::vector< std::vector<int> > neighborVec(Cluster, std::vector<int>());

	double PCA_KMeans_delta, KMeans_delta;
	struct timeval start, end;

	gettimeofday(&start, NULL);

	std::vector<int> recorder(Row);
	do
	{
		before = moving;
		/* preset cluster number recorder */
		memset(storage,0,sizeof(int)*Cluster);
		centerTemp = MatrixXf::Zero(Cluster, PC_Number);

	#pragma omp parallel for schedule(dynamic) num_threads(8)
		for (int i = 0; i < Cluster; ++i)
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
				for (int j = 0; j < Cluster; ++j)
				{
					temp = (cArray.row(i)-clusterCenter.row(j)).norm();
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
					recorder[i] = clusTemp;
					centerTemp.row(clusTemp)+=cArray.row(i);
				}
			}
		}

		moving = FLT_MIN;

	#pragma omp parallel for reduction(max:moving) num_threads(8)
		for (int i = 0; i < Cluster; ++i)
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
	}while(abs(moving-before)/before >= 1.0e-2 && tag < 20 && moving>0.01);

	gettimeofday(&end, NULL);
	
	float delta = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;

	tr.eventList.push_back("k-means iteration for PC takes ");
	tr.timeList.push_back(to_string(delta)+"s");

	std::multimap<int,int> groupMap;

	float entropy = 0.0;
	float probability;


	for (int i = 0; i < Cluster; ++i)
	{
		groupMap.insert(std::pair<int,int>(storage[i],i));
		if(storage[i]>0)
		{
			probability = float(storage[i])/float(Row);
			entropy += probability*log2f(probability);
		}
	}

	int groupNo = 0;
	int increasingOrder[Cluster];
	for (multimap<int,int>::iterator it = groupMap.begin(); it != groupMap.end(); ++it)
	{
		if(it->first>0)
		{
			increasingOrder[it->second] = (groupNo++);
		}
	}

	/* calculate the balanced entropy */
	entropy = -entropy/log2f(groupNo);


#pragma omp parallel for schedule(dynamic) num_threads(8)	
	for (int i = 0; i < Row; ++i)
	{
		group[i] = increasingOrder[recorder[i]];
		totalNum[i] = storage[recorder[i]];
	}

	float shortest, farDist, toCenter;
	int shortestIndex = 0, fartestIndex = 0, tempIndex = 0;
	std::vector<int> neighborTemp;

	for (int i = 0; i < Cluster; ++i)
	{
		if(storage[i]>0 && !neighborVec[i].empty())
		{
			neighborTemp = neighborVec[i];
			shortest = FLT_MAX;
			farDist = FLT_MIN;

			for (int j = 0; j < storage[i]; ++j)
			{
				tempIndex = neighborTemp[j];
				toCenter = (clusterCenter.row(i)-cArray.row(tempIndex)).norm();

				if(toCenter<shortest)
				{
					shortest = toCenter;
					shortestIndex = tempIndex;
				}
				if(toCenter>farDist)
				{
					farDist = toCenter;
					fartestIndex = tempIndex;
				}
			}
			closest.push_back(ExtractedLine(shortestIndex,increasingOrder[i]));
			furthest.push_back(ExtractedLine(fartestIndex,increasingOrder[i]));
		}
	}

	MatrixXf pcSing(PC_Number,Column);

#pragma omp parallel for schedule(dynamic) num_threads(8)
	for (int i = 0; i < PC_Number; ++i)
	{
		pcSing.row(i) = SingVec.row(i);
	}

	MatrixXf massPos = clusterCenter*pcSing;

	for (int i = 0; i < Cluster; ++i)
	{
		if(storage[i]>0)
		{
			massPos.row(i) += meanTrajectory.transpose();
			std::vector<float> vecTemp;
			for (int j = 0; j < Column; ++j)
			{
				vecTemp.push_back(massPos(i,j));
			}
			massCenter.push_back(MeanLine(vecTemp,increasingOrder[i]));
		}
	}

	std::cout << "Mean position re-mapped to trajectory space completed!" << std::endl;

/* Silhouette effect */
	gettimeofday(&start, NULL);

	sil.computeValue(cArray,group,groupNo,isPBF);

	std::cout << "Silhouette computation completed!" << std::endl;

	gettimeofday(&end, NULL);
	delta = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;

	tr.eventList.push_back("Clustering evaluation computing takes: ");
	tr.timeList.push_back(to_string(delta)+"s");

	ValidityMeasurement vm;
	vm.computeValue(cArray, group);
	tr.eventList.push_back("PCA Validity measure is: ");
	tr.timeList.push_back(to_string(vm.f_c));

	/* write value of the silhouette class */
	IOHandler::writeReadme(entropy, sil);
}



void PCA_Cluster::performAHC(const MatrixXf& cArray, 
							 const int& Row, 
							 const int& Column, 
							 const int& PC_Number, 
		      				 const MatrixXf& SingVec, 
		      				 const VectorXf& meanTrajectory, 
		      				 MatrixXf& clusterCenter, 
		      				 std::vector<MeanLine>& massCenter)
{
	return;
}


void PCA_Cluster::performDirectK_Means(const Eigen::MatrixXf& data, 
									   const int& Row, 
									   const int& Column, 
									   std::vector<MeanLine>& massCenter,
									   std::vector<int>& group, 
									   std::vector<int>& totalNum, 
									   std::vector<ExtractedLine>& closest,
									   std::vector<ExtractedLine>& furthest, 
									   const int& normOption,
									   TimeRecorder& tr,
									   Silhouette& sil)
{

	performFullK_MeansByClusters(data, Row, Column, massCenter, CLUSTER, group, 
								 totalNum, closest, furthest, normOption, tr, sil);
}


void PCA_Cluster::performPCA_Clustering(const Eigen::MatrixXf& data, 
										const int& Row, 
										const int& Column, 
										std::vector<MeanLine>& massCenter,
										std::vector<int>& group, 
										std::vector<int>& totalNum, 
										std::vector<ExtractedLine>& closest, 
										std::vector<ExtractedLine>& furthest, 
										const int& Cluster,
										TimeRecorder& tr,
										Silhouette& sil)
{
	MatrixXf cArray, SingVec;
	VectorXf meanTrajectory(Column);
	int PC_Number;

	performSVD(cArray, data, Row, Column, PC_Number, SingVec, meanTrajectory, tr);
	performPC_KMeans(cArray, Row, Column, PC_Number, SingVec, meanTrajectory, 
					 massCenter, Cluster, group, totalNum, closest, furthest, data, tr, sil);
}


void PCA_Cluster::performDirectK_Means(const Eigen::MatrixXf& data, 
									   const int& Row, 
									   const int& Column, 
									   std::vector<MeanLine>& massCenter,
									   std::vector<int>& group, 
									   std::vector<int>& totalNum, 
									   std::vector<ExtractedLine>& closest, 
									   std::vector<ExtractedLine>& furthest, 
									   const int& Cluster, 
									   const int& normOption,
									   TimeRecorder& tr,
									   Silhouette& sil)
{
	performFullK_MeansByClusters(data, Row, Column, massCenter, Cluster, group, 
								 totalNum, closest, furthest, normOption, tr, sil);
}


void PCA_Cluster::performFullK_MeansByClusters(const Eigen::MatrixXf& data, 
											   const int& Row, 
											   const int& Column, 
											   std::vector<MeanLine>& massCenter,
											   const int& Cluster, 
											   std::vector<int>& group, 
											   std::vector<int>& totalNum, 
											   std::vector<ExtractedLine>& closest, 
											   std::vector<ExtractedLine>& furthest, 
											   const int& normOption,
											   TimeRecorder& tr,
											   Silhouette& sil)
{	
	MetricPreparation object(Row, Column);
	object.preprocessing(data, Row, Column, normOption);

	MatrixXf clusterCenter;

	switch(initializationOption)
	{
	case 1:
		Initialization::generateRandomPos(clusterCenter, Column, data, Cluster);
		break;

	case 2:
		Initialization::generateFromSamples(clusterCenter, Column, data, Cluster);
		break;

	case 3:
		Initialization::generateFarSamples(clusterCenter, Column, data, Cluster, 
										   normOption, object);
		break;
	}

	std::cout << "Initialization completed!" << std::endl;

	float moving=1000, tempMoving,/* dist, tempDist, */before;
	int *storage = new int[Cluster]; // used to store number inside each cluster
	MatrixXf centerTemp;
	int tag = 0;
	std::vector< std::vector<int> > neighborVec(Cluster, std::vector<int>());

/* perform K-means with different metrics */
	std::cout << "K-means start!" << std::endl;	
	struct timeval start, end;
	gettimeofday(&start, NULL);
	std::vector<int> recorder(Row); //use to record which cluster the row belongs to

	do
	{
	/* reset storage number and weighted mean inside each cluster*/
		before=moving;
		memset(storage,0,sizeof(int)*Cluster);
		centerTemp = MatrixXf::Zero(Cluster,Column);

	/* clear streamline indices for each cluster */
	#pragma omp parallel for schedule(dynamic) num_threads(8)
		for (int i = 0; i < Cluster; ++i)
		{
			neighborVec[i].clear();
		}

	#pragma omp parallel num_threads(8)
		{
		#pragma omp for nowait
			for (int i = 0; i < Row; ++i)
			{
				int clusTemp;
				float dist = FLT_MAX;
				float tempDist;
				for (int j = 0; j < Cluster; ++j)
				{
					tempDist = getDisimilarity(clusterCenter.row(j),data,i,normOption,object);
					if(tempDist<dist)
					{
						dist = tempDist;
						clusTemp = j;
					}
				}
				recorder[i] = clusTemp;

			#pragma omp critical
				{
					storage[clusTemp]++;
					neighborVec[clusTemp].push_back(i);
					centerTemp.row(clusTemp)+=data.row(i);
				}
			}
		}
		moving = FLT_MIN;

	/* measure how much the current center moves from original center */	
	#pragma omp parallel for reduction(max:moving) num_threads(8)
		for (int i = 0; i < Cluster; ++i)
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
		std::cout << "K-means iteration " << ++tag << " completed, and moving is " << moving 
				  << "!" << std::endl;
	}while(abs(moving-before)/before >= 1.0e-2 && tag < 20 && moving > 0.01);
	
	gettimeofday(&end, NULL);
	double delta = ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;

	tr.eventList.push_back("k-means iteration takes ");
	tr.timeList.push_back(to_string(delta)+"s");

	std::multimap<int,int> groupMap;

	float entropy = 0.0, probability;
	int increasingOrder[Cluster];

	int nonZero = 0;
	for (int i = 0; i < Cluster; ++i)
	{
		groupMap.insert(std::pair<int,int>(storage[i],i));
		if(storage[i]>0)
		{
			probability=float(storage[i])/float(Row);
			entropy+=probability*log2f(probability);
			++nonZero;
		}
	}
	entropy = -entropy/log2f(nonZero);

	int groupNo = 0;
	for (std::multimap<int,int>::iterator it = groupMap.begin(); it != groupMap.end(); ++it)
	{
		if(it->first>0)
		{
			increasingOrder[it->second] = (groupNo++);
		}
	}
	/* finish tagging for each group */

	/* record labeling information */
	IOHandler::generateGroups(neighborVec);

	// set cluster group number and size number 
#pragma omp parallel for schedule(dynamic) num_threads(8)
	for (int i = 0; i < Row; ++i)
	{
		group[i] = increasingOrder[recorder[i]];
		totalNum[i] = storage[recorder[i]];
	}

	float shortest, toCenter, farDist;
	int shortestIndex = 0, tempIndex = 0, furthestIndex = 0;
	std::vector<int> neighborTemp;

	/* choose cloest and furthest streamlines to centroid streamlines */
	for (int i = 0; i < Cluster; ++i)
	{
		if(storage[i]>0)
		{

			neighborTemp = neighborVec[i];
			shortest = FLT_MAX;
			farDist = FLT_MIN;

			for (int j = 0; j < storage[i]; ++j)
			{
				// j-th internal streamlines 
				tempIndex = neighborTemp[j];
				toCenter = getDisimilarity(clusterCenter.row(i),data,tempIndex,normOption,object);

				/* update the closest index to centroid */
				if(toCenter<shortest)
				{
					shortest = toCenter;
					shortestIndex = tempIndex;
				}

				/* update the farthest index to centroid */
				if(toCenter>farDist)
				{
					farDist = toCenter;
					furthestIndex = tempIndex;
				}
			}
			closest.push_back(ExtractedLine(shortestIndex,increasingOrder[i]));
			furthest.push_back(ExtractedLine(furthestIndex,increasingOrder[i]));
			//distFile << std::endl;
		}
	}
	//distFile.close();

	std::vector<float> closeSubset;
	/* based on known cluster centroid, save them as vector for output */
	for (int i = 0; i < Cluster; ++i)
	{
		if(storage[i]>0)
		{
			for (int j = 0; j < Column; ++j)
			{
				closeSubset.push_back(clusterCenter(i,j));
			}
			massCenter.push_back(MeanLine(closeSubset,increasingOrder[i]));
			closeSubset.clear();
		}
	}
	delete[] storage;
	std::cout << "Has taken closest and furthest out!" << std::endl;


/* Silhouette computation started */

	//groupNo record group numbers */

	if(groupNo<=1)
		return;

	/* if the dataset is not PBF, then should record distance matrix for Gamma matrix compution */
	if(!isPBF)
	{
		deleteDistanceMatrix(data.rows());

		if(!getDistanceMatrix(data, normOption, object))
		{
			std::cout << "Failure to compute distance matrix!" << std::endl;
		}
	}

	gettimeofday(&start, NULL);
	std::cout << "Final cluster has " << groupNo << " groups!" << std::endl;

	sil.computeValue(normOption,data,Row,Column,group,object,groupNo,isPBF);

	std::cout << "Silhouette computation completed!" << std::endl;

	gettimeofday(&end, NULL);
	delta = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;

	tr.eventList.push_back("Clustering evaluation computing takes: ");
	tr.timeList.push_back(to_string(delta)+"s");

	ValidityMeasurement vm;
	vm.computeValue(normOption, data, group, object, isPBF);
	tr.eventList.push_back("K-means Validity measure is: ");
	tr.timeList.push_back(to_string(vm.f_c));

	/* write value of the silhouette class */
	IOHandler::writeReadme(entropy, sil);
}
