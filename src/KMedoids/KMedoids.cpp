#include "KMedoids.h"


extern bool isPBF;


KMedoids::KMedoids(const Parameter& pm, 
				   const Eigen::MatrixXf& data,
				   const int& numOfClusters)
		 :initialStates(pm.initialization), isSample(pm.isSample),
		  data(data), numOfClusters(numOfClusters)
{

}
	

KMedoids::~KMedoids()
{

}


void KMedoids::getMedoids(FeatureLine& fline,
						  const int& normOption,
						  Silhouette& sil,
						  EvaluationMeasure& measure,
						  TimeRecorder& tr) const
{
	MetricPreparation object(data.rows(), data.cols());
	object.preprocessing(data, data.rows(), data.cols(), normOption);

	MatrixXf clusterCenter;
	getInitCenter(clusterCenter, object, normOption);

	float moving=1000, tempMoving,/* dist, tempDist, */before;
	int *storage = new int[numOfClusters]; // used to store number inside each cluster
	MatrixXf centerTemp;
	int tag = 0;
	std::vector< std::vector<int> > neighborVec(numOfClusters, 
												std::vector<int>());

/* perform K-means with different metrics */
	std::cout << "K-medoids start!" << std::endl;
	const int& Row = data.rows();
	const int& Column = data.cols();	
	struct timeval start, end;
	gettimeofday(&start, NULL);
	std::vector<int> recorder(Row); //use to record which cluster the row belongs to

	do
	{
	/* reset storage number and weighted mean inside each cluster*/
		before=moving;
		memset(storage,0,sizeof(int)*numOfClusters);
		centerTemp = clusterCenter;

	/* clear streamline indices for each cluster */
	#pragma omp parallel for schedule(static) num_threads(8)
		for (int i = 0; i < numOfClusters; ++i)
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
				for (int j = 0; j < numOfClusters; ++j)
				{
					tempDist = getDisimilarity(clusterCenter.row(j),
								data,i,normOption,object);
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
				}
			}
		}

		computeMedoids(centerTemp, neighborVec, normOption, object);

		moving = FLT_MIN;

	/* measure how much the current center moves from original center */	
	#pragma omp parallel for reduction(max:moving) num_threads(8)
		for (int i = 0; i < numOfClusters; ++i)
		{
			if(storage[i]>0)
			{
				tempMoving = (centerTemp.row(i)-clusterCenter.row(i)).norm();
				clusterCenter.row(i) = centerTemp.row(i);
				if(moving<tempMoving)
					moving = tempMoving;
			}
		}
		std::cout << "K-means iteration " << ++tag << " completed, and moving is " << moving 
				  << "!" << std::endl;
	}while(abs(moving-before)/before >= 1.0e-2 && tag < 20/* && moving > 5.0*/);
	
	gettimeofday(&end, NULL);
	double delta = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;

	tr.eventList.push_back("K-medoids iteration takes: ");
	tr.timeList.push_back(to_string(delta)+"s");

	std::multimap<int,int> groupMap;
	float entropy = 0.0;
	float probability;
	vector<int> increasingOrder(numOfClusters);
	for (int i = 0; i < numOfClusters; ++i)
	{
		groupMap.insert(std::pair<int,int>(storage[i],i));
		if(storage[i]>0)
		{
			probability = float(storage[i])/float(Row);
			entropy += probability*log2f(probability);
		}
	}

	int groupNo = 0;
	for (std::multimap<int,int>::iterator it = groupMap.begin(); it != groupMap.end(); ++it)
	{
		if(it->first>0)
		{
			increasingOrder[it->second] = (groupNo++);
		}
	}

	entropy = -entropy/log2f(groupNo);
	/* finish tagging for each group */

	/* record labeling information */
	IOHandler::generateGroups(neighborVec);


	// set cluster group number and size number 
#pragma omp parallel for schedule(static) num_threads(8)
	for (int i = 0; i < Row; ++i)
	{
		fline.group[i] = increasingOrder[recorder[i]];
		fline.totalNum[i] = storage[recorder[i]];
	}

	float shortest, toCenter, farDist;
	int shortestIndex = 0, tempIndex = 0, furthestIndex = 0;
	std::vector<int> neighborTemp;

	/* choose cloest and furthest streamlines to centroid streamlines */
	for (int i = 0; i < numOfClusters; ++i)
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
				toCenter = getDisimilarity(clusterCenter.row(i),data,
							tempIndex,normOption,object);
				if(!isSample)
				{
					if(toCenter<shortest)
					{
						shortest = toCenter;
						shortestIndex = tempIndex;
					}
				}
				/* update the farthest index to centroid */
				if(toCenter>farDist)
				{
					farDist = toCenter;
					furthestIndex = tempIndex;
				}
			}
			if(!isSample)
				fline.closest.push_back(ExtractedLine(shortestIndex,increasingOrder[i]));
			fline.furthest.push_back(ExtractedLine(furthestIndex,increasingOrder[i]));
		}
	}

	std::vector<float> closeSubset;
	/* based on known cluster centroid, save them as vector for output */
	for (int i = 0; i < numOfClusters; ++i)
	{
		if(storage[i]>0)
		{
			for (int j = 0; j < Column; ++j)
			{
				closeSubset.push_back(clusterCenter(i,j));
			}
			fline.centerMass.push_back(MeanLine(closeSubset,increasingOrder[i]));
			closeSubset.clear();
		}
	}
	delete[] storage;
	std::cout << "Has taken closest and furthest out!" << std::endl;


/* Silhouette computation started */

	std::cout << "The finalized cluster size is: " << groupNo << std::endl;
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

	//groupNo record group numbers */
	gettimeofday(&start, NULL);

	sil.computeValue(normOption,data,Row,Column,fline.group,object,groupNo,isPBF);
	std::cout << "Silhouette computation completed!" << std::endl;

	gettimeofday(&end, NULL);
	delta = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;

	tr.eventList.push_back("Evaluation analysis would take: ");
	tr.timeList.push_back(to_string(delta)+"s");

	ValidityMeasurement vm;
	vm.computeValue(normOption, data, fline.group, object, isPBF);
	tr.eventList.push_back("K-medoids Validity measure is: ");
	tr.timeList.push_back(to_string(vm.f_c));

	/* store the evaluation value result */
	measure.silVec.push_back(sil.sAverage);
	measure.gammaVec.push_back(sil.gammaStatistic);
	measure.entropyVec.push_back(entropy);
	measure.dbIndexVec.push_back(sil.dbIndex);

}


void KMedoids::getInitCenter(MatrixXf& initialCenter,
							 const MetricPreparation& object,
							 const int& normOption) const
{
	switch(initialStates)
	{
	case 1:
		Initialization::generateRandomPos(initialCenter, 
				data.cols(), data, numOfClusters);
		break;

	default:
	case 2:
		Initialization::generateFromSamples(initialCenter, 
				data.cols(), data, numOfClusters);
		break;

	case 3:
		Initialization::generateFarSamples(initialCenter, 
				data.cols(), data, numOfClusters, normOption, object);
		break;
	}
	std::cout << "Initialization completed!" << std::endl;
}


void KMedoids::computeMedoids(MatrixXf& centerTemp, 
							  const vector<vector<int> >& neighborVec, 
							  const int& normOption, 
							  const MetricPreparation& object) const
{
	centerTemp = MatrixXf(numOfClusters,data.cols());
	if(isSample)//centroid is from samples with minimal L1 summation
				//use Voronoi iteration https://en.wikipedia.org/wiki/K-medoids
	{
	#pragma omp parallel for schedule(static) num_threads(8)
		for(int i=0;i<neighborVec.size();++i)
		{
			const vector<int>& clusMember = neighborVec[i];
			const int& clusSize = clusMember.size();
			MatrixXf mutualDist = MatrixXf::Zero(clusSize, clusSize);
			/*mutualDist to store mutual distance among lines of each cluster */
			for(int j=0;j<clusSize;++j)
			{
				for(int k=j+1;k<clusSize;++k)
				{
					mutualDist(j,k) = getDisimilarity(data,clusMember[j],
						clusMember[k],normOption,object);
					mutualDist(k,j) = mutualDist(j,k);
				}
			}

			float minL1_norm = FLT_MAX, rowSummation;
			int index = -1;
			for(int j=0;j<clusSize;++j)
			{
				rowSummation = mutualDist.row(j).sum();
				if(rowSummation<minL1_norm)
				{
					minL1_norm = rowSummation;
					index = j;
				}
			}
			centerTemp.row(i)=data.row(clusMember[index]);
		}
	}

	else//use Weiszfeld's algorithm to get geometric median
		//reference at https://en.wikipedia.org/wiki/Geometric_median
	{
		MatrixXf originCenter = centerTemp;
	#pragma omp parallel for schedule(static) num_threads(8)
		for(int i=0;i<numOfClusters;++i)
		{
			const vector<int>& clusMember = neighborVec[i];
			const int& clusSize = clusMember.size();
			float distToCenter, distInverse, percentage = 1.0;
			int tag = 0;
			while(tag<=10&&percentage>=0.02)
			{
				VectorXf numerator = VectorXf::Zero(data.cols());
				VectorXf previous = centerTemp.row(i);
				float denominator = 0;
				for(int j=0;j<clusSize;++j)
				{
					distToCenter = getDisimilarity(centerTemp.row(i),
							data,clusMember[j],normOption,object);
					distInverse = (distToCenter>1.0e-8)?1.0/distToCenter:1.0e8;
					numerator += data.row(clusMember[j])*distInverse;
					denominator += distInverse;
				}
				centerTemp.row(i) = numerator/denominator;
				percentage = (centerTemp.row(i)-previous).norm()/previous.norm();
				tag++;
			}
		}
	}
}

