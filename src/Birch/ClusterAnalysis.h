#ifndef _CLUSTERANALYSIS_H_
#define _CLUSTERANALYSIS_H_

#include "CFTree.h"
#include "IOHandler.h"
#include "Silhouette.h"
//#include "item_type.h"
#include <string>
#include <sstream>
#include <fstream>
#include <climits>
#include <string.h>
#include <time.h>
#include <sstream>


/* why use extern? Too many parameters passed in */
std::vector<string> activityList;
std::vector<double> timeList;

bool isPBF;

template<boost::uint32_t dim>
MetricPreparation CFTree<dim>::object = MetricPreparation();

template<boost::uint32_t dim>
int CFTree<dim>::normOption = -1;

template<boost::uint32_t dim>
int CFTree<dim>::totalNodes = 0;

typedef CFTree<1800u> cftree_type;

cftree_type::float_type birch_threshold;



struct FileIndex
{
	int vertexCount, maxElement;
	FileIndex()
	{}

	~FileIndex()
	{}
};


template<boost::uint32_t dim>
void getUserInput(const int& argc, 
				  char **argv,
				  std::vector<std::vector<float> >& trajectories,
				  Eigen::MatrixXf& equalArray,
				  std::vector<item_type<dim> >& items,
				  int& dimension,
				  FileIndex& fi)
{
	if( argc != 3 )
	{
		std::cout << "usage: birch (input-file) (dimension)" << std::endl;
		exit(1);
	}
	int samplingMethod;
	stringstream ss;
	ss << "../dataset/" << argv[1];


/* get the bool tag for isPBF */
	std::cout << "It is a PBF dataset? 1.Yes, 0.No" << std::endl;
	int PBFjudgement;
	std::cin >> PBFjudgement;
	assert(PBFjudgement==1||PBFjudgement==0);
	isPBF = (PBFjudgement==1);

	std::cout << "Please choose the sampling method? " << endl
	          << "1.filling, 2.uniform sampling." << std::endl;
	std::cin >> samplingMethod;
	assert(samplingMethod==1||samplingMethod==2);

	dimension = atoi(argv[2]);

	struct timeval start, end;
	double timeTemp;

	gettimeofday(&start, NULL);
	IOHandler::readFile(ss.str(), trajectories, fi.vertexCount, 
					    dimension, fi.maxElement);
	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec  - start.tv_sec) * 1000000u 
			   + end.tv_usec - start.tv_usec) / 1.e6;
	activityList.push_back("Read file takes: ");
	timeList.push_back(timeTemp);


	gettimeofday(&start, NULL);
	if(samplingMethod==1)
		IOHandler::expandArray(equalArray,trajectories,dimension,
							   fi.maxElement);
	else if(samplingMethod==2)
		IOHandler::sampleArray(equalArray,trajectories,dimension,
							   fi.maxElement);
	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec  - start.tv_sec) * 1000000u 
			   + end.tv_usec - start.tv_usec) / 1.e6;
	activityList.push_back("Pre-processing takes: ");
	timeList.push_back(timeTemp);


	load_items(equalArray, items);
	std::cout << items.size() << " items loaded" << std::endl;
}


template<typename T>
static void print_items( const std::string fname, T& items )
{
	struct _compare_item_id
	{
		bool operator()( const item_type<3>& lhs, const item_type<3>& rhs ) 
		const { return lhs.cid() < rhs.cid(); }
	};

	int maxGroup = INT_MIN;
	int belongGroup;
	for( std::size_t i = 0 ; i < items.size() ; i++ )
	{
//		for( std::size_t d = 0 ; d < cftree_type::fdim ; d++ )
//			fout << items[i].item[d] << " ";
		belongGroup = items[i].cid();
		if(belongGroup>maxGroup)
			maxGroup=belongGroup;
	}
}


template<boost::uint32_t dim>
static void load_items( const Eigen::MatrixXf& matrixData,
						std::vector<item_type<dim> >& items)
{
    items.resize(matrixData.rows());
#pragma omp parallel for schedule(dynamic) num_threads(8)
    for (int i = 0; i < items.size(); ++i)
    {
    	const Eigen::VectorXf& eachRow = matrixData.row(i);
    	cftree_type::item_vec_type item(eachRow.data(), eachRow.data()+eachRow.size());
		items[i] = &(item[0]);
    }
}


const float getMaxDist(const Eigen::MatrixXf& equalArray, 
					   const MetricPreparation& object, 
					   const int& normOption)
{
	const float& Percentage = 0.1;
	const int& chosen = int(Percentage*equalArray.rows());
	float result = -0.1;
#pragma omp parallel for reduction(max:result) num_threads(8)
	for (int i = 0; i < chosen; ++i)
	{
		for (int j = 0; j < equalArray.rows(); ++j)
		{
			if(i==j)
				continue;

			float dist;
			if(distanceMatrix)
				dist = distanceMatrix[i][j];
			else
				dist = getDisimilarity(equalArray.row(i),
						 equalArray.row(j),i,j,normOption,object);
			if(dist>result)
				result=dist;
		}
	}	
	return result;
}

template<boost::uint32_t dim>
void getBirchClusterTrial(const MetricPreparation& object,
						  const int& normOption,
						  std::vector<item_type<dim> >& items,
						  const float& distThreshold,
						  int& maxGroup,
						  std::vector<int>& item_cids)
{
	cftree_type tree(distThreshold, 0);
	tree.cftree_type::object = object;
	tree.cftree_type::normOption = normOption;
	tree.cftree_type::totalNodes = items.size();

	// phase 1 and 2: building, compacting when overflows memory limit
	for( std::size_t i = 0 ; i < items.size() ; i++ )
		tree.insert((float_type*)(&(items[i][0])));

	// phase 2 or 3: compacting? or clustering?
	// merging overlayed sub-clusters by rebuilding true
	std::cout << "Curve dimensionality is: " << cftree_type::fdim << std::endl;

	tree.rebuild(false);

	// phase 3: clustering sub-clusters using the existing clustering algorithm
	//cftree_type::cfentry_vec_type entries;
	std::vector<CFEntry<1800u> > entries;
	
	item_cids.clear();

	tree.cluster( entries );

	// phase 4: redistribution

	// @comment ts - it is also possible to another clustering algorithm hereafter
	//				for example, we have k initial points for k-means clustering algorithm
	//tree.redist_kmeans( items, entries, 0 );
    tree.redist(items.begin(), items.end(), entries, item_cids);
    maxGroup = INT_MIN;
    for (std::size_t i = 0; i < item_cids.size(); i++)
    {
    	int& itemCID = items[i].cid(); 
        itemCID = item_cids[i];
        if(maxGroup<itemCID)
        	maxGroup=itemCID;
    }
}
								


template<boost::uint32_t dim>
void getBirchClustering(std::vector<item_type<dim> >& items,
						char **argv,
						std::vector<std::vector<float> >& trajectories,
						const FileIndex& fi,
						const Eigen::MatrixXf& equalArray,
						const int& dimension,
						std::vector<int>& item_cids,
						int& maxGroup,
						int& normOption,
						string& fullName,
						MetricPreparation& object)
{

	std::cout << std::endl;
	std::cout << "Choose a norm from 0-11!" << std::endl;
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

	object = MetricPreparation(equalArray.rows(), equalArray.cols());
	object.preprocessing(equalArray, equalArray.rows(), equalArray.cols(), normOption);

	/* if the dataset is not PBF, then should record distance matrix for Gamma matrix compution */
	if(!isPBF)
	{
		deleteDistanceMatrix(equalArray.rows());

		if(!getDistanceMatrix(equalArray, normOption, object))
		{
			std::cout << "Failure to compute distance matrix!" << std::endl;
		}
	}

	const float distThreshold = getMaxDist(equalArray, object, normOption);

	std::cout << "Enter approximate number of clusters: " << std::endl;
	int requiredClusters;
	std::cin >> requiredClusters;
	const int& upperClusters = requiredClusters*1.2;
	const int& lowerClusters = requiredClusters*0.8;

	std::cout << "Sampled max distance is: " << distThreshold << std::endl;

	float left = 0, right = 0.5, middle;

	struct timeval start, end;
	double timeTemp;
	gettimeofday(&start, NULL);

	int iteration = 0;
	while(true&&iteration<15)
	{
		std::cout << "Iteration for birch clustering: " << ++iteration
				  << std::endl;
		middle = (left+right)/2.0;
		std::cout << "Weight is " << middle << std::endl;
		getBirchClusterTrial(object,normOption,items,middle*distThreshold,
							 maxGroup,item_cids);
		std::cout << maxGroup << std::endl;
		if(maxGroup<=upperClusters && maxGroup>=lowerClusters)
			break;
		else if(maxGroup>upperClusters)
			left = middle;
		else if(maxGroup<lowerClusters)
			right = middle;
	}
	birch_threshold = middle;

	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec  - start.tv_sec) * 1000000u 
			   + end.tv_usec - start.tv_usec) / 1.e6;
	activityList.push_back("Birch search required numbers takes: ");
	timeList.push_back(timeTemp);

    std::cout << "Max group is: " << maxGroup << std::endl;
    //print_items(argc >= 4 ? ss.str().c_str() : "item_cid.txt", items);
    stringstream ss;
    ss << "../dataset/" << argv[1] << "_full.vtk";
    fullName = ss.str();

    IOHandler::printVTK(fullName, trajectories, fi.vertexCount, dimension); 

}


void getClusterAnalysis(const vector<vector<float> >& trajectories,
						const FileIndex& fi,
						const MatrixXf& equalArray,
						const int& dimension, 
					    vector<int>& item_cids, 
					    const int& maxGroup, 
					    const int& normOption,
					    const string& fullName,
					    const MetricPreparation& object)
{
	int numClusters = maxGroup+1;
	std::vector<int> container(numClusters,0);
	for (int i = 0; i < item_cids.size(); ++i)
		++container[item_cids[i]];

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

	numClusters = groupNo;

/* compute balanced Entropy value for the clustering algorithm */
	const int& Row = equalArray.rows();
	float entropy = 0.0, probability;
	for(int i=0;i<container.size();++i)
	{
		probability = float(container[i])/float(Row);
		entropy+=probability*log2f(probability);
	}
	entropy = -entropy/log2f(numClusters);



#pragma omp parallel for schedule(dynamic) num_threads(8)
	for (int i = 0; i < item_cids.size(); ++i)
		item_cids[i]=increasingOrder[item_cids[i]];

	std::vector<std::vector<int> > storage(numClusters);
	for (int i = 0; i < item_cids.size(); ++i)
		storage[item_cids[i]].push_back(i);

	/* record labeling information */
	IOHandler::generateGroups(storage);


	IOHandler::printClusters(trajectories,item_cids,container, 
		 "norm"+to_string(normOption), fullName,dimension);

	struct timeval start, end;
	double timeTemp;

	gettimeofday(&start, NULL);
	Silhouette sil;
	sil.computeValue(normOption,equalArray,equalArray.rows(),
		equalArray.cols(),item_cids,object,numClusters,isPBF);
	gettimeofday(&end, NULL);
	timeTemp = ((end.tv_sec  - start.tv_sec) * 1000000u 
			   + end.tv_usec - start.tv_usec) / 1.e6;
	activityList.push_back("Silhouette calculation takes: ");
	timeList.push_back(timeTemp);


	/* compute the centroid coordinates of each clustered group */
	Eigen::MatrixXf centroid = MatrixXf::Zero(numClusters,equalArray.cols());
	vector<vector<float> > cenVec(numClusters);
#pragma omp parallel for schedule(dynamic) num_threads(8)
	for (int i=0;i<numClusters;++i)
	{
		const std::vector<int>& groupRow = storage[i];
		for (int j = 0; j < groupRow.size(); ++j)
		{
			centroid.row(i)+=equalArray.row(groupRow[j]);
		}		
		centroid.row(i)/=groupRow.size();
		const Eigen::VectorXf& vec = centroid.row(i);
		cenVec[i] = vector<float>(vec.data(), vec.data()+equalArray.cols());
	}

	vector<vector<float> > closest(numClusters);
	vector<vector<float> > furthest(numClusters);

#pragma omp parallel for schedule(dynamic) num_threads(8)
	for (int i=0;i<numClusters;++i)
	{
		float minDist = FLT_MAX;
		float maxDist = -10;
		int minIndex = -1, maxIndex = -1;
		const std::vector<int>& groupRow = storage[i];
		const Eigen::VectorXf& eachCentroid = centroid.row(i);
		for (int j = 0; j < groupRow.size(); ++j)
		{
			float distance = getDisimilarity(eachCentroid,equalArray,groupRow[j],normOption,object);
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
		closest[i] = trajectories[minIndex];
		furthest[i] = trajectories[maxIndex];
	}

	std::cout << "Finishing extracting features!" << std::endl;	
	IOHandler::printFeature("norm"+to_string(normOption)+"_closest.vtk", 
			closest, sil.sCluster, dimension);
	IOHandler::printFeature("norm"+to_string(normOption)+"_furthest.vtk",
			furthest, sil.sCluster, dimension);
	IOHandler::printFeature("norm"+to_string(normOption)+"_centroid.vtk", 
			cenVec, sil.sCluster,dimension);

	IOHandler::printToFull(trajectories, sil.sData, 
			"norm"+to_string(normOption)+"_SValueLine", fullName, dimension);
	IOHandler::printToFull(trajectories, item_cids, sil.sCluster, 
		      "norm"+to_string(normOption)+"_SValueCluster", fullName, dimension);

	IOHandler::generateReadme(activityList,timeList,normOption,
				numClusters,sil.sAverage,birch_threshold);

/* print entropy value for the clustering algorithm */
	IOHandler::writeReadme(entropy,sil);

/* measure closest and furthest rotation */
	std::vector<float> closestRot, furthestRot;
	const float& closestAverage = getRotation(closest, closestRot);
	const float& furthestAverage = getRotation(furthest, furthestRot);

	IOHandler::writeReadme(closestAverage, furthestAverage);

}


static float randf()
{
	return float(rand()/(double)RAND_MAX);
}


#endif
