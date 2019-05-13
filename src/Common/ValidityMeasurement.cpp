/*
 * ValidityMeasurement.cpp
 *
 *  Created on: Jun 24, 2018
 *      Author: lieyu
 */

#include "ValidityMeasurement.h"

ValidityMeasurement::ValidityMeasurement() {
	// TODO Auto-generated constructor stub

}

ValidityMeasurement::~ValidityMeasurement() {
	// TODO Auto-generated destructor stub
}


// function API for computing the validity measurement for general cases
void ValidityMeasurement::computeValue(const int& normOption, const MatrixXf& array,
		const std::vector<int>& group, const MetricPreparation& object, const bool& isPBF)
{
	std::cout << "Compute validity measurement..." << std::endl;
	// get how many different groups it totally has
	int max_group = -1;
	const int& num_node = group.size();
	for(int i=0; i<num_node; ++i)
	{
		if(group[i]==-1)
			continue;
		max_group = std::max(group[i], max_group);
	}
	max_group+=1;

	std::vector<std::vector<int> > storage(max_group);

	for(int i=0; i<num_node; ++i)
	{
		if(group[i]==-1)
			continue;
		storage[group[i]].push_back(i);
	}

	std::vector<std::tuple<float, float, float> > measureVec(max_group);

	for(int i=0; i<max_group; ++i)
	{
		getMST_Parent_Node(measureVec[i], storage[i], object, normOption, array, isPBF);
	}

	float minSc = 0, maxSc = 0, aver_sigma = 0, std_sigma = 0, std_variance;
	for(int i=0; i<max_group; ++i)
	{
		// get the min Sc by summation
		minSc+=std::get<1>(measureVec[i]);
		// get the max Sc by summation
		maxSc+=std::get<2>(measureVec[i]);
		std_variance = std::get<0>(measureVec[i]);
		// get the average variance and standard variation of variance
		aver_sigma+=std_variance;
		std_sigma+=std_variance*std_variance;
	}
	aver_sigma/=max_group;
	std_sigma = std_sigma/float(max_group-1)-float(max_group)/float(max_group-1)*aver_sigma*aver_sigma;

	if(std_sigma<1.0E-10)
	{
		std_sigma=1.0E-10;
	}
	std_sigma=sqrt(std_sigma);

	float h_DDc = aver_sigma+std_sigma;

	minSc/=float(max_group);
	maxSc/=float(max_group);

	// compute g1_Sc
	float g1_Sc = (1.0-minSc)*(1.0-maxSc);
	if(g1_Sc<0)
	{
		std::cout << "Negative number for g1_Sc computation!" << std::endl;
		exit(1);
	}
	g1_Sc = aver_sigma*sqrt(g1_Sc);

	// compute g2_Sc
	float g2_Sc = minSc*maxSc;
	if(g2_Sc<0)
	{
		std::cout << "Negative number for g2_Sc computation!" << std::endl;
		exit(1);
	}
	g2_Sc = aver_sigma/sqrt(g2_Sc);

	// compute g_Sc
	float g_Sc = (sqrt(g1_Sc*g2_Sc)+(g1_Sc+g2_Sc)/2.0)/2.0;

	// compoute f_c
	f_c = h_DDc*g_Sc;

	if(isnan(f_c) || isinf(f_c))
	{
		std::cout << "Error for f_c to have inf or nan values!" << std::endl;
		exit(1);
	}

	/* normalization of validity measurement */
	float min_dist = FLT_MAX, max_dist = -1.0;
	const int& row = array.rows();
#pragma omp parallel for reduction(min:min_dist) num_threads(8)
	for(int i=0; i<row; ++i)
	{
		for(int j=0; j<row; ++j)
		{
			if(i==j)
				continue;
			float dist;
			if(distanceMatrix)
				dist = distanceMatrix[i][j];
			else
				dist = getDisimilarity(array, i, j, normOption, object);
			min_dist = std::min(min_dist, dist);
		}
	}

#pragma omp parallel for reduction(max:max_dist) num_threads(8)
	for(int i=0; i<row; ++i)
	{
		for(int j=0; j<row; ++j)
		{
			if(i==j)
				continue;
			float dist;
			if(distanceMatrix)
				dist = distanceMatrix[i][j];
			else
				dist = getDisimilarity(array, i, j, normOption, object);
			max_dist = std::max(max_dist, dist);
		}
	}
	std::cout << "min dist is " << min_dist << ", and max is " << max_dist << std::endl;
	f_c/=(max_dist-min_dist)/**(max_dist-min_dist)*/;

	assert(!isnan(f_c) && !isinf(f_c));
	std::cout << "Validity measurement is " << f_c << std::endl;
}


// function API for computing the validity measurement for PCA case only
void ValidityMeasurement::computeValue(const MatrixXf& array, const std::vector<int>& group)
{
	std::cout << "Compute validity measurement..." << std::endl;
	// get how many different groups it totally has
	int max_group = -1;
	const int& num_node = group.size();
	for(int i=0; i<num_node; ++i)
	{
		if(group[i]==-1)
			continue;
		max_group = std::max(group[i], max_group);
	}
	max_group+=1;

	std::vector<std::vector<int> > storage(max_group);
	for(int i=0; i<num_node; ++i)
	{
		if(group[i]==-1)
			continue;
		storage[group[i]].push_back(i);
	}

	std::vector<std::tuple<float, float, float> > measureVec(max_group);

	for(int i=0; i<max_group; ++i)
	{
		getMST_Parent_Node(measureVec[i], storage[i], array);
	}

	float minSc = 0, maxSc = 0, aver_sigma = 0, std_sigma = 0, std_variance;
	for(int i=0; i<max_group; ++i)
	{
		// get the min Sc by summation
		minSc+=std::get<1>(measureVec[i]);
		// get the max Sc by summation
		maxSc+=std::get<2>(measureVec[i]);
		std_variance = std::get<0>(measureVec[i]);
		// get the average variance and standard variation of variance
		aver_sigma+=std_variance;
		std_sigma+=std_variance*std_variance;
	}
	aver_sigma/=float(max_group);
	std_sigma = std_sigma/float(max_group-1)-float(max_group)/float(max_group-1)*aver_sigma*aver_sigma;
	if(std_sigma<1.0E-10)
	{
		std_sigma = 1.0E-10;
	}
	else
		std_sigma=sqrt(std_sigma);

	float h_DDc = aver_sigma+std_sigma;

	minSc/=float(max_group);
	maxSc/=float(max_group);

	// compute g1_Sc
	float g1_Sc = (1.0-minSc)*(1.0-maxSc);
	if(g1_Sc<0)
	{
		std::cout << "Negative number for g1_Sc computation!" << std::endl;
		exit(1);
	}
	g1_Sc = aver_sigma*sqrt(g1_Sc);

	// compute g2_Sc
	float g2_Sc = minSc*maxSc;
	if(g2_Sc<0)
	{
		std::cout << "Negative number for g2_Sc computation!" << std::endl;
		exit(1);
	}
	g2_Sc = aver_sigma/sqrt(g2_Sc);

	// compute g_Sc
	float g_Sc = (sqrt(g1_Sc*g2_Sc)+(g1_Sc+g2_Sc)/2.0)/2.0;

	// compoute f_c
	f_c = h_DDc*g_Sc;

	if(isnan(f_c) || isinf(f_c))
	{
		std::cout << "Error for f_c to have inf or nan values!" << std::endl;
		exit(1);
	}

	/* normalization of validity measurement */
	float min_dist = FLT_MAX, max_dist = -1.0;
	const int& row = array.rows();
#pragma omp parallel for reduction(min:min_dist) num_threads(8)
	for(int i=0; i<row; ++i)
	{
		for(int j=0; j<row; ++j)
		{
			if(i==j)
				continue;
			min_dist = std::min(min_dist, (array.row(i)-array.row(j)).norm());
		}
	}

#pragma omp parallel for reduction(max:max_dist) num_threads(8)
	for(int i=0; i<row; ++i)
	{
		for(int j=0; j<row; ++j)
		{
			if(i==j)
				continue;
			max_dist = std::max(max_dist, (array.row(i)-array.row(j)).norm());
		}
	}
	std::cout << "min dist is " << min_dist << ", and max is " << max_dist << std::endl;
	f_c/=(max_dist-min_dist)/**(max_dist-min_dist)*/;

	assert(!isnan(f_c) && !isinf(f_c));
	std::cout << "Validity measurement is " << f_c << std::endl;
}



// get MST for each cluster given index and pair-wise distance for general cases
void ValidityMeasurement::getMST_Parent_Node(std::tuple<float, float, float>& values,
			const std::vector<int>& clusterNode, const MetricPreparation& object, const int& normOption,
			const MatrixXf& array, const bool& isPBF)
{
	using namespace boost;
	typedef adjacency_list < vecS, vecS, undirectedS, no_property, property < edge_weight_t, float > > Graph;
	typedef graph_traits < Graph >::edge_descriptor Edge;
	typedef graph_traits < Graph >::vertex_descriptor Vertex;
	typedef std::pair<int, int> E;

	const int& num_nodes = clusterNode.size();

	if(num_nodes<=1)
	{
		values = std::make_tuple(0.0,0.0,0.0);
		std::cout << "Find 1-candidate cluster!" << std::endl;
		return;
	}

	const int num_edges = num_nodes*(num_nodes-1)/2;
	Eigen::MatrixXf distM;
	// if distanceMatrix is not stored ahead of time
	if(isPBF)
		distM = Eigen::MatrixXf(num_nodes, num_nodes);
	// assign [source, destination] index pair and weight lists
	E *edge_array = new E[num_edges];
	float *weights = new float[num_edges], dist;
	int temp = 0;
	for(int i=0; i<num_nodes-1; ++i)
	{
		for(int j=i+1; j<num_nodes; ++j)
		{
			// assign index pair
			edge_array[temp] = std::make_pair(i, j);
			// assign weight list
			if(distanceMatrix)
				dist = distanceMatrix[clusterNode[i]][clusterNode[j]];
			else
				dist = getDisimilarity(array, clusterNode[i], clusterNode[j], normOption, object);

			if(isPBF)
			{
				distM(i,j) = dist;
				distM(j,i) = dist;
			}

			weights[temp] = dist;
			++temp;
		}
	}

#if defined(BOOST_MSVC) && BOOST_MSVC <= 1300
	Graph g(num_nodes);
	property_map<Graph, edge_weight_t>::type weightmap = get(edge_weight, g);
	for (std::size_t j = 0; j < num_edges; ++j) {
		Edge e; bool inserted;
		tie(e, inserted) = add_edge(edge_array[j].first, edge_array[j].second, g);
		weightmap[e] = weights[j];
	}
#else
	Graph g(edge_array, edge_array + num_edges, weights, num_nodes);
#endif
	property_map < Graph, edge_weight_t >::type weight = get(edge_weight, g);
	std::vector < Edge > spanning_tree;
	kruskal_minimum_spanning_tree(g, std::back_inserter(spanning_tree));

	if(edge_array!=NULL)
	{
		delete[] edge_array;
		edge_array = NULL;
	}

	if(weights!=NULL)
	{
		delete[] weights;
		weights = NULL;
	}

	// compute the standard deviation for the distance in MST
	double summation = 0.0, sq_summation = 0.0, average_mst_d, max_d_mst = -1.0;
	const int& MST_EDGE_NUM = num_nodes-1;

	for (std::vector < Edge >::iterator ei = spanning_tree.begin(); ei != spanning_tree.end(); ++ei)
	{
			dist = weight[*ei];
			max_d_mst = std::max(double(dist), max_d_mst);
			summation+=dist;
			sq_summation+=dist*dist;
	}

	float variance;

	if(MST_EDGE_NUM<=1)
	{
		variance = 0;
		average_mst_d = summation;
	}
	else
	{
		average_mst_d = summation/float(MST_EDGE_NUM);
		variance = sq_summation/float(MST_EDGE_NUM-1)-average_mst_d*summation/float(MST_EDGE_NUM-1);

		if(variance<1.0E-10)
		{
			variance = 1.0E-10;
		}
		variance = sqrt(variance);
	}

	// compute the inner Sc value for this cluster
	int min_index, max_index;
	float min_Sc = get_Sc_by_range(isPBF, distM, clusterNode, max_d_mst, object, normOption, array, min_index);
	float max_Sc = get_Sc_by_range(isPBF, distM, clusterNode, average_mst_d, object, normOption, array, max_index);

	min_Sc/=float(min_index);
	max_Sc/=float(max_index);

	// store the standard deviation, min Sc and max Sc in the tuple
	values = std::make_tuple(variance, min_Sc, max_Sc);
}


// get MST for each cluster given index and pair-wise distance, PCA case only
void ValidityMeasurement::getMST_Parent_Node(std::tuple<float, float, float>& values,
			const std::vector<int>& clusterNode, const MatrixXf& array)
{
	using namespace boost;
	typedef adjacency_list < vecS, vecS, undirectedS, no_property, property < edge_weight_t, float > > Graph;
	typedef graph_traits < Graph >::edge_descriptor Edge;
	typedef graph_traits < Graph >::vertex_descriptor Vertex;
	typedef std::pair<int, int> E;
	// get number of points in one cluster
	const int& num_nodes = clusterNode.size();
	if(num_nodes<=1)
	{
		values = std::make_tuple(0.0,0.0,0.0);
		return;
	}
	const int num_edges = num_nodes*(num_nodes-1)/2;

	Eigen::MatrixXf distM = Eigen::MatrixXf(num_nodes, num_nodes);
	// assign [source, destination] index pair and weight lists
	E *edge_array = new E[num_edges];
	float *weights = new float[num_edges], dist;
	int temp = 0;
	for(int i=0; i<num_nodes-1; ++i)
	{
		for(int j=i+1; j<num_nodes; ++j)
		{
			// assign index pair
			edge_array[temp] = std::make_pair(i, j);

			dist = (array.row(clusterNode[i])-array.row(clusterNode[j])).norm();
			distM(i,j) = dist;
			distM(j,i) = dist;

			weights[temp] = dist;
			++temp;
		}
	}


#if defined(BOOST_MSVC) && BOOST_MSVC <= 1300
	Graph g(num_nodes);
	property_map<Graph, edge_weight_t>::type weightmap = get(edge_weight, g);
	for (std::size_t j = 0; j < num_edges; ++j) {
		Edge e; bool inserted;
		tie(e, inserted) = add_edge(edge_array[j].first, edge_array[j].second, g);
		weightmap[e] = weights[j];
	}
#else
	Graph g(edge_array, edge_array + num_edges, weights, num_nodes);
#endif
	property_map < Graph, edge_weight_t >::type weight = get(edge_weight, g);
	std::vector < Edge > spanning_tree;

	kruskal_minimum_spanning_tree(g, std::back_inserter(spanning_tree));


	if(edge_array!=NULL)
	{
		delete[] edge_array;
		edge_array = NULL;
	}

	if(weights!=NULL)
	{
		delete[] weights;
		weights = NULL;
	}

	// compute the standard deviation for the distance in MST
	double summation = 0.0, sq_summation = 0.0, average_mst_d, max_d_mst = -1.0;
	const int& MST_EDGE_NUM = num_nodes-1;

	for (std::vector < Edge >::iterator ei = spanning_tree.begin(); ei != spanning_tree.end(); ++ei)
	{
		dist = weight[*ei];
		summation+=dist;
		sq_summation+=dist*dist;
		max_d_mst=std::max(max_d_mst, double(dist));
	}

	float variance;
	if(MST_EDGE_NUM<=1)
	{
		variance = 0;
		average_mst_d = summation;
	}
	else
	{
		average_mst_d=summation/float(MST_EDGE_NUM);
		variance = sq_summation/float(MST_EDGE_NUM-1)-average_mst_d*summation/float(MST_EDGE_NUM-1);

		if(variance<1.0E-10)
		{
			variance = 1.0E-10;
		}
		variance = sqrt(variance);
	}

	// compute the inner Sc value for this cluster
	int min_index, max_index;
	float min_Sc = get_Sc_by_range(distM, clusterNode, max_d_mst, array, min_index);
	float max_Sc = get_Sc_by_range(distM, clusterNode, average_mst_d, array, max_index);

	min_Sc/=float(min_index);
	max_Sc/=float(max_index);

	// store the standard deviation, min Sc and max Sc in the tuple
	values = std::make_tuple(variance, min_Sc, max_Sc);
}



// compute the Sc by input range value for general cases
const float ValidityMeasurement::get_Sc_by_range(const bool& isPBF, const Eigen::MatrixXf& distM,
		                    const std::vector<int>& clusterNode, const float& rangeValue,
							const MetricPreparation& object, const int& normOption, const MatrixXf& array,
							int& index)
{
	const int& node_number = clusterNode.size();
	float result = 0.0;

	index = 0;
	int inside_whole, inside_cluster;
	for(int i=0; i<node_number; ++i)
	{
		inside_whole = 0, inside_cluster = 0;
		// count how many points in N_epsi(P_i) for the whole dataset
	#pragma omp parallel num_threads(8)
		{
		#pragma omp for nowait
			for(int j=0; j<array.rows(); ++j)
			{
				// don't want to handle duplicates and itself
				if(clusterNode[i]==j)
					continue;
				float dist;
				if(distanceMatrix)
					dist = distanceMatrix[clusterNode[i]][j];
				else
					dist = getDisimilarity(array, clusterNode[i], j, normOption, object);

			#pragma omp critical
				{
					if(dist<=rangeValue)
					{
						++inside_whole;
					}
				}
			}

		}

		// count how many points in N_epsi(P_i) for current cluster
	#pragma omp parallel num_threads(8)
		{
		#pragma omp for nowait
			for(int j=0; j<node_number; ++j)
			{
				// don't want to handle duplicates and itself
				if(i==j)
					continue;
				float dist;
				if(isPBF)
					dist = distM(i,j);
				else
					dist = distanceMatrix[clusterNode[i]][clusterNode[j]];

			#pragma omp critical
				if(dist<=rangeValue)
					++inside_cluster;
			}

		}
		assert(inside_cluster<=inside_whole);
		if(inside_whole==0)
			continue;
		++index;
		result+=float(inside_cluster)/float(inside_whole);
	}
	return result;
}



// compute the Sc by input range value for PCA case only
const float ValidityMeasurement::get_Sc_by_range(const Eigen::MatrixXf& distM, const std::vector<int>& clusterNode,
												 const float& rangeValue, const MatrixXf& array, int& index)
{
	const int& node_number = clusterNode.size();
	float result = 0.0;

	index = 0;
	int inside_whole, inside_cluster;
	for(int i=0; i<node_number; ++i)
	{
		inside_whole = 0, inside_cluster = 0;
		// count how many points in N_epsi(P_i) for the whole dataset
	#pragma omp parallel num_threads(8)
		{
		#pragma omp for nowait
			for(int j=0; j<array.rows(); ++j)
			{
				// don't want to handle duplicates and itself
				if(i==j)
					continue;
				float dist = (array.row(clusterNode[i])-array.row(j)).norm();

			#pragma omp critical
				if(dist<=rangeValue)
					++inside_whole;
			}

		}

		// count how many points in N_epsi(P_i) for current cluster
	#pragma omp parallel num_threads(8)
		{
		#pragma omp for nowait
			for(int j=0; j<node_number; ++j)
			{
				// don't want to handle duplicates and itself
				if(i==j)
					continue;
				float dist = distM(i,j);

			#pragma omp critical
				if(dist<=rangeValue)
					++inside_cluster;
			}

		}
		assert(inside_cluster<=inside_whole);
		if(inside_whole==0)
			continue;
		result+=float(inside_cluster)/float(inside_whole);
		++index;
		assert(!std::isnan(result));
	}
	return result;
}

