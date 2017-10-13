
/*
 *  This file is part of birch-clustering-algorithm.
 *
 *  birch-clustering-algorithm is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  birch-clustering-algorithm is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with birch-clustering-algorithm.  If not, see <http://www.gnu.org/licenses/>.
 *
 *	Copyright (C) 2011 Taesik Yoon (otterrrr@gmail.com)
 */


/** Simple test code for birch-clustering algorithm
 *
 * BIRCH has 4 phases: building, compacting, clustering, redistribution.
 * 
 * building - building cftree inserting a new data-point
 * compacting - make cftree smaller enlarging the range of sub-clusters
 * clustering - clustering sub-clusters(summarized clusters) using the existing clustering algorithm
 * redistribution - labeling data-points to the closest center
 */

#include "ClusterAnalysis.h"

int main( int argc, char** argv)
{
	std::vector<std::vector<float> > trajectories;
	Eigen::MatrixXf equalArray;
	std::vector<item_type<750u> > items;
	int dimension, maxGroup, normOption;
	FileIndex fi;
	std::vector<int> item_cids;
	string fullName;
	MetricPreparation object;

	getUserInput(argc, argv, trajectories, equalArray, items, dimension, fi);

	getBirchClustering(items,argv,trajectories,fi,equalArray,
					   dimension, item_cids, maxGroup,normOption,
					   fullName, object);

	getClusterAnalysis(trajectories,fi,equalArray,dimension, 
					   item_cids, maxGroup, normOption, fullName,
					   object);
	return 0;
}
