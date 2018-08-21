#ifndef _PREDEFINED_H
#define _PREDEFINED_H

#include "IOHandler.h"
#include "Initialization.h"
#include "Silhouette.h"


// define a treeNode structure to store AHC clustering tree
struct AHC_node
{
	int index = -1;

	/* to alleviate the computational cost to traverse all node elements */
	std::vector<int> element;

	AHC_node(const int& index): index(index)
	{}

	AHC_node()
	{}
};

// remove two elements in template vector
template <class T>
void deleteVecElements(std::vector<T>& origine, const T& first, const T& second);


/* we will use a min-heap to perserve sorted distance for hirarchical clustering */
struct DistNode
{
	int first = -1, second = -1;
	float distance = -1.0;

	DistNode(const int& first, const int& second, const float& dist):first(first), second(second), distance(dist)
	{}

	DistNode()
	{}
};

#endif
