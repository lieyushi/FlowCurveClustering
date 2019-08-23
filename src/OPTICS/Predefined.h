/*
 * @brief The priliminary functions needed to perform the OTPICS clustering
 * @author Lieyu Shi
 */


#ifndef _PREDEFINED_H
#define _PREDEFINED_H

#include "IOHandler.h"
#include "Initialization.h"
#include "Silhouette.h"


/*
 * @brief The point type defined by density-based clustering algorithm, based on whether the number of neighborhood
 * reaches the minPts or not
 */
enum PointType
{
	CORE = 0,
	BORDER,
	NOISE
};


/*
 * @brief point node struct, it has point type (PointType), visited (bool), group (cluser label)
 */
struct PointNode
{
	int type;
	bool visited;
	int group;
	float reachabilityDist;
	float core_distance;
	vector<int> neighbor;
	PointNode():type(-1), visited(false), group(-1), reachabilityDist(-1.0), core_distance(-1.0)
	{}

	~PointNode()
	{}
};


/*
 * @brief used to store objects with index and distance in priority queue
 */
struct OrderedPoint
{
	int index;
	float reachabilityDist;

	OrderedPoint(const int& index, const float& reachabilityDist):
	index(index), reachabilityDist(reachabilityDist)
	{}

	OrderedPoint(): index(-1), reachabilityDist(-1.0)
	{}
};


/*
 * @brief The class to define the node of point, this is actually the linked list implementation
 */
struct pointNode
{
	OrderedPoint value;
	pointNode* next;

	pointNode(const OrderedPoint& value): value(value), next(NULL)
	{		
	}

	pointNode(): value(OrderedPoint()), next(NULL)
	{
	}

	~pointNode()
	{
		if(next)
		{
			delete next;
			next = NULL;
		}
	}
};


/*
 * @brief The linked list C++ implementation
 */
class LinkedList
{
public:
	pointNode* start;

	LinkedList(): start(NULL)
	{}

	LinkedList(const OrderedPoint& value): start(new pointNode(value))
	{}

	~LinkedList()
	{
		deleteNode(start);
	}

	void insertNode(pointNode *vNode)
	{
		if(start==NULL)
		{
			start = vNode;
			return;
		}
		else if(vNode->value.reachabilityDist<start->value.reachabilityDist)
		{
			vNode->next = start;
			start = vNode;
			return;
		}
		pointNode *temp = start;
		while(temp->next && temp->next->value.reachabilityDist<vNode->value.reachabilityDist)
		{
			temp=temp->next;
		}
		vNode->next = temp->next;
		temp->next = vNode;
	}


	void updateNode(const int& index, const float& newDist)
	{
		pointNode *temp;
		if(start==NULL)
		{
			std::cout << "Linked list is empty!" << std::endl;
			return;
		}
		if(start->value.index==index)
		{
			temp = start;
			start = start->next;
			temp->next = NULL;
			temp->value.reachabilityDist = newDist;
			insertNode(temp);
			return;
		}

		pointNode *before = start;
		temp = before->next;
		while(temp)
		{
			if(temp->value.index==index)
				break;
			temp=temp->next;
			before=before->next;
		}
		if(temp==NULL)
		{
			std::cout << "Indexed node is not inside the linked list!" << std::endl;
			return;
		}
		else
		{
			temp->value.reachabilityDist = newDist;
			before->next = temp->next;
			temp->next = NULL;
			insertNode(temp);
		}
	}

private:
	// a recursive deletion for all nodes of start */
	void deleteNode(pointNode*& pNode)
	{
		if(pNode==NULL)
			return;
		else if(pNode && pNode->next==NULL)
		{
			delete pNode;
			pNode = NULL;
			return;
		}
		deleteNode(pNode->next);
	}
};


/*
 * @brief The class object to define the data set
 */
struct DataSet
{
	vector<vector<float> > dataVec;
	Eigen::MatrixXf dataMatrix;
	int maxElements;
	int vertexCount;
	int dimension;

	string strName;
	string fullName;

	DataSet()
	{}

	~DataSet()
	{}
};

#endif
