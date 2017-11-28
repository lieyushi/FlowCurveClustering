#ifndef _PREDEFINED_H
#define _PREDEFINED_H

#include "../Common/IOHandler.h"
#include "../Common/Initialization.h"
#include "../Common/Silhouette.h"



struct DataSet
{
	vector<vector<float> > dataVec;	//original dataset
	Eigen::MatrixXf dataMatrix;	//sampled dataset
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

// define a treeNode structure to store AHC clustering tree
struct TreeNode
{
	TreeNode* left;
	TreeNode* right;
	int index;

	TreeNode():left(NULL), right(NULL), index(-1)
	{}

	TreeNode(const int& index): left(NULL), right(NULL), index(index)
	{}

	~TreeNode()
	{
	}

	bool operator==(TreeNode* others)
	{
		if(!this&&!others)
			return true;
		else if(!this||!others)
			return false;
		return (index==others->index)&&(left==others->left)&&(right==others->right);
	}

	bool operator<(TreeNode* others)
	{
		if(!this || !others)
			return false;
		return index<others->index;
	}
};

// dfs treeNode to find leaf nodes
void dfsTraversal(vector<int>& elementList, TreeNode* node);

// delete treeNode nodes
void deleteTreeNode(TreeNode* &node);

// remove two elements in template vector
template <class T>
void deleteVecElements(std::vector<T>& origine, const T& first, const T& second);

#endif
