/*
 * Predefined.cpp
 *
 *  Created on: Nov 26, 2017
 *      Author: lieyu
 */

#include "Predefined.h"


// dfs treeNode to find leaf nodes
void dfsTraversal(std::vector<int>& elementList, TreeNode* node)
{
	//empty node
	if(node==NULL)
		return;
	//leaf node
	if(node->left==NULL&&node->right==NULL)
		elementList.push_back(node->index);
	if(node->left)
		dfsTraversal(elementList, node->left);
	if(node->right)
		dfsTraversal(elementList, node->right);
}

// delete treeNode nodes
void deleteTreeNode(TreeNode* &node)
{
	TreeNode* temp = node;
	if(node)
	{
		delete node;
		node = NULL;
		deleteTreeNode(temp->left);
		deleteTreeNode(temp->right);
	}
}


// remove two elements in template vector
template <class T>
void deleteVecElements(std::vector<T>& original, const T& first, const T& second)
{
	std::size_t size = original.size();
	assert(size>2);
	vector<T> result(size-2);
	int tag = 0;
	for(int i=0;i<size;++i)
	{
		//meet with target elements, not copied
		if(original[i]==first || original[i]==second)
			continue;
		result[tag++]=original[i];
	}
	assert(tag==size-2);
	original = result;
}
