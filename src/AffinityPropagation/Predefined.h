/* Predefined.h file for spectral clustering
 * Needed to be decided which should be contained
 */


#ifndef _PREDEFINED_H
#define _PREDEFINED_H

#include "IOHandler.h"
#include "Initialization.h"
#include "Silhouette.h"



struct DataSet
{
	vector<vector<float> > dataVec;	//original dataset
	Eigen::MatrixXf dataMatrix;	//sampled dataset
	int maxElements = -1;
	int vertexCount = -1;
	int dimension = -1;

	string strName;
	string fullName;
	string dataName;

};


struct Ensemble
{
	int size;
	std::vector<int> element;
};



/* my template function for std::swap */
template <class T>
void mySwap(T& a, T& b)
{
    T temp = a;
    a = b;
    b = temp;
}


/* my template-based partition function */
template <class T>
int partition(std::vector<T>& array, const int& left, const int& right, const int& pivotIndex)
{
    T pivotValue = array[pivotIndex];
    mySwap(array[pivotIndex], array[right]);
    int storeIndex = left;
    for(int i=left; i<right;++i)
    {
        if(array[i]<pivotValue)
        {
            mySwap(array[storeIndex], array[i]);
            ++storeIndex;
        }
    }
    mySwap(array[right], array[storeIndex]);
    return storeIndex;
}


/* return k-th element in the unsorted array with quicks-selection algorithm, pseudocode referrenced at https://en.wikipedia.org/wiki/Quickselect */
template <class T>
/* left is left index while right is right index */
T select(std::vector<T>& array, int left, int right, const int& k)
{
    int pivotIndex;
    while(true)
    {
        if(left==right)
            return array[left];
        pivotIndex = (left+right)/2;
        pivotIndex = partition(array, left, right, pivotIndex);
        if(k==pivotIndex)
            return array[k];
        else if(k<pivotIndex)
            right = pivotIndex-1;
        else
            left = pivotIndex+1;
    }
}


#endif
