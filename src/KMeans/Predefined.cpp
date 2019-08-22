/*
 * Predefined.cpp
 *
 *  Created on: Nov 26, 2017
 *      Author: lieyu
 */

#include "Predefined.h"


/*
 * @brief remove two elements in template vector
 *
 * @param[out] origine The vector to be operated on
 * @param[in] first The first index
 * @param[in] second The second index
 */
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
