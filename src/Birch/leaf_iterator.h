/*
 * leaf_iterator.h
 *
 *  Created on: Feb 11, 2017
 *      Author: cotrik
 */

#ifndef LEAF_ITERATOR_H_
#define LEAF_ITERATOR_H_
#include "CFNodeLeaf.h"
//class CFNodeLeaf;
typedef std::random_access_iterator_tag iterator_category;
//template<boost::uint32_t dim>
//typedef CFNodeLeaf<dim>      T;
//typedef T               value_type;
//typedef T&              reference;
//typedef T*              pointer;
typedef std::ptrdiff_t  difference_type;
/** leaf iterator */
template<boost::uint32_t dim>
struct leaf_iterator : public std::forward_iterator_tag
{
    leaf_iterator( CFNodeLeaf<dim>* in_leaf ) : leaf( in_leaf ) {}
    leaf_iterator    operator++() { leaf = (CFNodeLeaf<dim>*)leaf->next.get(); return leaf_iterator(leaf); }
    bool             operator!=( const leaf_iterator rhs ) const { return !(leaf == rhs.leaf); }
    CFNodeLeaf<dim>& operator*() { return *leaf; }
    CFNodeLeaf<dim>* operator->() { return leaf; }

    CFNodeLeaf<dim>* leaf;
};



#endif /* LEAF_ITERATOR_H_ */
