/*
 * item_type.h
 *
 *  Created on: Feb 13, 2017
 *      Author: cotrik
 */

#ifndef ITEM_TYPE_H_
#define ITEM_TYPE_H_
template<boost::uint32_t dim>
class CFTree;

template<boost::uint32_t dim>
struct item_type
{
    item_type() : id(0) { std::fill( item, item + sizeof(item)/sizeof(item[0]), 0 ); }
    item_type( float* in_item ) : id(0) { std::copy(in_item, in_item+sizeof(item)/sizeof(item[0]), item); }
    float& operator[]( int i ) { return item[i]; }
    float operator[]( int i ) const { return item[i]; }
    std::size_t size() const { return sizeof(item)/sizeof(item[0]); }

    int& cid() { return id; }
    const int cid() const { return id; }

    //float item[cftree_type::fdim];
    float item[dim];
    int id;
};

#endif /* ITEM_TYPE_H_ */
