/*
 * CFNode.h
 *
 *  Created on: Feb 11, 2017
 *      Author: cotrik
 */
/** CFNode is composed of several CFEntries within page-size, and acts like B-tree node.
 *
 * CFNode should be page-sized for more efficient operation.
 * Like b-tree twist their node when removing and inserting node, CFTree perform similar operations on its own CFNodes.
 *
 * CFNode has two types: intermediate node leaf node, especially leaf node has additional pointers to neighbor leaves.
 */
#ifndef CFNODE_H_
#define CFNODE_H_
#include <vector>
#include <list>
#include <exception>
#include <assert.h>
#include <time.h>
#include <boost/cstdint.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/vector.hpp>
//#include "CFEntry.h"
#define PAGE_SIZE           (500*1024) /* assuming 4K page */
#define ARRAY_COUNT(a)      (sizeof(a)/sizeof(a[0]))
template<boost::uint32_t dim>
class CFEntry;
template<boost::uint32_t dim>
class CFNodeLeaf;
template<boost::uint32_t dim>
struct CFNode
{
    CFNode() : size(0) {
        //entries.resize((PAGE_SIZE - ( sizeof(CFNodeLeaf<dim> *)*2 /* 2 leaf node pointers */ + sizeof(std::size_t) /* size */ + sizeof(void*) /* vtptr */ )) / sizeof(CFEntry<dim> )/*max_entries*/);
    }
//    CFNode(const CFNode<dim>& other)
//    : size(other.size)
//    , entries(other.entries)
//    {
//    }
    virtual ~CFNode() {}
    virtual bool IsLeaf() const = 0;

    /** add new CFEntry to this CFNode */
    void Add( CFEntry<dim>& e )
    {
//        entries.push_back(e);
//        size++;
        assert( size < MaxEntrySize() );
        entries[size++] = e;
    }

    /** replace old CFEntry as new CFEntry */
    void Replace(CFEntry<dim>& old_entry, CFEntry<dim>& new_entry)
    {
        for (std::size_t i = 0; i < size; i++) {
            if (entries[i] == old_entry) {
                entries[i] = new_entry;
                return;
            }
        }
        // should not be reached here if replacement is successful
        assert(false);
    }

    /** Max # of CFEntries this CFNode could contain */
    std::size_t MaxEntrySize() const
    {
        return ARRAY_COUNT(entries);
    }

    /** CFNode is full, no more CFEntries can be in */
    bool IsFull() const
    {
        return size == MaxEntrySize();
    }

    /** CFNode has nothing  */
    bool IsEmpty() const
    {
        return size == 0;
    }

    std::size_t     size;   /** # CFEntries this CFNode contains */
    //std::vector<CFEntry<dim> > entries;
    CFEntry<dim>         entries[(PAGE_SIZE - ( sizeof(CFNodeLeaf<dim> *)*2 /* 2 leaf node pointers */ + sizeof(std::size_t) /* size */ + sizeof(void*) /* vtptr */ )) / sizeof(CFEntry<dim> )/*max_entries*/]; /** Array of CFEntries */
};

#endif /* CFNODE_H_ */
