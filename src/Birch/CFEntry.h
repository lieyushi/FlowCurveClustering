/*
 * CFEntry.h
 *
 *  Created on: Feb 11, 2017
 *      Author: cotrik
 */

#ifndef CFENTRY_H_
#define CFENTRY_H_

#include <vector>
#include <list>
#include <exception>
#include <assert.h>
#include <time.h>
#include <boost/cstdint.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/vector.hpp>
//#include "GeodesicDistance.h"
#include "CFNode.h"
//class CFEntry;
//class CFNode;
//class CFNodeItmd;
//class CFNodeLeaf;
/** this exception is produced when the current item size is not suitable. */
struct CFTreeInvalidItemSize : public std::exception {};

//enum { fdim = dim }; /** enum for recognizing dimension outside this class. */

typedef float float_type;  /** float type according to a precision - double, float, and so on. */
typedef std::vector<float_type> item_vec_type; /** vector of items. */
/** pointer type of CFNode.
 * shared_ptr is applied to preventing memory leakage.
 * This node pointer is deleted, having no referencers
 */
//typedef boost::shared_ptr<CFNode> cfnode_sptr_type;
//typedef std::pair<CFEntry*, CFEntry*> cfentry_pair_type; /** pair cfentry pointers. */
//typedef std::vector<CFEntry*> cfentry_ptr_vec_type; /** vector of cfentry pointers. */
//typedef float_type (*dist_func_type)(const CFEntry&, const CFEntry&); /** distance function pointer. */
//typedef std::vector<CFEntry> cfentry_vec_type; /** vector of cfentries. */

typedef boost::numeric::ublas::vector<float_type>           ublas_vec_type;         /* ublas vector in float_type. */
typedef boost::numeric::ublas::symmetric_matrix<float_type> ublas_sym_matrix_type;  /* ublas symmetric matrix in float_type. */

/** CFEntry is compact representation of group of data-points.
 * This entry contains linear sums each dimension and one square sum to represent data-points in this group
 */
template<boost::uint32_t dim>
class CFEntry
{
public:
    typedef boost::shared_ptr<CFNode<dim> > cfnode_sptr_type;
    typedef std::pair<CFEntry*, CFEntry*> cfentry_pair_type; /** pair cfentry pointers. */
    typedef std::vector<CFEntry*> cfentry_ptr_vec_type; /** vector of cfentry pointers. */
    typedef float_type (*dist_func_type)(const CFEntry&, const CFEntry&); /** distance function pointer. */
    typedef std::vector<CFEntry> cfentry_vec_type; /** vector of cfentries. */
    enum { fdim = dim }; /** enum for recognizing dimension outside this class. */
    /** Empty construct initialized with zeros */
    CFEntry() : n(0), sum_sq(0.0), id(0)
    {
        std::fill(sum, sum + dim, 0);
    }
    CFEntry(const CFEntry& other)
    : n(other.n)
    , sum_sq(other.n)
    , child (other.child)
    , id(0)
    {
        //std::fill(sum, sum + dim, 0);
        std::copy(other.sum, other.sum + dim, sum);
    }

//    float GeoDistance(const CFEntry& other) {
//        //gd.Init();
////        float d = gd.GetGeodesicDistance(id, other.id);
//        //gd.Destroy();
//        return d;
//    }

    /** Constructor when array of T type items come.
     * initialize CFEntry with one data-point
     */
    template<typename T>
    CFEntry( T* item ) : n(1), sum_sq(0.0), id(0)
    {
        std::copy( item, item + dim, sum );
        for( std::size_t i = 0 ; i < dim ; i++ )
            sum_sq += item[i] * item[i];
    }

    template<typename T>
    CFEntry( T* item, std::size_t id) : n(1), sum_sq(0.0), id(id)
    {
        std::copy( item, item + dim, sum );
        for( std::size_t i = 0 ; i < dim ; i++ )
            sum_sq += item[i] * item[i];
    }

    /** Constructor for root entry with children */
    CFEntry( const cfnode_sptr_type& in_child ) : n(0), sum_sq(0.0), child(in_child), id(0)
    {
        std::fill(sum, sum + dim, 0);
    }

    /** Operator returning a new CFEntry merging from two CFEntries */
    CFEntry operator+( const CFEntry& rhs )
    {
        CFEntry e;
        e.n = n + rhs.n;
        for( std::size_t i = 0 ; i < dim ; i++ )
            e.sum[i] = sum[i] + rhs.sum[i];
        e.sum_sq = sum_sq + rhs.sum_sq;

        return e;
    }

    /** Operator merging two CFEntries to the left-hand-side CFEntry */
    void operator+=( const CFEntry& e )
    {
        for( std::size_t i = 0 ; i < dim ; i++ )
        {
            float_type val = e.sum[i];
            sum[i] += val;
            sum_sq += val*val;
        }
        n += e.n;
    }

    /** Operator removing data-points from one CFEntry  */
    void operator-=( const CFEntry& e )
    {
        for( std::size_t i = 0 ; i < dim ; i++ )
        {
            float_type val = e.sum[i];
            sum[i] -= val;
            sum_sq -= val*val;
        }
        n -= e.n;
    }

    bool operator==( const CFEntry& e ) const
    {
        if (n != e.n || sum_sq != e.sum_sq || child != e.child) return false;
        for (int i = 0; i < dim; i++) if (sum[i] != e.sum[i]) return false;

        return true;
    }

    /** Does this CFEntry have children? */
    bool HasChild() const   { return child.get() != NULL; }

    std::size_t         n;          /* the number of data-points in */
    float_type          sum[dim];   /* linear sum of each dimension of n data-points */
    float_type          sum_sq;     /* square sum of n data-points */
    cfnode_sptr_type    child;      /* pointer to a child node */
    std::size_t         id;
   // static GeodesicDistance gd;
};

//template <boost::uint32_t dim>
//GeodesicDistance CFEntry<dim>::gd;

#endif /* CFENTRY_H_ */
