/*
 *  This file is part of birch-clustering-algorithm.
 *
 *  birch-clustering-algorithm is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  birch-clustering-algorithm is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with birch-clustering-algorithm.  If not, see <http://www.gnu.org/licenses/>.
 *
 *	Copyright (C) 2011 Taesik Yoon (otterrrr@gmail.com)
 */

#ifndef __CFTREE_H__
#define	__CFTREE_H__

#include "CFEntry.h"
#include "CFNode.h"
#include "CFNodeItmd.h"
#include "CFNodeLeaf.h"
#include "leaf_iterator.h"
#include "subcluster_summary.h"
#include "item_type.h"
#include "../Common/Distance.h"

#include <math.h>
#include <map>
#include <list>
#include <fstream>
#include <exception>
#include <assert.h>
#include <time.h>
#include <boost/cstdint.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/vector.hpp>
#ifndef FALSE
	#define FALSE 0
#endif

#ifndef TRUE
	#define TRUE 1
#endif

//struct item_type;

/** class CFTree ( clustering feature tree ).
 * 
 * according to the paper,
 * birch maintains btree-like data structure consisting of summarized clusters
 * 
 * @param dim  dimensions of item, this parameter should be fixed before compiling
 */

template<boost::uint32_t dim>
class CFTree
{
//	class CFEntry;
//	class CFNode;
//	class CFNodeItmd;
//	class CFNodeLeaf;

public:

	static MetricPreparation object;
	static int normOption;
	static int totalNodes;

	/** this exception is produced when the current item size is not suitable. */
	struct CFTreeInvalidItemSize : public std::exception {};

	enum { fdim = dim }; /** enum for recognizing dimension outside this class. */

	typedef	float float_type;	/** float type according to a precision - double, float, and so on. */
	typedef std::vector<float_type>	item_vec_type; /** vector of items. */
//	/** pointer type of CFNode.
//	 * shared_ptr is applied to preventing memory leakage.
//	 * This node pointer is deleted, having no referencers
//	 */
	typedef boost::shared_ptr<CFNode<dim> > cfnode_sptr_type;
	typedef std::pair<CFEntry<dim>*, CFEntry<dim>*> cfentry_pair_type; /** pair cfentry pointers. */
	typedef std::vector<CFEntry<dim>*> cfentry_ptr_vec_type; /** vector of cfentry pointers. */
	typedef float_type (*dist_func_type)(const CFEntry<dim>&, const CFEntry<dim>&); /** distance function pointer. */
	typedef std::vector<CFEntry<dim> > cfentry_vec_type; /** vector of cfentries. */

	typedef boost::numeric::ublas::vector<float_type>			ublas_vec_type;			/* ublas vector in float_type. */
	typedef boost::numeric::ublas::symmetric_matrix<float_type>	ublas_sym_matrix_type;	/* ublas symmetric matrix in float_type. */
	typedef std::vector< subcluster_summary > subsum_vec_type;

public:

	/** Euclidean Distance of Centroid */
	static float_type _DistD0( const CFEntry<dim>& lhs, const CFEntry<dim>& rhs )
	{
		float_type dist = 0.0;
		float_type tmp;
		for (std::size_t i = 0 ; i < dim ; i++) {
			tmp =  lhs.sum[i]/lhs.n - rhs.sum[i]/rhs.n;
			dist += tmp*tmp;
		}
		//assert(dist >= 0.0);
		dist = sqrt(dist);
		return (std::max)(dist, (float)0.0);
	}

	/** Manhattan Distance of Centroid */
	static float_type _DistD1( const CFEntry<dim>& lhs, const CFEntry<dim>& rhs )
	{
		float_type dist = 0.0;
		float_type tmp;
		for (std::size_t i = 0 ; i < dim ; i++) {
			tmp = std::abs(lhs.sum[i]/lhs.n - rhs.sum[i]/rhs.n);
			dist += tmp;
		}
		//assert(dist >= 0.0);
		return (std::max)(dist, (float)0.0);
	}

	/** Pairwise IntraCluster Distance */
	static float_type _DistD2( const CFEntry<dim>& lhs, const CFEntry<dim>& rhs )
	{
		float_type dot = 0.0;
		for(std::size_t i = 0 ; i < dim ; i++)
			dot += lhs.sum[i] * rhs.sum[i];

		float_type dist = ( rhs.n*lhs.sum_sq + lhs.n*rhs.sum_sq - 2*dot ) / (lhs.n*rhs.n);
		//assert(dist >= 0.0);
		return std::max(dist, (float)0.0);
	}

	/** Pairwise InterClusterDistance */
	static float_type _DistD3(const CFEntry<dim>& lhs, const CFEntry<dim>& rhs)
    {
        std::size_t tmpn = lhs.n + rhs.n;
        float_type tmp1, tmp2 = 0.0;
        for (std::size_t i = 0; i < dim; i++) {
            tmp1 = lhs.sum[i] + rhs.sum[i];
            tmp2 += tmp1 / tmpn * tmp1 / (tmpn - 1);
        }
        float_type dist = 2 * ((lhs.sum_sq + rhs.sum_sq) / (tmpn - 1) - tmp2);
        //assert(dist >= 0.0);
        return std::max(dist, (float)0.0);
    }

    /** Pairwise InterClusterDistance */
    static float_type _DistD4(const CFEntry<dim>& lhs, const CFEntry<dim>& rhs)
    {
        float_type dist = 0.0;
        float_type tmp;
        for (std::size_t i = 0 ; i < dim ; i++) {
            tmp =  lhs.sum[i]/lhs.n - rhs.sum[i]/rhs.n;
            dist += pow(tmp, 4);
        }
        dist = pow(dist, 1.0/4);
        //assert(dist >= 0.0);
        return (std::max)(dist, (float)0.0);
    }

    /** Pairwise InterClusterDistance */
    static float_type _DistD8(const CFEntry<dim>& lhs, const CFEntry<dim>& rhs)
    {
        float_type dist = 0.0;
        float_type tmp;
        for (std::size_t i = 0 ; i < dim ; i++) {
            tmp =  lhs.sum[i]/lhs.n - rhs.sum[i]/rhs.n;
            dist += pow(tmp, 8);
        }
        dist = pow(dist, 1.0/8);
        //assert(dist >= 0.0);
        return (std::max)(dist, (float)0.0);
    }

    /** Pairwise InterClusterDistance */
    static float_type _DistD16(const CFEntry<dim>& lhs, const CFEntry<dim>& rhs)
    {
        float_type dist = 0.0;
        float_type tmp;
        for (std::size_t i = 0 ; i < dim ; i++) {
            tmp =  lhs.sum[i]/lhs.n - rhs.sum[i]/rhs.n;
            dist += pow(tmp, 16);
        }
        dist = pow(dist, 1.0/16);
        //assert(dist >= 0.0);
        return (std::max)(dist, (float)0.0);
    }

    
    static float_type _DistLarry(const CFEntry<dim>& lhs, const CFEntry<dim>& rhs)
    {
    	Eigen::VectorXf rLeft(dim), rRight(dim);
    	rLeft = VectorXf::Map(&(lhs.sum[0]), dim);
    	rRight = VectorXf::Map(&(rhs.sum[0]), dim);
    	rLeft /= lhs.n;
    	rRight /= rhs.n;

    	//if(distanceMatrix)
    	//	return distanceMatrix[lhs.id][rhs.id];
    	//else
    		return getDisimilarity(rLeft,rRight,lhs.id,rhs.id,normOption,object);
    }

	static float_type _Diameter(const CFEntry<dim>& e)
    {
        if (e.n <= 1)
            return 0.0;

        float_type temp = 0.0;
        for (std::size_t i = 0; i < dim; i++)
            temp += e.sum[i] / e.n * e.sum[i] / (e.n - 1);

        float_type diameter = 2 * (e.sum_sq / (e.n - 1) - temp);

        //assert(diameter >= 0.0);
        return (std::max)(diameter, (float)0.0);
    }

	/** Radius of the CFEntry */
    static float_type _Radius(const CFEntry<dim>& e)
    {
        if (e.n <= 1)
            return 0.0;

        float_type tmp0, tmp1 = 0.0;
        for (std::size_t i = 0; i < dim; i++) {
            tmp0 = e.sum[i] / e.n;
            tmp1 += tmp0 * tmp0;
        }
        float_type radius = e.sum_sq / e.n - tmp1;

        //assert(radius >= 0.0);
        return std::max(radius, (float)0.0);
    }

private:
	/** Functor is used for choosing the closest one */
	struct CloseEntryLessThan
	{
		CloseEntryLessThan( const CFEntry<dim>& in_base_entry, const dist_func_type& in_dist_func ) : base_entry(in_base_entry), dist_func(in_dist_func) {}
		bool operator()(const CFEntry<dim>& lhs, const CFEntry<dim>& rhs)	{ return dist_func(lhs, base_entry) < dist_func(rhs, base_entry); }

		const CFEntry<dim>& base_entry;
		const dist_func_type& dist_func;
	};

public:


	/** CFTree construct with memory limit and designated distance functions
	 *
	 * @param in_dist_threshold range within a CFEntry
	 * @param in_mem_limit memory limit to which CFTree can utilize, if CFTree overflows this limit, then distance threshold become larger to rebuild more compact CFTree
	 * @param in_dist_func distance function between CFEntries
	 * @param in_dist_func distance function tests if a new data-point should be absorbed or not
	 **/
	CFTree( float_type in_dist_threshold, std::size_t in_mem_limit, dist_func_type in_dist_func = _DistLarry, dist_func_type in_absorb_dist_func = _DistLarry) :
		mem_limit(in_mem_limit), dist_threshold(in_dist_threshold), root( new CFNodeLeaf<dim>() ), dist_func(in_dist_func), absorb_dist_func(in_absorb_dist_func), node_cnt(1/* root node */),
		leaf_dummy( new CFNodeLeaf<dim>() )
	{
		((CFNodeLeaf<dim>*)leaf_dummy.get())->next = root;
	}

	void setDistThreshold(const float_type& in_dist_threshold)
	{
		dist_threshold = in_dist_threshold;
	}

	~CFTree(void) {}

	/** whether this CFTree is empty or not */
	bool empty() const { return root->IsEmpty(); }

	/** inserting one data-point */
	void insert( item_vec_type& item )
	{
		if( item.size() != dim )
			throw CFTreeInvalidItemSize();

		/*std::cout << "n is: " << CFEntry<dim>::n << std::endl;
		for (int i = 0; i < dim; ++i)
		{
		 	std::cout << cfentry_vec_type[0].sum[i] << std::endl;
		}*/

		insert(&item[0]);
	}

	/** inserting one data-point with T typed */
	template<typename T>
	void insert(T* item)
	{
	    static std::size_t id = 0;
		//CFEntry<dim> e(item);
	    //CFEntry<dim> e(item, id++);
	    CFEntry<dim> e(item, id%totalNodes);
	    ++id;
		insert(e);
		//if (id % 5448 == 0)
		//std::cout << id << " item\n";
	}

	/** inserting a new entry */
	void insert(CFEntry<dim>& e)
    {
        bool bsplit;
        insert(root.get(), e, bsplit);

        // there's no exception for the root as regard to splitting, indeed
        if (bsplit) {
            split_root(e);
        }

        std::size_t curr_mem = node_cnt * sizeof(CFNode<dim> );
        if (mem_limit > 0 && node_cnt * sizeof(CFNode<dim> ) > mem_limit) {
            rebuild();
        }
    }

	/** get the beginning of leaf iterators */
	leaf_iterator<dim> leaf_begin() { return leaf_iterator<dim>( (CFNodeLeaf<dim>*)((CFNodeLeaf<dim>*)leaf_dummy.get())->next.get()); }
	/** get the end of leaf iterators  */
	leaf_iterator<dim> leaf_end() { return leaf_iterator<dim>(NULL); }

	/** get leaf entries */
	void get_entries( cfentry_vec_type& out_entries )
	{
		std::size_t n_leaf_entries = 0;
		//leaf_iterator<dim> it = leaf_begin();
		for( leaf_iterator<dim> it = leaf_begin() ; it != leaf_end() ; ++it)
			n_leaf_entries += it->size;

		out_entries.clear();
		out_entries.reserve(n_leaf_entries);
		for( leaf_iterator<dim> it = leaf_begin() ; it != leaf_end() ; ++it )
			std::copy( it->entries, it->entries + it->size, std::back_inserter(out_entries) );
	}

private:

	void insert(CFNode<dim>* node, CFEntry<dim>& new_entry, bool &bsplit)
    {
        // empty node, it might be root node at first insertion
        if (node->IsEmpty()) {
            node->Add(new_entry);
            bsplit = false;
            return;
        }

        CFEntry<dim>& close_entry = *find_close(node, new_entry);

        // non-leaf
        if (close_entry.HasChild()) {
            insert(close_entry.child.get(), new_entry, bsplit);

            // no more split
            if (!bsplit)
                close_entry += (new_entry);
            // split here
            else
                split(*node, close_entry, new_entry, bsplit);
        }
        //leaf
        else {
            // absorb
            if (absorb_dist_func(close_entry, new_entry) < dist_threshold) {
                close_entry += (new_entry);
                bsplit = false;
            }
            // add new_entry
            else if (node->size < node->MaxEntrySize()) {
                node->Add(new_entry);
                bsplit = false;
            }
            // handle with the split cond. at parent-level
            else {
                bsplit = true;
            }
        }
    }

	CFEntry<dim>* find_close( CFNode<dim>* node, CFEntry<dim>& new_entry )
	{
		CFEntry<dim>* begin = &node->entries[0];
		CFEntry<dim>* end = begin + node->size;
		CFEntry<dim>* e = std::min_element( begin, end, CloseEntryLessThan(new_entry, dist_func) );
		return e != end ? e : NULL;
	}

	void split( CFNode<dim>& node, CFEntry<dim>& close_entry, CFEntry<dim>& new_entry, bool& bsplit )
	{
		CFNode<dim>* old_node = close_entry.child.get();
		assert( old_node != NULL );

		// make the list of entries, old entries
		cfentry_ptr_vec_type entries;
		entries.reserve( old_node->MaxEntrySize() + 1 );
		for( std::size_t i = 0 ; i < root->MaxEntrySize() ; i++ )
			entries.push_back(&old_node->entries[i]);
		entries.push_back(&new_entry);

		// find the farthest entry pair
		cfentry_pair_type far_pair;
		find_farthest_pair( entries, far_pair );

		bool node_is_leaf = old_node->IsLeaf();

		// make two split nodes
		cfnode_sptr_type node_lhs( node_is_leaf ? (CFNode<dim>*) new CFNodeLeaf<dim>() : (CFNode<dim>*) new CFNodeItmd<dim>() );
		cfnode_sptr_type node_rhs( node_is_leaf ? (CFNode<dim>*) new CFNodeLeaf<dim>() : (CFNode<dim>*) new CFNodeItmd<dim>() );

		// two entries for new root node
		// and connect child node to the entries
		CFEntry<dim> entry_lhs( node_lhs );
		CFEntry<dim> entry_rhs( node_rhs );

		if( node_is_leaf )
		{
			assert( node_lhs->IsLeaf() && node_rhs->IsLeaf() );
			
			CFNodeLeaf<dim>* leaf_node = (CFNodeLeaf<dim>*)old_node;

			cfnode_sptr_type prev = leaf_node->prev;
			cfnode_sptr_type next = leaf_node->next;

			if( prev != NULL )
				((CFNodeLeaf<dim>*)prev.get())->next = node_lhs;
			if( next != NULL )
				((CFNodeLeaf<dim>*)next.get())->prev = node_rhs;

			((CFNodeLeaf<dim>*)node_lhs.get())->prev = prev;
			((CFNodeLeaf<dim>*)node_lhs.get())->next = node_rhs;
			((CFNodeLeaf<dim>*)node_rhs.get())->prev = node_lhs;
			((CFNodeLeaf<dim>*)node_rhs.get())->next = next;
		}

		// rearrange old entries to new entries
		rearrange(entries, far_pair, entry_lhs, entry_rhs);

		// one old entry is divided into to new entries
		// so the first one is included instead of old ones
		node.Replace(close_entry, entry_lhs);

		// the full node indicates that this node have to be split as well
		bsplit = node.IsFull();

		// copy the second entry newly created into return variable 'new_entry'
		if( bsplit )
			new_entry = entry_rhs;
		// if affordable, not split, add the second entry to the node
		else
			node.Add(entry_rhs);

		// for statistics and mornitoring memory usage
		node_cnt++;
	}

	void split_root( CFEntry<dim>& e )
	{
		// make the list of entries, old entries
		cfentry_ptr_vec_type entries;
		entries.reserve(root->MaxEntrySize() + 1);
		for( std::size_t i = 0 ; i < root->MaxEntrySize() ; i++ )
			entries.push_back(&root->entries[i]);
		entries.push_back(&e);

		// find the farthest entry pair
		cfentry_pair_type far_pair;
		find_farthest_pair( entries, far_pair );

		bool root_is_leaf = root->IsLeaf();

		// make two split nodes
		cfnode_sptr_type node_lhs( root_is_leaf ? (CFNode<dim>*) new CFNodeLeaf<dim>() : (CFNode<dim>*) new CFNodeItmd<dim>() );
		cfnode_sptr_type node_rhs( root_is_leaf ? (CFNode<dim>*) new CFNodeLeaf<dim>() : (CFNode<dim>*) new CFNodeItmd<dim>() );

		// two entries for new root node
		// and connect child node to the entries
		CFEntry<dim> entry_lhs( node_lhs );
		CFEntry<dim> entry_rhs( node_rhs );

		// new root node result in two entries each of which has split node respectively
		cfnode_sptr_type new_root( new CFNodeItmd<dim>() );

		// update prev/next links of newly created leaves
		if( root_is_leaf )
		{
			assert( node_lhs->IsLeaf() && node_rhs->IsLeaf() );
			((CFNodeLeaf<dim>*)leaf_dummy.get())->next = node_lhs;
			((CFNodeLeaf<dim>*)node_lhs.get())->prev = leaf_dummy;
			((CFNodeLeaf<dim>*)node_lhs.get())->next = node_rhs;
			((CFNodeLeaf<dim>*)node_rhs.get())->prev = node_lhs;
		}

		// rearrange old entries to new entries
		rearrange( entries, far_pair, entry_lhs, entry_rhs );

		// substitute new_root to 'root' variable
		new_root->Add(entry_lhs);
		new_root->Add(entry_rhs);
		root = new_root;

		// for statistics and mornitoring memory usage
		node_cnt++;
	}

	void rearrange( cfentry_ptr_vec_type& entries, cfentry_pair_type& far_pair, CFEntry<dim>& entry_lhs, CFEntry<dim>& entry_rhs )
	{
		entry_lhs.child->Add(*far_pair.first);
		entry_lhs += *far_pair.first;
		entry_rhs.child->Add(*far_pair.second);
		entry_rhs += *far_pair.second;

		for (std::size_t i = 0; i < entries.size(); i++) {
            CFEntry<dim>& e = *entries[i];
            if (&e == far_pair.first || &e == far_pair.second)
                continue;

            float_type dist_first = dist_func(*far_pair.first, e);
            float_type dist_second = dist_func(*far_pair.second, e);

            CFEntry<dim>& e_update = dist_first < dist_second ? entry_lhs : entry_rhs;
            e_update.child->Add(e);
            e_update += e;
        }
	}

	void find_farthest_pair( std::vector<CFEntry<dim>*>& entries, /* out */cfentry_pair_type& far_pair )
	{
		//assert( entries.size() >= 2 );

		float_type max_dist = -1.0;
        for (std::size_t i = 0; i < entries.size() - 1; i++) {
            for (std::size_t j = i + 1; j < entries.size(); j++) {
                CFEntry<dim>& e1 = *entries[i];
                CFEntry<dim>& e2 = *entries[j];

                float_type dist = dist_func(e1, e2);
                if (max_dist < dist) {
                    max_dist = dist;
                    far_pair.first = &e1;
                    far_pair.second = &e2;
                }
            }
        }
	}

	float_type average_dist_closest_pair_leaf_entries()
	{
		std::size_t total_n = 0;
		float_type	total_d = 0.0;
		float_type	dist;

		// determine new threshold
		CFNodeLeaf<dim>* leaf = (CFNodeLeaf<dim>*)leaf_dummy.get();
		while (leaf != NULL) {
            if (leaf->size >= 2) {
                std::vector<float_type> min_dists(leaf->size, (std::numeric_limits<float_type>::max)());
                for (std::size_t i = 0; i < leaf->size - 1; i++) {
                    for (std::size_t j = i + 1; j < leaf->size; j++) {
                        dist = dist_func(leaf->entries[i], leaf->entries[j]);
                        dist = dist >= 0.0 ? sqrt(dist) : 0.0;
                        if (min_dists[i] > dist)
                            min_dists[i] = dist;
                        if (min_dists[j] > dist)
                            min_dists[j] = dist;
                    }
                }
                for (std::size_t i = 0; i < leaf->size; i++)
                    total_d += min_dists[i];
                total_n += leaf->size;
            }

            // next leaf
            leaf = (CFNodeLeaf<dim>*) leaf->next.get();
        }
		return total_d / total_n;
	}
public:
	/** rebuild tree from the existing leaf entries.
	 *
	 * rebuilding cftree is regarded as clustering, because there could be overlapped cfentries.
	 * birch guarantees datapoints in cfentries within a range, but two data-points within a range can be separated to different cfentries
	 *
	 * @param extend	if true, the size of tree reaches to memory limit, so distance threshold enlarges.
	 *					in case of both true and false, rebuilding CFTree from the existing leaves.
	 */
	void rebuild( bool extend = true )
	{
		if( extend )
		{
			// decide the next threshold
			float_type new_threshold = std::pow(average_dist_closest_pair_leaf_entries(),2);
			dist_threshold = dist_threshold > new_threshold ? dist_threshold*2 : new_threshold;
		}

		// construct a new tree by inserting all the node from the previous tree
		CFTree<dim> new_tree( dist_threshold, mem_limit );
		
		CFNodeLeaf<dim>* leaf = (CFNodeLeaf<dim>*) leaf_dummy.get();
		//std::cout << "leaf size is " << leaf->size << std::endl;
		assert(leaf!=NULL);
        while (leaf != NULL) {
            for (std::size_t i = 0; i < leaf->size; i++)
                new_tree.insert(leaf->entries[i]);

            // next leaf
            leaf = (CFNodeLeaf<dim>*) leaf->next.get();
        }

		// really I'd like to replace the previous tree to the new one by
		// stating " *this = new_tree; ", but it doesn't work because 'this' is const pointer
		// copy root and dummy_node
		// copy statistics variable

		root = new_tree.root;
		leaf_dummy = new_tree.leaf_dummy;
		node_cnt = new_tree.node_cnt;
	}

private:

	// data structure
	cfnode_sptr_type	root;
	cfnode_sptr_type	leaf_dummy;	/* start node of leaves */
	
	// parameters
	std::size_t			mem_limit;
	float_type			dist_threshold;
	dist_func_type		dist_func;
	dist_func_type		absorb_dist_func;

	// statistics
	std::size_t			node_cnt;

/* phase 3 - applying a global clustering algorithm to subclusters */
//#include "CFTree_CFCluster.h"
private:
    bool _has_differences(const std::vector<ublas_vec_type>& lhs, const std::vector<ublas_vec_type>& rhs) const
    {
        assert(lhs.size() == rhs.size());

        for (std::size_t i = 0; i < lhs.size(); i++) {
            const ublas_vec_type diff = lhs[i] - rhs[i];
            if (norm_2(diff) > std::numeric_limits<float_type>::epsilon())
                return true;
        }
        return false;
    }

    bool _has_differences(const std::vector<ublas_vec_type>& lhs, const std::vector<ublas_vec_type>& rhs, std::vector<bool>& active )
    {
        assert( lhs.size() == rhs.size() );

        for (std::size_t i = 0; i < lhs.size(); i++){
            const ublas_vec_type diff = lhs[i] - rhs[i];
            active[i] = norm_2(diff) > std::numeric_limits<float_type>::epsilon();
        }
        return std::count( active.begin(), active.end(), true ) > 0;
    }

public:
//    template<typename item_list_type>
//    void redist_kmeans(item_list_type& items, cfentry_vec_type& entries, std::size_t iteration = 2)
//    {
//        using namespace boost::numeric::ublas;
//
//        if (items.empty())
//            return;
//
//        assert(items[0].size() == dim);
//
//        // start from k means from k entries
//        std::vector<ublas_vec_type> prev_means(entries.size());
//        std::vector<ublas_vec_type> means(entries.size());
//        for (std::size_t i = 0; i < means.size(); i++) {
//            prev_means[i].resize(dim);
//            prev_means[i].clear();
//
//            CFEntry<dim>& e = entries[i];
//            ublas_vec_type& mean = means[i];
//
//            mean.resize(dim);
//            std::copy(e.sum, e.sum + dim, mean.begin());
//            mean /= e.n;
//        }
//
//        // until it is converged
//        std::size_t iteration_count = 0;
//        if( iteration == 0 )
//            iteration = (std::numeric_limits<std::size_t>::max)();
//
//        bool active = true;
//
//        //while( iteration_count < iteration && _has_differences( prev_means, means ) )
//        while (iteration_count < iteration && active) {
//            active = false;
//            for (item_list_type::iterator item_it = items.begin(); item_it != items.end(); ++item_it) {
//                item_list_type::value_type & item = *item_it;
//                float_type min_dist = (std::numeric_limits<float_type>::max)();
//
//                int prev_cid = item.cid();
//
//                for (std::size_t cid = 0; cid < means.size(); ++cid) {
//                    ublas_vec_type diff(dim);
//                    std::transform(&item[0], &item[0] + dim, means[cid].begin(), diff.begin(),
//                            std::minus<float_type>());
//                    float_type dist = norm_2(diff);
//
//                    if (min_dist > dist) {
//                        min_dist = dist;
//                        item_it->cid() = cid;
//                    }
//                }
//
//                if (prev_cid != item.cid())
//                    active = true;
//            }
//
//            std::stringstream ss;
//            ss << "k-means_iteration" << iteration_count << ".txt";
//            std::ofstream fout(ss.str().c_str());
//
//            for (std::size_t c = 0; c < means.size(); c++) {
//                fout << "(" << c << ") ";
//                for (std::size_t d = 0; d < dim; d++)
//                    fout << means[c][d] << (d == dim - 1 ? "" : ",");
//                fout << std::endl;
//            }
//            fout << std::endl;
//
//            for( std::size_t i = 0 ; i < items.size() ; i++ )
//                fout << i << ":" << items[i].cid() << std::endl;
//            fout.close();
//
//            // store means to prev_means and zeoring means
//            prev_means = means;
//            for( std::size_t i = 0 ; i < means.size() ; i++ )
//                means[i].clear();
//
//            // rearrange means and count # items for each cluster
//            std::vector<std::size_t> mean_counts(means.size(), 0);
//            for (item_list_type::iterator item_it = items.begin(); item_it != items.end(); ++item_it) {
//                item_list_type::value_type & item = *item_it;
//                std::transform(&item[0], &item[0] + dim, means[item.cid()].begin(), means[item.cid()].begin(), std::plus<float_type>());
//                ++mean_counts[item.cid()];
//            }
//
//            // averaging means to calculate centroids
//            for( std::size_t i = 0 ; i < means.size() ; i++ )
//                means[i] /= mean_counts[i];
//
//            iteration_count++;
//        }
//        //std::cout << "iteration count = " << iteration_count << std::endl;
//    }



    //typedef std::vector< subcluster_summary > subsum_vec_type;

//    struct subcluster_lessthan_norm
//    {
//        bool operator() ( const subcluster_summary& lhs, const subcluster_summary& rhs ) const { return (lhs.norm) < (rhs.norm); }
//    };

    //template<typename _iter>
    //template<cfentry_vec_type::iterator _iter>
   // void redist(std::vector<item_type<4824u> >::iterator begin, std::vector<item_type<4824u> >::iterator end, /*cfentry_vec_type*/std::vector<CFEntry<4824u> >& entries, std::vector<int>& out_cid)
	void redist(std::vector<item_type<4824u> >::iterator begin, std::vector<item_type<4824u> >::iterator end, /*cfentry_vec_type*/std::vector<CFEntry<4824u> >& entries, std::vector<int>& out_cid)
    {
        using namespace boost::numeric::ublas;

        // prepare summaries for each subcluster
        // summaries = ( center, radius, norm )
        subsum_vec_type subclusters;
        subclusters.reserve(entries.size());
        for (std::size_t i = 0; i < entries.size(); i++) {
            const CFEntry<dim>& e = entries[i];
            ublas_vec_type center(dim);
            std::copy(e.sum, e.sum + dim, center.begin());
            center /= e.n;
            subclusters.push_back(subcluster_summary(center, _Radius(e), std::sqrt(inner_prod(center, center))));
        }

        std::sort(subclusters.begin(), subclusters.end(), subcluster_lessthan_norm());

        // in addition to an individual summary for each subcluster
        // calculate pairwise euclidean distances of subclusters
        std::size_t n = subclusters.size();
        ublas_sym_matrix_type dist_mat(n, n);
        for (std::size_t i = 0; i < n - 1; i++) {
            for (std::size_t j = i + 1; j < n; j++) {
                ublas_vec_type diff = subclusters[i].center - subclusters[j].center;
                dist_mat(i, j) = inner_prod(diff, diff);
            }
        }

        out_cid.clear();
        out_cid.reserve(end - begin);
        for (std::vector<item_type<4824u> >::iterator it = begin; it != end; it++) {
            ublas_vec_type v(dim);
            std::copy(&(*it)[0], &(*it)[0] + dim, v.begin());
            out_cid.push_back(_redist(v, subclusters, dist_mat));
        }
    }

private:

    int _redist(ublas_vec_type& tmpv, subsum_vec_type& subsums, ublas_sym_matrix_type& dist_mat)
    {
        int imin, imax, i, k, n, start, end, median;
        float d, tmpnorm, idist, tmpdist;
        ublas_vec_type diff;

        i = 0;
        n = (int) subsums.size();
        tmpnorm = std::sqrt(inner_prod(tmpv, tmpv));

        // i=ClosestNorm(tmpnorm,norms,0,n-1);
        // for efficiency, replace recursion by iteration
        start = 0;
        end = n - 1;
        while (start < end) {
            if (end - start == 1) {
                float_type norm_end = subsums[end].norm;
                float_type norm_start = subsums[start].norm;

                i = tmpnorm > norm_end ? end : tmpnorm < norm_start ? start :
                    tmpnorm - norm_start < norm_end - tmpnorm ? start : end;
                start = end = i;
            } else {
                median = (start + end) / 2;
                float_type norm_med = subsums[median].norm;
                if (tmpnorm > norm_med)
                    start = median;
                else
                    end = median;
            }
        }

        diff = tmpv - subsums[i].center;
        idist= inner_prod(diff, diff);

        // imin=MinLargerThan(tmpnorm-sqrt(idist),norms,0,n-1);
        // imax=MaxSmallerThan(tmpnorm+sqrt(idist),norms,0,n-1);
        // for efficiency, replace recursion by iteration

        tmpdist = tmpnorm - sqrt(idist);
        start = 0;
        end = n - 1;
        while (start < end) {
            median = (start + end) / 2;
            float_type norm_med = subsums[median].norm;
            if (tmpdist > norm_med)
                start = median + 1;
            else
                end = median;
        }
        imin = start;

        tmpdist = tmpnorm + sqrt(idist);
        start = 0;
        end = n - 1;
        while (start < end) {
            median = (start + end + 1) / 2;
            float_type norm_med = subsums[median].norm;
            if (tmpdist < norm_med)
                end = median - 1;
            else
                start = median;
        }
        imax = start;

        // ClosestCenter(i,idist,tmpv,centers,imin,imax,matrix,n);
        // for efficiency, replace procedure by inline
        k = imin;
        while (k <= imax) {
            if (dist_mat(k, i) <= 4 * idist) {
                diff = tmpv - subsums[k].center;
                d = inner_prod(diff, diff);
                if (d < idist) {
                    idist = d;
                    i = k;
                }
            }
            k++;
        }
        return i;
    }
/* phase 4 - redistribute actual data points to subclusters */
//#include "CFTree_Redist.h"
public:
    void cluster(cfentry_vec_type& entries)
    {
        get_entries(entries);
        _cluster(entries);
    }

    private:

    struct HierarchicalClustering
    {
        typedef boost::numeric::ublas::symmetric_matrix<float_type> dist_matrix_type;

        HierarchicalClustering(int n, dist_func_type& in_dist_func)
                : size(n), step(-1), ii(n), jj(n), cf(n), dd(n), dist_func(in_dist_func), chain(n + 1), chainptr(-1), stopchain(FALSE)
        {
        }

        void merge(cfentry_vec_type& entries)
        {
            int nentry = (int) entries.size();
            int i, j, n1, n2;

            int CurI, PrevI, NextI;
            int uncheckcnt = nentry;

            std::vector<int> checked(nentry);
            for (i = 0; i < nentry; i++)
                checked[i] = i + 1;
            // 0: invalid after being merged to other entries
            // positive 1..nentry+1 :     original entries
            // negative -1..-(nentry-1) : merged entries

            // get initial distances
            // std::vector<float_type> dist(nentry*(nentry-1)/2);
            dist_matrix_type dist(nentry, nentry);

            for (i = 0; i < nentry - 1; i++)
                for (j = i + 1; j < nentry; j++)
                    dist(i, j) = dist_func(entries[i], entries[j]);

            CurI = rand() % nentry;         // step1
            chain[++chainptr] = CurI;

            while (uncheckcnt > 1) {
                // step4
                if (chainptr == -1) {
                    chainptr++;
                    chain[chainptr] = pick_one_unchecked(nentry, &checked[0]);
                }
                PrevI = chainptr > 0 ? chain[chainptr - 1] : -1;
                stopchain = FALSE;

                // step2
                while (stopchain == FALSE) {
                    CurI = chain[chainptr];
                    NextI = nearest_neighbor(CurI, nentry, &checked[0], dist);

                    // it is impossible NextI be -1 because uncheckcnt>1
                    if (NextI == PrevI)
                        stopchain = TRUE;
                    else {
                        chain[++chainptr] = NextI;
                        PrevI = CurI;
                    }
                } // end of while for step 2

                step++;

                // step3
                ii[step] = checked[CurI];
                jj[step] = checked[NextI];

                dd[step] = dist(CurI, NextI);

                bool curr_org = checked[CurI] > 0;
                bool next_org = checked[NextI] > 0;

                CFEntry<dim>& curr_entry = curr_org ? entries[CurI] : cf[-checked[CurI] - 1];
                CFEntry<dim>& next_entry = next_org ? entries[NextI] : cf[-checked[NextI] - 1];
                n1 = curr_entry.n;
                n2 = next_entry.n;
                cf[step] = curr_entry + next_entry;

                update_distance(n1, n2, CurI, NextI, nentry, &checked[0], dist);
                uncheckcnt--;
                checked[CurI] = -(step + 1);
                checked[NextI] = 0;
                chainptr--;
                chainptr--;
            } //end of while (uncheckcnt>1)

            // prepare for SplitHierarchy
            stopchain = FALSE;
            chainptr = 0;
            chain[chainptr] = -(step + 1);
        }

        void split(float_type ft)
        {
            while (_split(ft))
                /* nothing to do */;
        }

        short _split(float_type ft)
        {
            int i, j;

            if (chainptr == size)
                return FALSE;

            i = largest_merge(chainptr, ft);
            if (i != -1) {
                j = -chain[i] - 1;
                chain[i] = ii[j];
                chain[++chainptr] = jj[j];
                return TRUE;
            }

            stopchain = TRUE;
            return FALSE;
        }

        void result( /* inout */cfentry_vec_type& entries)
        {
            int j;
            std::vector<CFEntry<dim> > tmpentries(chainptr + 1);
            for (j = 0; j <= chainptr; j++) {
                tmpentries[j] = chain[j] < 0 ? cf[-chain[j] - 1] : entries[chain[j] - 1];
            }
            entries = tmpentries;
        }

        /* for SplitHierarchy use only */
        int farthest_merge(int chainptr)
        {
            if (chainptr <= 0)
                return chainptr;

            float d, dmax = 0;
            int i, imax = -1;
            for (i = 0; i <= chainptr; i++) {
                if (chain[i] < 0) {
                    d = dd[-chain[i] - 1];
                    if (d > dmax) {
                        imax = i;
                        dmax = d;
                    }
                }
            }
            return imax;
        }

        /* for SplitHierarchy use only */
        int largest_merge(int chainptr, float_type dist_threshold)
        {
            for (int i = 0; i <= chainptr; i++) {
                if (chain[i] < 0) {
                    if (_Diameter(cf[-chain[i] - 1]) > dist_threshold)
                        return i;
                }
            }
            return -1;
        }

        /* for MergeHierarchy use only */
        int nearest_neighbor(int CurI, int n, int *checked, dist_matrix_type& dist)
        {
            int imin = 0;
            float d, dmin = (std::numeric_limits<float_type>::max)();
            for (int i = 0; i < n; i++) {
                if (i == CurI || checked[i] == 0)
                    continue;

                d = dist(i, CurI);
                if (d < dmin) {
                    dmin = d;
                    imin = i;
                }
            }

            return dmin < (std::numeric_limits<float_type>::max)() ? imin : -1;
        }

        /* for MergeHierarchy use only */
        int pick_one_unchecked(int n, int *checked)
        {
            int i, j = rand() % n;
            for (i = 0; i < n; i++)
                if (checked[(i + j) % n] != 0)
                    return (i + j) % n;
            return -1;
        }

        /* for MergeHierarchy use only */
        void update_distance(int n1, int n2, int CurI, int NextI, int n, int *checked, dist_matrix_type& dist)
        {
            for (int i = 0; i < n; i++) {
                if (i == CurI || i == NextI)
                    continue;

                if (checked[i] != 0)
                    dist(i, CurI) = (n1 * dist(i, CurI) + n2 * dist(i, NextI)) / (n1 + n2);
            }
        }

        int size;
        int step;
        std::vector<int> ii;
        std::vector<int> jj;
        std::vector<CFEntry<dim> > cf;
        std::vector<float_type> dd;

        std::vector<int> chain;
        int chainptr;
        short stopchain;
        dist_func_type& dist_func;
    };

    void refine_cluster(cfentry_vec_type& entries)
    {
        std::vector<bool> merged(entries.size(), false);

        std::vector<std::size_t> not_visited;
        not_visited.reserve(entries.size());

        for (std::size_t i = 0; i < entries.size(); i++)
            not_visited.push_back(i);

        int temp = 0;
        while (!not_visited.empty()) {
            // pick any entry
            // std::cout << "iterations " << temp++ << std::endl;
            CFEntry<dim>& ref_entry = entries[not_visited.back()];
            CFEntry<dim> curr_entry = ref_entry;
            not_visited.pop_back();

            bool something_merged = false;
            for (std::size_t i = 0; i < not_visited.size(); i++) {
                // index of next item
                std::size_t v = not_visited[i];
                CFEntry<dim>& e = entries[v];
                if (dist_func(ref_entry, e) <= dist_threshold / 2) {
                    curr_entry += e;
                    merged[v] = true;
                    something_merged = true;
                }
            }

            if (something_merged)
                ref_entry = curr_entry;

            // remove if visited
            struct _remove_if_merged
            {
                _remove_if_merged(std::vector<bool>& in_visited)
                        : visited(in_visited)
                {
                }
                bool operator()(const std::size_t i) const
                {
                    return visited[i];
                }
            private:
                std::vector<bool>& visited;
            };

            // prepare for next loop, removing visited indices
            if (something_merged)
                not_visited.erase(std::remove_if(not_visited.begin(), not_visited.end(), _remove_if_merged(merged)),
                        not_visited.end());
        }

        struct _remove_if_merged_by_item
        {
            typedef CFEntry<dim> argument_type;
            _remove_if_merged_by_item(cfentry_vec_type& in_entries, std::vector<bool>& in_merged)
                    : entries(in_entries), merged(in_merged)
            {
            }
            bool operator()(const CFEntry<dim>& e) const
            {
                std::ptrdiff_t i = (&e - &entries[0]);
                return merged[i];
            }

        private:
            std::vector<bool>& merged;
            cfentry_vec_type& entries;
        };
        entries.erase(std::remove_if(entries.begin(), entries.end(), _remove_if_merged_by_item(entries, merged)),
                entries.end());
    }

    void _cluster(cfentry_vec_type& entries)
    {
        int n = (int) entries.size();

        if (n <= 1)
            return;

        if (dist_func == _DistD0 || dist_func == _DistD1
            || dist_func == _DistLarry) {
            refine_cluster(entries);
        } else {
            HierarchicalClustering h(n - 1, dist_func);
            h.merge(entries);
            h.split(dist_threshold);
            h.result(entries);
        }
    }
public:
    //static GeodesicDistance gd;
};





#endif
