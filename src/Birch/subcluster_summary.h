/*
 * subcluster_summary.h
 *
 *  Created on: Feb 11, 2017
 *      Author: cotrik
 */

#ifndef SUBCLUSTER_SUMMARY_H_
#define SUBCLUSTER_SUMMARY_H_

/************************************************************************/
/* The original redistribution code of birch
/* In my view point, it could be burdensome due to O(n^2) cost
/************************************************************************/
struct subcluster_summary
{
    subcluster_summary( ) : radius(0.0), norm(0.0) {}
    subcluster_summary( const ublas_vec_type& in_center, 
    					const float_type& in_radius, 
    					const float_type& in_norm ) : 
    					center( in_center ), 
    					radius(in_radius), 
    					norm(in_norm) {}
    //subcluster_summary( float_type* in_center, float_type& in_radius, float_type& in_norm ) : center( in_center, in_center + dim ), radius(in_radius), norm(in_norm) {}

    ublas_vec_type  center;
    float_type      radius;
    float_type      norm;
};

struct subcluster_lessthan_norm
{
    bool operator() ( const subcluster_summary& lhs, const subcluster_summary& rhs ) const { return (lhs.norm) < (rhs.norm); }
};

#endif /* SUBCLUSTER_SUMMARY_H_ */
