/*
 * CFNodeLeaf.h
 *
 *  Created on: Feb 11, 2017
 *      Author: cotrik
 */

#ifndef CFNODELEAF_H_
#define CFNODELEAF_H_
#include "CFNode.h"
/** CFNode which is leaf */
template<boost::uint32_t dim>
struct CFNodeLeaf : public CFNode<dim>
{
    typedef boost::shared_ptr<CFNode<dim> > cfnode_sptr_type;
    CFNodeLeaf() : CFNode<dim>() {}
    virtual ~CFNodeLeaf(){}
    virtual bool IsLeaf() const { return true; }

    cfnode_sptr_type prev;  /** previous CFNode */
    cfnode_sptr_type next;  /** next CFNode */
};

#endif /* CFNODELEAF_H_ */
