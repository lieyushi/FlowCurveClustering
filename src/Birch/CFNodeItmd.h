/*
 * CFNodeItmd.h
 *
 *  Created on: Feb 11, 2017
 *      Author: cotrik
 */

#ifndef CFNODEITMD_H_
#define CFNODEITMD_H_

#include "CFNode.h"
/** CFNode which is intermediate */
template<boost::uint32_t dim>
struct CFNodeItmd : public CFNode<dim>
{
    CFNodeItmd() : CFNode<dim>() {}
    virtual ~CFNodeItmd() {};
    virtual bool                IsLeaf() const { return false; }
};


#endif /* CFNODEITMD_H_ */
