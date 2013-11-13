/**
 * vof.h
 *
 * by Mohammad Elmi
 * Created on 10 Nov 2013
 */

#ifndef _VOF_H_
#define _VOF_H_

#include "domain.h"

namespace aban2
{

class vof
{
    domain *d;
public:
    vof(domain *_d);
    ~vof();

    void calc_nbs(); 
    void calc_dists(); 
    void advect();
};

}

#endif