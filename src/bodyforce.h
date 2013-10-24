/**
 * solver.h
 *
 * by Mohammad Elmi
 * Created on 23 Oct 2013
 */

#ifndef _BODYFORCE_H_
#define _BODYFORCE_H_

#include <algorithm>

#include "domain.h"

namespace aban2
{
class bodyforce
{
public:
    static void add_g(domain *d)
    {
        for (int i = 0; i < d->n; ++i)
            d->ustar[1][i] += -7 * d->dt;
    }
};

}

#endif