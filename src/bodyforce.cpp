/**
 * solver.cpp
 *
 * by Mohammad Elmi
 * Created on 24 Oct 2013
 */

#include "bodyforce.h"

namespace aban2
{
void bodyforce::add_g(domain *d)
{
    for (int i = 0; i < d->n; ++i)
        d->ustar[1][i] += -7 * d->dt;
}
}
