#include "fsm.h"
#include <algorithm>
#include <cmath>

namespace aban2
{

void fsm::detect_on_interface_cells()
{
    std::fill_n(on_interface, d->n, false);

    //calculating fullnesses
    for (size_t i = 0; i < d->ndir[0]; ++i)
        for (size_t j = 0; j < d->ndir[1]; ++j)
            for (size_t k = 0; k < d->ndir[2]; ++k)
            {
                size_t ix;
                if (!d->exists(i, j, k, ix)) continue;
                fullnesses[ix] = get_fullness(ix);
            }

    //calculating on_interfaces
    for (size_t i = 0; i < d->ndir[0]; ++i)
        for (size_t j = 0; j < d->ndir[1]; ++j)
            for (size_t k = 0; k < d->ndir[2]; ++k)
            {
                size_t ix;
                if (!d->exists(i, j, k, ix)) continue;
                on_interface[ix] = is_on_interface(i, j, k, ix);
            }
}

bool fsm::is_on_interface(size_t i, size_t j, size_t k, size_t ix)
{
    if (fullnesses[ix] == fullness::half) return true;
    if (fullnesses[ix] == fullness::empty) return false;

    // if full
    size_t neighb_ix;
    if (d->exists_and_inside(i + 1, j, k, neighb_ix))if (fullnesses[neighb_ix] == fullness::empty) return true;
    if (d->exists_and_inside(i - 1, j, k, neighb_ix))if (fullnesses[neighb_ix] == fullness::empty) return true;
    if (d->exists_and_inside(i, j + 1, k, neighb_ix))if (fullnesses[neighb_ix] == fullness::empty) return true;
    if (d->exists_and_inside(i, j - 1, k, neighb_ix))if (fullnesses[neighb_ix] == fullness::empty) return true;
    if (d->exists_and_inside(i, j, k + 1, neighb_ix))if (fullnesses[neighb_ix] == fullness::empty) return true;
    if (d->exists_and_inside(i, j, k - 1, neighb_ix))if (fullnesses[neighb_ix] == fullness::empty) return true;
    return false;
}

fsm::fullness fsm::get_fullness(size_t ix)
{
    if (d->vof[ix] < 1e-6)return fullness::empty;
    if (d->vof[ix] > 1 - 1e-6)return fullness::full;
    return fullness::half;
}

void fsm::init_ls()
{
    for (size_t i = 0; i < d->n; ++i)
        if (!on_interface[i])
            if (fullnesses[i] == fullness::full)
                d->ls[i] = neginf;
            else
                d->ls[i] = posinf;
        else //for debug
            d->ls[i] = 0;
}

void fsm::set_new_phi(size_t i, size_t j, size_t k, size_t ix)
{
    double coeff = d->ls[ix] > 0 ? 1 : -1;
    double around[][] = {{posinf, posinf}, {posinf, posinf}, {posinf, posinf}};

    // around values
    size_t neighb_ix;
    if (d->exists_and_inside(i + 1, j, k, neighb_ix))around[0][0] = coeff * d->ls[neighb_ix];
    if (d->exists_and_inside(i - 1, j, k, neighb_ix))around[0][1] = coeff * d->ls[neighb_ix];
    if (d->exists_and_inside(i, j + 1, k, neighb_ix))around[1][0] = coeff * d->ls[neighb_ix];
    if (d->exists_and_inside(i, j - 1, k, neighb_ix))around[1][1] = coeff * d->ls[neighb_ix];
    if (d->exists_and_inside(i, j, k + 1, neighb_ix))around[2][0] = coeff * d->ls[neighb_ix];
    if (d->exists_and_inside(i, j, k - 1, neighb_ix))around[2][1] = coeff * d->ls[neighb_ix];

    double xmin =
    {
        std::min(around[0][0], around[0][1]),
        std::min(around[1][0], around[1][1]),
        std::min(around[2][0], around[2][1]),
        posinf
    };

    std::sort(xmin, xmin + 3); //we already know that posinf is the biggest value

    double xbar;
    for (int p = 0; p < NDIRS; ++p)
    {
        double a = 1, b = 0, c = -h2;
        for (int i = 0; i < p + 1; ++i)
        {
            b -= 2.0 * xmin[i];
            c += xmin[i] * xmin[i];
        }
        xbar = (-b + std::sqrt(b * b - 4.0 * a * c)) / (2.0 * a);
        if (xbar < xbar[p + 1]) break;
    }

    if (coeff * d->ls[ix] > xbar)d->ls[ix] = xbar;
}

void fsm::set_for_bounds(size_t dir, bool asc, size_t &start, int &limit, int &step)
{
    if (asc)
    {
        start = 0;
        limit = d->ndir[dir];
        step = +1;
    }
    else
    {
        start = d->ndir[dir] - 1;
        limit = -1;
        step = -1;
    }
}

double fsm::calculate_phis(bool iasc, bool jasc, bool kasc)
{
    size_t istart, jstart, kstart;
    int ilimit, jlimit, klimit;
    int istep, jstep, kstep;

    set_for_bounds(0, iasc, istart, ilimit, istep);
    set_for_bounds(1, jasc, jstart, jlimit, jstep);
    set_for_bounds(2, kasc, kstart, klimit, kstep);

    for (int i = istart; i != ilimit; i += istep)
        for (int j = jstart; j != jlimit; j += jstep)
            for (int k = kstart; k != klimit; k += kstep)
            {
                size_t ix;
                if (!d->exists(i, j, k, ix))continue;
                if (on_interface[ix])continue;
                set_new_phi(i, j, k, ix);
            }
}

fsm::fsm(domain *_d): d(_d)
{
    h2 = d->delta * d->delta;
    max_dist = max_dist_coeff * d->delta;
    posinf = std::sqrt((d->ndir[0] * d->ndir[0] + d->ndir[1] * d->ndir[1] + d->ndir[2] * d->ndir[2] + 1) * d->delta * d->delta);
    neginf = -posinf;
    on_interface = new bool[d->n];
    fullnesses = new fullness[d->n];
}

fsm::~fsm()
{
    delete[] on_interface;
    delete[] fullnesses;
}

void fsm::redist()
{
    detect_on_interface_cells();

    init_ls();

    calculate_phis(true, true, true);
    calculate_phis(false, true, true);

    calculate_phis(false, false, true);
    calculate_phis(true, false, true);

#ifdef THREE_D
    calculate_phis(true, false, false);
    calculate_phis(false, false, false);

    calculate_phis(false, true, false);
    calculate_phis(true, true, false);
#endif
}

} // aban2