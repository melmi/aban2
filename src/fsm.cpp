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
                size_t no;
                if (!d->exists(i, j, k, no)) continue;
                fullnesses[no] = get_fullness(no);
            }

    //calculating on_interfaces
    for (size_t i = 0; i < d->ndir[0]; ++i)
        for (size_t j = 0; j < d->ndir[1]; ++j)
            for (size_t k = 0; k < d->ndir[2]; ++k)
            {
                size_t no;
                if (!d->exists(i, j, k, no)) continue;
                on_interface[no] = is_on_interface(i, j, k, no);
            }
}

bool fsm::is_on_interface(size_t i, size_t j, size_t k, size_t no)
{
    if (fullnesses[no] == fullness::half) return true;
    if (fullnesses[no] == fullness::empty) return false;

    // if full
    size_t neighb_no;
    if (d->exists_and_inside(i + 1, j, k, neighb_no))if (fullnesses[neighb_no] == fullness::empty) return true;
    if (d->exists_and_inside(i - 1, j, k, neighb_no))if (fullnesses[neighb_no] == fullness::empty) return true;
    if (d->exists_and_inside(i, j + 1, k, neighb_no))if (fullnesses[neighb_no] == fullness::empty) return true;
    if (d->exists_and_inside(i, j - 1, k, neighb_no))if (fullnesses[neighb_no] == fullness::empty) return true;
    if (d->exists_and_inside(i, j, k + 1, neighb_no))if (fullnesses[neighb_no] == fullness::empty) return true;
    if (d->exists_and_inside(i, j, k - 1, neighb_no))if (fullnesses[neighb_no] == fullness::empty) return true;
    return false;
}

fsm::fullness fsm::get_fullness(size_t no)
{
    if (d->vof[no] < 1e-6)return fullness::empty;
    if (d->vof[no] > 1.0 - 1e-6)return fullness::full;
    return fullness::half;
}

void fsm::init_ls()
{
    for (size_t i = 0; i < d->n; ++i)
        if (!on_interface[i])
            if (fullnesses[i] == fullness::full)
                d->ls[i] = -inf;
            else
                d->ls[i] = inf;
}

void fsm::set_new_phi(size_t i, size_t j, size_t k, size_t no)
{
    double coeff = d->ls[no] > 0 ? 1 : -1;
    double around[3][2] = {{inf, inf}, {inf, inf}, {inf, inf}};

    // around values
    size_t neighb_no;
    if (d->exists_and_inside(i + 1, j, k, neighb_no))around[0][0] = coeff * d->ls[neighb_no];
    if (d->exists_and_inside(i - 1, j, k, neighb_no))around[0][1] = coeff * d->ls[neighb_no];
    if (d->exists_and_inside(i, j + 1, k, neighb_no))around[1][0] = coeff * d->ls[neighb_no];
    if (d->exists_and_inside(i, j - 1, k, neighb_no))around[1][1] = coeff * d->ls[neighb_no];
    if (d->exists_and_inside(i, j, k + 1, neighb_no))around[2][0] = coeff * d->ls[neighb_no];
    if (d->exists_and_inside(i, j, k - 1, neighb_no))around[2][1] = coeff * d->ls[neighb_no];

    double xmin[] =
    {
        std::min(around[0][0], around[0][1]),
        std::min(around[1][0], around[1][1]),
        std::min(around[2][0], around[2][1]),
        inf
    };

    std::sort(xmin, xmin + NDIRS); //we already know that inf is the biggest value

    double xbar;
    for (int p = 0; p < NDIRS; ++p)
    {
        double a = 0, b = 0, c = -h2;
        for (int i = 0; i < p + 1; ++i)
        {
            a += 1.0;
            b += -2.0 * xmin[i];
            c += xmin[i] * xmin[i];
        }
        xbar = (-b + std::sqrt(b * b - 4.0 * a * c)) / (2.0 * a);

        if (xbar < xmin[p + 1]) break;
    }

    if (coeff * d->ls[no] > xbar)d->ls[no] = coeff * xbar;
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
                size_t no;
                if (!d->exists(i, j, k, no))continue;
                if (on_interface[no])continue;
                set_new_phi(i, j, k, no);
            }
}

fsm::fsm(domain *_d): d(_d)
{
    h2 = d->delta * d->delta;
    max_dist = max_dist_coeff * d->delta;
    inf = std::sqrt((d->ndir[0] * d->ndir[0] + d->ndir[1] * d->ndir[1] + d->ndir[2] * d->ndir[2] + 1) * d->delta * d->delta);
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