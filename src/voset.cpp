#include "voset.h"
#include "fsm.h"
#include "ireconst.h"
#include "gradient.h"

#include <algorithm>
#include <cmath>

namespace aban2
{
void voset::detect_on_interface_cells()
{
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

bool voset::is_on_interface(size_t i, size_t j, size_t k, size_t no)
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

fullness voset::get_fullness(size_t no)
{
    if (d->vof[no] < 1e-6)return fullness::empty;
    if (d->vof[no] > 1.0 - 1e-6)return fullness::full;
    return fullness::half;
}

double voset::err()
{
    double result = 0;
    for (size_t i = 0; i < d->n; ++i)
        if (on_interface[i])
            result = std::max(result, std::abs(d->ls[i] - old_ls[i]));
    return result;
}

double voset::calculate_normals(double *phi)
{
    for (size_t dir = 0; dir < NDIRS; ++dir)
        for (size_t irow = 0; irow < d->nrows[dir]; ++irow)
        {
            mesh_row *row = d->rows[dir] + irow;
            double *line = d->extract_scalars(row, phi);
            double *grad = gradient::get_1d_row(d, row, line, &bcondition::vof, dir);
            d->insert_scalars(row, d->nb[dir], grad);
            delete[] line;
            delete[] grad;
        }

    for (size_t i = 0; i < d->n; ++i)
    {
        auto n = vector::from_data(d->nb, i);
        n.normalize(1e-6);

        d->nb[0][i] = n.x;
        d->nb[1][i] = n.y;
        d->nb[2][i] = n.z;
    }
}

void voset::calculate_interface_distances()
{
    for (size_t i = 0; i < d->n; ++i)
        if (on_interface[i])
        {
            auto n = vector::from_data(d->nb, i);
            reconsts[i] = ireconst::from_volume(celldims, n, h3 * d->vof[i]);
            d->ls[i] = 0.5 * (celldims * reconsts[i].m) - reconsts[i].alpha;
        }
}

voset::voset(domain *_d): d(_d)
{
    h3 = d->delta * d->delta * d->delta;
    celldims = {d->delta, d->delta, d->delta};
    old_ls = new double[d->n];
    on_interface = new bool[d->n];
    fullnesses = new fullness[d->n];
    reconsts = new ireconst[d->n];
    redistancer = new fsm(d, on_interface, fullnesses);
}

voset::~voset()
{
    delete[] old_ls;
    delete[] on_interface;
    delete[] fullnesses;
    delete[] reconsts;
    delete redistancer;
}

void voset::make_ls()
{
    detect_on_interface_cells();

    std::fill_n(old_ls, d->n, 0);
    calculate_normals(d->vof);
    calculate_interface_distances();
    redistancer->redist();
    while (err() > 1e-6)
    {
        std::copy_n(d->ls, d->n, old_ls);
        calculate_normals(d->ls);
        calculate_interface_distances();
        redistancer->redist();
    }
}

} // aban2