#include "vof.h"

#include "diffusion.h"
#include "gradient.h"

aban2::vof::vof(aban2::domain *_d): d(_d)
{
}

aban2::vof::~vof()
{
}

void aban2::vof::calc_nbs()
{
    std::copy_n(d->vof, d->n, d->smooth_vof);
    for (size_t dir = 0; dir < NDIRS; ++dir)
        for (size_t irow = 0; irow < d->nrows[dir]; ++irow)
        {
            mesh_row *row = d->rows[dir] + irow;
            double *line = d->extract_scalars(row, d->smooth_vof);
            diffusion::diffuse(d, row, line, 1, &bcondition::voftype, &bcondition::vof, dir);
            d->insert_scalars(row, d->smooth_vof, line);
            delete[] line;
        }

    for (size_t dir = 0; dir < NDIRS; ++dir)
        for (size_t irow = 0; irow < d->nrows[dir]; ++irow)
        {
            mesh_row *row = d->rows[dir] + irow;
            double *line = d->extract_scalars(row, d->smooth_vof);
            double *grad = gradient::get_1d_row(d, row, line, &bcondition::vof, dir);
            d->insert_scalars(row, d->nb[dir], line);
            delete[] line;
            delete[] grad;
        }
}

void aban2::vof::calc_dists()
{
    // x^3-(x-a)^3+(x-b)^3 = a^3-3 a^2 x+3 a x^2-b^3+3 b^2 x-3 b x^2+x^3
    //                     = x^3 + (3a-3b) x^2 + (-3a^2+3b^2) x + (a^3-b^3)

    for (int i = 0; i < d->n; ++i)
    {
        vector n = vector::from_data(d->nb, d->cellnos[i]);
        double dmax = d->delta * (n.x + n.y + n.z);

    }
}

void aban2::vof::advect()
{

}
