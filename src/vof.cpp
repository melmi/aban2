#include "vof.h"

#include "diffusion.h"
#include "gradient.h"
#include "ireconst.h"
#include "voset.h"

namespace aban2
{

double vof::get_vof_flux(mesh_row *row, size_t i, double udt)
{
    size_t no = d->cellno(row, i);
    // if (vos->fullnesses[no] != fullness::half)return std::abs(udt) * d->vof[no];
    auto remaining = vos->reconsts[no].get_remaining(row->dir, udt);
    auto result = vos->reconsts[no].volume - remaining.volume;
    vos->reconsts[no] = remaining;
    return result;
}

void vof::advect_row(mesh_row *row)
{
    double flux_vof;

    double *uf = d->extract_scalars(row, d->uf[row->dir]);
    double *row_mass = d->extract_scalars(row, vof_mass);

    for (size_t face = 0; face < row->n - 1; ++face)
    {
        double udt = uf[face] * d->dt;
        size_t from, to;
        if (udt > 0)to = (from = face) + 1; else from = (to = face) + 1;

        flux_vof = get_vof_flux(row, from, udt);

        // if (row_mass[from] - flux_vof < 0)flux_vof = h3 - row_mass[from];
        // if (row_mass[to] + flux_vof > h3)flux_vof = h3 - row_mass[to];

        row_mass[from] -= flux_vof;
        row_mass[to] += flux_vof;
    }

    // auto startbc = d->boundaries[row->start_code];
    // auto endbc   = d->boundaries[row->end_code  ];

    // double uf_start =
    //     startbc->face_val(&bcondition::q, d->q[row->dir], row, bcside::start, row->dir) / d->_rho;
    // double uf_end   =
    //     endbc  ->face_val(&bcondition::q, d->q[row->dir], row, bcside::end  , row->dir) / d->_rho;

    // row_mass[0         ] +=
    //     startbc->face_val(&bcondition::vof, d->vof, row, bcside::start, row->dir) * uf_start * d->dt;
    // row_mass[row->n - 1] -=
    //     endbc  ->face_val(&bcondition::vof, d->vof, row, bcside::end  , row->dir) * uf_end   * d->dt;

    d->insert_scalars(row, vof_mass, row_mass);

    delete[] row_mass;
    delete[] uf;
}

vof::vof(aban2::domain *_d): d(_d)
{
    h3 = d->delta * d->delta * d->delta;
    vos = new voset(d);
    vof_mass = new double[d->n];
}

vof::~vof()
{
    delete vos;
    delete[] vof_mass;
}

void vof::reconstruct_full_and_empty_cells()
{
    vector c {d->delta, d->delta, d->delta};
    // vector dummy_n {1, 0, 0};
    for (size_t i = 0; i < d->n; ++i)
        if (!vos->on_interface[i])
            if (vos->fullnesses[i] == fullness::full)
                vos->reconsts[i] = ireconst::from_full(c, vector::from_data(d->nb, i));
            else if (vos->fullnesses[i] == fullness::empty)
                vos->reconsts[i] = ireconst::from_empty(c, vector::from_data(d->nb, i));
}

void vof::advect()
{
    for (size_t i = 0; i < d->n; ++i) vof_mass[i] = d->vof[i] * h3;

    vos->make_ls();
    reconstruct_full_and_empty_cells();
    xfirst = !xfirst;

    if (xfirst)
        for (size_t dir = 0; dir < NDIRS; ++dir)
            for (size_t irow = 0; irow < d->nrows[dir]; ++irow)
            {
                mesh_row *row = d->rows[dir] + irow;
                advect_row(row);
            }
    else
        for (int dir = NDIRS - 1; dir > -1; --dir)
            for (size_t irow = 0; irow < d->nrows[dir]; ++irow)
            {
                mesh_row *row = d->rows[dir] + irow;
                advect_row(row);
            }
    vos->make_ls();

    for (size_t i = 0; i < d->n; ++i) d->vof[i] = vof_mass[i] / h3;
}

} // aban2

