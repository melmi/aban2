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
    if (vos->fullnesses[no] != fullness::half)return udt * d->vof[no];
    return vos->reconsts[i].get_flux(row->dir, udt, vector::from_data(d->nb, no));
}

void vof::advect_row(mesh_row *row)
{
    double flux_vof, flux_q[3];

    double *uf = d->extract_scalars(row, d->uf[row->dir]);
    double *row_mass = d->extract_scalars(row, vof_mass);

    for (size_t face = 0; face < d->n - 1; ++face)
    {
        size_t i = face, ii = face + 1;
        double udt = uf[i] * d->dt;

        if (udt > 0)
            flux_vof = get_vof_flux(row, i, udt);
        else
            flux_vof = get_vof_flux(row, ii, udt);

        row_mass[i ] -= flux_vof;
        row_mass[ii] += flux_vof;
    }

    auto startbc = d->boundaries[row->start_code];
    auto endbc   = d->boundaries[row->end_code  ];

    double uf_start = startbc->face_val(&bcondition::q, d->q[row->dir], row, bcside::start, row->dir) / d->_rho;
    double uf_end   = endbc  ->face_val(&bcondition::q, d->q[row->dir], row, bcside::end  , row->dir) / d->_rho;

    row_mass[0       ] += startbc->face_val(&bcondition::vof, d->vof, row, bcside::start, row->dir) * uf_start * d->dt;
    row_mass[d->n - 1] -= endbc  ->face_val(&bcondition::vof, d->vof, row, bcside::end  , row->dir) * uf_end   * d->dt;

    d->insert_scalars(row, vof_mass, row_mass);

    delete[] row_mass;
    delete[] uf;
}

vof::vof(aban2::domain *_d): d(_d)
{
    h3 = d->delta * d->delta * d->delta;
    vos = new voset(d);
    double *vof_mass = new double[d->n];
}

vof::~vof()
{
    delete vos;
    delete[] vof_mass;
}

void vof::advect()
{
    for (size_t i = 0; i < d->n; ++i) vof_mass[i] = d->vof[i] * h3;

    for (size_t dir = 0; dir < NDIRS; ++dir)
        for (size_t irow = 0; irow < d->nrows[dir]; ++irow)
            mesh_row *row = d->rows[dir] + irow;
}

} // aban2

