#include "vof.h"

#include <cmath>
#include <memory>
#include "gradient.h"
#include "volreconst.h"

namespace aban2
{

vof::vof(aban2::domain *_d): d(_d), _fsnorm(_d)
{
    vcell = d->delta * d->delta * d->delta;
    mass = new double[d->n];
    rhou[0] = new double[d->n];
    rhou[1] = new double[d->n];
    rhou[2] = new double[d->n];
    fullnesses = new fullness[d->n];
    on_interface = new bool[d->n];
    reconsts = new volreconst *[d->n];
    std::fill_n(reconsts, d->n, nullptr);
}

vof::~vof()
{
    delete[] mass;
    delete[] rhou[0];
    delete[] rhou[1];
    delete[] rhou[2];
    delete[] fullnesses;
    delete[] on_interface;
    delete_reconsts();
    delete[] reconsts;
}

void vof::set_fullnesses()
{
    for (size_t i = 0; i < d->ndir[0]; ++i)
        for (size_t j = 0; j < d->ndir[1]; ++j)
            for (size_t k = 0; k < d->ndir[2]; ++k)
            {
                size_t no;
                if (!d->exists(i, j, k, no)) continue;

                if (d->vof[no] < 1e-6)
                    fullnesses[no] = fullness::empty;
                else if (d->vof[no] > 1.0 - 1e-6)
                    fullnesses[no] = fullness::full;
                else
                    fullnesses[no] = fullness::half;
            }
}

bool vof::is_empty(size_t i, size_t j, size_t k)
{
    size_t no;
    if (d->exists_and_inside(i, j, k, no))
        return fullnesses[no] == fullness::empty;
    return false;
}

bool vof::is_on_interface(size_t i, size_t j, size_t k, size_t no)
{
    if (fullnesses[no] == fullness::half) return true;
    if (fullnesses[no] == fullness::empty) return false;

    // if full
    return
        is_empty(i + 1, j, k) ||
        is_empty(i - 1, j, k) ||
        is_empty(i, j + 1, k) ||
        is_empty(i, j - 1, k) ||
        is_empty(i, j, k + 1) ||
        is_empty(i, j, k - 1);
}

void vof::detect_interfacial_cells()
{
    for (size_t i = 0; i < d->ndir[0]; ++i)
        for (size_t j = 0; j < d->ndir[1]; ++j)
            for (size_t k = 0; k < d->ndir[2]; ++k)
            {
                size_t no;
                if (!d->exists(i, j, k, no)) continue;
                on_interface[no] = is_on_interface(i, j, k, no);
            }
}

void vof::create_reconsts()
{
    vector c {d->delta, d->delta, d->delta};

    for (int i = 0; i < d->n; ++i)
        if (on_interface[i])
            reconsts[i] = volreconst::from_volume(c, -vector::from_data(d->nb, i), mass[i]);
}

void vof::delete_reconsts()
{
    for (int i = 0; i < d->n; ++i)
        if (reconsts[i] != nullptr)
        {
            delete reconsts[i];
            reconsts[i] = nullptr;
        }
}

std::tuple<double, vector> vof::get_flux(mesh_row *row, size_t i, double udt, double ***grad_ustar)
{
    size_t no = d->cellno(row, i);
    double v0, v1, vt, c;
    vector q0, q1, qt;

    vt = d->delta * d->delta * std::abs(udt);
    c = ((udt > 0 ? 1.0 : -1.0) * d->delta - udt) / 2.0;
    qt.cmpnt[row->dir] = vt * c;

    if (on_interface[no])
        std::tie(v1, q1) = reconsts[no]->get_flux(row->dir, udt);
    else
    {
        if (fullnesses[no] == fullness::full)
        {
            v1 = vt;
            q1 = qt;
        }
        else
        {
            v1 = 0;
            // q1={0,0,0}; // default
        }
    }

    v0 = vt - v1;
    q0 = qt - q1;

    vector ustar = vector::from_data(d->ustar, no);
    vector grad[]
    {
        vector::from_data(grad_ustar[0], no),
        vector::from_data(grad_ustar[1], no),
#ifdef THREE_D
        vector::from_data(grad_ustar[2], no)
#else
        vector()
#endif
    };

    vector flux_ustar0 = d->rho0 * (ustar * v0 /*+ vector(q0 * grad[0], q0 * grad[1], q0 * grad[2])*/);
    vector flux_ustar1 = d->rho1 * (ustar * v1 /*+ vector(q1 * grad[0], q1 * grad[1], q1 * grad[2])*/);

    return std::make_tuple(v1, flux_ustar0 + flux_ustar1);
}

void vof::advect_row(mesh_row *row, double ***grad_ustar)
{
    std::tuple<double, vector> flux;
    size_t n = row->n;
    double a = d->delta * d->delta;
    double *row_uf = d->extract_scalars(row, d->uf[row->dir]);
    double *row_mass = d->extract_scalars(row, mass);
    double *row_rhou[3]
    {
        d->extract_scalars(row, rhou[0]),
        d->extract_scalars(row, rhou[1]),
        d->extract_scalars(row, rhou[2])
    };

    for (size_t face = 0; face < n - 1; ++face)
    {
        double udt = row_uf[face] * d->dt;
        size_t from, to;
        if (udt > 0)to = (from = face) + 1;
        else from = (to = face) + 1;

        flux = get_flux(row, from, udt, grad_ustar);

        row_mass[from] -= std::get<0>(flux);
        row_mass[to]   += std::get<0>(flux);

        for (int i = 0; i < 3; ++i)
        {
            row_rhou[i][from] -= std::get<1>(flux).cmpnt[i];
            row_rhou[i][to]   += std::get<1>(flux).cmpnt[i];
        }
    }

    double uf_start = flowbc::bc_u_getter[row->dir](d, row, d->u[row->dir], bcside::start);
    double uf_end   = flowbc::bc_u_getter[row->dir](d, row, d->u[row->dir], bcside::end);

    double vof_start = flowbc::bc_vof_getter(d, row, d->vof, bcside::start);
    double vof_end   = flowbc::bc_vof_getter(d, row, d->vof, bcside::end);

    mass[0    ] += vof_start * a * uf_start * d->dt;
    mass[n - 1] -= vof_end   * a * uf_end   * d->dt;

    double rhov_start = d->rho_bar(vof_start) * a * uf_start * d->dt;
    double rhov_end   = d->rho_bar(vof_end  ) * a * uf_end   * d->dt;
    for (int i = 0; i < 3; ++i)
    {
        rhou[i][0    ] += rhov_start * flowbc::bc_u_getter[i](d, row, d->u[i], bcside::start);
        rhou[i][n - 1] -= rhov_end   * flowbc::bc_u_getter[i](d, row, d->u[i], bcside::end);
    }

    d->insert_scalars(row, mass, row_mass);
    d->insert_scalars(row, rhou[0], row_rhou[0]);
    d->insert_scalars(row, rhou[1], row_rhou[1]);
    d->insert_scalars(row, rhou[2], row_rhou[2]);

    delete[] row_uf;
    delete[] row_mass;
    delete[] row_rhou[0];
    delete[] row_rhou[1];
    delete[] row_rhou[2];
}

void vof::correct_vofs(double *grad_uf_dir)
{
    // applying correction terms
    for (size_t i = 0; i < d->n; ++i)
        if (d->vof[i] > 0.5)
            mass[i] += grad_uf_dir[i] * vcell * d->dt;
}

void vof::calculate_normals()
{
    set_fullnesses();
    detect_interfacial_cells();

    size_t no;

    for (size_t i = 0; i < d->ndir[0]; ++i)
        for (size_t j = 0; j < d->ndir[1]; ++j)
            for (size_t k = 0; k < d->ndir[2]; ++k)
                if (d->exists(i, j, k, no))
                    if (on_interface[no])
                        _fsnorm.get_normal(i, j, k, no).to_data(d->nb, no);
}

void vof::calculate_masses_from_vars()
{
    for (size_t i = 0; i < d->n; ++i)
    {
        mass[i] = d->vof[i] * vcell;
        double rhov = d->rho[i] * vcell;
        rhou[0][i] = rhov * d->ustar[0][i];
        rhou[1][i] = rhov * d->ustar[1][i];
        rhou[2][i] = rhov * d->ustar[2][i];
    }
}

void vof::calculate_vars_from_masses()
{
    for (size_t i = 0; i < d->n; ++i)
    {
        d->vof[i] = mass[i] / vcell;
        d->rho[i] = d->rho_bar(d->vof[i]);
        double rhov = d->rho[i] * vcell;
        d->ustar[0][i] = rhou[0][i] / rhov;
        d->ustar[1][i] = rhou[1][i] / rhov;
        d->ustar[2][i] = rhou[2][i] / rhov;
    }
}

void vof::advect()
{
    start_dir = (start_dir + 1) % NDIRS;
    auto grad_uf = gradient::of_uf(d);

    for (size_t idir = 0; idir < NDIRS; ++idir)
    {
        calculate_normals();
        calculate_masses_from_vars();
        delete_reconsts();
        create_reconsts();

        size_t dir = (start_dir + idir) % NDIRS;

        double ***grad_ustar = gradient::of_vec(d, d->ustar, flowbc::bc_u_getter);

        for (size_t irow = 0; irow < d->nrows[dir]; ++irow)
            advect_row(d->rows[dir] + irow, grad_ustar);

        correct_vofs(grad_uf[dir][dir]);

        calculate_vars_from_masses();

        domain::delete_var(3, grad_ustar);
    }

    domain::delete_var(3, grad_uf);
}

int vof_err::is_inside(vector n, double alpha, vector x)
{
    return n * x > alpha;
}

double vof_err::err(vector n1, double alpha1, vector n2, double alpha2, double h)
{
    const int n = 101;
    double h_n = h / (n - 1);
    double result = 0;

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
        {
            vector p {h_n * i, h_n * j, 0};
            result += std::abs(is_inside(n1, alpha1, p) - is_inside(n2, alpha2, p));
        }
    result *= h * h / (n * n);
    return result;
}

double vof_err::compare(vof &v)
{
    double sum = 0;
    double n = 0;
    vector c {v.d->delta, v.d->delta, v.d->delta};
    double h3 = c.x * c.y * c.z;

    for (size_t i = 0; i < v.d->n; ++i)
        if (v.on_interface[i])
        {
            n += 1.0;

            vector n1 = -vector::from_data(v.d->nb, i);
            for (int dir = 0; dir < 3; ++dir)n1.cmpnt[dir] = std::abs(n1.cmpnt[dir]);
            auto vol = std::unique_ptr<volreconst>(volreconst::from_volume(c, n1, v.d->vof[i] * h3));
            double alpha1 = vol->alpha;

            vector n2 = vector::from_data(v.d->u, i);
            for (int dir = 0; dir < 3; ++dir)n2.cmpnt[dir] = std::abs(n2.cmpnt[dir]);
            double alpha2 = v.d->p[i];

            sum += err(n1, alpha1, n2, alpha2, v.d->delta);
        }

    return sum;
}

} // aban2




