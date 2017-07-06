#include "vof.h"

#include <cmath>
#include <memory>
#include "gradient.h"
#include "volreconst.h"

namespace aban2
{

    vof_t::vof_t(aban2::domain_t *_d): d(_d), _fsnorm(_d)
    {
        mass = new double[d->n];
        rhou0[0] = new double[d->n];
        rhou0[1] = new double[d->n];
        rhou0[2] = new double[d->n];
        rhou1[0] = new double[d->n];
        rhou1[1] = new double[d->n];
        rhou1[2] = new double[d->n];
        original_vof = new double[d->n];
        fullnesses = new fullness[d->n];
        on_interface = new bool[d->n];
        reconsts = new volreconst *[d->n];
        std::fill_n(reconsts, d->n, nullptr);
    }

    vof_t::~vof_t()
    {
        delete[] mass;
        delete[] rhou0[0];
        delete[] rhou0[1];
        delete[] rhou0[2];
        delete[] rhou1[0];
        delete[] rhou1[1];
        delete[] rhou1[2];
        delete[] original_vof;
        delete[] fullnesses;
        delete[] on_interface;
        delete_reconsts();
        delete[] reconsts;
    }

    void vof_t::set_fullnesses()
    {
        for (size_t i = 0; i < d->ndir[0]; ++i)
            for (size_t j = 0; j < d->ndir[1]; ++j)
                for (size_t k = 0; k < d->ndir[2]; ++k)
                {
                    size_t no;
                    if (!d->is_inside(i, j, k, no)) continue;

                    if (d->vof[no] < 1e-6)
                        fullnesses[no] = fullness::empty;
                    else if (d->vof[no] > 1.0 - 1e-6)
                        fullnesses[no] = fullness::full;
                    else
                        fullnesses[no] = fullness::half;
                }
    }

    bool vof_t::is_empty(size_t i, size_t j, size_t k, int di, int dj, int dk)
    {
        size_t no;
        if (d->is_inside(i, j, k, di, dj, dk, no))
            return fullnesses[no] == fullness::empty;
        return false;
    }

    bool vof_t::is_on_interface(size_t i, size_t j, size_t k, size_t no)
    {
        if (fullnesses[no] == fullness::half) return true;
        if (fullnesses[no] == fullness::empty) return false;

        // if full
        return
            is_empty(i, j, k, +1, 0, 0) ||
            is_empty(i, j, k, -1, 0, 0) ||
            is_empty(i, j, k, 0, +1, 0) ||
            is_empty(i, j, k, 0, -1, 0) ||
            is_empty(i, j, k, 0, 0, +1) ||
            is_empty(i, j, k, 0, 0, -1);
    }

    void vof_t::detect_interfacial_cells()
    {
        for (size_t i = 0; i < d->ndir[0]; ++i)
            for (size_t j = 0; j < d->ndir[1]; ++j)
                for (size_t k = 0; k < d->ndir[2]; ++k)
                {
                    size_t no;
                    if (!d->is_inside(i, j, k, no)) continue;
                    on_interface[no] = is_on_interface(i, j, k, no);
                }
    }

    void vof_t::create_reconsts()
    {
        vector c {d->delta, d->delta, d->delta};

        for (size_t i = 0; i < d->n; ++i)
            if (on_interface[i])
                reconsts[i] = volreconst::from_volume(c, -vector::from_data(d->nb, i), mass[i]);
    }

    void vof_t::delete_reconsts()
    {
        for (size_t i = 0; i < d->n; ++i)
            if (reconsts[i] != nullptr)
            {
                delete reconsts[i];
                reconsts[i] = nullptr;
            }
    }

    std::tuple<double, vector, vector> vof_t::get_flux(row_t *row, size_t i, double udt, double ***grad_ustar)
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

        vector flux_ustar0 = d->rho0 * (ustar * v0 + vector(q0 * grad[0], q0 * grad[1], q0 * grad[2]));
        vector flux_ustar1 = d->rho1 * (ustar * v1 + vector(q1 * grad[0], q1 * grad[1], q1 * grad[2]));

        // vector flux_ustar0 = d->rho0 * ustar * v0;
        // vector flux_ustar1 = d->rho1 * ustar * v1;

        return std::make_tuple(v1, flux_ustar0, flux_ustar1);
    }

    std::tuple<double, vector, vector> vof_t::get_bc_flux(row_t *row, double ***grad_ustar, bcside side)
    {
        size_t i = side == bcside::start ? 0 : row->n - 1;
        double uf_b = flowbc::bc_u_getter[row->dir](d, row, d->u[row->dir], side);
        if ((uf_b < 0 && side == bcside::start) || (uf_b > 0 && side == bcside::end))
            return get_flux(row, i, uf_b * d->dt, grad_ustar);

        double vof_b = flowbc::bc_vof_getter(d, row, d->vof, side);
        vector ustar_b
        {
            flowbc::bc_u_getter[0](d, row, d->ustar[0], side),
                flowbc::bc_u_getter[1](d, row, d->ustar[1], side),
                flowbc::bc_u_getter[2](d, row, d->ustar[2], side)
        };

        double v = d->aface * uf_b * d->dt;
        return std::make_tuple(
                vof_b * v,
                d->rho0 * vof_b * v * ustar_b,
                d->rho1 * (1.0 - vof_b) * v * ustar_b);
    }

    void vof_t::advect_row(row_t *row, double ***grad_ustar)
    {
        size_t n = row->n;
        double *row_uf = d->extract_scalars(row, d->uf[row->dir]);
        double *row_mass = d->extract_scalars(row, mass);
        double *row_rhou0[3]
        {
            d->extract_scalars(row, rhou0[0]),
                d->extract_scalars(row, rhou0[1]),
                d->extract_scalars(row, rhou0[2])
        };
        double *row_rhou1[3]
        {
            d->extract_scalars(row, rhou1[0]),
                d->extract_scalars(row, rhou1[1]),
                d->extract_scalars(row, rhou1[2])
        };

        for (size_t face = 0; face < n - 1; ++face)
        {
            double udt = row_uf[face] * d->dt;
            size_t from, to;
            if (udt > 0)to = (from = face) + 1;
            else from = (to = face) + 1;

            auto flux = get_flux(row, from, udt, grad_ustar);

            row_mass[from] -= std::get<0>(flux);
            row_mass[to]   += std::get<0>(flux);

            for (int i = 0; i < 3; ++i)
            {
                row_rhou0[i][from] -= std::get<1>(flux).cmpnt[i];
                row_rhou0[i][to]   += std::get<1>(flux).cmpnt[i];

                row_rhou1[i][from] -= std::get<2>(flux).cmpnt[i];
                row_rhou1[i][to]   += std::get<2>(flux).cmpnt[i];
            }
        }

        auto flux_start = get_bc_flux(row, grad_ustar, bcside::start);
        auto flux_end   = get_bc_flux(row, grad_ustar, bcside::end);

        mass[0    ] += std::get<0>(flux_start);
        mass[n - 1] -= std::get<0>(flux_end  );

        for (int i = 0; i < 3; ++i)
        {
            rhou0[i][0    ] += std::get<1>(flux_start).cmpnt[i];
            rhou0[i][n - 1] -= std::get<1>(flux_end  ).cmpnt[i];

            rhou1[i][0    ] += std::get<2>(flux_start).cmpnt[i];
            rhou1[i][n - 1] -= std::get<2>(flux_end  ).cmpnt[i];
        }

        d->insert_scalars(row, mass, row_mass);
        d->insert_scalars(row, rhou0[0], row_rhou0[0]);
        d->insert_scalars(row, rhou0[1], row_rhou0[1]);
        d->insert_scalars(row, rhou0[2], row_rhou0[2]);
        d->insert_scalars(row, rhou1[0], row_rhou1[0]);
        d->insert_scalars(row, rhou1[1], row_rhou1[1]);
        d->insert_scalars(row, rhou1[2], row_rhou1[2]);

        delete[] row_uf;
        delete[] row_mass;
        delete[] row_rhou0[0];
        delete[] row_rhou0[1];
        delete[] row_rhou0[2];
        delete[] row_rhou1[0];
        delete[] row_rhou1[1];
        delete[] row_rhou1[2];
    }

    void vof_t::correct_vofs(size_t dir)
    {
        auto grad_uf = gradient::of_uf_dir(d, dir);

        // applying correction terms
        for (size_t i = 0; i < d->n; ++i)
            if (original_vof[i] > 0.5)
                mass[i] += grad_uf[i] * d->vcell * d->dt;

        delete[] grad_uf;
    }

    void vof_t::calculate_normals()
    {
        size_t no;

        for (size_t i = 0; i < d->ndir[0]; ++i)
            for (size_t j = 0; j < d->ndir[1]; ++j)
                for (size_t k = 0; k < d->ndir[2]; ++k)
                    if (d->is_inside(i, j, k, no))
                        if (on_interface[no])
                            _fsnorm.get_normal(i, j, k, no).to_data(d->nb, no);
    }

    void vof_t::calculate_vof_masses_from_vars()
    {
        for (size_t i = 0; i < d->n; ++i)
            mass[i] = d->vof[i] * d->vcell;
    }

    void vof_t::calculate_ustar_masses_from_vars(double ***grad_ustar)
    {
        for (size_t i = 0; i < d->n; ++i)
        {
            vector rhou0p, rhou1p;
            vector ustar = vector::from_data(d->ustar, i);
            if (fullnesses[i] == fullness::half)
            {
                double v1 = reconsts[i]->get_volume();
                vector q1 = reconsts[i]->get_moments();
                double v0 = d->vcell - v1;
                vector q0 = -q1;

                vector grad[]
                {
                    vector::from_data(grad_ustar[0], i),
                        vector::from_data(grad_ustar[1], i),
#ifdef THREE_D
                        vector::from_data(grad_ustar[2], i)
#else
                            vector()
#endif
                };

                rhou0p = d->rho0 * (ustar * v0 + vector(q0 * grad[0], q0 * grad[1], q0 * grad[2]));
                rhou1p = d->rho1 * (ustar * v1 + vector(q1 * grad[0], q1 * grad[1], q1 * grad[2]));
            }
            else
            {
                double _vof = d->vof[i];
                rhou0p = (1.0 - _vof) * d->rho0 * ustar * d->vcell;
                rhou1p = _vof * d->rho1 * ustar * d->vcell;
            }

            rhou0p.to_data(rhou0, i);
            rhou1p.to_data(rhou1, i);
        }
    }

    void vof_t::calculate_vars_from_masses()
    {
        for (size_t i = 0; i < d->n; ++i)
        {
            d->vof[i] = mass[i] / d->vcell;
            d->rho[i] = d->rho_bar(d->vof[i]);
            d->ustar[0][i] = (rhou0[0][i] / d->rho0 + rhou1[0][i] / d->rho1) / d->vcell;
            d->ustar[1][i] = (rhou0[1][i] / d->rho0 + rhou1[1][i] / d->rho1) / d->vcell;
            d->ustar[2][i] = (rhou0[2][i] / d->rho0 + rhou1[2][i] / d->rho1) / d->vcell;
        }
    }

    void vof_t::advect()
    {
        std::copy_n(d->vof, d->n, original_vof);
        start_dir = (start_dir + 1) % NDIRS;

        for (size_t idir = 0; idir < NDIRS; ++idir)
        {
            set_fullnesses();
            detect_interfacial_cells();
            calculate_normals();

            // reconstruct
            calculate_vof_masses_from_vars();
            create_reconsts();
            auto grad_ustar = gradient::of_vec(d, d->ustar, flowbc::bc_u_getter);
            calculate_ustar_masses_from_vars(grad_ustar);

            // evolve
            size_t dir = (start_dir + idir) % NDIRS;

            for (size_t irow = 0; irow < d->nrows[dir]; ++irow)
                advect_row(d->rows[dir] + irow, grad_ustar);

            correct_vofs(dir);

            // average
            calculate_vars_from_masses();

            delete_reconsts();
            domain_t::delete_var(3, grad_ustar);
        }
    }

    int vof_err_t::is_inside(vector n, double alpha, vector x)
    {
        return n * x > alpha;
    }

    double vof_err_t::err(vector n1, double alpha1, vector n2, double alpha2, double h)
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

    double vof_err_t::compare(vof_t &v)
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
                double alpha1 = vol->get_alpha();

                vector n2 = vector::from_data(v.d->u, i);
                for (int dir = 0; dir < 3; ++dir)n2.cmpnt[dir] = std::abs(n2.cmpnt[dir]);
                double alpha2 = v.d->p[i];

                sum += err(n1, alpha1, n2, alpha2, v.d->delta);
            }

        return sum;
    }

} // aban2
