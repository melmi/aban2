#include "vof.h"

#include <cmath>
#include <memory>
#include "gradient.h"
#include "volreconst.h"

namespace aban2
{

vof::vof(aban2::domain *_d): d(_d)
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

vof::neighbs_t vof::get_nighb_vals(size_t i, size_t j, size_t k)
{
    neighbs_t result;
    size_t no;

    for (int ii = 0; ii < 3; ++ii)
        for (int jj = 0; jj < 3; ++jj)
        #ifdef THREE_D
            for (int kk = 0; kk < 3; ++kk)
                if (d->exists_and_inside(i + ii - 1, j + jj - 1, k + kk - 1, no))
                    result[ii][jj][kk] = d->vof[no];
                else
                    result[ii][jj][kk] = -1;
        #else
            if (d->exists_and_inside(i + ii - 1, j + jj - 1, k, no))
                result[ii][jj][0] = result[ii][jj][1] = result[ii][jj][2] = d->vof[no];
            else
                result[ii][jj][0] = result[ii][jj][1] = result[ii][jj][2] = -1;
        #endif

    return result;
}

double vof::col_sum(double delta, double v1, double v2, double v3)
{
    if (v1 < -0.5 || v2 < -0.5 || v3 < -0.5)return -1;
    return delta * (v1 + v2 + v3);
}

double vof::dir_sign(double v1, double v2, double v3)
{
    if (v1 > -epsilon && v3 > -epsilon) return (v3 - v1) / (std::abs(v3 - v1) + epsilon);
    if (v1 > -epsilon && v2 > -epsilon) return (v2 - v1) / (std::abs(v2 - v1) + epsilon);
    if (v3 > -epsilon && v2 > -epsilon) return (v3 - v2) / (std::abs(v3 - v2) + epsilon);
    return 0;
}

double vof::get_column_grad(double p, double f, double b, double delta)
{
    if (f < -0.5) return (p - b) / delta;
    if (b < -0.5) return (f - p) / delta;

    double gf = (f - p);// / delta;
    double gb = (p - b);// / delta;

    // return std::abs(gb) > std::abs(gf) ? gb : gf; // Zalski
    // return 0.5 * (gf + bf); // CC

    double yf = p + 1.5 * delta * gb;
    double yb = p - 1.5 * delta * gf;

    bool cf = 0 <= yf && yf <= 3;
    bool cb = 0 <= yb && yb <= 3;

    if (cf && !cb) return gf;
    if (cb && !cf) return gb;

    return 0.5 * (gf + gb);
}

bool vof::create_column_candidate(vector &v, size_t dir, double sign, molecule_t molecule)
{
    if (molecule.p < -0.5 || std::abs(sign) < 0.5 ||
            (molecule.f1 < -0.5 && molecule.b1 < -0.5) ||
            (molecule.f2 < -0.5 && molecule.b2 < -0.5)) return false;

    size_t dir1 = (dir + 1) % 3, dir2 = (dir + 2) % 3;

    v.cmpnt[dir] = sign;
    v.cmpnt[dir1] = get_column_grad(molecule.p, molecule.f1, molecule.b1, d->delta);
    v.cmpnt[dir2] = get_column_grad(molecule.p, molecule.f2, molecule.b2, d->delta);

    return true;
}

void vof::create_column_candidates(const neighbs_t &n, vector *candidates, bool *ok)
{
    double center[3], dir1f[3], dir1b[3], dir2f[3], dir2b[3], sign[3];

    sign[0] = dir_sign(n[0][1][1], n[1][1][1], n[2][1][1]);
    sign[1] = dir_sign(n[1][0][1], n[1][1][1], n[1][2][1]);
    sign[2] = dir_sign(n[1][1][0], n[1][1][1], n[1][1][2]);

    center[0] = col_sum(d->delta, n[0][1][1], n[1][1][1], n[2][1][1]);
    center[1] = col_sum(d->delta, n[1][0][1], n[1][1][1], n[1][2][1]);
    center[2] = col_sum(d->delta, n[1][1][0], n[1][1][1], n[1][1][2]);

    dir1f[0] = col_sum(d->delta, n[0][2][1], n[1][2][1], n[2][2][1]);
    dir1f[1] = col_sum(d->delta, n[1][0][2], n[1][1][2], n[1][2][2]);
    dir1f[2] = col_sum(d->delta, n[2][1][0], n[2][1][1], n[2][1][2]);

    dir1b[0] = col_sum(d->delta, n[0][0][1], n[1][0][1], n[2][0][1]);
    dir1b[1] = col_sum(d->delta, n[1][0][0], n[1][1][0], n[1][2][0]);
    dir1b[2] = col_sum(d->delta, n[0][1][0], n[0][1][1], n[0][1][2]);

    dir2f[0] = col_sum(d->delta, n[0][1][2], n[1][1][2], n[2][1][2]);
    dir2f[1] = col_sum(d->delta, n[2][0][1], n[2][1][1], n[2][2][1]);
    dir2f[2] = col_sum(d->delta, n[1][2][0], n[1][2][1], n[1][2][2]);

    dir2b[0] = col_sum(d->delta, n[0][1][0], n[1][1][0], n[2][1][0]);
    dir2b[1] = col_sum(d->delta, n[0][0][1], n[0][1][1], n[0][2][1]);
    dir2b[2] = col_sum(d->delta, n[1][0][0], n[1][0][1], n[1][0][2]);

    for (size_t i = 0; i < 3; ++i)
        ok[i] = create_column_candidate(candidates[i], i, sign[i], {center[i], dir1f[i], dir1b[i], dir2f[i], dir2b[i]});
}

bool vof::create_young_candidate(const neighbs_t &n, vector &candidate)
{
    for (int ii = 0; ii < 3; ++ii)
        for (int jj = 0; jj < 3; ++jj)
            for (int kk = 0; kk < 3; ++kk)
                if (n[ii][jj][kk] < -0.5)return false;

    candidate.x = (
                      (1.0 * (n[2][0][0] + n[2][0][2] + n[2][2][0] + n[2][2][2]) +
                       2.0 * (n[2][1][0] + n[2][1][2] + n[2][0][1] + n[2][2][1]) +
                       4.0 * (n[2][1][1])) -
                      (1.0 * (n[0][0][0] + n[0][0][2] + n[0][2][0] + n[0][2][2]) +
                       2.0 * (n[0][1][0] + n[0][1][2] + n[0][0][1] + n[0][2][1]) +
                       4.0 * (n[0][1][1]))
                  ) / 2. / d->delta;

    candidate.y = (
                      (1.0 * (n[0][2][0] + n[0][2][2] + n[2][2][0] + n[2][2][2]) +
                       2.0 * (n[1][2][0] + n[1][2][2] + n[0][2][1] + n[2][2][1]) +
                       4.0 * (n[1][2][1])) -
                      (1.0 * (n[0][0][0] + n[0][0][2] + n[2][0][0] + n[2][0][2]) +
                       2.0 * (n[1][0][0] + n[1][0][2] + n[0][0][1] + n[2][0][1]) +
                       4.0 * (n[1][0][1]))
                  ) / 2. / d->delta;

    candidate.z = (
                      (1.0 * (n[0][0][2] + n[0][2][2] + n[2][0][2] + n[2][2][2]) +
                       2.0 * (n[1][0][2] + n[1][2][2] + n[0][1][2] + n[2][1][2]) +
                       4.0 * (n[1][1][2])) -
                      (1.0 * (n[0][0][0] + n[0][2][0] + n[2][0][0] + n[2][2][0]) +
                       2.0 * (n[1][0][0] + n[1][2][0] + n[0][1][0] + n[2][1][0]) +
                       4.0 * (n[1][1][0]))
                  ) / 2. / d->delta;

    return true;
}

void vof::relax_center_val(int i, int j, int k, neighbs_t &n)
{
    if (n[i][j][k] < -0.5)n[i][j][k] = n[1][1][1];
}

void vof::relax_edge_val(int i, int j, int k, neighbs_t &n)
{
    if (n[i][j][k] >= -0.5)return;

    int i1 = i, i2 = i, j1 = j, j2 = j, k1 = k, k2 = k;

    if (i == 1) {j1 = (j + 1) % 3; k2 = (k + 1) % 3;}
    if (j == 1) {i1 = (i + 1) % 3; k2 = (k + 1) % 3;}
    if (k == 1) {i1 = (i + 1) % 3; j2 = (j + 1) % 3;}

    n[i][j][k] = 0.5 * (n[i1][j1][k1] + n[i2][j2][k2]);
}

void vof::relax_corner_val(int i, int j, int k, neighbs_t &n)
{
    if (n[i][j][k] >= -0.5)return;

    int ii = (i + 1) % 3;
    int jj = (j + 1) % 3;
    int kk = (k + 1) % 3;

    n[i][j][k] = (n[ii][j][k] + n[i][jj][k] + n[i][j][kk]) / 3.0;
}

void vof::relax_neighb_vals(neighbs_t &n)
{
    relax_center_val(0, 1, 1, n);
    relax_center_val(2, 1, 1, n);
    relax_center_val(1, 0, 1, n);
    relax_center_val(1, 2, 1, n);
    relax_center_val(1, 1, 0, n);
    relax_center_val(1, 1, 2, n);

    relax_edge_val(1, 0, 0, n);
    relax_edge_val(1, 0, 2, n);
    relax_edge_val(1, 2, 0, n);
    relax_edge_val(1, 2, 2, n);
    relax_edge_val(0, 1, 0, n);
    relax_edge_val(0, 1, 2, n);
    relax_edge_val(2, 1, 0, n);
    relax_edge_val(2, 1, 2, n);
    relax_edge_val(0, 0, 1, n);
    relax_edge_val(0, 2, 1, n);
    relax_edge_val(2, 0, 1, n);
    relax_edge_val(2, 2, 1, n);

    relax_corner_val(0, 0, 0, n);
    relax_corner_val(0, 0, 2, n);
    relax_corner_val(0, 2, 0, n);
    relax_corner_val(0, 2, 2, n);
    relax_corner_val(2, 0, 0, n);
    relax_corner_val(2, 0, 2, n);
    relax_corner_val(2, 2, 0, n);
    relax_corner_val(2, 2, 2, n);

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                std::cout << i << " " << j << " " << k << " " << n[i][j][k] << std::endl;
}

inline vector taxicab_normalized(const vector &v)
{
    double size = std::abs(v.x) + std::abs(v.y) + std::abs(v.z);
    return v * (1.0 / size);
}

void vof::set_normal(size_t i, size_t j, size_t k, size_t no)
{
    if (!on_interface[no]) return;

    vector y; // young's normal vector
    // bool y_ok; // young's normal ok
    vector c[3]; // column normal vector
    bool c_ok[3];
    vector *final;

    auto neighb_vals = get_nighb_vals(i, j, k);
    create_column_candidates(neighb_vals, c, c_ok);
    // y_ok = create_young_candidate(neighb_vals, y);

    if (c_ok[0] || c_ok[1] || c_ok[2])
    {
        double m0[] {0, 0, 0};
        size_t dir = 0;
        for (size_t i = 0; i < 3; ++i)
            if (c_ok[i])
            {
                m0[i] = std::abs(taxicab_normalized(c[i]).cmpnt[i]);
                if (m0[i] > m0[dir])dir = i;
            }

        // if (y_ok)
        //     if (std::abs(taxicab_normalized(y).cmpnt[dir]) < m0[dir])
        //         final = &y;
        //     else
        //         final = c + dir;
        // else
        final = c + dir;
    }
    // else if (y_ok) // Young helper
    //     final = &y;
    else
    {
        relax_neighb_vals(neighb_vals);
        create_young_candidate(neighb_vals, y);
        final = &y;
    }

    final->normalize();
    final->to_data(d->nb, no);
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
        if (udt > 0)to = (from = face) + 1; else from = (to = face) + 1;

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
                    set_normal(i, j, k, no);
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
        delete_reconsts();
        create_reconsts();
        calculate_masses_from_vars();

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

