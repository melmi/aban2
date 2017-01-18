#include "vof.h"

#include <cmath>
#include <memory>
#include "gradient.h"
#include "volreconst.h"

namespace aban2
{

fsnorm::fsnorm(aban2::domain_t *_d): d(_d)
{
}

fsnorm::~fsnorm()
{
}

fsnorm::neighbs_t fsnorm::get_nighb_vals(size_t i, size_t j, size_t k)
{
    neighbs_t result;
    size_t no;

    for (int ii = 0; ii < 3; ++ii)
        for (int jj = 0; jj < 3; ++jj)
#ifdef THREE_D
            for (int kk = 0; kk < 3; ++kk)
                if (d->exists_and_inside(i, j, k, ii - 1, jj - 1, kk - 1, no))
                    result[ii][jj][kk] = d->vof[no];
                else
                    result[ii][jj][kk] = -1;
#else
            if (d->is_inside(i, j, k, ii - 1, jj - 1, k, no))
                result[ii][jj][0] = result[ii][jj][1] = result[ii][jj][2] = d->vof[no];
            else
                result[ii][jj][0] = result[ii][jj][1] = result[ii][jj][2] = -1;
#endif

    return result;
}

double fsnorm::col_sum(double delta, double v1, double v2, double v3)
{
    if (v1 < -0.5 || v2 < -0.5 || v3 < -0.5)return -1;
    return delta * (v1 + v2 + v3);
}

double fsnorm::dir_sign(double v1, double v2, double v3)
{
    if (v1 > -epsilon && v3 > -epsilon) return (v3 - v1) / (std::abs(v3 - v1) + epsilon);
    if (v1 > -epsilon && v2 > -epsilon) return (v2 - v1) / (std::abs(v2 - v1) + epsilon);
    if (v3 > -epsilon && v2 > -epsilon) return (v3 - v2) / (std::abs(v3 - v2) + epsilon);
    return 0;
}

double fsnorm::get_column_grad(double p, double f, double b, double delta)
{
    if (f < -0.5) return (p - b) / delta;
    if (b < -0.5) return (f - p) / delta;

    double gf = (f - p) / delta;
    double gb = (p - b) / delta;

    // return std::abs(gb) > std::abs(gf) ? gb : gf; // Zalski
    // return 0.5 * (gf + bf); // CC

    double yf = p + 1.5 * delta * gb;
    double yb = p - 1.5 * delta * gf;

    bool cf = 0 <= yf && yf <= 3 * delta;
    bool cb = 0 <= yb && yb <= 3 * delta;

    if (cf && !cb) return gf;
    if (cb && !cf) return gb;

    return 0.5 * (gf + gb);
}

bool fsnorm::create_column_candidate(vector &v, size_t dir, double sign, molecule_t molecule)
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

void fsnorm::create_column_candidates(const neighbs_t &n, vector *candidates, bool *ok)
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

bool fsnorm::create_young_candidate(const neighbs_t &n, vector &candidate)
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

void fsnorm::relax_center_val(int i, int j, int k, neighbs_t &n)
{
    if (n[i][j][k] < -0.5)n[i][j][k] = n[1][1][1];
}

void fsnorm::relax_edge_val(int i, int j, int k, neighbs_t &n)
{
    if (n[i][j][k] >= -0.5)return;

    int i1 = i, i2 = i, j1 = j, j2 = j, k1 = k, k2 = k;

    if (i == 1)
    {
        j1 = (j + 1) % 3;
        k2 = (k + 1) % 3;
    }
    if (j == 1)
    {
        i1 = (i + 1) % 3;
        k2 = (k + 1) % 3;
    }
    if (k == 1)
    {
        i1 = (i + 1) % 3;
        j2 = (j + 1) % 3;
    }

    n[i][j][k] = 0.5 * (n[i1][j1][k1] + n[i2][j2][k2]);
}

void fsnorm::relax_corner_val(int i, int j, int k, neighbs_t &n)
{
    if (n[i][j][k] >= -0.5)return;

    int ii = (i + 1) % 3;
    int jj = (j + 1) % 3;
    int kk = (k + 1) % 3;

    n[i][j][k] = (n[ii][j][k] + n[i][jj][k] + n[i][j][kk]) / 3.0;
}

void fsnorm::relax_neighb_vals(neighbs_t &n)
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
}

inline vector taxicab_normalized(const vector &v)
{
    double size = std::abs(v.x) + std::abs(v.y) + std::abs(v.z);
    return v * (1.0 / size);
}

vector fsnorm::get_normal(size_t i, size_t j, size_t k, size_t no)
{
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
        int selected_dir = -1;
        for (size_t i = 0; i < 3; ++i)
            if (c_ok[i])
            {
                m0[i] = std::abs(taxicab_normalized(c[i]).cmpnt[i]);
                if (selected_dir == -1 || m0[i] > m0[selected_dir])selected_dir = i;
            }

        // if (y_ok)
        //     if (std::abs(taxicab_normalized(y).cmpnt[selected_dir]) < m0[selected_dir])
        //         final = &y;
        //     else
        //         final = c + selected_dir;
        // else
        final = c + selected_dir;
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
    return *final;
}
}

