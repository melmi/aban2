/**
 * vof.h
 *
 * by Mohammad Elmi
 * Created on 10 Nov 2013
 */

#ifndef _VOF_H_
#define _VOF_H_

#include "common.h"
#include <array>
#include <tuple>
#include "domain.h"
#include "volreconst.h"

namespace aban2
{

enum class fullness
{
    full, empty, half
};

class vof
{
    typedef std::array<std::array<std::array<double, 3>, 3>, 3 > neighbs_t;
    struct molecule_t
    {
        double p, f1, b1, f2, b2;
    };
    static const double epsilon;

    domain *d;
    double vcell; //cell volume
    size_t start_dir = 0;
    double *mass;
    double *rhou[3];
    fullness *fullnesses;
    bool *on_interface;
    volreconst **reconsts;

    void set_fullnesses();
    inline bool is_empty(size_t i, size_t j, size_t k);
    bool is_on_interface(size_t i, size_t j, size_t k, size_t no);
    void detect_interfacial_cells();
    neighbs_t get_nighb_vals(size_t i, size_t j, size_t k);
    void relax_neighb_vals(neighbs_t &n);
    inline static double col_sum(double delta, double v1, double v2, double v3);
    inline static double dir_sign(double v1, double v2, double v3);
    inline static double get_column_grad(double p, double f, double b, double delta);
    bool create_column_candidate(vector &v, size_t dir, double sign, molecule_t molecule);
    void create_column_candidates(const neighbs_t &n, vector *candidates, bool *ok);
    bool create_young_candidate(const neighbs_t &n, vector &v);
    void set_normal(size_t i, size_t j, size_t k, size_t no);
    void calculate_normals();
    void create_reconsts();
    void delete_reconsts();
    inline std::tuple<double, vector> get_flux(mesh_row *row, size_t i, double udt, double ***grad_ustar);
    void advect_row(mesh_row *row, double ***grad_ustar);
    void correct_vofs(double *grad_uf_dir);
    void calculate_masses_from_vars();
    void calculate_vars_from_masses();
public:
    vof(domain *_d);
    ~vof();

    void advect();

    friend class vof_err;
};

const double vof::epsilon = 1e-8;

class vof_err
{
    inline static int is_inside(vector n, double alpha, vector x);
    inline static double err(vector n1, double alpha1, vector n2, double alpha2, double h);
public:
    static double compare(vof &v);
};
}

#endif